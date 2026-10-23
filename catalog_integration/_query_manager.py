"""Resumable catalogue-query runner for OCVRC_Github.

All statuses are per SOURCE_ID. Failed queries never masquerade as no matches.
Checkpoint writes use atomic replacement; the source fingerprint blocks stale reuse.
"""
from pathlib import Path
import hashlib
import json
import os
import time
import uuid
import re
import pandas as pd

DONE = {'MATCH', 'NO_MATCH'}
PENDING = {'QUERY_FAILED', 'NOT_QUERIED_SERVICE_UNAVAILABLE'}
STATUS_COLUMNS = ['SOURCE_ID', 'STATUS', 'ATTEMPTS', 'ERROR_MESSAGE', 'RA', 'DEC']


def atomic_pickle(frame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f'.{path.name}.{uuid.uuid4().hex}.tmp')
    try:
        frame.to_pickle(temp)
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def signature(df, survey, radius):
    # Identifies the entire input catalogue, including source order and positions.
    rows = [(str(r.SOURCE_ID), float(r.RA), float(r.DEC)) for r in df.itertuples()]
    data = json.dumps([survey, float(radius), rows], separators=(',', ':'))
    return hashlib.sha256(data.encode()).hexdigest()


def run_queries(survey, sources, query_one, output_dir, *, radius, retries=2,
                delay=3., failure_limit=5, save_every=10,
                start=0, end=None, resume=True):
    """query_one(source_row, radius) -> (status, result_dict, list_of_measurement_dicts).

    status is MATCH or NO_MATCH; exceptions trigger retry/failure logging.
    On resume, completed rows are reused and failed/skipped rows retried.
    """
    if sources.SOURCE_ID.isna().any() or sources.SOURCE_ID.duplicated().any():
        raise ValueError('Input SOURCE_ID must be unique and non-null')
    if not {'RA', 'DEC'}.issubset(sources.columns):
        raise ValueError('Sources must contain RA and DEC')
    if failure_limit < 1 or save_every < 1 or retries < 0:
        raise ValueError('Invalid query-manager settings')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    meta_path = output_dir / f'{survey}_query_meta.json'
    status_path = output_dir / f'{survey}_query_status.pkl'
    results_path = output_dir / f'{survey}.pkl'
    measurements_path = output_dir / f'{survey}_measurements.pkl'
    failures_path = output_dir / f'{survey}_failures.pkl'
    fingerprint = signature(sources, survey, radius)
    status, results, measurements = {}, {}, {}
    if resume and meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding='utf-8'))
        if meta.get('fingerprint') != fingerprint:
            raise ValueError(f'{survey}: input sources/coordinates or search radius changed. '
                             'Move previous checkpoint files aside or set RESUME = FALSE.')
        if not all(p.exists() for p in (status_path, results_path, measurements_path)):
            raise ValueError(f'{survey}: incomplete checkpoint files; do not resume')
        for _, row in pd.read_pickle(status_path).iterrows():
            status[row['SOURCE_ID']] = row.to_dict()
        for _, row in pd.read_pickle(results_path).iterrows():
            results[row['SOURCE_ID']] = row.to_dict()
        frame = pd.read_pickle(measurements_path)
        if 'EPOCH_CAT' in frame.columns and 'EPOCH' not in frame.columns:
            frame = frame.rename(columns={'EPOCH_CAT': 'EPOCH'})
        for sid, sub in frame.groupby('SOURCE_ID'):
            measurements[sid] = sub.to_dict('records')
    elif resume and any(p.exists() for p in (meta_path, status_path)):
        raise ValueError(f'{survey}: existing checkpoint has no complete metadata; '
                         'move old checkpoint files aside or set RESUME = FALSE')

    source_ids = sources.SOURCE_ID.tolist()
    # Build explicit record for every source, including those outside selected bounds.
    for i, row in sources.iterrows():
        sid = row.SOURCE_ID
        status.setdefault(sid, {'SOURCE_ID': sid, 'STATUS': 'NOT_QUERIED',
                                'ATTEMPTS': 0, 'ERROR_MESSAGE': '',
                                'RA': float(row.RA), 'DEC': float(row.DEC)})
    stop = len(sources) if end is None or end <= 0 else min(len(sources), end)
    if not 0 <= start <= stop <= len(sources):
        raise ValueError('Invalid query row bounds')

    def save():
        status_df = pd.DataFrame([status[s] for s in source_ids], columns=STATUS_COLUMNS)
        result_df = pd.DataFrame([results[s] for s in source_ids if s in results])
        records = [m for s in source_ids for m in measurements.get(s, [])]
        cols = ['SOURCE_ID', 'OBSERVATORY', 'EPOCH', 'BAND', 'UNITS',
                'MEASUREMENT', 'ERROR', 'POS_ERROR', 'NEG_ERROR', 'NEIGHBOR']
        measure_df = pd.DataFrame(records, columns=cols)
        failure_df = status_df[status_df.STATUS.isin(PENDING)].copy()
        # Metadata is written last, so its fingerprint only points to full datasets.
        atomic_pickle(status_df, status_path)
        atomic_pickle(result_df, results_path)
        atomic_pickle(measure_df, measurements_path)
        atomic_pickle(failure_df, failures_path)
        temp = meta_path.with_name(meta_path.name + '.tmp')
        temp.write_text(json.dumps({'fingerprint': fingerprint, 'survey': survey,
                                    'radius_arcsec': radius, 'source_count': len(sources)}, indent=2),
                        encoding='utf-8')
        os.replace(temp, meta_path)
        return result_df, measure_df, failure_df, status_df

    consecutive_failures = 0
    service_unavailable = False
    processed = 0
    for pos in range(start, stop):
        row = sources.iloc[pos]
        sid = row.SOURCE_ID
        if resume and status[sid]['STATUS'] in DONE:
            continue
        if service_unavailable or consecutive_failures >= failure_limit:
            status[sid].update(STATUS='NOT_QUERIED_SERVICE_UNAVAILABLE',
                               ERROR_MESSAGE=('HTTP 503: Service Unavailable' if service_unavailable else
                                              f'Circuit breaker: {failure_limit} consecutive failures'))
            continue
        error = ''
        for attempt in range(retries + 1):
            status[sid]['ATTEMPTS'] += 1
            try:
                outcome, record, records = query_one(row, radius)
                if outcome not in DONE:
                    raise ValueError(f'Unexpected query outcome {outcome!r}')
                record = dict(record)
                record['SOURCE_ID'] = sid
                results[sid] = record
                measurements[sid] = [dict(m, SOURCE_ID=sid) for m in records]
                status[sid].update(STATUS=outcome, ERROR_MESSAGE='')
                consecutive_failures = 0
                print(f'{survey} row {pos} SOURCE_ID={sid}: {outcome}')
                break
            except Exception as exc:
                raw = str(exc)
                # A definite 503 means the service is offline, not a source-level failure.
                # Stop retrying this source AND stop submitting requests for other sources.
                if re.search(r'\b503\b', raw):
                    service_unavailable = True
                    error = 'HTTP 503: Service Unavailable'
                    break
                error = f'{type(exc).__name__}: {raw.splitlines()[0][:160] if raw else type(exc).__name__}'
                if attempt < retries:
                    time.sleep(delay * (2 ** attempt))
        else:
            pass
        # This also handles an explicit 503 that broke out of the retry loop.
        if status[sid]['STATUS'] not in DONE:
            results.pop(sid, None)
            measurements.pop(sid, None)
            status[sid].update(STATUS='QUERY_FAILED', ERROR_MESSAGE=error)
            consecutive_failures += 1
            print(f'{survey} row {pos} SOURCE_ID={sid}: QUERY_FAILED ({error})')
            if service_unavailable:
                print(f'{survey}: HTTP 503 received; skipping remaining pending queries.')
        processed += 1
        if processed % save_every == 0:
            save()
    output = save()
    counts = output[3]['STATUS'].value_counts().to_dict()
    print(f'{survey} QUERY SUMMARY: {counts}; measurement rows={len(output[1])}')
    return output
