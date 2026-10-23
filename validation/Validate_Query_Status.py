"""OCVRC_Github: read-only validation of external catalogue query status.

Run directly in Spyder. Uses only OCVRC_Github/output files; does not
query remote services and does not modify data or SQLite.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / 'output'
FINAL = OUTPUT / 'final'
SURVEYS = ('Gaia', '2MASS', 'GALEX', 'LoTSS', 'ROSAT', 'XMM', 'Chandra', 'eROSITA')
DONE = {'MATCH', 'NO_MATCH'}
PENDING = {'QUERY_FAILED', 'NOT_QUERIED_SERVICE_UNAVAILABLE', 'NOT_QUERIED'}
VALID = DONE | PENDING
FIELDS = ['SOURCE_ID', 'OBSERVATORY', 'EPOCH', 'BAND', 'UNITS',
          'MEASUREMENT', 'ERROR', 'POS_ERROR', 'NEG_ERROR', 'NEIGHBOR']

checks = []
summary = []


def record(result, name, detail=''):
    checks.append({'RESULT': result, 'CHECK': name, 'DETAIL': str(detail)})
    print(f'{result:<7} {name}' + (f' — {detail}' if detail != '' else ''))


def check(ok, name, detail=''):
    record('PASS' if bool(ok) else 'FAIL', name, detail)


def read_frame(path, label):
    if not path.is_file():
        record('FAIL', f'{label}: file exists', str(path))
        return None
    try:
        data = pd.read_pickle(path)
    except Exception as exc:
        record('FAIL', f'{label}: readable pickle', f'{type(exc).__name__}: {exc}')
        return None
    if not isinstance(data, pd.DataFrame):
        record('FAIL', f'{label}: DataFrame', type(data).__name__)
        return None
    return data


def ids_equal(frame, source_ids):
    return (frame['SOURCE_ID'].notna().all()
            and not frame['SOURCE_ID'].duplicated().any()
            and set(frame['SOURCE_ID']) == source_ids)


print('\nOCVRC_Github QUERY STATUS VALIDATION')
print('=' * 70)

sources = read_frame(FINAL / 'source.pkl', 'Final source')
final_meas = read_frame(FINAL / 'measurement.pkl', 'Final measurements')

if sources is not None:
    required = {'SOURCE_ID', 'RA', 'DEC'}
    check(required <= set(sources.columns), 'Final source: required identifier/position columns')
    if required <= set(sources.columns):
        check(sources.SOURCE_ID.notna().all() and sources.SOURCE_ID.is_unique,
              'Final source: IDs nonmissing and unique')
        source_ids = set(sources['SOURCE_ID'].dropna())
    else:
        source_ids = None
else:
    source_ids = None

if final_meas is not None:
    check({'SOURCE_ID', 'OBSERVATORY', 'BAND', 'MEASUREMENT'} <= set(final_meas),
          'Final measurements: required comparison columns')

for survey in SURVEYS:
    print(f'\n{survey.upper()}')
    print('-' * 70)
    status = read_frame(OUTPUT / f'{survey}_query_status.pkl', f'{survey} status')
    measurements = read_frame(OUTPUT / f'{survey}_measurements.pkl', f'{survey} measurements')
    failures_path = OUTPUT / f'{survey}_failures.pkl'
    failures = read_frame(failures_path, f'{survey} failures')
    result = {'SURVEY': survey, 'MATCH': 0, 'NO_MATCH': 0, 'QUERY_FAILED': 0,
              'NOT_QUERIED_SERVICE_UNAVAILABLE': 0, 'NOT_QUERIED': 0,
              'MEASUREMENT_ROWS': 0, 'FAILURE_ROWS': 0}

    if status is None:
        summary.append(result)
        continue
    if not {'SOURCE_ID', 'STATUS'} <= set(status):
        record('FAIL', f'{survey}: status has SOURCE_ID and STATUS')
        summary.append(result)
        continue

    check(status['SOURCE_ID'].notna().all() and status['SOURCE_ID'].is_unique,
          f'{survey}: exactly one status record per SOURCE_ID')
    if source_ids is not None:
        check(ids_equal(status, source_ids), f'{survey}: covers all final sources',
              f'{len(status)} records vs {len(source_ids)} sources')
    bad_status = status.loc[~status.STATUS.isin(VALID), 'STATUS']
    check(bad_status.empty, f'{survey}: recognized status codes',
          bad_status.astype(str).unique().tolist() if not bad_status.empty else '')
    counts = status.STATUS.value_counts()
    for name in VALID:
        result[name] = int(counts.get(name, 0))

    # Each MATCH/NO_MATCH represents a successful query. In particular,
    # a NO_MATCH is NOT the same as QUERY_FAILED or service unavailable.
    incomplete = int(status.STATUS.isin(PENDING).sum())
    record('WARNING' if incomplete else 'PASS', f'{survey}: all queries completed',
           f'{incomplete} incomplete/outage sources (not nondetections)' if incomplete else '')

    # The checkpoint stores the input sky position for each source. Reconcile
    # with current final-source positions to catch stale checkpoint reuse.
    if sources is not None and {'RA', 'DEC'} <= set(status):
        pos = status[['SOURCE_ID', 'RA', 'DEC']].merge(
            sources[['SOURCE_ID', 'RA', 'DEC']], on='SOURCE_ID',
            how='inner', suffixes=('_checkpoint', '_final'))
        valid_pos = (np.isclose(pd.to_numeric(pos.RA_checkpoint, errors='coerce'),
                                pd.to_numeric(pos.RA_final, errors='coerce'), rtol=0, atol=1e-8)
                     & np.isclose(pd.to_numeric(pos.DEC_checkpoint, errors='coerce'),
                                  pd.to_numeric(pos.DEC_final, errors='coerce'), rtol=0, atol=1e-8))
        check(bool(valid_pos.all()), f'{survey}: status coordinates match final source',
              f'{int((~valid_pos).sum())} differing positions')

    if 'ATTEMPTS' in status:
        tries = pd.to_numeric(status.ATTEMPTS, errors='coerce')
        check(tries.notna().all() and (tries >= 0).all()
              and np.isclose(tries, np.round(tries)).all(),
              f'{survey}: valid query attempt counts')
        attempted = status.STATUS.isin(DONE | {'QUERY_FAILED'})
        check(bool((tries.loc[attempted] >= 1).all()),
              f'{survey}: attempted queries have positive attempt count')

    if failures is not None:
        result['FAILURE_ROWS'] = len(failures)
        if {'SOURCE_ID', 'STATUS'} <= set(failures):
            check(failures.SOURCE_ID.notna().all() and failures.SOURCE_ID.is_unique,
                  f'{survey}: failures have unique source IDs')
            expected = set(status.loc[status.STATUS.isin(
                {'QUERY_FAILED', 'NOT_QUERIED_SERVICE_UNAVAILABLE'}), 'SOURCE_ID'])
            actual = set(failures.SOURCE_ID)
            check(actual == expected, f'{survey}: failure file agrees with statuses',
                  f'{len(actual)} recorded vs {len(expected)} failed/skipped')
            check(failures.STATUS.isin(PENDING).all(),
                  f'{survey}: failure file has no successful queries')
        else:
            record('FAIL', f'{survey}: failures file has SOURCE_ID and STATUS')

    if measurements is not None:
        result['MEASUREMENT_ROWS'] = len(measurements)
        missing = [c for c in FIELDS if c not in measurements]
        check(not missing, f'{survey}: measurement output schema', missing if missing else '')
        if 'SOURCE_ID' in measurements:
            completed = set(status.loc[status.STATUS.eq('MATCH'), 'SOURCE_ID'])
            measured = set(measurements.SOURCE_ID.dropna())
            check(measurements.SOURCE_ID.notna().all() and measured <= completed,
                  f'{survey}: measurements belong only to MATCH sources',
                  f'{len(measured)} measured sources, {len(completed)} matched sources')
            no_rows = int(status.STATUS.eq('NO_MATCH').sum())
            record('PASS', f'{survey}: NO_MATCH is a valid empty measurement result',
                   f'{no_rows} nondetections; {len(measurements)} measurement rows total')
            # A MATCH can legitimately yield no measured flux in usable bands.
            no_flux_matches = completed - measured
            if no_flux_matches:
                record('WARNING', f'{survey}: matches without usable measurement rows',
                       f'{len(no_flux_matches)} sources (may be legitimate)')
        if 'MEASUREMENT' in measurements:
            v = pd.to_numeric(measurements.MEASUREMENT, errors='coerce')
            check(v.notna().all() and bool(np.isfinite(v).all()),
                  f'{survey}: cached measurements numeric and finite')

    # Optional metadata produced by the shared query manager.
    meta_path = OUTPUT / f'{survey}_query_meta.json'
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding='utf-8'))
            check(meta.get('survey') == survey and meta.get('source_count') == len(status),
                  f'{survey}: checkpoint metadata survey and source count')
            check(isinstance(meta.get('fingerprint'), str)
                  and len(meta['fingerprint']) == 64,
                  f'{survey}: checkpoint fingerprint recorded')
        except Exception as exc:
            record('FAIL', f'{survey}: checkpoint metadata readable',
                   f'{type(exc).__name__}: {exc}')
    else:
        record('WARNING', f'{survey}: checkpoint metadata missing',
               'Resume safety fingerprint not independently verifiable')

    # Independent output cross-check: number of rows in survey output should
    # match the number included for that survey in the final table.
    # Compare by normalized OBSERVATORY; this is distinct per survey here.
    labels = {'Gaia': {'gaia'}, '2MASS': {'2mass'}, 'GALEX': {'galex'},
              'LoTSS': {'lofar'}, 'ROSAT': {'rosat'}, 'XMM': {'xmm'},
              'Chandra': {'chandra'}, 'eROSITA': {'erosita'}}
    if (final_meas is not None and {'OBSERVATORY'} <= set(final_meas)
            and measurements is not None):
        final_count = int(final_meas.OBSERVATORY.astype(str).str.casefold()
                          .isin(labels[survey]).sum())
        check(final_count == len(measurements),
              f'{survey}: final measurement row count matches survey output',
              f'final={final_count}, cached={len(measurements)}')

    summary.append(result)


df_query_status_validation = pd.DataFrame(checks, columns=['RESULT', 'CHECK', 'DETAIL'])
df_query_survey_summary = pd.DataFrame(summary)
print('\nQUERY STATUS BY SURVEY')
print('=' * 70)
print(df_query_survey_summary.to_string(index=False))
print('\nQUERY STATUS VALIDATION SUMMARY')
print('=' * 70)
print(df_query_status_validation.RESULT.value_counts().reindex(
    ['PASS', 'WARNING', 'FAIL'], fill_value=0).to_string())
print('Read-only: no data files modified; no network queries performed.')
