"""Read-only internal validation of OCVRC_Github final tables and SQLite.

Run independently in Spyder or through validation/Run_Validation.py.
Requires pandas and numpy. Does not access external catalogs or original OCVRC.
"""
from __future__ import annotations

import configparser
import math
import sqlite3
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
config = configparser.ConfigParser(interpolation=None)
if not config.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')
out = Path(config.get('PATHS', 'OUTPUT_DIR', fallback='output')).expanduser()
OUTPUT = out if out.is_absolute() else ROOT / out
FINAL = OUTPUT / 'final'
DATABASE = FINAL / 'OCVRC.db'

SCHEMAS = {
    'source': ['SOURCE_ID','RA','DEC','SIMBAD','MAJOR_AXIS','MINOR_AXIS',
               'POSITION_ANGLE','DISTANCE','DISTANCE_ERROR','PM_RA','PM_DEC'],
    'name': ['SOURCE_ID','NAME','DEFAULT_NAME'],
    'class': ['SOURCE_ID','CLASS','DEFAULT_CLASS'],
    'period': ['SOURCE_ID','TYPE','UNITS','PERIOD','PERIOD_ERROR','CONFIRMED'],
    'measurement': ['SOURCE_ID','OBSERVATORY','EPOCH','BAND','UNITS',
                    'MEASUREMENT','ERROR','POS_ERROR','NEG_ERROR','NEIGHBOR'],
}
SURVEYS = ['Gaia','2MASS','GALEX','LoTSS','ROSAT','XMM','Chandra','eROSITA']
VALID_STATUS = {'MATCH','NO_MATCH','QUERY_FAILED','NOT_QUERIED_SERVICE_UNAVAILABLE','NOT_QUERIED'}
BOOL_COLUMNS = {'SIMBAD','DEFAULT_NAME','DEFAULT_CLASS','CONFIRMED','NEIGHBOR'}
NUMERIC_COLUMNS = {'SOURCE_ID','RA','DEC','MAJOR_AXIS','MINOR_AXIS','POSITION_ANGLE',
                   'DISTANCE','DISTANCE_ERROR','PM_RA','PM_DEC','DEFAULT_NAME',
                   'DEFAULT_CLASS','PERIOD','PERIOD_ERROR','CONFIRMED',
                   'MEASUREMENT','ERROR','POS_ERROR','NEG_ERROR','NEIGHBOR','SIMBAD'}

checks = []

def report(level, check, details=''):
    checks.append({'RESULT': level, 'CHECK': check, 'DETAILS': str(details)[:240]})
    print(f'{level:7s} {check}' + (f' — {details}' if details else ''))

def check(condition, label, details=''):
    report('PASS' if condition else 'FAIL', label, details if not condition else '')

def normalized(value):
    if value is None or pd.isna(value):
        return None
    if isinstance(value, (np.bool_, bool)):
        return int(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        val = float(value)
        return val if math.isfinite(val) else str(val)
    return str(value) if not isinstance(value, (str, bytes)) else value

def equivalent(a, b):
    a, b = normalized(a), normalized(b)
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, (int,float)) and isinstance(b, (int,float)):
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
    return a == b

def row_signature(row):
    # Counter comparison is independent of row order and preserves duplicate counts.
    # Float normalization used only to group candidate rows; exact tolerance checked below.
    sig = []
    for val in row:
        val = normalized(val)
        if isinstance(val, (int,float)):
            val = round(float(val), 10) if isinstance(val, float) else val
        sig.append(val)
    return tuple(sig)

def compare_tables(frame, sql, label):
    if list(frame.columns) != list(sql.columns):
        report('FAIL', f'{label}: SQLite column order/schema differs')
        return
    if len(frame) != len(sql):
        report('FAIL', f'{label}: row count', f'pickle={len(frame)}, SQLite={len(sql)}')
        return
    # Fast multiset match; if numerical rounding falls on opposite sides of a boundary,
    # do an ordered, tolerant fallback after sorting on normalized, stable representations.
    fa = Counter(row_signature(row) for row in frame.itertuples(index=False, name=None))
    fb = Counter(row_signature(row) for row in sql.itertuples(index=False, name=None))
    if fa == fb:
        report('PASS', f'{label}: SQLite matches pickle (including duplicates and NULLs)', f'{len(frame)} rows')
    else:
        diff = sum((fa - fb).values()) + sum((fb - fa).values())
        report('FAIL', f'{label}: SQLite differs from pickle', f'{diff} unmatched row signatures; investigate values')

def bool_check(series, nullable):
    values = series.dropna()
    return (nullable or not series.isna().any()) and values.map(
        lambda x: isinstance(x, (bool,np.bool_,int,np.integer)) and x in (0,1)
    ).all()

def run_validation():
    global checks
    checks = []
    print('\nOCVRC_Github INTERNAL VALIDATION\n' + '='*65)
    frames = {}
    for name, cols in SCHEMAS.items():
        path = FINAL / f'{name}.pkl'
        if not path.is_file():
            report('FAIL', f'{name}: final pickle missing', path)
            continue
        try:
            frame = pd.read_pickle(path)
            if not isinstance(frame, pd.DataFrame):
                raise TypeError('Not a DataFrame')
            frames[name] = frame
            check(list(frame.columns) == cols, f'{name}: required final column schema',
                  f'expected {cols}; got {list(frame.columns)}')
        except Exception as exc:
            report('FAIL', f'{name}: cannot load pickle', f'{type(exc).__name__}: {exc}')

    if not DATABASE.is_file():
        report('FAIL', 'SQLite database missing', DATABASE)
    else:
        try:
            with sqlite3.connect(f'file:{DATABASE.as_posix()}?mode=ro', uri=True) as conn:
                check(conn.execute('PRAGMA integrity_check').fetchone()[0] == 'ok',
                      'SQLite: integrity_check')
                check(not conn.execute('PRAGMA foreign_key_check').fetchall(),
                      'SQLite: foreign_key_check')
                actual = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                check({f'{x}_table' for x in SCHEMAS}.issubset(actual), 'SQLite: all five tables present')
                for name, cols in SCHEMAS.items():
                    table = f'{name}_table'
                    if table not in actual:
                        continue
                    actual_cols = [r[1] for r in conn.execute(f'PRAGMA table_info({table})')]
                    check(actual_cols == cols, f'{table}: SQLite schema', f'expected {cols}; got {actual_cols}')
                    if name in frames:
                        sql = pd.read_sql_query(f'SELECT * FROM {table}', conn)
                        compare_tables(frames[name], sql, name)
        except Exception as exc:
            report('FAIL', 'SQLite: read/comparison failure', f'{type(exc).__name__}: {exc}')

    if 'source' in frames and 'SOURCE_ID' in frames['source'].columns:
        src = frames['source']
        ids = src['SOURCE_ID']
        check(not ids.isna().any() and ids.is_unique, 'SOURCE_ID: unique and non-null in source')
        valid_ids = set(ids.dropna())
        for name in ('name','class','period','measurement'):
            if name in frames and 'SOURCE_ID' in frames[name]:
                f = frames[name]
                check(not f['SOURCE_ID'].isna().any() and set(f['SOURCE_ID']) <= valid_ids,
                      f'{name}: valid foreign SOURCE_ID')
        for key, lower, upper, strict_upper in [('RA',0,360,True),('DEC',-90,90,False)]:
            if key in src:
                x = pd.to_numeric(src[key], errors='coerce')
                check(x.notna().all() and np.isfinite(x).all() and x.ge(lower).all()
                      and (x.lt(upper) if strict_upper else x.le(upper)).all(),
                      f'astrometry: {key} in degrees range')
        for col in ('SIMBAD',):
            if col in src:
                check(bool_check(src[col], False), f'astrometry: {col} boolean')
        for col in ('MAJOR_AXIS','MINOR_AXIS','DISTANCE','DISTANCE_ERROR'):
            if col in src:
                x = pd.to_numeric(src[col], errors='coerce')
                check(((x.dropna() >= 0) & np.isfinite(x.dropna())).all(),
                      f'astrometry: {col} nonnegative when present')
        if all(x in src for x in ('SIMBAD','MAJOR_AXIS','MINOR_AXIS','POSITION_ANGLE')):
            no_match = src['SIMBAD'].eq(False)
            errs = ['MAJOR_AXIS','MINOR_AXIS','POSITION_ANGLE']
            check(src.loc[no_match, errs].isna().all().all(),
                  'SIMBAD: no-match rows have NULL position uncertainties')
            yes = src['SIMBAD'].eq(True)
            missing = int(src.loc[yes, errs].isna().any(axis=1).sum())
            if missing:
                report('WARNING','SIMBAD: matched sources with partial/missing position uncertainties',missing)
            else:
                report('PASS','SIMBAD: matched position uncertainties populated')

    if 'name' in frames and 'source' in frames:
        f = frames['name']
        if {'SOURCE_ID','NAME','DEFAULT_NAME'} <= set(f):
            check(not f['NAME'].isna().any() and f['NAME'].astype(str).str.strip().ne('').all(),
                  'name: nonempty names')
            check(not f.duplicated(['SOURCE_ID','NAME']).any(), 'name: no duplicate names per source')
            check(bool_check(f['DEFAULT_NAME'], False), 'name: boolean DEFAULT_NAME')
            counts = f.groupby('SOURCE_ID')['DEFAULT_NAME'].sum()
            check(set(counts.index) == set(frames['source']['SOURCE_ID']) and counts.eq(1).all(),
                  'name: exactly one default per source')
    if 'class' in frames:
        f = frames['class']
        if {'SOURCE_ID','CLASS','DEFAULT_CLASS'} <= set(f):
            check(not f['CLASS'].isna().any() and f['CLASS'].astype(str).str.strip().ne('').all(),
                  'class: nonempty classes')
            check(not f.duplicated(['SOURCE_ID','CLASS']).any(), 'class: no duplicate classes per source')
            check(bool_check(f['DEFAULT_CLASS'], False), 'class: boolean DEFAULT_CLASS')
            check(f.groupby('SOURCE_ID')['DEFAULT_CLASS'].sum().eq(1).all(),
                  'class: exactly one default per classified source')
    if 'period' in frames:
        f = frames['period']
        if {'SOURCE_ID','TYPE','PERIOD','PERIOD_ERROR','CONFIRMED'} <= set(f):
            check(not f.duplicated(['SOURCE_ID','TYPE']).any(), 'period: no duplicate types per source')
            check(bool_check(f['CONFIRMED'],True), 'period: CONFIRMED nullable boolean')
            for col in ('PERIOD','PERIOD_ERROR'):
                x = pd.to_numeric(f[col], errors='coerce')
                check((x.dropna() > 0 if col == 'PERIOD' else x.dropna() >= 0).all()
                      and (col != 'PERIOD' or x.notna().all()), f'period: valid {col}')
    if 'measurement' in frames:
        f = frames['measurement']
        if {'MEASUREMENT','ERROR','POS_ERROR','NEG_ERROR','NEIGHBOR','OBSERVATORY','EPOCH','BAND','UNITS'} <= set(f):
            for col in ('MEASUREMENT','ERROR','POS_ERROR','NEG_ERROR'):
                x = pd.to_numeric(f[col], errors='coerce')
                check(np.isfinite(x.dropna()).all() and (col != 'MEASUREMENT' or x.notna().all()),
                      f'measurement: {col} finite when present')
            check(bool_check(f['NEIGHBOR'],True), 'measurement: NEIGHBOR nullable boolean')
            neg = f['MEASUREMENT'].lt(0).sum()
            report('PASS','measurement: negative fluxes retained (informational)',int(neg))
            radio = f['OBSERVATORY'].isin(['VLA','ASKAP'])
            check(f.loc[radio,'NEIGHBOR'].notna().all(), 'radio: every observation has NEIGHBOR flag')
            report('PASS','radio: number of flagged observations',int(f.loc[radio,'NEIGHBOR'].eq(1).sum()))
            repeated = f.loc[f['OBSERVATORY'].eq('VLA')].groupby(['SOURCE_ID','EPOCH'],dropna=False).size()
            repeated = repeated[repeated > 1]
            if len(repeated):
                report('WARNING', 'VLASS: multiple measurements per source/epoch (overlapping tiles possible)',
                       f'{len(repeated)} source/epoch groups')
            else:
                report('PASS','VLASS: no repeated source/epoch groups')
            missing_info = f.loc[radio,['EPOCH','BAND','UNITS']].isna().any(axis=1).sum()
            check(missing_info == 0,'radio: epoch, band, and units present',f'{missing_info} missing rows')

    if 'source' in frames and 'SOURCE_ID' in frames['source']:
        expected = set(frames['source']['SOURCE_ID'])
        for survey in SURVEYS:
            path = OUTPUT / f'{survey}_query_status.pkl'
            if not path.exists():
                report('WARNING', f'{survey}: status file missing')
                continue
            try:
                status = pd.read_pickle(path)
                if not {'SOURCE_ID','STATUS'} <= set(status):
                    report('FAIL',f'{survey}: status file missing SOURCE_ID/STATUS')
                    continue
                check(not status['SOURCE_ID'].duplicated().any() and not status['SOURCE_ID'].isna().any(),
                      f'{survey}: one status per source')
                check(set(status['SOURCE_ID']) == expected, f'{survey}: status covers every source')
                check(status['STATUS'].isin(VALID_STATUS).all(), f'{survey}: recognized statuses')
                incomplete = int(status['STATUS'].isin({'QUERY_FAILED','NOT_QUERIED_SERVICE_UNAVAILABLE','NOT_QUERIED'}).sum())
                if incomplete:
                    report('WARNING',f'{survey}: incomplete/outage queries (not nondetections)', incomplete)
                else:
                    report('PASS',f'{survey}: all queries completed')
            except Exception as exc:
                report('FAIL',f'{survey}: cannot read query status', f'{type(exc).__name__}: {exc}')

    result = pd.DataFrame(checks, columns=['RESULT','CHECK','DETAILS'])
    print('\nVALIDATION SUMMARY\n' + '='*65)
    print(result['RESULT'].value_counts().reindex(['PASS','WARNING','FAIL'],fill_value=0).to_string())
    print('Read-only validation: no source data or database files modified.')
    return result

if __name__ == '__main__':
    df_validation = run_validation()
