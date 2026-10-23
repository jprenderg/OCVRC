"""Build the five final OCVRC tables from separately tested stage outputs.

Run in Spyder current console to inspect df_source, df_name, df_class,
df_period, df_measurement, and df_survey_summary.

Writes output/final/*.pkl; preserves intermediate input files.
"""
from pathlib import Path
import configparser
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
config = configparser.ConfigParser(interpolation=None)
if not config.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')

def project_path(option, fallback):
    value = Path(config.get('PATHS', option, fallback=fallback))
    return value if value.is_absolute() else ROOT / value

OUTPUT = project_path('OUTPUT_DIR', 'output')
FINAL = OUTPUT / 'final'
SOURCE_COLS = ['SOURCE_ID', 'RA', 'DEC', 'SIMBAD',
               'MAJOR_AXIS', 'MINOR_AXIS', 'POSITION_ANGLE',
               'DISTANCE', 'DISTANCE_ERROR', 'PM_RA', 'PM_DEC']
NAME_COLS = ['SOURCE_ID', 'NAME', 'DEFAULT_NAME']
CLASS_COLS = ['SOURCE_ID', 'CLASS', 'DEFAULT_CLASS']
PERIOD_COLS = ['SOURCE_ID', 'TYPE', 'UNITS', 'PERIOD', 'PERIOD_ERROR', 'CONFIRMED']
MEAS_COLS = ['SOURCE_ID', 'OBSERVATORY', 'EPOCH', 'BAND', 'UNITS',
             'MEASUREMENT', 'ERROR', 'POS_ERROR', 'NEG_ERROR', 'NEIGHBOR']
SURVEYS = ['Gaia', '2MASS', 'GALEX', 'LoTSS', 'ROSAT', 'XMM', 'Chandra', 'eROSITA']


def read_table(path):
    if not path.is_file():
        raise FileNotFoundError(f'Required input missing: {path}')
    value = pd.read_pickle(path)
    if not isinstance(value, pd.DataFrame):
        raise TypeError(f'Expected pandas DataFrame: {path}')
    value = value.copy()
    # Read cached survey outputs produced under the previous schema.
    if 'EPOCH_CAT' in value.columns and 'EPOCH' not in value.columns:
        value = value.rename(columns={'EPOCH_CAT': 'EPOCH'})
    return value


def validate_source_ids(table, label, known_ids, allow_subset=True):
    if 'SOURCE_ID' not in table.columns:
        raise ValueError(f'{label}: missing SOURCE_ID')
    if table['SOURCE_ID'].isna().any():
        raise ValueError(f'{label}: null SOURCE_ID')
    if not pd.to_numeric(table['SOURCE_ID'], errors='coerce').notna().all():
        raise ValueError(f'{label}: nonnumeric SOURCE_ID')
    if not (table['SOURCE_ID'].astype(float) % 1 == 0).all():
        raise ValueError(f'{label}: noninteger SOURCE_ID')
    table['SOURCE_ID'] = table['SOURCE_ID'].astype('int64')
    ids = set(table['SOURCE_ID'])
    if ids - known_ids:
        raise ValueError(f'{label}: SOURCE_ID values absent from source table: {sorted(ids-known_ids)[:20]}')
    if not allow_subset and ids != known_ids:
        raise ValueError(f'{label}: SOURCE_ID coverage mismatch; missing {sorted(known_ids-ids)[:20]}')


# Base tables are the authoritative result of the most recent ingestion.
df_source_initial = read_table(OUTPUT / 'source.pkl')
# Support ingestion outputs generated before the public column renaming.
df_source_initial = df_source_initial.rename(columns={
    'COO_ERR_MAJ': 'MAJOR_AXIS', 'COO_ERR_MIN': 'MINOR_AXIS',
    'COO_ERR_ANGLE': 'POSITION_ANGLE'})
for column in ['SOURCE_ID', 'RA', 'DEC']:
    if column not in df_source_initial:
        raise ValueError(f'source.pkl: missing {column}')
if df_source_initial['SOURCE_ID'].isna().any() or df_source_initial['SOURCE_ID'].duplicated().any():
    raise ValueError('source.pkl: null or duplicate SOURCE_ID')
df_source_initial['SOURCE_ID'] = df_source_initial['SOURCE_ID'].astype('int64')
known_ids = set(df_source_initial['SOURCE_ID'])
if not known_ids:
    raise ValueError('source.pkl is empty')
for key, low, high in [('RA', 0, 360), ('DEC', -90, 90)]:
    values = pd.to_numeric(df_source_initial[key], errors='coerce')
    if not values.notna().all() or not np.isfinite(values).all():
        raise ValueError(f'source.pkl has invalid {key}')
    good = values.ge(low) & (values.lt(high) if key == 'RA' else values.le(high))
    if not good.all():
        raise ValueError(f'source.pkl has {key} outside expected degrees range')

df_name = read_table(OUTPUT / 'name.pkl')
df_class = read_table(OUTPUT / 'class.pkl')
df_period = read_table(OUTPUT / 'period.pkl')
for label, table in [('name', df_name), ('class', df_class), ('period', df_period)]:
    validate_source_ids(table, label, known_ids)

if df_name.duplicated(['SOURCE_ID', 'NAME']).any():
    raise ValueError('Duplicate SOURCE_ID/NAME records')
if set(df_name['SOURCE_ID']) != known_ids:
    raise ValueError('Every source needs a name')
if not df_name.groupby('SOURCE_ID')['DEFAULT_NAME'].sum().eq(1).all():
    raise ValueError('Each source must have exactly one DEFAULT_NAME')
if df_name['NAME'].isna().any() or df_name['NAME'].astype(str).str.strip().eq('').any():
    raise ValueError('Blank NAME value')
if df_class.duplicated(['SOURCE_ID', 'CLASS']).any():
    raise ValueError('Duplicate SOURCE_ID/CLASS records')
if not df_class.empty and not df_class.groupby('SOURCE_ID')['DEFAULT_CLASS'].sum().eq(1).all():
    raise ValueError('Every classified source must have exactly one DEFAULT_CLASS')
if df_period.duplicated(['SOURCE_ID', 'TYPE']).any():
    raise ValueError('Duplicate SOURCE_ID/period TYPE records')

# Gaia's enriched source table must belong to the current ingestion run.
df_gaia_source = read_table(OUTPUT / 'source_gaia.pkl')
validate_source_ids(df_gaia_source, 'source_gaia', known_ids, allow_subset=False)
if df_gaia_source['SOURCE_ID'].duplicated().any():
    raise ValueError('source_gaia.pkl contains duplicate SOURCE_ID')
for field in ['DISTANCE', 'DISTANCE_ERROR', 'PM_RA', 'PM_DEC']:
    if field not in df_gaia_source:
        raise ValueError(f'source_gaia.pkl missing {field}')
# Guard against stale Gaia outputs after a new CSV was ingested.
coords = df_source_initial[['SOURCE_ID', 'RA', 'DEC']].merge(
    df_gaia_source[['SOURCE_ID', 'RA', 'DEC']], on='SOURCE_ID', suffixes=('_INPUT', '_GAIA'),
    validate='one_to_one')
for key in ['RA', 'DEC']:
    if not np.allclose(coords[f'{key}_INPUT'], coords[f'{key}_GAIA'], rtol=0, atol=1e-9):
        raise ValueError('Gaia source file has stale source coordinates; rerun Get_Gaia_DR3.py')

# Construct the published source schema, never exporting ingestion metadata.
# SIMBAD_FLAG records whether a 2-arcsec SIMBAD match was found.
if 'SIMBAD_FLAG' not in df_source_initial.columns:
    raise ValueError('source.pkl: missing SIMBAD_FLAG')
flags = df_source_initial['SIMBAD_FLAG']
if flags.isna().any():
    raise ValueError('source.pkl: SIMBAD_FLAG contains null values')
# Do not use astype(bool) on arbitrary strings; e.g. bool('False') is True.
if not flags.map(lambda x: isinstance(x, (bool, np.bool_))).all():
    raise ValueError('source.pkl: SIMBAD_FLAG must contain True/False booleans')
# These columns must be collected during SIMBAD ingestion; never invent values.
POS_FIELDS = ['MAJOR_AXIS', 'MINOR_AXIS', 'POSITION_ANGLE']
missing_pos = [c for c in POS_FIELDS if c not in df_source_initial.columns]
if missing_pos:
    raise ValueError(
        f'source.pkl missing SIMBAD coordinate error fields {missing_pos}. '
        'Replace Build_Source_Name_Class.py with the updated version and '
        'rerun source ingestion before building the final database.')
df_source = df_source_initial[['SOURCE_ID', 'RA', 'DEC', *POS_FIELDS]].copy()
for col in POS_FIELDS:
    df_source[col] = pd.to_numeric(df_source[col], errors='raise')
    if not np.isfinite(df_source[col].dropna()).all():
        raise ValueError(f'source.pkl contains infinite {col} values')
# No SIMBAD counterpart means no SIMBAD coordinate error ellipse.
for col in POS_FIELDS:
    df_source.loc[~flags.astype(bool), col] = np.nan
df_source['SIMBAD'] = flags.astype(bool).to_numpy()
df_source = df_source.merge(
    df_gaia_source[['SOURCE_ID', 'DISTANCE', 'DISTANCE_ERROR', 'PM_RA', 'PM_DEC']],
    on='SOURCE_ID', how='left', validate='one_to_one')
for field in ['DISTANCE', 'DISTANCE_ERROR', 'PM_RA', 'PM_DEC']:
    df_source[field] = pd.to_numeric(df_source[field], errors='raise')
df_source = df_source[SOURCE_COLS].copy()

# Publish only the defined columns for the other four database tables.
for label, table, columns in [('name', df_name, NAME_COLS),
                              ('class', df_class, CLASS_COLS),
                              ('period', df_period, PERIOD_COLS)]:
    missing = sorted(set(columns) - set(table.columns))
    if missing:
        raise ValueError(f'{label}.pkl: missing required columns {missing}')
df_name = df_name[NAME_COLS].copy()
df_class = df_class[CLASS_COLS].copy()
df_period = df_period[PERIOD_COLS].copy()

# Preserve all finite radio fits, including low-S/N and NEIGHBOR=1.
# This is the radio-stage measurement.pkl, not this builder's final output.
df_radio = read_table(OUTPUT / 'measurement.pkl')
frames = []
summary = []

def add_frame(label, table):
    if not set(MEAS_COLS).issubset(table.columns):
        raise ValueError(f'{label}: missing columns {sorted(set(MEAS_COLS)-set(table.columns))}')
    frame = table[MEAS_COLS].copy()
    validate_source_ids(frame, label, known_ids)
    frame['MEASUREMENT'] = pd.to_numeric(frame['MEASUREMENT'], errors='coerce')
    for col in ['ERROR', 'POS_ERROR', 'NEG_ERROR', 'NEIGHBOR']:
        frame[col] = pd.to_numeric(frame[col], errors='coerce')
    if (~np.isfinite(frame['MEASUREMENT'].dropna())).any():
        raise ValueError(f'{label}: nonfinite measurements')
    frame = frame.dropna(subset=['MEASUREMENT']).reset_index(drop=True)
    frames.append(frame)
    summary.append({'SURVEY': label, 'STATUS': 'Included', 'MEASUREMENT_ROWS': len(frame)})

add_frame('VLASS/RACS', df_radio)
for survey in SURVEYS:
    path = OUTPUT / f'{survey}_measurements.pkl'
    status_path = OUTPUT / f'{survey}_query_status.pkl'
    if not path.exists():
        summary.append({'SURVEY': survey, 'STATUS': 'NOT RUN', 'MEASUREMENT_ROWS': 0})
        continue
    frame = read_table(path)
    add_frame(survey, frame)
    if status_path.exists():
        statuses = read_table(status_path)
        validate_source_ids(statuses, survey + '_status', known_ids)
        counts = statuses['STATUS'].value_counts()
        failures = int(counts.get('QUERY_FAILED', 0))
        skipped = int(counts.get('NOT_QUERIED_SERVICE_UNAVAILABLE', 0))
        pending = int(counts.get('NOT_QUERIED', 0))
        # Replace 'Included' with a status reflecting actual query completeness.
        summary[-1]['STATUS'] = ('INCOMPLETE / QUERY FAILURES' if failures or skipped
                                 else 'INCOMPLETE / NOT QUERIED' if pending
                                 else 'COMPLETE')
        if failures or skipped or pending:
            print(f'WARNING: {survey}: {failures} failed, {skipped} service-unavailable, '
                  f'{pending} not queried. These are NOT nondetections.')
    else:
        failures_path = OUTPUT / f'{survey}_failures.pkl'
        if failures_path.exists() and len(read_table(failures_path)):
            summary[-1]['STATUS'] = 'INCOMPLETE / LEGACY QUERY FAILURES'
            print(f'WARNING: {survey} has failed queries (legacy format).')
        else:
            summary[-1]['STATUS'] = 'STATUS UNVERIFIED (LEGACY OUTPUT)'

df_measurement = pd.concat(frames, ignore_index=True)[MEAS_COLS]
df_survey_summary = pd.DataFrame(summary)

# Validate and convert NEIGHBOR to nullable integer
valid = (
    df_measurement["NEIGHBOR"].isna()
    | df_measurement["NEIGHBOR"].isin([0, 1])
)

if not valid.all():
    raise ValueError("Invalid NEIGHBOR values.")

df_measurement["NEIGHBOR"] = (
    df_measurement["NEIGHBOR"].astype("Int64")
)

# No partial output is written until all validations have passed.
FINAL.mkdir(parents=True, exist_ok=True)
for filename, df in [('source.pkl', df_source), ('name.pkl', df_name),
                     ('class.pkl', df_class), ('period.pkl', df_period),
                     ('measurement.pkl', df_measurement)]:
    df.to_pickle(FINAL / filename)

print('\nFINAL DATABASE SUMMARY')
print('=' * 62)
for name, df in [('source', df_source), ('name', df_name), ('class', df_class),
                 ('period', df_period), ('measurement', df_measurement)]:
    print(f'{name:16s} {len(df):8,d} rows')
print('\nSURVEY INPUT STATUS')
print(df_survey_summary.to_string(index=False))
print(f'\nFinal tables saved to: {FINAL}')
print('Intermediate stage files were not overwritten.')
