"""OCVRC_Github: read-only, standalone measurement validation.

Run directly in Spyder or with Python. Validates final measurement.pkl,
independently of the original OCVRC database or external services.
"""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PATH = ROOT / 'output' / 'final' / 'measurement.pkl'

MEAS_COLS = [
    'SOURCE_ID', 'OBSERVATORY', 'EPOCH', 'BAND', 'UNITS',
    'MEASUREMENT', 'ERROR', 'POS_ERROR', 'NEG_ERROR', 'NEIGHBOR'
]

checks = []


def record(result, description, detail=''):
    checks.append({'RESULT': result, 'CHECK': description, 'DETAIL': str(detail)})
    print(f'{result:<7} {description}' + (f' — {detail}' if detail else ''))


def check(condition, description, detail=''):
    record('PASS' if bool(condition) else 'FAIL', description, detail)


def populated(series):
    return series.notna() & series.astype('string').str.strip().ne('').fillna(False)


print('\nOCVRC_Github MEASUREMENT VALIDATION')
print('=' * 70)

if not PATH.exists():
    record('FAIL', 'measurement.pkl exists', str(PATH))
    df_measurement_validation = pd.DataFrame(checks)
    raise FileNotFoundError(PATH)

frame = pd.read_pickle(PATH)
check(isinstance(frame, pd.DataFrame), 'measurement.pkl contains a DataFrame')
missing = [c for c in MEAS_COLS if c not in frame.columns]
check(not missing, 'all required measurement columns present', missing if missing else '')
if missing:
    df_measurement_validation = pd.DataFrame(checks)
    raise ValueError(f'Missing measurement columns: {missing}')

extra = [c for c in frame.columns if c not in MEAS_COLS]
check(not extra, 'no unexpected measurement columns', extra if extra else '')
check(len(frame) > 0, 'measurement table is not empty', f'{len(frame)} rows')
check(frame['SOURCE_ID'].notna().all(), 'SOURCE_ID nonmissing')

for col in ('OBSERVATORY', 'EPOCH', 'BAND', 'UNITS'):
    blank = ~populated(frame[col])
    check(not blank.any(), f'{col}: populated for each observation', f'{int(blank.sum())} missing/blank')

numeric = {}
for col in ('MEASUREMENT', 'ERROR', 'POS_ERROR', 'NEG_ERROR'):
    values = pd.to_numeric(frame[col], errors='coerce')
    numeric[col] = values
    bad = frame[col].notna() & values.isna()
    check(not bad.any(), f'{col}: no nonnumeric entries', f'{int(bad.sum())} invalid')
    infinite = values.notna() & ~np.isfinite(values.to_numpy(dtype=float, na_value=np.nan))
    check(not infinite.any(), f'{col}: finite where present', f'{int(infinite.sum())} nonfinite')

measurement = numeric['MEASUREMENT']
error = numeric['ERROR']
check(measurement.notna().all(), 'MEASUREMENT: no missing measured values', f'{int(measurement.isna().sum())} missing')
negative_errors = error.notna() & (error < 0)
check(not negative_errors.any(), 'ERROR: nonnegative where provided', f'{int(negative_errors.sum())} negative')
zero_errors = error.notna() & (error == 0)
record('WARNING' if zero_errors.any() else 'PASS', 'ERROR: zero uncertainties', f'{int(zero_errors.sum())} rows' if zero_errors.any() else '')
# Accept a symmetric ERROR or a complete pair of asymmetric uncertainties.
# An incomplete asymmetric pair is not sufficient on its own.
pos_error = numeric['POS_ERROR']
neg_error = numeric['NEG_ERROR']
missing_uncertainty = error.isna() & (pos_error.isna() | neg_error.isna())
record('WARNING' if missing_uncertainty.any() else 'PASS',
       'Uncertainty: symmetric ERROR or both asymmetric errors present',
       f'{int(missing_uncertainty.sum())} rows lack usable uncertainty' if missing_uncertainty.any() else '')

neighbor = pd.to_numeric(frame['NEIGHBOR'], errors='coerce')
invalid_neighbor = (frame['NEIGHBOR'].notna() & neighbor.isna()) | (neighbor.notna() & ~neighbor.isin([0, 1]))
check(not invalid_neighbor.any(), 'NEIGHBOR: only 0, 1, or missing', f'{int(invalid_neighbor.sum())} invalid')

observatory = frame['OBSERVATORY'].astype('string').str.strip().str.upper()
radio = observatory.isin(['VLA', 'ASKAP']).fillna(False)
other = ~radio
check(frame.loc[radio, 'NEIGHBOR'].notna().all(), 'VLASS/RACS: NEIGHBOR always populated', f'{int(frame.loc[radio, "NEIGHBOR"].isna().sum())} missing')
nonradio_flags = frame.loc[other, 'NEIGHBOR'].notna()
check(not nonradio_flags.any(), 'Other surveys: NEIGHBOR not applicable (NULL)', f'{int(nonradio_flags.sum())} unexpected flags')

expected_radio = {
    'VLA': ('2-4 GHz', 'mJy'),
    'ASKAP': ('887.5 MHz', 'mJy'),
}
for telescope, (band, units) in expected_radio.items():
    mask = observatory.eq(telescope).fillna(False)
    if not mask.any():
        record('WARNING', f'{telescope}: observations present', 'no rows; coverage may be absent')
        continue
    actual_band = frame.loc[mask, 'BAND'].astype('string').str.strip()
    actual_units = frame.loc[mask, 'UNITS'].astype('string').str.strip()
    check(actual_band.eq(band).all(), f'{telescope}: expected frequency band',
          f'expected {band}; observed {sorted(actual_band.dropna().unique().tolist())}')
    check(actual_units.eq(units).all(), f'{telescope}: flux units',
          f'expected {units}; observed {sorted(actual_units.dropna().unique().tolist())}')

negative_radio = radio & measurement.lt(0)
record('PASS', 'Radio: negative measurements preserved (informational)', f'{int(negative_radio.sum())} rows')
record('PASS', 'Radio: neighbor-flagged measurements (informational)',
       f'{int((radio & neighbor.eq(1).fillna(False)).sum())} rows')

# Exact duplicates are suspicious; source+epoch alone is intentionally NOT
# considered a duplicate (overlapping VLASS image tiles can be legitimate).
duplicate = frame.duplicated(subset=MEAS_COLS, keep=False)
record('WARNING' if duplicate.any() else 'PASS', 'Exact duplicate measurement rows',
       f'{int(duplicate.sum())} rows' if duplicate.any() else '')

vlass = observatory.eq('VLA').fillna(False)
if vlass.any():
    groups = frame.loc[vlass].groupby(['SOURCE_ID', 'EPOCH'], dropna=False).size()
    overlaps = groups[groups > 1]
    record('WARNING' if len(overlaps) else 'PASS',
           'VLASS: repeated source/epoch combinations (possible overlapping tiles)',
           f'{len(overlaps)} groups' if len(overlaps) else '')

# Summarize survey-level coverage without assuming all sources are in a survey.
df_measurement_survey_summary = (
    frame.assign(_OBSERVATORY=observatory, _NEGATIVE=measurement.lt(0),
                 _NEIGHBOR_ONE=neighbor.eq(1).fillna(False))
    .groupby('_OBSERVATORY', dropna=False)
    .agg(ROWS=('SOURCE_ID', 'size'), SOURCES=('SOURCE_ID', 'nunique'),
         NEGATIVE_VALUES=('_NEGATIVE', 'sum'), NEIGHBOR_FLAGGED=('_NEIGHBOR_ONE', 'sum'))
    .reset_index().rename(columns={'_OBSERVATORY': 'OBSERVATORY'})
)

print('\nMEASUREMENTS BY OBSERVATORY')
print('-' * 70)
print(df_measurement_survey_summary.to_string(index=False))

df_measurement_validation = pd.DataFrame(checks)
print('\nMEASUREMENT VALIDATION SUMMARY')
print('=' * 70)
print(df_measurement_validation['RESULT'].value_counts().reindex(['PASS', 'WARNING', 'FAIL'], fill_value=0).to_string())
print('Read-only: no files modified; no network queries performed.')
