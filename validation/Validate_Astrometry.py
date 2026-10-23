"""Independent, read-only astrometric validation for OCVRC_Github.

Run in Spyder with %runfile .../validation/Validate_Astrometry.py --wdir
Uses only local output files; makes no internet requests or database changes.
"""
from pathlib import Path
import configparser
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
cfg = configparser.ConfigParser(interpolation=None)
if not cfg.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')
p = Path(cfg.get('PATHS', 'OUTPUT_DIR', fallback='output')).expanduser()
OUTPUT = p if p.is_absolute() else ROOT / p
FINAL = OUTPUT / 'final'

results = []
def record(result, label, details=''):
    results.append({'RESULT': result, 'CHECK': label, 'DETAILS': str(details)})
    print(f'{result:7} {label}' + (f' — {details}' if details else ''))

def check(ok, label, details=''):
    record('PASS' if ok else 'FAIL', label, '' if ok else details)

def numeric(col):
    original = source[col]
    n = pd.to_numeric(original, errors='coerce')
    malformed = original.notna() & n.isna()
    check(not malformed.any(), f'{col}: numeric or missing', f'{malformed.sum()} malformed')
    invalid_inf = n.notna() & ~np.isfinite(n)
    check(not invalid_inf.any(), f'{col}: finite when present', f'{invalid_inf.sum()} nonfinite')
    return n

print('\nOCVRC_Github ASTROMETRY VALIDATION')
print('=' * 66)
source = pd.read_pickle(FINAL / 'source.pkl')
required = ['SOURCE_ID','RA','DEC','SIMBAD','MAJOR_AXIS','MINOR_AXIS',
            'POSITION_ANGLE','DISTANCE','DISTANCE_ERROR','PM_RA','PM_DEC']
missing = [c for c in required if c not in source]
if missing:
    raise ValueError(f'Source table missing columns: {missing}')

check(not source['SOURCE_ID'].isna().any() and source['SOURCE_ID'].is_unique,
      'SOURCE_ID: unique and present')
ra = numeric('RA')
dec = numeric('DEC')
check(bool((ra.notna() & ra.ge(0) & ra.lt(360)).all()), 'RA: 0 <= RA < 360 deg')
check(bool((dec.notna() & dec.ge(-90) & dec.le(90)).all()), 'DEC: -90 <= DEC <= 90 deg')

sim = source['SIMBAD']
valid_sim = sim.notna() & sim.map(lambda x: isinstance(x, (bool, np.bool_)) or x in (0,1))
check(bool(valid_sim.all()), 'SIMBAD: valid nonmissing match flags')
matched = sim.eq(True).fillna(False)
unmatched = sim.eq(False).fillna(False)
record('PASS', 'SIMBAD: match counts', f'{matched.sum()} matched; {unmatched.sum()} unmatched')

# ============================================================
# SIMBAD POSITIONAL UNCERTAINTIES
# ============================================================

maj = numeric('MAJOR_AXIS')
minor = numeric('MINOR_AXIS')
angle = numeric('POSITION_ANGLE')

# Check that positional uncertainties are nonnegative
check(
    bool((maj.dropna() >= 0).all()),
    'SIMBAD ellipse: major uncertainty >= 0 mas'
)

check(
    bool((minor.dropna() >= 0).all()),
    'SIMBAD ellipse: minor uncertainty >= 0 mas'
)

# SIMBAD can report coordinate uncertainties with PA = 90 degrees,
# in which case the first uncertainty may be smaller than the second.
invalid_ellipse = (
    maj.notna()
    & minor.notna()
    & (maj < minor)
    & ~angle.eq(90).fillna(False)
)

check(
    not invalid_ellipse.any(),
    'SIMBAD: positional uncertainty ordering consistent with angle convention',
    f'{invalid_ellipse.sum()} violations'
)

# Sources without SIMBAD matches should have no SIMBAD uncertainties
check(
    bool((~unmatched | (maj.isna() & minor.isna() & angle.isna())).all()),
    'SIMBAD ellipse: unmatched sources have no SIMBAD uncertainties'
)

# Missing uncertainty values may be legitimate for matched SIMBAD objects
incomplete = matched & (maj.isna() | minor.isna() | angle.isna())

record(
    'WARNING' if incomplete.any() else 'PASS',
    'SIMBAD ellipse: uncertainties available for matched sources',
    f'{incomplete.sum()} matched sources missing at least one value'
    if incomplete.any() else ''
)

# Position angle plausibility
bad_angle = angle.notna() & ~angle.between(-360, 360)

record(
    'WARNING' if bad_angle.any() else 'PASS',
    'SIMBAD ellipse: position angle plausible',
    f'{bad_angle.sum()} values outside [-360, 360] deg'
    if bad_angle.any() else ''
)

par_dist = numeric('DISTANCE')
par_err = numeric('DISTANCE_ERROR')
pm_ra = numeric('PM_RA')
pm_dec = numeric('PM_DEC')
check(bool((par_dist.dropna() > 0).all()), 'Gaia: distance > 0 kpc when present')
check(bool((par_err.dropna() >= 0).all()), 'Gaia: distance uncertainty >= 0 kpc')
check(bool((~par_err.notna() | par_dist.notna()).all()),
      'Gaia: distance error not present without distance')
partial_pm = pm_ra.notna() ^ pm_dec.notna()
record('WARNING' if partial_pm.any() else 'PASS', 'Gaia: PM_RA and PM_DEC populated together',
       f'{partial_pm.sum()} sources have only one component' if partial_pm.any() else '')
record('PASS', 'Gaia: astrometry coverage',
       f'{par_dist.notna().sum()} distances; {pm_ra.notna().sum()} PM_RA; {pm_dec.notna().sum()} PM_DEC')

# Compare final astrometry to the existing local intermediate files. This tests
# data transfer through the pipeline, not scientific correctness of SIMBAD/Gaia.
for filename, columns, label in [
    ('source.pkl', ['RA','DEC'], 'ingestion'),
    ('source_gaia.pkl', ['DISTANCE','DISTANCE_ERROR','PM_RA','PM_DEC'], 'Gaia')
]:
    path = OUTPUT / filename
    if not path.exists():
        record('WARNING', f'{label}: intermediate file available', f'{path} not found')
        continue
    stage = pd.read_pickle(path)
    if 'SOURCE_ID' not in stage or any(c not in stage for c in columns):
        record('FAIL', f'{label}: intermediate schema', f'missing SOURCE_ID or {columns}')
        continue
    if stage['SOURCE_ID'].isna().any() or stage['SOURCE_ID'].duplicated().any():
        record('FAIL', f'{label}: intermediate SOURCE_ID', 'missing/duplicate SOURCE_ID')
        continue
    joined = source[['SOURCE_ID'] + columns].merge(
        stage[['SOURCE_ID'] + columns], on='SOURCE_ID', how='outer',
        suffixes=('_final', '_stage'), indicator=True, validate='one_to_one'
    )
    if not joined['_merge'].eq('both').all():
        record('FAIL', f'{label}: source coverage', f'{joined["_merge"].ne("both").sum()} unmatched IDs')
        continue
    failed = []
    for col in columns:
        left = pd.to_numeric(joined[f'{col}_final'], errors='coerce').to_numpy(dtype=float)
        right = pd.to_numeric(joined[f'{col}_stage'], errors='coerce').to_numpy(dtype=float)
        wrong = ~np.isclose(left, right, rtol=1e-12, atol=1e-12, equal_nan=True)
        if wrong.any():
            failed.append(f'{col}: {wrong.sum()}')
    check(not failed, f'{label}: final astrometry matches intermediate values', ', '.join(failed))

# Check the SIMBAD flag against ingestion without making a SIMBAD network query.
path = OUTPUT / 'source.pkl'
if path.exists():
    stage = pd.read_pickle(path)
    if {'SOURCE_ID','SIMBAD_FLAG'}.issubset(stage.columns) and stage['SOURCE_ID'].is_unique:
        comp = source[['SOURCE_ID','SIMBAD']].merge(
            stage[['SOURCE_ID','SIMBAD_FLAG']], on='SOURCE_ID', how='outer', indicator=True
        )
        ok = comp['_merge'].eq('both').all() and comp['SIMBAD'].eq(comp['SIMBAD_FLAG']).fillna(False).all()
        check(bool(ok), 'SIMBAD: final flags match ingestion flags')
    else:
        record('WARNING', 'SIMBAD: ingestion flag comparison', 'Intermediate flags unavailable')

# Final result remains available in Spyder Variable Explorer.
df_astrometry_validation = pd.DataFrame(results)
print('\nASTROMETRY VALIDATION SUMMARY')
print('=' * 66)
print(df_astrometry_validation['RESULT'].value_counts().reindex(['PASS','WARNING','FAIL'], fill_value=0).to_string())
print('Read-only: no files modified; no network queries performed.')
