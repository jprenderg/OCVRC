"""Run all four internal OCVRC_Github validation scripts.

Independent read-only checks. Warnings (e.g. XMM outage or overlapping
VLASS tiles) are reported but do not fail the pipeline. Validation failures
raise an exception AFTER printing a consolidated summary.
"""
from pathlib import Path
import runpy
import pandas as pd

HERE = Path(__file__).resolve().parent

SUITES = [
    ('Database', 'Validate_Database.py', 'df_validation'),
    ('Astrometry', 'Validate_Astrometry.py', 'df_astrometry_validation'),
    ('Measurements', 'Validate_Measurements.py', 'df_measurement_validation'),
    ('Query Status', 'Validate_Query_Status.py', 'df_query_status_validation'),
]

all_results = []
summary_rows = []

print('\nOCVRC_Github COMPREHENSIVE INTERNAL VALIDATION')
print('=' * 70)

for suite, filename, variable in SUITES:
    script = HERE / filename
    if not script.is_file():
        raise FileNotFoundError(f'Missing {suite} validation script: {script}')

    print(f'\n[VALIDATION: {suite.upper()}] {filename}')
    namespace = runpy.run_path(str(script), run_name='__main__')
    if variable not in namespace:
        raise RuntimeError(f'{filename} did not expose {variable}')

    frame = namespace[variable].copy()
    if 'RESULT' not in frame.columns:
        raise RuntimeError(f'{filename}: result DataFrame lacks RESULT')
    invalid = ~frame['RESULT'].isin(['PASS', 'WARNING', 'FAIL'])
    if invalid.any():
        raise RuntimeError(f'{filename}: unexpected RESULT values')
    counts = frame['RESULT'].value_counts()
    summary_rows.append({
        'SUITE': suite,
        'PASS': int(counts.get('PASS', 0)),
        'WARNING': int(counts.get('WARNING', 0)),
        'FAIL': int(counts.get('FAIL', 0)),
    })
    frame.insert(0, 'SUITE', suite)
    all_results.append(frame.rename(columns={'DETAILS': 'DETAIL'}))

# Visible in Spyder's Variable Explorer, without overwriting the old
# df_validation variable used by the database-only validation script.
df_validation_summary = pd.DataFrame(summary_rows)
df_validation_all = pd.concat(all_results, ignore_index=True)
df_validation = all_results[0].drop(columns='SUITE').copy()

print('\nCONSOLIDATED VALIDATION SUMMARY')
print('=' * 70)
print(df_validation_summary.to_string(index=False))
print('-' * 70)
totals = df_validation_summary[['PASS', 'WARNING', 'FAIL']].sum()
print(f"TOTAL: {totals['PASS']} PASS, {totals['WARNING']} WARNING, {totals['FAIL']} FAIL")
print('Read-only validation: no astronomical queries or data modifications.')

if totals['FAIL']:
    raise RuntimeError(f"Internal validation failed: {totals['FAIL']} failed checks.")
