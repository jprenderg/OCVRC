"""Gaia DR3 astrometry and photometry with resumable, per-source queries."""
from pathlib import Path
import configparser
import numpy as np
import pandas as pd
from astroquery.gaia import Gaia
from _catalog_query_common import managed
ROOT = Path(__file__).resolve().parents[1]
config = configparser.ConfigParser(interpolation=None)
if not config.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')
def config_path(section, key, default):
    p = Path(config.get(section, key, fallback=default)).expanduser()
    return p if p.is_absolute() else ROOT / p
SOURCE_PATH = config_path('PATHS', 'SOURCE_PKL', 'output/source.pkl')
OUTPUT_DIR = config_path('PATHS', 'OUTPUT_DIR', 'output')
RADIUS = config.getfloat('GAIA', 'MATCH_RADIUS_ARCSEC', fallback=2.0)
START_ROW = config.getint('GAIA', 'START_ROW', fallback=0)
END_ROW = config.getint('GAIA', 'END_ROW', fallback=-1)
if not 0 < RADIUS <= 60:
    raise ValueError('GAIA MATCH_RADIUS_ARCSEC must be in (0, 60].')
GAIA_COLS = [
    'source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec',
    'phot_g_mean_mag', 'phot_g_mean_flux', 'phot_g_mean_flux_error',
    'phot_bp_mean_mag', 'phot_bp_mean_flux', 'phot_bp_mean_flux_error',
    'phot_rp_mean_mag', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error',
]
PHOT_COLUMNS = ['GAIA_G_MAG', 'GAIA_G_MAG_ERROR', 'GAIA_BP_MAG',
                'GAIA_BP_MAG_ERROR', 'GAIA_RP_MAG', 'GAIA_RP_MAG_ERROR']
RAD_DEG = RADIUS / 3600.0
FACTOR = 2.5 / np.log(10.0)


def finite_number(v):
    """Convert null and masked Gaia values to NaN without accidental zeros."""
    if v is None or np.ma.is_masked(v):
        return np.nan
    try:
        v = float(v)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def retrieve(ra, dec):
    adql = f'''SELECT TOP 100 {', '.join(GAIA_COLS)},
      DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra:.12f}, {dec:.12f})) AS sep_deg
      FROM gaiadr3.gaia_source
      WHERE 1=CONTAINS(POINT('ICRS', ra, dec),
                      CIRCLE('ICRS', {ra:.12f}, {dec:.12f}, {RAD_DEG:.12f}))
      ORDER BY sep_deg ASC'''
    result = Gaia.launch_job_async(adql, dump_to_file=False).get_results()
    if len(result) == 0:
        return None
    match = result[0]
    return match if finite_number(match['sep_deg']) * 3600 <= RADIUS else None


def flux_mag_error(match, band):
    f = finite_number(match[f'phot_{band}_mean_flux'])
    e = finite_number(match[f'phot_{band}_mean_flux_error'])
    return FACTOR * e / f if np.isfinite(f) and f > 0 and np.isfinite(e) and e >= 0 else np.nan


def make_photometry(df):
    parts = []
    for band, prefix in [('G', 'GAIA_G'), ('BP', 'GAIA_BP'), ('RP', 'GAIA_RP')]:
        x = pd.DataFrame({
            'SOURCE_ID': df['SOURCE_ID'], 'OBSERVATORY': 'Gaia',
            'EPOCH': 'DR3', 'BAND': band, 'UNITS': 'mag',
            'MEASUREMENT': df[f'{prefix}_MAG'], 'ERROR': df[f'{prefix}_MAG_ERROR'],
            'POS_ERROR': np.nan, 'NEG_ERROR': np.nan, 'NEIGHBOR': np.nan,
        })
        parts.append(x)
    return pd.concat(parts, ignore_index=True).dropna(subset=['MEASUREMENT']).reset_index(drop=True)



def query_one(src, radius):
    match = retrieve(float(src.RA), float(src.DEC))
    record = {'GAIA_FLAG': match is not None, 'GAIA_SOURCE_ID': None,
              'GAIA_SEPARATION_ARCSEC': np.nan, 'DISTANCE': np.nan,
              'DISTANCE_ERROR': np.nan, 'PM_RA': np.nan, 'PM_DEC': np.nan}
    record.update({col: np.nan for col in PHOT_COLUMNS})
    if match is not None:
        record['GAIA_SOURCE_ID'] = str(match['source_id'])
        record['GAIA_SEPARATION_ARCSEC'] = finite_number(match['sep_deg']) * 3600
        plx, eplx = finite_number(match['parallax']), finite_number(match['parallax_error'])
        if np.isfinite(plx) and plx > 0:
            record['DISTANCE'] = 1 / plx  # kpc; simple inverse parallax
            if np.isfinite(eplx) and eplx >= 0:
                record['DISTANCE_ERROR'] = eplx / plx**2
        record['PM_RA'] = finite_number(match['pmra'])
        record['PM_DEC'] = finite_number(match['pmdec'])
        for band, prefix in [('g','GAIA_G'),('bp','GAIA_BP'),('rp','GAIA_RP')]:
            record[f'{prefix}_MAG'] = finite_number(match[f'phot_{band}_mean_mag'])
            record[f'{prefix}_MAG_ERROR'] = flux_mag_error(match, band)
    rows = []
    for band, prefix in [('G','GAIA_G'),('BP','GAIA_BP'),('RP','GAIA_RP')]:
        value = record[f'{prefix}_MAG']
        if np.isfinite(value):
            rows.append({'OBSERVATORY': 'Gaia','EPOCH': 'DR3', 'BAND': band,
                         'UNITS': 'mag', 'MEASUREMENT': value,
                         'ERROR': record[f'{prefix}_MAG_ERROR'],
                         'POS_ERROR': np.nan, 'NEG_ERROR': np.nan, 'NEIGHBOR': np.nan})
    return ('MATCH' if match is not None else 'NO_MATCH'), record, rows

df_source_initial = pd.read_pickle(SOURCE_PATH).copy().reset_index(drop=True)
df_gaia, df_gaia_measurement, df_gaia_failures, df_query_status = managed(
    'Gaia', df_source_initial, query_one, OUTPUT_DIR, RADIUS, 'GAIA',
    config, START_ROW, END_ROW)
# Keep ALL source rows, even sources with failed/pending Gaia queries.
astro_cols = ['GAIA_FLAG','GAIA_SOURCE_ID','GAIA_SEPARATION_ARCSEC',
              'DISTANCE','DISTANCE_ERROR','PM_RA','PM_DEC']
for col in astro_cols:
    if col not in df_gaia.columns:
        df_gaia[col] = np.nan
df_source = df_source_initial.merge(df_gaia[['SOURCE_ID'] + astro_cols],
    on='SOURCE_ID', how='left', validate='one_to_one')
# Ensure the builder's coordinate freshness test has a complete source table.
df_source.to_pickle(OUTPUT_DIR / 'source_gaia.pkl')
print(f'GAIA DR3: matches={int(df_query_status.STATUS.eq("MATCH").sum())}, '
      f'measurements={len(df_gaia_measurement)}, '
      f'query failures={len(df_gaia_failures)}')
