"""Query ASTRON's LoTSS DR2 SOURCE catalogue by position.

Input: output/source.pkl, from source ingestion.
Outputs (pandas pickle only):
    output/LoTSS.pkl
    output/LoTSS_measurements.pkl
    output/LoTSS_failures.pkl

The source catalogue contains catalogue *detections* only; a non-match is
NOT a zero-flux or a forced-photometry measurement.
Run in Spyder's current IPython console to inspect the df_* DataFrames.
"""

from pathlib import Path
from configparser import ConfigParser
from io import BytesIO
from time import sleep
from _catalog_query_common import managed
import numpy as np
import pandas as pd
import requests
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ConfigParser(inline_comment_prefixes=('#', ';'))
if not CONFIG.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(f'Missing {ROOT / "config.txt"}')

def project_path(section, key, default):
    p = Path(CONFIG.get(section, key, fallback=default)).expanduser()
    return p if p.is_absolute() else ROOT / p

SOURCE_PATH = project_path('PATHS', 'SOURCE_PKL', 'output/source.pkl')
OUTPUT_DIR = project_path('PATHS', 'OUTPUT_DIR', 'output')
RADIUS = CONFIG.getfloat('LOTSS', 'MATCH_RADIUS_ARCSEC', fallback=2.0)
START_ROW = CONFIG.getint('LOTSS', 'START_ROW', fallback=0)
END_ROW = CONFIG.getint('LOTSS', 'END_ROW', fallback=-1)
TIMEOUT = CONFIG.getfloat('LOTSS', 'TIMEOUT_SECONDS', fallback=60)
RETRIES = CONFIG.getint('LOTSS', 'MAX_RETRIES', fallback=2)
DELAY = CONFIG.getfloat('LOTSS', 'RETRY_DELAY_SECONDS', fallback=3)
SCS_URL = 'https://vo.astron.nl/lotss_dr2/q/src_cone/scs.xml'

if not 0 < RADIUS <= 60:
    raise ValueError('LOTSS MATCH_RADIUS_ARCSEC must be between 0 and 60')
if RETRIES < 0:
    raise ValueError('MAX_RETRIES must be >= 0')

# Mirror the original OCVRC LoTSS fields, retaining catalogue provenance.
COLUMNS = ['SOURCE_ID', 'LOTSS_MATCH', 'LOTSS_SEP_ARCSEC',
           'Source_Name', 'RA', 'DEC', 'E_RA', 'E_DEC', 'Peak_flux',
           'E_Peak_flux', 'Total_flux', 'E_Total_flux', 'Isl_rms',
           'S_Code', 'Mosaic_ID']
CATALOG_FIELDS = COLUMNS[3:]


def number(x):
    try:
        if np.ma.is_masked(x) or x is None:
            return np.nan
        z = float(x)
        return z if np.isfinite(z) else np.nan
    except (ValueError, TypeError):
        return np.nan


def plain(x):
    if x is None or np.ma.is_masked(x):
        return None
    if isinstance(x, bytes):
        return x.decode('utf-8', errors='replace').strip()
    return str(x).strip()


def catalog_query(ra, dec):
    params = {'RA':ra, 'DEC':dec, 'SR':RADIUS/3600.,'VERB':3,'MAXREC':1000}
    response = requests.get(SCS_URL, params=params, timeout=TIMEOUT)
    response.raise_for_status()
    if b'QUERY_STATUS' in response.content and b'ERROR' in response.content[:10000]:
        raise RuntimeError('LoTSS service returned QUERY_STATUS=ERROR')
    return Table.read(BytesIO(response.content), format='votable')

def choose_nearest(table, ra, dec):
    if len(table) == 0:
        return None, np.nan
    if not {'RA', 'DEC'}.issubset(table.colnames):
        raise ValueError(f'LoTSS response missing RA or DEC: {table.colnames}')
    ras = np.array([number(x) for x in table['RA']])
    decs = np.array([number(x) for x in table['DEC']])
    valid = np.isfinite(ras) & np.isfinite(decs)
    if not valid.any():
        return None, np.nan
    positions = SkyCoord(ras[valid] * u.deg, decs[valid] * u.deg)
    target = SkyCoord(ra * u.deg, dec * u.deg)
    sep = target.separation(positions).arcsec
    j = int(np.argmin(sep))
    if sep[j] > RADIUS:
        return None, np.nan
    index = np.flatnonzero(valid)[j]
    return table[index], float(sep[j])



def query_one(src, radius):
    ra,dec=float(src.RA),float(src.DEC)
    if not (0<=ra<360 and -90<=dec<=90):
        raise ValueError(f'Invalid RA/DEC: {ra}, {dec}')
    table=catalog_query(ra,dec)
    closest,sep=choose_nearest(table,ra,dec)
    record={'LOTSS_MATCH':closest is not None,'LOTSS_SEP_ARCSEC':sep,
            **{field:np.nan for field in CATALOG_FIELDS}}
    for field in ('Source_Name','S_Code','Mosaic_ID'):
        record[field]=None
    if closest is not None:
        for field in CATALOG_FIELDS:
            if field in table.colnames:
                record[field]=(plain(closest[field]) if field in
                               ('Source_Name','S_Code','Mosaic_ID') else number(closest[field]))
    rows=[]
    if closest is not None and np.isfinite(record['Peak_flux']):
        rows.append({'OBSERVATORY':'LOFAR','EPOCH':'LoTSS DR2',
                     'BAND':'120-168 MHz','UNITS':'mJy',
                     'MEASUREMENT':record['Peak_flux'],
                     'ERROR':record['E_Peak_flux'], 'POS_ERROR':np.nan,
                     'NEG_ERROR':np.nan,'NEIGHBOR':np.nan})
    return ('MATCH' if closest is not None else 'NO_MATCH'),record,rows

df_source=pd.read_pickle(SOURCE_PATH).copy().reset_index(drop=True)
df_lotss,df_lotss_measurement,df_lotss_failures,df_query_status=managed(
    'LoTSS',df_source,query_one,OUTPUT_DIR,RADIUS,'LOTSS',CONFIG,START_ROW,END_ROW)
print(f'LoTSS DR2: matches={int(df_query_status.STATUS.eq("MATCH").sum())}, '
      f'measurement rows={len(df_lotss_measurement)}, failed={len(df_lotss_failures)}')
