"""Shared I/O and catalog-matching helpers for OCVRC_Github photometry."""
from pathlib import Path
import configparser
import time
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u

ROOT = Path(__file__).resolve().parents[1]
CONFIG = configparser.ConfigParser(interpolation=None)
if not CONFIG.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')

def path_from_config(section, key, default):
    p = Path(CONFIG.get(section, key, fallback=default)).expanduser()
    return p if p.is_absolute() else ROOT / p

def number(x):
    if x is None or np.ma.is_masked(x):
        return np.nan
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan

def load_source():
    path = path_from_config('PATHS', 'SOURCE_PKL', 'output/source.pkl')
    df = pd.read_pickle(path).copy().reset_index(drop=True)
    missing = {'SOURCE_ID','RA','DEC'} - set(df.columns)
    if missing:
        raise ValueError(f'Missing source columns: {sorted(missing)}')
    if df.SOURCE_ID.isna().any() or df.SOURCE_ID.duplicated().any():
        raise ValueError('SOURCE_ID is missing or duplicated')
    return df

def bounds(df, section):
    start = CONFIG.getint(section, 'START_ROW', fallback=0)
    end = CONFIG.getint(section, 'END_ROW', fallback=-1)
    end = len(df) if end <= 0 else min(end, len(df))
    if not (0 <= start <= end):
        raise ValueError(f'Invalid {section} bounds: {start}:{end}, size={len(df)}')
    return start, end

def query_with_retry(fn, section):
    retries = CONFIG.getint(section,'MAX_RETRIES',fallback=2)
    delay = CONFIG.getfloat(section,'RETRY_DELAY_SECONDS',fallback=3)
    for attempt in range(retries + 1):
        try:
            return fn()
        except Exception:
            if attempt == retries:
                raise
            time.sleep(delay)

def nearest(table, coord, ra_key='ra', dec_key='dec', radius_arcsec=2):
    if len(table) == 0:
        return None, np.nan
    ra = np.array([number(v) for v in table[ra_key]], dtype=float)
    dec = np.array([number(v) for v in table[dec_key]], dtype=float)
    valid = np.isfinite(ra) & np.isfinite(dec)
    if not valid.any():
        return None, np.nan
    idx = np.flatnonzero(valid)
    coords = SkyCoord(ra[idx]*u.deg, dec[idx]*u.deg)
    sep = coord.separation(coords).arcsec
    pos = int(np.argmin(sep))
    return (table[int(idx[pos])], float(sep[pos])) if sep[pos] <= radius_arcsec else (None, np.nan)

def pivot(df, mappings, observatory, epoch):
    records=[]
    for band, value_col, error_col in mappings:
        for row in df.itertuples(index=False):
            values = row._asdict()
            flux = number(values[value_col])
            if not np.isfinite(flux):
                continue
            records.append({'SOURCE_ID':values['SOURCE_ID'], 'OBSERVATORY':observatory,
                            'EPOCH':epoch, 'BAND':band, 'UNITS':'mag',
                            'MEASUREMENT':flux, 'ERROR':number(values[error_col]),
                            'POS_ERROR':np.nan, 'NEG_ERROR':np.nan, 'NEIGHBOR':np.nan})
    columns=['SOURCE_ID','OBSERVATORY','EPOCH','BAND','UNITS','MEASUREMENT',
             'ERROR','POS_ERROR','NEG_ERROR','NEIGHBOR']
    return pd.DataFrame(records,columns=columns)
