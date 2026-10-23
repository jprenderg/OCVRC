"""Shared X-ray querying utilities; no data are silently discarded on request failures."""
from pathlib import Path
from configparser import ConfigParser
from io import BytesIO
from time import sleep
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse_single_table
from astropy import units as u

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ConfigParser(inline_comment_prefixes=('#',';'))
if not CONFIG.read(ROOT/'config.txt', encoding='utf-8'):
    raise FileNotFoundError(f'Missing config: {ROOT / "config.txt"}')

def path_option(key, default):
    p = Path(CONFIG.get('PATHS', key, fallback=default)).expanduser()
    return p if p.is_absolute() else ROOT/p

OUTPUT = path_option('OUTPUT_DIR', 'output')
SOURCE = path_option('SOURCE_PKL', 'output/source.pkl')
OUTPUT.mkdir(parents=True, exist_ok=True)

def settings(section):
    return dict(radius=CONFIG.getfloat(section,'MATCH_RADIUS_ARCSEC', fallback=5),
                start=CONFIG.getint(section,'START_ROW', fallback=0),
                end=CONFIG.getint(section,'END_ROW', fallback=-1),
                timeout=CONFIG.getfloat(section,'TIMEOUT_SECONDS', fallback=60),
                retries=CONFIG.getint(section,'MAX_RETRIES', fallback=2),
                delay=CONFIG.getfloat(section,'RETRY_DELAY_SECONDS', fallback=3))

def number(v):
    try:
        if np.ma.is_masked(v): return np.nan
        f=float(v)
        return f if np.isfinite(f) else np.nan
    except (TypeError, ValueError, OverflowError): return np.nan

def table_frame(data):
    if isinstance(data, pd.DataFrame): return data.copy()
    return data.to_pandas() if hasattr(data,'to_pandas') else pd.DataFrame(data)

def parse_votable(content):
    return parse_single_table(BytesIO(content)).to_table().to_pandas()

def cone_get(url, params, timeout):
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    if not r.content: raise RuntimeError('Empty response')
    return parse_votable(r.content)

def column(df, *alternatives):
    names = {str(c).lower():c for c in df.columns}
    for c in alternatives:
        if c.lower() in names: return names[c.lower()]
    return None

def closest(df, ra, dec, radius):
    ra_key=column(df,'ra','raj2000')
    dec_key=column(df,'dec','dej2000')
    if df.empty or ra_key is None or dec_key is None: return None, np.nan
    rr=pd.to_numeric(df[ra_key],errors='coerce'); dd=pd.to_numeric(df[dec_key],errors='coerce')
    ok=rr.between(0,360,inclusive='left') & dd.between(-90,90)
    if not ok.any(): return None,np.nan
    df=df.loc[ok].reset_index(drop=True)
    target=SkyCoord(ra*u.deg,dec*u.deg)
    points=SkyCoord(pd.to_numeric(df[ra_key]).to_numpy()*u.deg,pd.to_numeric(df[dec_key]).to_numpy()*u.deg)
    separations=target.separation(points).arcsec
    idx=int(np.argmin(separations)); sep=float(separations[idx])
    return (df.iloc[idx],sep) if sep<=radius else (None,sep)

def get_value(row, *keys):
    if row is None: return np.nan
    names={str(k).lower():k for k in row.index}
    for k in keys:
        if k.lower() in names: return number(row[names[k.lower()]])
    return np.nan

def get_text(row, *keys):
    if row is None: return None
    names={str(k).lower():k for k in row.index}
    for k in keys:
        if k.lower() in names:
            v=row[names[k.lower()]]
            if v is not None and not pd.isna(v):
                return v.decode(errors='replace') if isinstance(v,bytes) else str(v)
    return None

def run(survey, query, fields, bands, observatory, epoch, units, default_radius):
    """Query with checkpoint/resume and return traditional three DataFrames."""
    from _query_manager import run_queries
    opts = settings(survey)
    inp = pd.read_pickle(SOURCE)
    cfg = CONFIG
    failure_limit = cfg.getint('QUERY_MANAGEMENT', 'CONSECUTIVE_FAILURE_LIMIT', fallback=5)
    save_every = cfg.getint('QUERY_MANAGEMENT', 'SAVE_EVERY', fallback=10)
    resume = cfg.getboolean('QUERY_MANAGEMENT', 'RESUME', fallback=True)

    def query_one(item, radius):
        ra, dec = number(item['RA']), number(item['DEC'])
        if not (np.isfinite(ra) and np.isfinite(dec) and 0 <= ra < 360 and -90 <= dec <= 90):
            raise ValueError('Invalid coordinates')
        df = table_frame(query(ra, dec, radius, opts['timeout']))
        # An empty response is a genuine NO_MATCH. A malformed nonempty response is not.
        if not df.empty and (column(df, 'ra', 'raj2000') is None or
                             column(df, 'dec', 'dej2000') is None):
            raise ValueError('Query response lacks RA/DEC columns')
        match, sep = closest(df, ra, dec, radius)
        record = {f'{survey}_MATCH': match is not None,
                  'SEP_ARCSEC': sep, 'MATCH_NAME': get_text(match, 'name', 'iauname', 'src_name')}
        record.update({key: get_value(match, *options) for key, options in fields.items()})
        rows = []
        if match is not None:
            for band, flux_col, error_col, pos_col, neg_col in bands:
                flux = record.get(flux_col, np.nan)
                if not np.isfinite(flux):
                    continue
                rows.append({'OBSERVATORY': observatory, 'EPOCH': epoch,
                             'BAND': band, 'UNITS': units, 'MEASUREMENT': flux,
                             'ERROR': record.get(error_col, np.nan) if error_col else np.nan,
                             'POS_ERROR': record.get(pos_col, np.nan) if pos_col else np.nan,
                             'NEG_ERROR': record.get(neg_col, np.nan) if neg_col else np.nan,
                             'NEIGHBOR': np.nan})
        return ('MATCH' if match is not None else 'NO_MATCH'), record, rows

    df_survey, df_measurement, df_failures, df_status = run_queries(
        survey, inp, query_one, OUTPUT, radius=opts['radius'],
        retries=opts['retries'], delay=opts['delay'],
        failure_limit=failure_limit, save_every=save_every,
        start=max(0, opts['start']), end=opts['end'], resume=resume)
    # Expose query status in Spyder's Variable Explorer as well.
    globals()['df_query_status'] = df_status
    return df_survey, df_measurement, df_failures
