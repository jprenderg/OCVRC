"""GALEX NUV/FUV photometry with per-source checkpoint and resume."""
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Catalogs
from _photometry_common import CONFIG, number, load_source, bounds, path_from_config
from _catalog_query_common import managed
SECTION = 'GALEX'
RADIUS = CONFIG.getfloat(SECTION, 'MATCH_RADIUS_ARCSEC', fallback=2.)
if not 0 < RADIUS <= 60: raise ValueError('Invalid GALEX radius')
output_dir = path_from_config('PATHS','OUTPUT_DIR','output')
df_source=load_source()
start,end=bounds(df_source,SECTION)
FIELDS=['GALEX_NUV','GALEX_NUV_ERR','GALEX_FUV','GALEX_FUV_ERR','GALEX_SEP_ARCSEC']

def query_one(src, radius):
    coord=SkyCoord(float(src.RA)*u.deg,float(src.DEC)*u.deg)
    table=Catalogs.query_region(coord,radius=radius*u.arcsec,catalog='GALEX')
    rec={k:np.nan for k in FIELDS}
    if len(table):
        if 'distance_arcmin' not in table.colnames:
            raise ValueError('GALEX response missing distance_arcmin')
        sep=np.array([number(v) for v in table['distance_arcmin']])*60
        good=np.flatnonzero(np.isfinite(sep)&(sep<=radius))
        if len(good):
            ix=int(good[np.argmin(sep[good])]); found=table[ix]
            for out,key in [('GALEX_NUV','nuv_mag'),('GALEX_NUV_ERR','nuv_magerr'),
                            ('GALEX_FUV','fuv_mag'),('GALEX_FUV_ERR','fuv_magerr')]:
                rec[out]=number(found[key])
            rec['GALEX_SEP_ARCSEC']=float(sep[ix])
    rows=[]
    for band,val,err in [('NUV','GALEX_NUV','GALEX_NUV_ERR'),
                         ('FUV','GALEX_FUV','GALEX_FUV_ERR')]:
        if np.isfinite(rec[val]):
            rows.append({'OBSERVATORY':'GALEX','EPOCH':'GALEX','BAND':band,
                         'UNITS':'mag','MEASUREMENT':rec[val],'ERROR':rec[err],
                         'POS_ERROR':np.nan,'NEG_ERROR':np.nan,'NEIGHBOR':np.nan})
    found=np.isfinite(rec['GALEX_SEP_ARCSEC'])
    return ('MATCH' if found else 'NO_MATCH'),rec,rows

df_galex, df_galex_measurement, df_galex_failures, df_query_status = managed(
    'GALEX',df_source,query_one,output_dir,RADIUS,SECTION,CONFIG,start,end)
print(f'GALEX: matches={int(df_query_status.STATUS.eq("MATCH").sum())}, '
      f'photometry rows={len(df_galex_measurement)}, failed={len(df_galex_failures)}')
