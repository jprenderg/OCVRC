"""2MASS point-source photometry with per-source checkpoint and resume."""
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.ipac.irsa import Irsa
from _photometry_common import CONFIG, number, load_source, bounds, nearest, path_from_config
from _catalog_query_common import managed
SECTION = '2MASS'
RADIUS = CONFIG.getfloat(SECTION, 'MATCH_RADIUS_ARCSEC', fallback=2.)
if not 0 < RADIUS <= 60: raise ValueError('Invalid 2MASS radius')
output_dir = path_from_config('PATHS', 'OUTPUT_DIR', 'output')
df_source = load_source()
start, end = bounds(df_source, SECTION)
FIELDS = ['TMASS_J_MAG','TMASS_H_MAG','TMASS_KS_MAG','TMASS_J_ERR',
          'TMASS_H_ERR','TMASS_KS_ERR','TMASS_SEP_ARCSEC']

def query_one(src, radius):
    coord = SkyCoord(float(src.RA)*u.deg, float(src.DEC)*u.deg)
    table = Irsa.query_region(coord, catalog='fp_psc', spatial='Cone', radius=radius*u.arcsec)
    found, sep = nearest(table, coord, radius_arcsec=radius)
    rec = {k: np.nan for k in FIELDS}
    if found is not None:
        for out,key in [('TMASS_J_MAG','j_m'),('TMASS_H_MAG','h_m'),
                        ('TMASS_KS_MAG','k_m'),('TMASS_J_ERR','j_cmsig'),
                        ('TMASS_H_ERR','h_cmsig'),('TMASS_KS_ERR','k_cmsig')]:
            rec[out] = number(found[key])
        rec['TMASS_SEP_ARCSEC'] = sep
    rows=[]
    for band, val, err in [('J','TMASS_J_MAG','TMASS_J_ERR'),
                            ('H','TMASS_H_MAG','TMASS_H_ERR'),
                            ('KS','TMASS_KS_MAG','TMASS_KS_ERR')]:
        if np.isfinite(rec[val]):
            rows.append({'OBSERVATORY':'2MASS','EPOCH':'fp_psc', 'BAND':band,
                         'UNITS':'mag','MEASUREMENT':rec[val],'ERROR':rec[err],
                         'POS_ERROR':np.nan,'NEG_ERROR':np.nan,'NEIGHBOR':np.nan})
    return ('MATCH' if found is not None else 'NO_MATCH'), rec, rows

df_2mass, df_2mass_measurement, df_2mass_failures, df_query_status = managed(
    '2MASS', df_source, query_one, output_dir, RADIUS, SECTION, CONFIG, start, end)
print(f'2MASS: matches={int(df_query_status.STATUS.eq("MATCH").sum())}, '
      f'photometry rows={len(df_2mass_measurement)}, failed={len(df_2mass_failures)}')
