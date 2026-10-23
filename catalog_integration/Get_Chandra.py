"""Chandra CSC cone search via CDA VO. Run in current Spyder console."""
from _xray_common import *
def query(ra,dec,radius,timeout):
    return cone_get('https://cda.cfa.harvard.edu/cscvo/coneSearch',
                    {'RA':ra,'DEC':dec,'SR':radius/3600,'VERB':3},timeout)
fields={'flux_aper_b':['flux_aper_b','fluxaper_b'],
        'flux_aper_lolim_b':['flux_aper_lolim_b','fluxaper_lolim_b'],
        'flux_aper_hilim_b':['flux_aper_hilim_b','fluxaper_hilim_b']}
# Derive asymmetric errors after obtaining the confidence bounds; these are not direct catalogue errors.
def with_errors(ra,dec,radius,timeout):
    df=query(ra,dec,radius,timeout)
    if not df.empty:
        lower=column(df,'flux_aper_lolim_b','fluxaper_lolim_b')
        upper=column(df,'flux_aper_hilim_b','fluxaper_hilim_b')
        flux=column(df,'flux_aper_b','fluxaper_b')
        if all(x is not None for x in (lower,upper,flux)):
            df['plus_flux_err']=pd.to_numeric(df[upper],errors='coerce')-pd.to_numeric(df[flux],errors='coerce')
            df['min_flux_err']=pd.to_numeric(df[flux],errors='coerce')-pd.to_numeric(df[lower],errors='coerce')
    return df
fields.update({'plus_flux_err':['plus_flux_err'],'min_flux_err':['min_flux_err']})
df_chandra,df_chandra_measurement,df_chandra_failures=run('Chandra',with_errors,fields,
    [('0.08–10 keV','flux_aper_b',None,'plus_flux_err','min_flux_err')],
    'Chandra','CSC 2.1','erg cm-2 s-1',2)

# Normalize the catalog identifier for both new and RESUMED cached results.
# The previous version used an empty EPOCH; resume can reuse those rows.
if not df_chandra_measurement.empty:
    mask = df_chandra_measurement['EPOCH'].isna() | (
        df_chandra_measurement['EPOCH'].astype('string').str.strip().eq('').fillna(False)
    )
    if mask.any():
        df_chandra_measurement.loc[mask, 'EPOCH'] = 'CSC 2.1'
        df_chandra_measurement.to_pickle(OUTPUT / 'Chandra_measurements.pkl')
        print(f'Chandra: labeled {int(mask.sum())} cached measurements CSC 2.1')

# Expose status DataFrame in Spyder Variable Explorer.
df_query_status = pd.read_pickle(OUTPUT / 'Chandra_query_status.pkl')
