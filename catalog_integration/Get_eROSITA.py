"""eROSITA DR2_Main (eRASS:3) via MPE cone-search service. Run in Spyder."""
from _xray_common import *
BANDS={'1':'0.2-2.3 keV','P1':'0.2-0.5 keV','P2':'0.5-1.0 keV',
       'P3':'1.0-2.0 keV','P4':'2.0-5.0 keV','P5':'5.0-8.0 keV',
       'P6':'4.0-10.0 keV','P7':'5.1-6.1 keV','P8':'6.2-7.1 keV','P9':'7.2-8.2 keV'}
def query(ra,dec,radius,timeout):
    return cone_get('https://erosita.mpe.mpg.de/erodat/catalogue/SCS',
            {'CAT':'DR2_Main','RA':ra,'DEC':dec,'SR':radius/3600,'VERB':3},timeout)
fields={'pos_err':['pos_err'],'det_like':['det_like']}
for band in BANDS:
    for out,original in [('flux','ML_FLUX'),('flux_err','ML_FLUX_ERR'),
                         ('plus_flux_err','ML_FLUX_UPERR'),('min_flux_err','ML_FLUX_LOWERR')]:
        fields[f'{out}_{band}']=[f'{original}_{band}']
bands=[(energy,f'flux_{band}',f'flux_err_{band}',
        f'plus_flux_err_{band}',f'min_flux_err_{band}') for band,energy in BANDS.items()]
df_erosita,df_erosita_measurement,df_erosita_failures=run('eROSITA',query,fields,bands,
    'eRosita','eRASS:3 DR2','erg cm-2 s-1',10)

# Expose status DataFrame in Spyder Variable Explorer.
df_query_status = pd.read_pickle(OUTPUT / 'eROSITA_query_status.pkl')
