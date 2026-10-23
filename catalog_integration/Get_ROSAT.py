"""ROSAT 2RXS cone search via HEASARC. Run in current Spyder console."""
from _xray_common import *
from astroquery.heasarc import Heasarc
import warnings
from astroquery.exceptions import NoResultsWarning
warnings.filterwarnings('ignore',category=NoResultsWarning)
heasarc=Heasarc()
def query(ra,dec,radius,timeout):
    return heasarc.query_region(SkyCoord(ra*u.deg,dec*u.deg),catalog='rass2rxs',radius=radius*u.arcsec)
fields={'ROSAT_COUNT_RATE':['Count_Rate'], 'ROSAT_COUNT_RATE_ERR':['Count_Rate_Error'],
        'ROSAT_EXPOSURE':['Exposure'],'ROSAT_DET_LIKE':['Detection_Likelihood'],
        'ROSAT_HR1':['Hardness_Ratio_1'],'ROSAT_HR1_ERR':['Hardness_Ratio_1_Error'],
        'ROSAT_HR2':['Hardness_Ratio_2'],'ROSAT_HR2_ERR':['Hardness_Ratio_2_Error']}
df_rosat,df_rosat_measurement,df_rosat_failures=run('ROSAT',query,fields,
    [('0.1-2 keV','ROSAT_COUNT_RATE','ROSAT_COUNT_RATE_ERR',None,None)],
    'ROSAT','rass2rxs','counts s-1',30)

# Expose status DataFrame in Spyder Variable Explorer.
df_query_status = pd.read_pickle(OUTPUT / 'ROSAT_query_status.pkl')
