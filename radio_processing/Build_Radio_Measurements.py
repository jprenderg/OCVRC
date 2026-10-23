"""Combine VLASS and RACS fitted fluxes. All finite observations kept."""
from pathlib import Path
import configparser
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
CFG = configparser.ConfigParser(interpolation=None)
CFG.read(ROOT/'config.txt',encoding='utf-8')
image_dir = Path(CFG.get('PATHS','IMAGE_OUTPUT_DIR',fallback='output/images'))
image_dir = image_dir if image_dir.is_absolute() else ROOT/image_dir
output_dir = Path(CFG.get('PATHS','OUTPUT_DIR',fallback='output'))
output_dir = output_dir if output_dir.is_absolute() else ROOT/output_dir
output_dir.mkdir(parents=True,exist_ok=True)
cols=['SOURCE_ID','OBSERVATORY','EPOCH','BAND','UNITS','MEASUREMENT',
      'ERROR','POS_ERROR','NEG_ERROR','NEIGHBOR']
frames=[]
for survey in ('VLASS','RACS'):
    path=image_dir/f'{survey}_flux_measurements.pkl'
    if not path.exists():
        raise FileNotFoundError(path)
    frame=pd.read_pickle(path).copy()
    required={'SOURCE_ID','MEASUREMENT','ERROR','NEIGHBOR_DETECTION'}
    if not required.issubset(frame.columns):
        raise ValueError(f'{survey}: missing {required-set(frame.columns)}')
    frame['OBSERVATORY'] = 'VLA' if survey=='VLASS' else 'ASKAP'
    frame['EPOCH'] = frame['EPOCH'] if survey=='VLASS' else 'RACS'
    frame['BAND'] = '2-4 GHz' if survey=='VLASS' else '887.5 MHz'
    frame['UNITS']='mJy'
    frame['NEIGHBOR']=frame['NEIGHBOR_DETECTION']
    frame['POS_ERROR']=np.nan
    frame['NEG_ERROR']=np.nan
    frames.append(frame[cols])
df_measurement = pd.concat(frames,ignore_index=True)
df_measurement['MEASUREMENT']=pd.to_numeric(df_measurement['MEASUREMENT'],errors='coerce')
df_measurement=df_measurement.dropna(subset=['MEASUREMENT']).reset_index(drop=True)
path_out=output_dir/'measurement.pkl'
df_measurement.to_pickle(path_out)
print(f'Saved {len(df_measurement)} VLASS/RACS measurements to {path_out}')
print('Flagged observations kept:',int((df_measurement['NEIGHBOR']==1).sum()))
