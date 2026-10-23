# -*- coding: utf-8 -*-
"""Download RACS DR1 FITS cutouts through CASDA; preserve OCVRC PSF input columns.

Requires astroquery >= 0.4.7 and a free CSIRO OPAL account.
Run in Spyder's current console to inspect df_* variables.
"""
from pathlib import Path
import configparser
import os
import time
import warnings

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.casda import Casda

# ============================================================
# CONFIGURATION
# ============================================================
ROOT = Path(__file__).resolve().parents[1]
CFG = configparser.ConfigParser(interpolation=None)
if not CFG.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')


def cfg_path(section, option):
    p = Path(CFG.get(section, option)).expanduser()
    return p if p.is_absolute() else ROOT / p


path_catalog = cfg_path('PATHS', 'SOURCE_PKL')
output_dir = cfg_path('PATHS', 'IMAGE_OUTPUT_DIR')
output_dir.mkdir(parents=True, exist_ok=True)
cache_dir = cfg_path('PATHS', 'RACS_CACHE_DIR') if CFG.has_option('PATHS', 'RACS_CACHE_DIR') else output_dir / 'racs_cache'
cache_dir.mkdir(parents=True, exist_ok=True)
BOX_SIZE = CFG.getint('RACS', 'BOX_SIZE', fallback=47)
CUTOUT_RADIUS_ARCMIN = CFG.getfloat('RACS', 'CUTOUT_RADIUS_ARCMIN', fallback=2.0)
USE_CACHE = CFG.getboolean('RACS', 'USE_CUTOUT_CACHE', fallback=True)
CASDA_USERNAME = os.environ.get('CASDA_USERNAME', CFG.get('RACS', 'CASDA_USERNAME', fallback='')).strip()
START_ROW = CFG.getint('IMAGE_RETRIEVAL', 'START_ROW', fallback=0)
end_txt = CFG.get('IMAGE_RETRIEVAL', 'END_ROW', fallback='-1').strip()
END_ROW = -1 if end_txt.upper() in ('NONE', '') else int(end_txt)
INDEX = CFG.getint('IMAGE_RETRIEVAL', 'DISPLAY_INDEX', fallback=-1)
if BOX_SIZE < 3 or BOX_SIZE % 2 != 1:
    raise ValueError('RACS BOX_SIZE must be an odd integer >= 3')
if CUTOUT_RADIUS_ARCMIN <= 0:
    raise ValueError('CUTOUT_RADIUS_ARCMIN must be positive')

# ============================================================
# READ CATALOG
# ============================================================
df_catalog = pd.read_pickle(path_catalog).reset_index(drop=True)
required = {'SOURCE_ID', 'RA', 'DEC'}
if not required.issubset(df_catalog.columns):
    raise ValueError(f'Missing catalog columns: {sorted(required - set(df_catalog.columns))}')
n_catalog = len(df_catalog)
START_ROW = max(0, START_ROW)
END_ROW = n_catalog if END_ROW <= 0 else min(END_ROW, n_catalog)
if START_ROW >= END_ROW:
    raise ValueError(f'Invalid loop bounds: START_ROW={START_ROW}, END_ROW={END_ROW}, catalog size={n_catalog}')
df_processed = df_catalog.iloc[START_ROW:END_ROW].copy().reset_index(drop=True)
if not np.isfinite(df_processed[['RA', 'DEC']].to_numpy(dtype=float)).all():
    raise ValueError('Non-finite source coordinates')
if not df_processed['DEC'].between(-90, 90).all():
    raise ValueError('Declination outside -90 to +90 degrees')
partial_run = START_ROW != 0 or END_ROW != n_catalog
suffix = '_test' if partial_run else ''
path_out = output_dir / f'RACS_images{suffix}.pkl'
path_bad = output_dir / f'RACS_bad_images{suffix}.pkl'
path_missing = output_dir / f'RACS_missing_sources{suffix}.pkl'
print(f'\nCATALOG\n========================================\nTotal sources: {n_catalog}\nStart row: {START_ROW}\nEnd row: {END_ROW}\nSources to process: {len(df_processed)}')

# ============================================================
# CASDA QUERY AND DOWNLOAD
# CASDA staging/cutouts require free OPAL authentication.
# Do not store OPAL passwords in config.txt.
# ============================================================
casda = Casda()
logged_in = False


def ensure_login():
    global logged_in, CASDA_USERNAME
    if logged_in:
        return
    if not CASDA_USERNAME:
        CASDA_USERNAME = input('CASDA OPAL username (email): ').strip()
    if not CASDA_USERNAME:
        raise RuntimeError('CASDA OPAL username is required to request cutouts.')
    casda.login(username=CASDA_USERNAME)  # secure interactive password/keyring handling
    logged_in = True


def available_racs_products(centre):
    # Query is positional; avoid selecting the nearest tile based solely on its centre.
    found = casda.query_region(centre, radius=1*u.arcmin)
    if found is None or len(found) == 0:
        return found
    public = casda.filter_out_unreleased(found)
    if len(public) == 0:
        return public
    names = np.asarray(public['filename']).astype(str)
    collections = np.asarray(public['obs_collection']).astype(str)
    mask = (collections == 'The Rapid ASKAP Continuum Survey') & np.char.startswith(names, 'RACS-DR1_') & np.char.endswith(names, 'A.fits')
    return public[mask]


def obtain_fits(source_id, ra, dec):
    # Coordinate in cache key prevents stale reuse when SOURCE_ID is reassigned.
    cache_file = cache_dir / f'RACS_{source_id}_{ra:.7f}_{dec:.7f}.fits'
    if USE_CACHE and cache_file.is_file():
        return cache_file, 'CACHE'
    centre = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    selected = available_racs_products(centre)
    if selected is None or len(selected) == 0:
        raise LookupError('No released RACS DR1 convolved tile at source coordinates')
    # Prefer nearest tile centre if multiple overlapping public convolved images exist.
    if len(selected) > 1 and {'s_ra', 's_dec'}.issubset(selected.colnames):
        coords = SkyCoord(np.asarray(selected['s_ra'], dtype=float)*u.deg,
                          np.asarray(selected['s_dec'], dtype=float)*u.deg)
        chosen = int(np.argmin(centre.separation(coords).deg))
        selected = selected[chosen:chosen+1]
    ensure_login()
    urls = casda.cutout(selected, coordinates=centre, radius=CUTOUT_RADIUS_ARCMIN*u.arcmin)
    if not urls:
        raise RuntimeError('CASDA returned no cutout download URLs')
    download_dir = cache_dir / f'_download_{source_id}'
    download_dir.mkdir(exist_ok=True)
    files = casda.download_files(urls, savedir=str(download_dir))
    fits_paths = [Path(p) for p in files if str(p).lower().endswith(('.fits', '.fits.gz'))]
    if not fits_paths:
        # Some CASDA versions return a file listing or include checksums.
        fits_paths = [p for p in download_dir.rglob('*') if p.is_file() and p.name.lower().endswith(('.fits', '.fits.gz'))]
    if not fits_paths:
        raise RuntimeError(f'CASDA returned no FITS files: {files}')
    src = fits_paths[0]
    if src.suffix.lower() == '.gz':
        with fits.open(src) as hdul:
            hdul.writeto(cache_file, overwrite=True)
    else:
        src.replace(cache_file)
    return cache_file, str(selected['filename'][0])


def extract_cutout(fits_path, source_id, ra, dec, filename):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        with fits.open(fits_path, memmap=False) as hdul:
            hdu = next((h for h in hdul if h.data is not None), None)
            if hdu is None:
                raise ValueError('FITS has no image data')
            header = hdu.header
            image = np.squeeze(np.asarray(hdu.data))
            if image.ndim != 2:
                raise ValueError(f'Unexpected FITS dimensions: {image.shape}')
            wcs = WCS(header).celestial
            bmaj, bmin, bpa = (float(header[k]) for k in ('BMAJ', 'BMIN', 'BPA'))
            px, py = np.abs(proj_plane_pixel_scales(wcs))
            x, y = wcs.world_to_pixel_values(ra, dec)
            if not np.isfinite([x, y]).all():
                raise ValueError('Non-finite WCS source location')
            xi, yi = int(np.round(x)), int(np.round(y))
            half = BOX_SIZE // 2
            ny, nx = image.shape
            if not (half <= xi < nx-half and half <= yi < ny-half):
                raise ValueError('Source outside complete cutout boundaries')
            cutout = image[yi-half:yi+half+1, xi-half:xi+half+1].copy()
            if cutout.shape != (BOX_SIZE, BOX_SIZE) or not np.isfinite(cutout).all():
                raise ValueError('Incomplete or non-finite cutout')
            return dict(SOURCE_ID=source_id, RA=ra, DEC=dec, FILENAME=filename,
                        SOURCE_X=half+x-xi, SOURCE_Y=half+y-yi,
                        BMAJ=bmaj, BMIN=bmin, BPA=bpa,
                        PIXSCALE_X=float(px), PIXSCALE_Y=float(py), IMAGE=cutout)

# ============================================================
# PROCESS SOURCES
# ============================================================
columns = ['SOURCE_ID', 'RA', 'DEC', 'FILENAME', 'SOURCE_X', 'SOURCE_Y',
           'BMAJ', 'BMIN', 'BPA', 'PIXSCALE_X', 'PIXSCALE_Y', 'IMAGE']
results, rejected, tile_rows = [], [], []
start_time = time.monotonic()
for idx, row in df_processed.iterrows():
    sid, ra, dec = row['SOURCE_ID'], float(row['RA']), float(row['DEC'])
    print(f'Row {idx}/{len(df_processed)} SOURCE_ID {sid}', flush=True)
    try:
        fits_path, filename = obtain_fits(sid, ra, dec)
        results.append(extract_cutout(fits_path, sid, ra, dec, filename))
        tile_rows.append(dict(SOURCE_ID=sid, FITS_FILE=filename, FITS_PATH=str(fits_path)))
    except Exception as exc:
        print(f'RACS exception SOURCE_ID {sid}: {exc}', flush=True)
        rejected.append(dict(INDEX=idx, SOURCE_ID=sid, FITS_FILE='', REASON=str(exc)))

df_final = pd.DataFrame(results, columns=columns).sort_values('SOURCE_ID').reset_index(drop=True)
df_bad = pd.DataFrame(rejected, columns=['INDEX', 'SOURCE_ID', 'FITS_FILE', 'REASON'])
df_missing = df_processed[~df_processed['SOURCE_ID'].isin(df_final['SOURCE_ID'])].copy().reset_index(drop=True)
df_tiles = pd.DataFrame(tile_rows, columns=['SOURCE_ID', 'FITS_FILE', 'FITS_PATH'])
df_final.to_pickle(path_out)
df_bad.to_pickle(path_bad)
df_missing.to_pickle(path_missing)

df_racs_images = df_final
df_racs_bad_images = df_bad
df_racs_missing = df_missing
df_racs_tiles = df_tiles
print(f'\nRACS IMAGE SUMMARY\n========================================\nCatalog sources: {n_catalog}\nProcessed sources: {len(df_processed)}\nStored cutouts: {len(df_final)}\nUnique sources with images: {df_final.SOURCE_ID.nunique()}\nSources without images: {len(df_missing)}\nRejected cutouts: {len(df_bad)}\nElapsed time: {(time.monotonic()-start_time)/60:.1f} min')
if not df_bad.empty:
    print('\nREJECTION REASONS\n', df_bad.REASON.value_counts().head(20))
print('\nSAVED FILES\n' + '\n'.join(map(str, (path_out, path_bad, path_missing))))
if INDEX >= 0 and INDEX < len(df_final):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.imshow(df_final.iloc[INDEX]['IMAGE'], origin='lower', cmap='gray')
    plt.title(f"RACS SOURCE_ID {df_final.iloc[INDEX]['SOURCE_ID']}")
    plt.colorbar(label='Jy/beam')
    plt.show()
else:
    print('Image display disabled.' if INDEX < 0 else 'DISPLAY_INDEX outside available image range.')
