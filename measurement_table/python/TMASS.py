# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 15:28:40 2025

@author: Joe
"""

from astroquery.irsa import Irsa
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import numpy as np


def tmass(source_id, ra_raw, dec_raw):
    
    print("\nStart tmass\n")

    radius = 15 * u.arcsec
    
    Irsa.list_catalogs()
    
    # --- Create coordinate ---
    coord = SkyCoord(ra=ra_raw * u.deg, dec=dec_raw * u.deg, frame='icrs')
    
    # --- Query IRSA ---
    res = Irsa.query_region(coord, catalog="fp_psc", spatial="Cone", radius=radius)
    
    if len(res) == 0:

        df_J = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': '2MASS',
            'CAT': 'fp_psc',
            'BAND': 'J',
            'UNITS': 'mag',
            'MEASUREMENT': None,
            'ERROR': None
        }])
    
        df_H = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': '2MASS',
            'CAT': 'fp_psc',
            'BAND': 'H',
            'UNITS': 'mag',
            'MEASUREMENT': None,
            'ERROR': None
        }])
    
        df_KS = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': '2MASS',
            'CAT': 'fp_psc',
            'BAND': 'KS',
            'UNITS': 'mag',
            'MEASUREMENT': None,
            'ERROR': None
        }])
    
        df_out = pd.concat([df_J, df_H, df_KS], ignore_index=True)
    
        return df_out
    
    
    # --- Get result coordinates (strip to floats first) ---
    res_ra = np.array(res["ra"], dtype=float)
    res_dec = np.array(res["dec"], dtype=float)
    src_coords = SkyCoord(res_ra * u.deg, res_dec * u.deg)
    
    # --- Find closest match ---
    sep = coord.separation(src_coords)
    idx = np.argmin(sep.arcsec)
    
    # --- Save results to df_out ---
    TMASS_J_MAG  = float(res["j_m"][idx])
    TMASS_H_MAG  = float(res["h_m"][idx])
    TMASS_KS_MAG = float(res["k_m"][idx])
    
    TMASS_J_ERR  = float(res["j_cmsig"][idx])
    TMASS_H_ERR  = float(res["h_cmsig"][idx])
    TMASS_KS_ERR = float(res["k_cmsig"][idx])
    
    df_J = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': '2MASS',
        'CAT': 'fp_psc',
        'BAND': 'J',
        'UNITS': 'mag',
        'MEASUREMENT': TMASS_J_MAG,
        'ERROR': TMASS_J_ERR
    }])
    
    df_H = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': '2MASS',
        'CAT': 'fp_psc',
        'BAND': 'H',
        'UNITS': 'mag',
        'MEASUREMENT': TMASS_H_MAG,
        'ERROR': TMASS_H_ERR
    }])
    
    df_KS = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': '2MASS',
        'CAT': 'fp_psc',
        'BAND': 'KS',
        'UNITS': 'mag',
        'MEASUREMENT': TMASS_KS_MAG,
        'ERROR': TMASS_KS_ERR
    }])

    df_out = pd.concat([df_J, df_H, df_KS], ignore_index=True)
    
    return df_out


unique_id = 1
ra = 261.62038
dec = -39.71886

# temp = tmass(unique_id, ra, dec)