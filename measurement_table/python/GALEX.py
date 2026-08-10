# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 07:11:30 2024

@author: Joe
"""

from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import numpy as np

def galex(source_id, ra_deg, dec_deg):
    
    print("\nStart galex\n")

    search_radius = 15 * u.arcsec

    coord = SkyCoord(
        ra=ra_deg * u.deg,
        dec=dec_deg * u.deg,
        frame='icrs'
    )

    galex_catalog = Catalogs.query_region(
        coord,
        radius=search_radius,
        catalog="GALEX"
    )

    df_GALEX = galex_catalog.to_pandas()

    # ------------------------------------------------------------
    # If no GALEX source found, return rows filled with None
    # ------------------------------------------------------------

    if df_GALEX.empty:

        df_NUV = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': 'GALEX',
            'CAT': 'GALEX',
            'BAND': 'NUV',
            'UNITS': 'mag',
            'MEASUREMENT': None,
            'ERROR': None
        }])

        df_FUV = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': 'GALEX',
            'CAT': 'GALEX',
            'BAND': 'FUV',
            'UNITS': 'mag',
            'MEASUREMENT': None,
            'ERROR': None
        }])

        df_out = pd.concat([df_NUV, df_FUV], ignore_index=True)

        return df_out

    # ------------------------------------------------------------
    # Use first returned source
    # ------------------------------------------------------------

    GALEX_NUV = df_GALEX.at[0, 'nuv_mag']
    GALEX_NUV_ERR = df_GALEX.at[0, 'nuv_magerr']

    GALEX_FUV = df_GALEX.at[0, 'fuv_mag']
    GALEX_FUV_ERR = df_GALEX.at[0, 'fuv_magerr']

    df_NUV = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'GALEX',
        'CAT': 'GALEX',
        'BAND': 'NUV',
        'UNITS': 'mag',
        'MEASUREMENT': GALEX_NUV,
        'ERROR': GALEX_NUV_ERR
    }])

    df_FUV = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'GALEX',
        'CAT': 'GALEX',
        'BAND': 'FUV',
        'UNITS': 'mag',
        'MEASUREMENT': GALEX_FUV,
        'ERROR': GALEX_FUV_ERR
    }])

    df_out = pd.concat([df_NUV, df_FUV], ignore_index=True)

    return df_out













