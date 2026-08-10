
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 16:56:28 2024

@author: josephprendergast
"""

from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import numpy as np


def gaia(source_id, ra, dec):
    
    print("\nStart gaia\n")

    coord = SkyCoord(
        ra=ra * u.degree,
        dec=dec * u.degree,
        frame='icrs'
    )

    radius = 15 * u.arcsec

    query = f"""
    SELECT 
        phot_g_mean_mag,
        phot_g_mean_flux,
        phot_g_mean_flux_error,

        phot_bp_mean_mag,
        phot_bp_mean_flux,
        phot_bp_mean_flux_error,

        phot_rp_mean_mag,
        phot_rp_mean_flux,
        phot_rp_mean_flux_error,

        DISTANCE(
            POINT('ICRS', ra, dec),
            POINT('ICRS', {coord.ra.deg}, {coord.dec.deg})
        ) AS dist_deg,

        ra,
        dec
    FROM 
        gaiadr3.gaia_source
    WHERE 
        CONTAINS(
            POINT('ICRS', ra, dec), 
            CIRCLE('ICRS', {coord.ra.deg}, {coord.dec.deg}, {radius.to(u.deg).value})
        ) = 1
    ORDER BY dist_deg ASC
    """

    job = Gaia.launch_job(query)
    results = job.get_results()

    df_results = results.to_pandas()

    # ------------------------------------------------------------
    # If Gaia query returns no rows, return placeholder rows
    # ------------------------------------------------------------

    if df_results.empty:

        df_g = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': 'Gaia',
            'CAT': 'DR3',
            'BAND': 'G',
            'UNITS': 'mag',
            'MEASUREMENT': None,
            'ERROR': None,
        }])

        df_bp = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': 'Gaia',
            'CAT': 'DR3',
            'BAND': 'BP',
            'UNITS': 'mag',
            'MEASUREMENT': None,
            'ERROR': None,
        }])

        df_rp = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': 'Gaia',
            'CAT': 'DR3',
            'BAND': 'RP',
            'UNITS': 'mag',
            'MEASUREMENT': None,
            'ERROR': None,
        }])

        df_out = pd.concat([df_g, df_bp, df_rp], ignore_index=True)

        return df_out

    # ------------------------------------------------------------
    # If Gaia query returns rows, use closest source
    # ------------------------------------------------------------

    closest_source = df_results.iloc[0]

    GAIA_G_MAG = closest_source['phot_g_mean_mag']
    GAIA_BP_MAG = closest_source['phot_bp_mean_mag']
    GAIA_RP_MAG = closest_source['phot_rp_mean_mag']

    g_flux = closest_source['phot_g_mean_flux']
    bp_flux = closest_source['phot_bp_mean_flux']
    rp_flux = closest_source['phot_rp_mean_flux']

    g_flux_err = closest_source['phot_g_mean_flux_error']
    bp_flux_err = closest_source['phot_bp_mean_flux_error']
    rp_flux_err = closest_source['phot_rp_mean_flux_error']

    factor = 2.5 / np.log(10)

    GAIA_G_MAG_ERROR = factor * g_flux_err / g_flux
    GAIA_BP_MAG_ERROR = factor * bp_flux_err / bp_flux
    GAIA_RP_MAG_ERROR = factor * rp_flux_err / rp_flux

    df_g = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'Gaia',
        'CAT': 'DR3',
        'BAND': 'G',
        'UNITS': 'mag',
        'MEASUREMENT': GAIA_G_MAG,
        'ERROR': GAIA_G_MAG_ERROR,
    }])

    df_bp = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'Gaia',
        'CAT': 'DR3',
        'BAND': 'BP',
        'UNITS': 'mag',
        'MEASUREMENT': GAIA_BP_MAG,
        'ERROR': GAIA_BP_MAG_ERROR,
    }])

    df_rp = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'Gaia',
        'CAT': 'DR3',
        'BAND': 'RP',
        'UNITS': 'mag',
        'MEASUREMENT': GAIA_RP_MAG,
        'ERROR': GAIA_RP_MAG_ERROR,
    }])

    df_out = pd.concat([df_g, df_bp, df_rp], ignore_index=True)

    return df_out

# ra = 270.12583
# dec = -34.60364
# unique_id = 1

# temp = gaia(unique_id, ra, dec)








    
