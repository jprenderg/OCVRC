# -*- coding: utf-8 -*-
"""
Created on Mon May 25 21:23:57 2026

@author: Joe
"""

import requests
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io.votable import parse_single_table
from io import BytesIO


def chandra(source_id, ra_deg, dec_deg):
    
    print("\nStart chandra\n")

    search_radius_arcsec = 15

    # ------------------------------------------------------------
    # Chandra CSC cone search
    # ------------------------------------------------------------

    base_url = "https://cda.cfa.harvard.edu/cscvo/coneSearch"

    radius_deg = search_radius_arcsec / 3600.0

    coord = SkyCoord(
        ra=ra_deg * u.deg,
        dec=dec_deg * u.deg
    )

    params = {
        "RA": ra_deg,
        "DEC": dec_deg,
        "SR": radius_deg,
        "VERB": 3
    }

    response = requests.get(
        base_url,
        params=params,
        timeout=60
    )

    response.raise_for_status()

    vot = parse_single_table(BytesIO(response.content))

    df_results = vot.to_table().to_pandas()

    # ------------------------------------------------------------
    # Handle empty query result
    # ------------------------------------------------------------

    if df_results.empty:

        df_out = pd.DataFrame([{
            'SOURCE_ID': source_id,
            'OBSERVATORY': 'Chandra',
            'CAT': '',
            'BAND': '0.08–10 keV',
            'UNITS': 'erg cm-2 s-1',
            'MEASUREMENT': None,
            'POS_ERROR': None,
            'NEG_ERROR': None
            
        }])

        return df_out

    # ------------------------------------------------------------
    # Find closest source
    # ------------------------------------------------------------

    coords_results = SkyCoord(
        ra=pd.to_numeric(
            df_results["ra"],
            errors="coerce"
        ).to_numpy() * u.deg,

        dec=pd.to_numeric(
            df_results["dec"],
            errors="coerce"
        ).to_numpy() * u.deg
    )

    df_results["sep_arcsec"] = (
        coord.separation(coords_results).arcsec
    )

    df_closest = (
        df_results
        .sort_values("sep_arcsec")
        .reset_index(drop=True)
        .iloc[0]
    )

    # ------------------------------------------------------------
    # Extract flux and asymmetric errors
    # ------------------------------------------------------------

    flux_aper_b = df_closest["flux_aper_b"]

    flux_aper_lolim_b = df_closest["flux_aper_lolim_b"]

    flux_aper_hilim_b = df_closest["flux_aper_hilim_b"]

    plus_flux_err = (
        flux_aper_hilim_b - flux_aper_b
    )

    minus_flux_err = (
        flux_aper_b - flux_aper_lolim_b
    )


    # ------------------------------------------------------------
    # Output dataframe
    # ------------------------------------------------------------

    df_out = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'Chandra',
        'CAT': '',
        'BAND': '0.08–10 keV',
        'UNITS': 'erg cm-2 s-1',
        'MEASUREMENT': flux_aper_b,
        'POS_ERROR': plus_flux_err,
        'NEG_ERROR': minus_flux_err,
    }])

    return df_out


