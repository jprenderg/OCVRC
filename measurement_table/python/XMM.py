# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 12:47:18 2026

@author: Joe
"""

from astroquery.esa.xmm_newton import XMMNewton
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd


def xmm(source_id, ra_deg, dec_deg):
    
    print("\nStart xmm\n")

    # ------------------------------------------------------------
    # Settings
    # ------------------------------------------------------------
    
    table_name = "xsa.v_epic_source_cat"
    search_radius_arcsec = 15
    radius_deg = search_radius_arcsec / 3600.0
    
    # ------------------------------------------------------------
    # Defaults
    # ------------------------------------------------------------
    
    ep_1_flux = None
    ep_1_flux_err = None
    ep_2_flux = None
    ep_2_flux_err = None
    ep_3_flux = None
    ep_3_flux_err = None
    ep_4_flux = None
    ep_4_flux_err = None
    ep_5_flux = None
    ep_5_flux_err = None
    
    # ------------------------------------------------------------
    # Query XMM
    # ------------------------------------------------------------
    
    try:
        coord = SkyCoord(
            ra=ra_deg * u.deg,
            dec=dec_deg * u.deg,
            frame="icrs"
        )
    
        query = f"""
        SELECT ra, dec,
               ep_1_flux, ep_1_flux_err,
               ep_2_flux, ep_2_flux_err,
               ep_3_flux, ep_3_flux_err,
               ep_4_flux, ep_4_flux_err,
               ep_5_flux, ep_5_flux_err
        FROM {table_name}
        WHERE 1 = CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_deg})
        )
        """
    
        results = XMMNewton.query_xsa_tap(query, output_format="votable")
        df_results = results.to_pandas()
    
        if not df_results.empty:
    
            coords_results = SkyCoord(
                ra=df_results["ra"].values * u.deg,
                dec=df_results["dec"].values * u.deg,
                frame="icrs"
            )
    
            df_results["sep_arcsec"] = coord.separation(coords_results).arcsec
            closest = df_results.loc[df_results["sep_arcsec"].idxmin()]
    
            ep_1_flux = None if pd.isna(closest["ep_1_flux"]) else float(closest["ep_1_flux"])
            ep_1_flux_err = None if pd.isna(closest["ep_1_flux_err"]) else float(closest["ep_1_flux_err"])
    
            ep_2_flux = None if pd.isna(closest["ep_2_flux"]) else float(closest["ep_2_flux"])
            ep_2_flux_err = None if pd.isna(closest["ep_2_flux_err"]) else float(closest["ep_2_flux_err"])
    
            ep_3_flux = None if pd.isna(closest["ep_3_flux"]) else float(closest["ep_3_flux"])
            ep_3_flux_err = None if pd.isna(closest["ep_3_flux_err"]) else float(closest["ep_3_flux_err"])
    
            ep_4_flux = None if pd.isna(closest["ep_4_flux"]) else float(closest["ep_4_flux"])
            ep_4_flux_err = None if pd.isna(closest["ep_4_flux_err"]) else float(closest["ep_4_flux_err"])
    
            ep_5_flux = None if pd.isna(closest["ep_5_flux"]) else float(closest["ep_5_flux"])
            ep_5_flux_err = None if pd.isna(closest["ep_5_flux_err"]) else float(closest["ep_5_flux_err"])
    
    except Exception as e:
        print(f"XMM query failed for UNIQUE_ID {unique_id}: {e}")
    
    # ------------------------------------------------------------
    # Output dataframe
    # ------------------------------------------------------------
    
    df_1 = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'XMM',
        'CAT': 'xsa.v_epic_source_cat',
        'BAND': '0.2–0.5 keV',
        'UNITS': 'erg cm-2 s-1',
        'MEASUREMENT': ep_1_flux,
        'ERROR': ep_1_flux_err
    }])
    
    df_2 = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'XMM',
        'CAT': 'xsa.v_epic_source_cat',
        'BAND': '0.5–1.0 keV',
        'UNITS': 'erg cm-2 s-1',
        'MEASUREMENT': ep_2_flux,
        'ERROR': ep_2_flux_err
    }])
    
    
    df_3 = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'XMM',
        'CAT': 'xsa.v_epic_source_cat',
        'BAND': '1.0–2.0 keV',
        'UNITS': 'erg cm-2 s-1',
        'MEASUREMENT': ep_3_flux,
        'ERROR': ep_3_flux_err
    }])
    
    
    df_4 = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'XMM',
        'CAT': 'xsa.v_epic_source_cat',
        'BAND': '2.0–4.5 keV',
        'UNITS': 'erg cm-2 s-1',
        'MEASUREMENT': ep_4_flux,
        'ERROR': ep_4_flux_err
    }])
    
    
    df_5 = pd.DataFrame([{
        'SOURCE_ID': source_id,
        'OBSERVATORY': 'XMM',
        'CAT': 'xsa.v_epic_source_cat',
        'BAND': '4.5–12.0 keV',
        'UNITS': 'erg cm-2 s-1',
        'MEASUREMENT': ep_5_flux,
        'ERROR': ep_5_flux_err
    }])
    
    df_out = pd.concat([df_1, df_2, df_3, df_4, df_5], ignore_index=True)
    
    return df_out

unique_id = 1
ra_deg = 267.881
dec_deg = -29.49896
# temp = xmm(unique_id, ra_deg, dec_deg)

