# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 08:47:04 2026

@author: Joe
"""

from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd
import logging


logging.getLogger("astroquery").setLevel(logging.ERROR)


def add_gaia_dist(df_source_interm):
    
    print("\nStart add_gaia_dist\n")
    """
    Adds DISTANCE and DISTANCE_ERROR in kpc using the nearest
    Gaia source within 2 arcseconds.
    """

    # Make sure positional indexing works correctly.
    df_source_interm = df_source_interm.reset_index(drop=True)

    for col in ["DISTANCE", "DISTANCE_ERROR"]:
        if col not in df_source_interm.columns:
            df_source_interm[col] = np.nan

    search_radius_arcsec = 2.0

    Gaia.ROW_LIMIT = 1

    for i in range(len(df_source_interm)):

        # Print only every 100 rows.
        if i % 100 == 0:
            print(f"Processing {i} of {len(df_source_interm)}")

        ra = df_source_interm.at[i, "RA"]
        dec = df_source_interm.at[i, "DEC"]

        if pd.isna(ra) or pd.isna(dec):
            continue

        coord = SkyCoord(
            ra=float(ra) * u.deg,
            dec=float(dec) * u.deg,
            frame="icrs"
        )

        try:

            result = Gaia.cone_search(
                coord,
                radius=search_radius_arcsec * u.arcsec
            ).get_results()

            if result is None or len(result) == 0:
                continue

            parallax = result[0]["parallax"]
            parallax_error = result[0]["parallax_error"]

            if np.ma.is_masked(parallax):
                continue

            parallax = float(parallax)

            if not np.isfinite(parallax) or parallax <= 0:
                continue

            # parallax in mas gives distance in kpc.
            df_source_interm.at[i, "DISTANCE"] = 1.0 / parallax

            if not np.ma.is_masked(parallax_error):

                parallax_error = float(parallax_error)

                if np.isfinite(parallax_error):

                    df_source_interm.at[i, "DISTANCE_ERROR"] = (
                        parallax_error / parallax**2
                    )

        except Exception as e:
            print(f"Row {i}: {e}")

    return df_source_interm