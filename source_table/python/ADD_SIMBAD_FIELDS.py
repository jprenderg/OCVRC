# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 19:30:38 2026

@author: Joe
"""

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd


def add_simbad_fields(df_source_interm):
    """
    Adds the following columns to df_source_interm:

        SOURCE_NAME
        MAJOR_AXIS
        MINOR_AXIS
        POSITION_ANGLE

    using the nearest SIMBAD source within search_radius_arcsec.
    """
    
    print("\nStart add_simbad_fields\n")

    # Create output columns
    new_cols = [
        "SOURCE_NAME",
        "MAJOR_AXIS",
        "MINOR_AXIS",
        "POSITION_ANGLE"
    ]

    for col in new_cols:
        df_source_interm[col] = np.nan

    search_radius_arcsec = 2.0
    search_radius_deg = search_radius_arcsec / 3600.0

    for i in range(len(df_source_interm)):

        print(f"Processing {i+1} of {len(df_source_interm)}")

        ra0 = df_source_interm.at[i, "RA"]
        dec0 = df_source_interm.at[i, "DEC"]

        try:

            query = f"""
            SELECT
                basic.main_id,
                basic.ra,
                basic.dec,
                basic.galdim_majaxis,
                basic.galdim_minaxis,
                basic.galdim_angle
            FROM basic
            WHERE 1 = CONTAINS(
                POINT('ICRS', basic.ra, basic.dec),
                CIRCLE('ICRS', {ra0}, {dec0}, {search_radius_deg})
            )
            """

            table = Simbad.query_tap(query)

            if table is None or len(table) == 0:
                continue

            df_query = table.to_pandas()

            source_coord = SkyCoord(
                ra=ra0 * u.deg,
                dec=dec0 * u.deg,
                frame="icrs"
            )

            simbad_coords = SkyCoord(
                ra=df_query["ra"].values * u.deg,
                dec=df_query["dec"].values * u.deg,
                frame="icrs"
            )

            df_query["SEPARATION"] = (
                source_coord.separation(simbad_coords).arcsec
            )

            j = df_query["SEPARATION"].idxmin()

            df_source_interm.at[i, "SOURCE_NAME"] = df_query.at[j, "main_id"]
            df_source_interm.at[i, "MAJOR_AXIS"] = df_query.at[j, "galdim_majaxis"]
            df_source_interm.at[i, "MINOR_AXIS"] = df_query.at[j, "galdim_minaxis"]
            df_source_interm.at[i, "POSITION_ANGLE"] = df_query.at[j, "galdim_angle"]

        except Exception as e:
            print(f"Row {i}: {e}")

    return df_source_interm