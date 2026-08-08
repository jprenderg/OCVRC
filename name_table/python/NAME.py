# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 07:00:23 2026

@author: Joe
"""

import pandas as pd
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u


def name(df_source):
    """
    Create a name table from a source table.

    Input dataframe columns:
        SOURCE_ID
        RA
        DEC

    Returns dataframe columns:
        SOURCE_ID
        NAME
        DEFAULT_NAME
    """

    search_radius_arcsec = 15.0
    search_radius_deg = search_radius_arcsec / 3600.0

    rows = []

    for i in range(len(df_source)):

        print(f"Processing row {i+1} of {len(df_source)}")

        source_id = df_source.loc[i, "SOURCE_ID"]
        ra0 = df_source.loc[i, "RA"]
        dec0 = df_source.loc[i, "DEC"]

        try:

            # ----------------------------------------------------
            # Find all SIMBAD objects within search radius
            # ----------------------------------------------------

            query = f"""
            SELECT
                basic.oid,
                basic.ra,
                basic.dec,
                basic.main_id
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

            # ----------------------------------------------------
            # Find nearest SIMBAD object
            # ----------------------------------------------------

            source_coord = SkyCoord(
                ra=ra0 * u.deg,
                dec=dec0 * u.deg,
                frame="icrs"
            )

            query_coords = SkyCoord(
                ra=df_query["ra"].values * u.deg,
                dec=df_query["dec"].values * u.deg,
                frame="icrs"
            )

            df_query["SEPARATION_ARCSEC"] = (
                source_coord.separation(query_coords).arcsec
            )

            jmin = df_query["SEPARATION_ARCSEC"].idxmin()

            oid = int(df_query.loc[jmin, "oid"])
            main_id = str(df_query.loc[jmin, "main_id"])

            # ----------------------------------------------------
            # Add default (main) name
            # ----------------------------------------------------

            rows.append({
                "SOURCE_ID": source_id,
                "NAME": main_id,
                "DEFAULT_NAME": True
            })

            # ----------------------------------------------------
            # Query aliases
            # ----------------------------------------------------

            query_alias = f"""
            SELECT id
            FROM ident
            WHERE oidref = {oid}
            """

            table_alias = Simbad.query_tap(query_alias)

            if table_alias is None or len(table_alias) == 0:
                continue

            df_alias = table_alias.to_pandas()

            aliases = (
                df_alias["id"]
                .astype(str)
                .dropna()
                .unique()
                .tolist()
            )

            aliases = sorted(
                alias for alias in aliases
                if alias != main_id
            )

            # ----------------------------------------------------
            # Add aliases
            # ----------------------------------------------------

            for alias in aliases:

                rows.append({
                    "SOURCE_ID": source_id,
                    "NAME": alias,
                    "DEFAULT_NAME": False
                })

        except Exception as e:

            print(f"ERROR at row {i}: {e}")
            continue

    # ------------------------------------------------------------
    # Build output dataframe
    # ------------------------------------------------------------

    df_names = pd.DataFrame(
        rows,
        columns=[
            "SOURCE_ID",
            "NAME",
            "DEFAULT_NAME"
        ]
    )
    
    df_names["NAME"] = df_names["NAME"].str.replace(r"^V\* ", "", regex=True)

    return df_names