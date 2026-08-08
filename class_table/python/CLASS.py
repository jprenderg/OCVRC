# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 18:52:06 2026

@author: Joe
"""

import pandas as pd
import numpy as np
from astroquery.simbad import Simbad


def source_class(source_id, ra, dec):
    """
    Query SIMBAD for classifications within 5 arcsec of one source.

    Parameters
    ----------
    source_id : int
        Source ID to attach to returned classes.

    ra : float
        RA in degrees.

    dec : float
        DEC in degrees.

    Returns
    -------
    df_out : pandas.DataFrame
        Columns: SOURCE_ID, CLASS, DEFAULT_CLASS
    """

    search_radius_arcsec = 5.0
    search_radius_deg = search_radius_arcsec / 3600.0

    empty = pd.DataFrame(
        columns=[
            "SOURCE_ID",
            "CLASS",
            "DEFAULT_CLASS"
        ]
    )

    if pd.isna(ra) or pd.isna(dec):
        return empty

    if not np.isfinite(ra) or not np.isfinite(dec):
        return empty

    query = f"""
        SELECT
            basic.oid,
            basic.ra,
            basic.dec,
            basic.main_id,
            basic.otype AS primary_class,
            otypes.otype AS other_class
        FROM basic
        LEFT JOIN otypes
            ON basic.oid = otypes.oidref
        WHERE
            1 = CONTAINS(
                POINT('ICRS', basic.ra, basic.dec),
                CIRCLE('ICRS', {float(ra)}, {float(dec)}, {search_radius_deg})
            )
    """

    rows = []

    try:
        result = Simbad.query_tap(query)

        if result is None or len(result) == 0:
            return empty

        df_result = result.to_pandas()

        df_temp = df_result.copy()

        df_temp["RA_DIFF"] = (
            (df_temp["ra"].astype(float) - float(ra))
            * np.cos(np.radians(float(dec)))
        )

        df_temp["DEC_DIFF"] = df_temp["dec"].astype(float) - float(dec)

        df_temp["SEP_DEG"] = np.sqrt(
            df_temp["RA_DIFF"] ** 2 +
            df_temp["DEC_DIFF"] ** 2
        )

        df_temp = df_temp[
            df_temp["SEP_DEG"] <= search_radius_deg
        ]

        if len(df_temp) == 0:
            return empty

        closest_oid = df_temp.sort_values("SEP_DEG").iloc[0]["oid"]

        df_match = df_temp[
            df_temp["oid"] == closest_oid
        ].copy()

        primary_class = df_match["primary_class"].iloc[0]

        if pd.notna(primary_class):
            rows.append({
                "SOURCE_ID": int(source_id),
                "CLASS": primary_class,
                "DEFAULT_CLASS": 1
            })

        for cls in df_match["other_class"].dropna().unique():

            if cls == primary_class:
                continue

            rows.append({
                "SOURCE_ID": int(source_id),
                "CLASS": cls,
                "DEFAULT_CLASS": 0
            })

    except Exception as e:
        print("SIMBAD class query error:", e)
        return empty

    df_out = pd.DataFrame(rows)

    if len(df_out) == 0:
        return empty

    df_out = (
        df_out
        .drop_duplicates()
        .sort_values(
            ["SOURCE_ID", "DEFAULT_CLASS", "CLASS"],
            ascending=[True, False, True]
        )
        .reset_index(drop=True)
    )
    
    df_out["CLASS"] = df_out["CLASS"].replace("CV*", "CV")
    df_out["CLASS"] = df_out["CLASS"].replace("*", "star")

    return df_out