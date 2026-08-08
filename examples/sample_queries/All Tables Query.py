# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 18:51:16 2026

@author: Joe
"""

import sqlite3
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repo_root))

from config import OUTPUT_DB, MASTER_DB


# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

# Reads from database update
db_path = OUTPUT_DB

# Read from master database
# db_path = MASTER_DB

ra = 304.39133
dec = -3.66411

search_radius_arcsec = 300.0


# ------------------------------------------------------------
# Open database
# ------------------------------------------------------------

conn = sqlite3.connect(db_path)


# ------------------------------------------------------------
# Load source_table and find closest SOURCE_ID
# ------------------------------------------------------------

df_source_all = pd.read_sql_query(
    """
    SELECT
        SOURCE_ID,
        PRIMARY_FLAG,
        RA,
        DEC,
        MAJOR_AXIS,
        MINOR_AXIS,
        POSITION_ANGLE,
        DISTANCE,
        DISTANCE_ERROR
    FROM source_table
    WHERE RA IS NOT NULL
      AND DEC IS NOT NULL
    """,
    conn
)

if df_source_all.empty:
    conn.close()
    raise ValueError("source_table contains no sources with valid coordinates.")

target = SkyCoord(
    ra=ra * u.deg,
    dec=dec * u.deg,
    frame="icrs"
)

coords = SkyCoord(
    ra=df_source_all["RA"].astype(float).values * u.deg,
    dec=df_source_all["DEC"].astype(float).values * u.deg,
    frame="icrs"
)

sep_arcsec = target.separation(coords).arcsec
idx = int(np.argmin(sep_arcsec))

source_id = int(df_source_all.iloc[idx]["SOURCE_ID"])
closest_sep = float(sep_arcsec[idx])

print("\nClosest primary source:")
print("SOURCE_ID =", source_id)
print("Separation arcsec =", closest_sep)


# ------------------------------------------------------------
# Load primary source row
# ------------------------------------------------------------

df_source = pd.read_sql_query(
    """
    SELECT *
    FROM source_table
    WHERE SOURCE_ID = ?
    """,
    conn,
    params=(source_id,)
)

if df_source.empty:
    conn.close()
    raise ValueError(f"SOURCE_ID {source_id} was not found.")

print("\nPrimary Source_Table row:")
print(df_source.to_string(index=False))


# ------------------------------------------------------------
# Find nearby non-primary sources within 5 arcminutes
# ------------------------------------------------------------

search_radius_deg = search_radius_arcsec / 3600.0

ra0 = float(df_source.iloc[0]["RA"])
dec0 = float(df_source.iloc[0]["DEC"])

cos_dec = np.cos(np.radians(dec0))

ra_min = ra0 - search_radius_deg / cos_dec
ra_max = ra0 + search_radius_deg / cos_dec
dec_min = dec0 - search_radius_deg
dec_max = dec0 + search_radius_deg

df_neighbors = pd.read_sql_query(
    """
    SELECT *
    FROM source_table
    WHERE PRIMARY_FLAG = 0
      AND RA BETWEEN ? AND ?
      AND DEC BETWEEN ? AND ?
    """,
    conn,
    params=(ra_min, ra_max, dec_min, dec_max)
)

if not df_neighbors.empty:

    c_primary = SkyCoord(
        ra=ra0 * u.deg,
        dec=dec0 * u.deg,
        frame="icrs"
    )

    c_neighbors = SkyCoord(
        ra=df_neighbors["RA"].astype(float).values * u.deg,
        dec=df_neighbors["DEC"].astype(float).values * u.deg,
        frame="icrs"
    )

    df_neighbors["SEPARATION_ARCSEC"] = (
        c_primary.separation(c_neighbors).arcsec
    )

    df_neighbors = df_neighbors[
        (df_neighbors["SEPARATION_ARCSEC"] <= search_radius_arcsec)
        & (df_neighbors["SOURCE_ID"] != source_id)
    ]

    df_neighbors = (
        df_neighbors
        .sort_values("SEPARATION_ARCSEC")
        .reset_index(drop=True)
    )


print("\nNearby PRIMARY_FLAG = 0 sources within 5 arcminutes:")

if df_neighbors.empty:
    print("None found.")
else:
    print(df_neighbors.to_string(index=False))


# ------------------------------------------------------------
# SOURCE_ID values for primary and neighbors
# ------------------------------------------------------------

related_source_ids = [source_id]

if not df_neighbors.empty:
    related_source_ids.extend(
        df_neighbors["SOURCE_ID"].astype(int).tolist()
    )

# Remove duplicate IDs while preserving order.
related_source_ids = list(dict.fromkeys(related_source_ids))

print("\nPrimary and neighbor SOURCE_ID values:")
print(related_source_ids)

placeholders = ",".join("?" for _ in related_source_ids)


# ------------------------------------------------------------
# Name Table: primary and neighbors
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("Name_Table — primary and neighbors")
print("-" * 80)

df_name = pd.read_sql_query(
    f"""
    SELECT
        s.PRIMARY_FLAG,
        n.*
    FROM name_table AS n
    JOIN source_table AS s
        ON n.SOURCE_ID = s.SOURCE_ID
    WHERE n.SOURCE_ID IN ({placeholders})
    ORDER BY
        s.PRIMARY_FLAG DESC,
        n.SOURCE_ID,
        n.DEFAULT_NAME DESC,
        n.NAME
    """,
    conn,
    params=related_source_ids
)

if df_name.empty:
    print("No rows found.")
else:
    print(df_name.to_string(index=False))


# ------------------------------------------------------------
# Measurement Table: primary and neighbors
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("Measurement_Table — primary and neighbors")
print("-" * 80)

df_measurement = pd.read_sql_query(
    f"""
    SELECT
        s.PRIMARY_FLAG,
        m.*
    FROM measurement_table AS m
    JOIN source_table AS s
        ON m.SOURCE_ID = s.SOURCE_ID
    WHERE m.SOURCE_ID IN ({placeholders})
    ORDER BY
        s.PRIMARY_FLAG DESC,
        m.SOURCE_ID,
        m.OBSERVATORY,
        m.EPOCH,
        m.BAND
    """,
    conn,
    params=related_source_ids
)

if df_measurement.empty:
    print("No rows found.")
else:
    print(df_measurement.to_string(index=False))


# ------------------------------------------------------------
# Class Table: primary and neighbors
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("Class_Table — primary and neighbors")
print("-" * 80)

df_class = pd.read_sql_query(
    f"""
    SELECT
        s.PRIMARY_FLAG,
        c.*
    FROM class_table AS c
    JOIN source_table AS s
        ON c.SOURCE_ID = s.SOURCE_ID
    WHERE c.SOURCE_ID IN ({placeholders})
    ORDER BY
        s.PRIMARY_FLAG DESC,
        c.SOURCE_ID,
        c.DEFAULT_CLASS DESC,
        c.CLASS
    """,
    conn,
    params=related_source_ids
)

if df_class.empty:
    print("No rows found.")
else:
    print(df_class.to_string(index=False))


# ------------------------------------------------------------
# Period Table: primary source only
# ------------------------------------------------------------

print("\n" + "-" * 80)
print("Period_Table — primary source only")
print("-" * 80)

df_period = pd.read_sql_query(
    """
    SELECT *
    FROM period_table
    WHERE SOURCE_ID = ?
    ORDER BY TYPE, PERIOD
    """,
    conn,
    params=(source_id,)
)

if df_period.empty:
    print("No rows found.")
else:
    print(df_period.to_string(index=False))


# ------------------------------------------------------------
# Close database
# ------------------------------------------------------------

conn.close()