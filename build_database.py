from pathlib import Path
import sqlite3
import numpy as np
import pandas as pd
import gc
import sys

import astropy.units as u
from astropy.coordinates import SkyCoord

from config import MASTER_DB, OUTPUT_DB, OUTPUT_SUMMARY

from source_table.python.VLASS_NEIGHBORS import vlass_neighbors
from source_table.python.RACS_NEIGHBORS import racs_neighbors
from source_table.python.ADD_SIMBAD_FIELDS import add_simbad_fields
from source_table.python.ADD_GAIA_DIST import add_gaia_dist

from name_table.python.NAME import name

from measurement_table.python.VLASS import vlass
from measurement_table.python.RACS import racs
from measurement_table.python.TMASS import tmass
from measurement_table.python.GAIA import gaia
from measurement_table.python.GALEX import galex
from measurement_table.python.CHANDRA import chandra
from measurement_table.python.ROSAT import rosat
from measurement_table.python.XMM import xmm

from class_table.python.CLASS import source_class


# ------------------------------------------------------------
# User inputs
# ------------------------------------------------------------

# ra = 304.39133
# dec = -3.66411

ra = 245.4470133276900
dec = -22.8862182469700

# Random
# ra = 245.
# dec = -21.


# ------------------------------------------------------------
# User-supplied values
#
# Used if no SIMBAD match is found for the primary source.
# ------------------------------------------------------------

user_source_name = "MySourceName"
user_class_name = "MyClassName"

PORB = "myOrbit"
PORB_ERR = None

PSPIN = "mySpin"
PSPIN_ERR = None


# ------------------------------------------------------------
# Neighbor duplicate tolerance
#
# If a RACS neighbor is within this distance of a VLASS
# neighbor, the RACS neighbor is removed and the VLASS
# neighbor is retained.
# ------------------------------------------------------------

neighbor_match_radius_arcsec = 2.0


# ------------------------------------------------------------
# Output mode
#
# True:
#     Delete the current output database and summary CSV
#     before beginning.
#
# False:
#     Keep the current output database and summary CSV
#     and append the new source(s).
# ------------------------------------------------------------

start_fresh = True


# ------------------------------------------------------------
# Check input database
# ------------------------------------------------------------

if not MASTER_DB.exists():

    raise FileNotFoundError(
        f"Input database not found: {MASTER_DB}"
    )


# ------------------------------------------------------------
# Create output directories
# ------------------------------------------------------------

OUTPUT_DB.parent.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT_SUMMARY.parent.mkdir(
    parents=True,
    exist_ok=True
)


# ------------------------------------------------------------
# Prepare output files
# ------------------------------------------------------------

if start_fresh:

    print()
    print("Starting fresh.")
    print()

    if OUTPUT_DB.exists():

        OUTPUT_DB.unlink()

        print("Deleted existing database:")
        print()
        print(f"  {OUTPUT_DB}")
        print()

    if OUTPUT_SUMMARY.exists():

        OUTPUT_SUMMARY.unlink()

        print("Deleted existing summary CSV:")
        print()
        print(f"  {OUTPUT_SUMMARY}")
        print()

else:

    print()
    print("Append mode.")
    print()

    if OUTPUT_DB.exists():

        print("Existing output database will be retained:")
        print()
        print(f"  {OUTPUT_DB}")
        print()

    else:

        print("No existing output database found.")
        print()
        print("A new output database will be created:")
        print()
        print(f"  {OUTPUT_DB}")
        print()

    if OUTPUT_SUMMARY.exists():

        print("Existing summary CSV will be retained:")
        print()
        print(f"  {OUTPUT_SUMMARY}")
        print()

    else:

        print("No existing summary CSV found.")
        print()
        print("A new summary CSV will be created:")
        print()
        print(f"  {OUTPUT_SUMMARY}")
        print()


db_path = OUTPUT_DB


# ------------------------------------------------------------
# Determine next SOURCE_ID
#
# start_fresh = True:
#     Begin SOURCE_ID numbering at 1.
#
# start_fresh = False:
#     Begin at maximum SOURCE_ID currently in OUTPUT_DB + 1.
# ------------------------------------------------------------

if start_fresh:

    primary_source_id = 1

else:

    output_max_id = 0

    if (
        OUTPUT_DB.exists()
        and OUTPUT_DB.stat().st_size > 0
    ):

        conn = sqlite3.connect(
            OUTPUT_DB
        )

        try:

            table_exists = conn.execute(
                """
                SELECT name
                FROM sqlite_master
                WHERE type = 'table'
                  AND name = 'source_table'
                """
            ).fetchone()

            if table_exists is not None:

                result = conn.execute(
                    """
                    SELECT MAX(SOURCE_ID)
                    FROM source_table
                    """
                ).fetchone()[0]

                if result is not None:

                    output_max_id = int(
                        result
                    )

        finally:

            conn.close()

    primary_source_id = (
        output_max_id + 1
    )


print()
print(
    f"New primary SOURCE_ID: "
    f"{primary_source_id}"
)
print()


# ============================================================
# Source Table
# ============================================================


# ------------------------------------------------------------
# Primary source
# ------------------------------------------------------------

df_primary = pd.DataFrame({
    "PRIMARY_FLAG": [True],
    "RA": [ra],
    "DEC": [dec],
})


# ------------------------------------------------------------
# VLASS neighbors
# ------------------------------------------------------------

df_vlass_neighbors = vlass_neighbors(
    ra,
    dec
)

df_vlass_neighbors[
    "PRIMARY_FLAG"
] = False


print(
    f"VLASS neighbors found: "
    f"{len(df_vlass_neighbors)}"
)


# ------------------------------------------------------------
# RACS neighbors
# ------------------------------------------------------------

df_racs_neighbors = racs_neighbors(
    ra,
    dec
)

df_racs_neighbors[
    "PRIMARY_FLAG"
] = False


print(
    f"RACS neighbors found: "
    f"{len(df_racs_neighbors)}"
)


# ------------------------------------------------------------
# Remove RACS neighbors that duplicate VLASS neighbors
#
# A RACS neighbor is considered the same source as a VLASS
# neighbor when the positions differ by no more than
# neighbor_match_radius_arcsec.
#
# VLASS has priority:
#     VLASS neighbor -> retained
#     matching RACS neighbor -> deleted
#
# No other duplicate checking is performed.
# ------------------------------------------------------------

if (
    len(df_vlass_neighbors) > 0
    and len(df_racs_neighbors) > 0
):

    vlass_coords = SkyCoord(
        ra=(
            df_vlass_neighbors["RA"]
            .astype(float)
            .to_numpy()
            * u.deg
        ),
        dec=(
            df_vlass_neighbors["DEC"]
            .astype(float)
            .to_numpy()
            * u.deg
        ),
        frame="icrs"
    )

    racs_coords = SkyCoord(
        ra=(
            df_racs_neighbors["RA"]
            .astype(float)
            .to_numpy()
            * u.deg
        ),
        dec=(
            df_racs_neighbors["DEC"]
            .astype(float)
            .to_numpy()
            * u.deg
        ),
        frame="icrs"
    )

    # For every RACS neighbor, find the nearest VLASS neighbor.
    _, separation, _ = (
        racs_coords.match_to_catalog_sky(
            vlass_coords
        )
    )

    duplicate_racs = (
        separation.arcsec
        <= neighbor_match_radius_arcsec
    )

    number_removed = int(
        np.sum(duplicate_racs)
    )

    df_racs_neighbors = (
        df_racs_neighbors.loc[
            ~duplicate_racs
        ]
        .reset_index(drop=True)
    )

    print(
        f"RACS neighbors matching VLASS "
        f"within "
        f"{neighbor_match_radius_arcsec:.1f} arcsec: "
        f"{number_removed}"
    )

    print(
        f"RACS neighbors retained: "
        f"{len(df_racs_neighbors)}"
    )


# ------------------------------------------------------------
# Combine primary source and neighbors
# ------------------------------------------------------------

df_source = pd.concat(
    [
        df_primary,
        df_vlass_neighbors,
        df_racs_neighbors
    ],
    ignore_index=True,
)


# ------------------------------------------------------------
# Add SIMBAD fields
# ------------------------------------------------------------

df_source = add_simbad_fields(
    df_source
)


# ------------------------------------------------------------
# Add Gaia distance
# ------------------------------------------------------------

df_source = add_gaia_dist(
    df_source
)


# ------------------------------------------------------------
# Keep primary source even if it has no SIMBAD match.
#
# Remove neighbors that have no SOURCE_NAME.
# ------------------------------------------------------------

df_source = pd.concat(
    [
        df_source.iloc[[0]],

        df_source.iloc[1:].dropna(
            subset=["SOURCE_NAME"]
        )
    ],
    ignore_index=True
)


# ------------------------------------------------------------
# Determine whether primary source has SIMBAD match
# ------------------------------------------------------------

if pd.isna(
    df_source.at[
        0,
        "SOURCE_NAME"
    ]
):

    no_simbad = True

else:

    no_simbad = False


# ------------------------------------------------------------
# Use user-supplied name if primary has no SIMBAD match
# ------------------------------------------------------------

if no_simbad:

    df_source.at[
        0,
        "SOURCE_NAME"
    ] = user_source_name


# ------------------------------------------------------------
# Assign SOURCE_ID values
# ------------------------------------------------------------

df_source[
    "SOURCE_ID"
] = np.arange(
    primary_source_id,
    primary_source_id
    + len(df_source),
)


# ------------------------------------------------------------
# Save copy containing SOURCE_NAME for name-table generation
# ------------------------------------------------------------

df_source_for_name = (
    df_source.copy()
)


# ------------------------------------------------------------
# Remove temporary SOURCE_NAME column from source table
# ------------------------------------------------------------

df_source = df_source.drop(
    columns=[
        "SOURCE_NAME"
    ]
)


print(
    f"New source rows: "
    f"{len(df_source)}"
)


# ============================================================
# Name Table
# ============================================================

if no_simbad:

    df_name = pd.DataFrame({
        "SOURCE_ID": [
            primary_source_id
        ],
        "NAME": [
            user_source_name
        ],
        "DEFAULT_NAME": [
            True
        ]
    })

else:

    df_name = name(
        df_source_for_name
    )


print(
    f"New name rows: "
    f"{len(df_name)}"
)


# ============================================================
# Measurement Table
# ============================================================

df_vlass_1 = pd.DataFrame()
df_vlass_23 = pd.DataFrame()
df_racs = pd.DataFrame()


# ------------------------------------------------------------
# Radio measurements for primary source and neighbors
# ------------------------------------------------------------

for _, row in df_source.iterrows():

    source_id = row[
        "SOURCE_ID"
    ]

    ra_source = row[
        "RA"
    ]

    dec_source = row[
        "DEC"
    ]


    # --------------------------------------------------------
    # VLASS epoch 1
    # --------------------------------------------------------

    df_temp = vlass(
        1,
        source_id,
        ra_source,
        dec_source
    )

    df_vlass_1 = pd.concat(
        [
            df_vlass_1,
            df_temp
        ],
        ignore_index=True
    )


    # --------------------------------------------------------
    # VLASS epochs 2/3
    # --------------------------------------------------------

    df_temp = vlass(
        23,
        source_id,
        ra_source,
        dec_source
    )

    df_vlass_23 = pd.concat(
        [
            df_vlass_23,
            df_temp
        ],
        ignore_index=True
    )


    # --------------------------------------------------------
    # RACS
    # --------------------------------------------------------

    df_temp = racs(
        source_id,
        ra_source,
        dec_source
    )

    df_racs = pd.concat(
        [
            df_racs,
            df_temp
        ],
        ignore_index=True
    )


# ------------------------------------------------------------
# Combine VLASS measurements
# ------------------------------------------------------------

df_vlass = pd.concat(
    [
        df_vlass_1,
        df_vlass_23
    ],
    ignore_index=True
)


# ------------------------------------------------------------
# Identify primary source
# ------------------------------------------------------------

df_primary_new = df_source[
    df_source[
        "PRIMARY_FLAG"
    ] == True
].copy()


primary_row = (
    df_primary_new.iloc[0]
)


primary_id = (
    primary_row[
        "SOURCE_ID"
    ]
)


primary_ra = (
    primary_row[
        "RA"
    ]
)


primary_dec = (
    primary_row[
        "DEC"
    ]
)


# ------------------------------------------------------------
# Multi-wavelength measurements for primary source
# ------------------------------------------------------------

df_tmass = tmass(
    primary_id,
    primary_ra,
    primary_dec
)


df_gaia = gaia(
    primary_id,
    primary_ra,
    primary_dec
)


df_galex = galex(
    primary_id,
    primary_ra,
    primary_dec
)


df_chandra = chandra(
    primary_id,
    primary_ra,
    primary_dec
)


df_rosat = rosat(
    primary_id,
    primary_ra,
    primary_dec
)


df_xmm = xmm(
    primary_id,
    primary_ra,
    primary_dec
)


# ------------------------------------------------------------
# Combine measurements
# ------------------------------------------------------------

df_measurement = pd.concat(
    [
        df_vlass,
        df_racs,
        df_tmass,
        df_gaia,
        df_galex,
        df_chandra,
        df_rosat,
        df_xmm,
    ],
    ignore_index=True,
)


df_measurement = (
    df_measurement.rename(
        columns={
            "CAT": "EPOCH"
        }
    )
)


print(
    f"New measurement rows: "
    f"{len(df_measurement)}"
)


# ============================================================
# Class Table
# ============================================================

df_class = pd.DataFrame()


for _, row in df_source.iterrows():

    source_id = row[
        "SOURCE_ID"
    ]

    ra_source = row[
        "RA"
    ]

    dec_source = row[
        "DEC"
    ]


    df_temp = source_class(
        source_id,
        ra_source,
        dec_source
    )


    df_class = pd.concat(
        [
            df_class,
            df_temp
        ],
        ignore_index=True
    )


# ------------------------------------------------------------
# Use user-supplied class if primary has no SIMBAD match
# ------------------------------------------------------------

if no_simbad:

    df_class = (
        df_class[
            df_class[
                "SOURCE_ID"
            ]
            != primary_source_id
        ]
        .reset_index(
            drop=True
        )
    )


    df_primary_class = pd.DataFrame({
        "SOURCE_ID": [
            primary_source_id
        ],
        "CLASS": [
            user_class_name
        ],
        "DEFAULT_CLASS": [
            True
        ],
    })


    df_class = pd.concat(
        [
            df_primary_class,
            df_class
        ],
        ignore_index=True
    )


print(
    f"New class rows: "
    f"{len(df_class)}"
)


# ============================================================
# Period Table
# ============================================================

df_period = pd.DataFrame({
    "SOURCE_ID": [
        primary_id,
        primary_id
    ],

    "TYPE": [
        "Orbital",
        "Spin"
    ],

    "UNITS": [
        "Hours",
        "Seconds"
    ],

    "PERIOD": [
        PORB,
        PSPIN
    ],

    "PERIOD_ERROR": [
        PORB_ERR,
        PSPIN_ERR
    ],
})


print(
    f"New period rows: "
    f"{len(df_period)}"
)


# ============================================================
# Summary dataframe
# ============================================================


# ------------------------------------------------------------
# Default name
# ------------------------------------------------------------

df_default_name = (
    df_name.loc[
        df_name[
            "DEFAULT_NAME"
        ].isin(
            [1, True]
        ),
        [
            "SOURCE_ID",
            "NAME"
        ]
    ]
    .drop_duplicates(
        subset="SOURCE_ID",
        keep="first"
    )
)


# ------------------------------------------------------------
# Default class
# ------------------------------------------------------------

df_default_class = (
    df_class.loc[
        df_class[
            "DEFAULT_CLASS"
        ].isin(
            [1, True]
        ),
        [
            "SOURCE_ID",
            "CLASS"
        ]
    ]
    .drop_duplicates(
        subset="SOURCE_ID",
        keep="first"
    )
)


# ------------------------------------------------------------
# Orbital period
# ------------------------------------------------------------

df_porb = (
    df_period.loc[
        df_period[
            "TYPE"
        ]
        .astype(str)
        .str.casefold()
        .eq(
            "orbital"
        ),
        [
            "SOURCE_ID",
            "PERIOD"
        ]
    ]
    .drop_duplicates(
        subset="SOURCE_ID",
        keep="first"
    )
    .rename(
        columns={
            "PERIOD": "PORB"
        }
    )
)


# ------------------------------------------------------------
# Ensure measurement values are numeric
#
# Non-numeric values are converted to NaN.
# The original df_measurement is not modified.
# ------------------------------------------------------------

df_measurement_numeric = (
    df_measurement.copy()
)


df_measurement_numeric[
    "MEASUREMENT_NUMERIC"
] = pd.to_numeric(
    df_measurement_numeric[
        "MEASUREMENT"
    ],
    errors="coerce"
)


# ------------------------------------------------------------
# Gaia G magnitude
# ------------------------------------------------------------

gaia_mask = (
    df_measurement_numeric[
        "OBSERVATORY"
    ]
    .astype(str)
    .str.casefold()
    .eq(
        "gaia"
    )
    &
    df_measurement_numeric[
        "BAND"
    ]
    .astype(str)
    .str.casefold()
    .eq(
        "g"
    )
)


df_gaia_summary = (
    df_measurement_numeric.loc[
        gaia_mask,
        [
            "SOURCE_ID",
            "MEASUREMENT_NUMERIC"
        ]
    ]
    .drop_duplicates(
        subset="SOURCE_ID",
        keep="first"
    )
    .rename(
        columns={
            "MEASUREMENT_NUMERIC":
                "GAIA_G_MAG"
        }
    )
)


# ------------------------------------------------------------
# RACS flux
# ------------------------------------------------------------

racs_mask = (
    df_measurement_numeric[
        "OBSERVATORY"
    ]
    .astype(str)
    .str.casefold()
    .eq(
        "askap"
    )
)


df_racs_summary = (
    df_measurement_numeric.loc[
        racs_mask,
        [
            "SOURCE_ID",
            "MEASUREMENT_NUMERIC"
        ]
    ]
    .drop_duplicates(
        subset="SOURCE_ID",
        keep="first"
    )
    .rename(
        columns={
            "MEASUREMENT_NUMERIC":
                "RACS_FLUX"
        }
    )
)


# ------------------------------------------------------------
# VLASS maximum flux and number of numerical detections
# ------------------------------------------------------------

vlass_mask = (
    df_measurement_numeric[
        "OBSERVATORY"
    ]
    .astype(str)
    .str.casefold()
    .eq(
        "vla"
    )
)


df_vlass_summary = (
    df_measurement_numeric.loc[
        vlass_mask
    ]
    .groupby(
        "SOURCE_ID",
        as_index=False
    )
    .agg(
        VLASS_MAX_FLUX=(
            "MEASUREMENT_NUMERIC",
            "max"
        ),

        VLASS_NUM_DETECTIONS=(
            "MEASUREMENT_NUMERIC",
            "count"
        )
    )
)


# ------------------------------------------------------------
# Convert decimal degrees to sexagesimal
# ------------------------------------------------------------

def ra_deg_to_hms(ra_deg):
    """
    Convert RA from decimal degrees to HH:MM:SS.SS
    """
    if pd.isna(ra_deg):
        return None

    total_hours = float(ra_deg) / 15.0

    hours = int(total_hours)
    total_minutes = (total_hours - hours) * 60.0
    minutes = int(total_minutes)
    seconds = (total_minutes - minutes) * 60.0

    return f"{hours:02d}:{minutes:02d}:{seconds:05.2f}"


def dec_deg_to_dms(dec_deg):
    """
    Convert DEC from decimal degrees to DD:MM:SS.S
    """
    if pd.isna(dec_deg):
        return None

    sign = "+" if float(dec_deg) >= 0 else "-"
    dec_abs = abs(float(dec_deg))

    degrees = int(dec_abs)
    total_minutes = (dec_abs - degrees) * 60.0
    minutes = int(total_minutes)
    seconds = (total_minutes - minutes) * 60.0

    return f"{sign}{degrees:02d}:{minutes:02d}:{seconds:04.1f}"


# ------------------------------------------------------------
# Construct summary dataframe
#
# The SQLite database contains the primary source and neighbors.
# The summary CSV contains ONLY the primary source.
# ------------------------------------------------------------

df_summary = (
    df_source.loc[
        df_source["PRIMARY_FLAG"].isin([1, True]),
        [
            "SOURCE_ID",
            "PRIMARY_FLAG",
            "RA",
            "DEC"
        ]
    ]

    .merge(
        df_default_name,
        on="SOURCE_ID",
        how="left"
    )

    .merge(
        df_default_class,
        on="SOURCE_ID",
        how="left"
    )

    .merge(
        df_porb,
        on="SOURCE_ID",
        how="left"
    )

    .merge(
        df_gaia_summary,
        on="SOURCE_ID",
        how="left"
    )

    .merge(
        df_racs_summary,
        on="SOURCE_ID",
        how="left"
    )

    .merge(
        df_vlass_summary,
        on="SOURCE_ID",
        how="left"
    )
)


# ------------------------------------------------------------
# Sources with no VLA measurement rows should have zero
# detections rather than NaN.
# ------------------------------------------------------------

df_summary[
    "VLASS_NUM_DETECTIONS"
] = (
    df_summary[
        "VLASS_NUM_DETECTIONS"
    ]
    .fillna(0)
    .astype(int)
)

 
# ------------------------------------------------------------
# Explicit final column order
# ------------------------------------------------------------

df_summary = df_summary[
    [
        "SOURCE_ID",
        "NAME",
        "RA",
        "DEC",
        "CLASS",
        "PORB",
        "GAIA_G_MAG",
        "RACS_FLUX",
        "VLASS_MAX_FLUX",
        "VLASS_NUM_DETECTIONS",
    ]
]

df_summary["RA"] = df_summary["RA"].apply(ra_deg_to_hms)
df_summary["DEC"] = df_summary["DEC"].apply(dec_deg_to_dms)

df_summary["RA"] = df_summary["RA"].str.replace(":", " ", regex=False)
df_summary["DEC"] = df_summary["DEC"].str.replace(":", " ", regex=False)

# ------------------------------------------------------------
# Replace NaN with None
# ------------------------------------------------------------

df_summary = (
    df_summary.where(
        pd.notna(
            df_summary
        ),
        None
    )
)


# ============================================================
# Write / append summary CSV
# ============================================================

if (
    not start_fresh
    and OUTPUT_SUMMARY.exists()
    and OUTPUT_SUMMARY.stat().st_size > 0
):

    # Append rows only.
    # Do not write column names again.

    df_summary.to_csv(
        OUTPUT_SUMMARY,
        mode="a",
        header=False,
        index=False
    )


    print(
        f"Appended "
        f"{len(df_summary)} rows "
        f"to summary CSV."
    )

else:

    # New CSV.
    # Include column names.

    df_summary.to_csv(
        OUTPUT_SUMMARY,
        mode="w",
        header=True,
        index=False
    )


    print(
        f"Created summary CSV with "
        f"{len(df_summary)} rows."
    )


print()
print(
    df_summary
)
print()


# ============================================================
# Add new rows to output database
# ============================================================

tables = {
    "source_table":
        df_source,

    "name_table":
        df_name,

    "measurement_table":
        df_measurement,

    "class_table":
        df_class,

    "period_table":
        df_period,
}


conn = sqlite3.connect(
    db_path
)


try:

    for table_name, df in tables.items():

        df_sql = df.replace({
            np.nan: None
        })


        if len(df_sql) == 0:

            print(
                f"Skipping "
                f"{table_name}: "
                f"no rows"
            )

            continue


        df_sql.to_sql(
            table_name,
            conn,
            if_exists="append",
            index=False,
        )


        print(
            f"Appended "
            f"{len(df_sql)} rows "
            f"to {table_name}"
        )


    conn.commit()


finally:

    conn.close()


# ------------------------------------------------------------
# Final cleanup
# ------------------------------------------------------------

gc.collect()


print()
print(
    "Update complete."
)
print()


print(
    "Updated database:"
)
print(
    f"  {OUTPUT_DB}"
)
print()


print(
    "Summary CSV:"
)
print(
    f"  {OUTPUT_SUMMARY}"
)
print()


# sys.exit()
