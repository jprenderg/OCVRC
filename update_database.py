from pathlib import Path
import sqlite3
import numpy as np
import pandas as pd
import gc
import sys
import shutil

from config import MASTER_DB, OUTPUT_DB, OUTPUT_SUMMARY

from source_table.python.MAXID import maxid
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

# AR Sco
ra = 304.39133
dec = -3.66411

# Random
# ra = 245.
# dec = -21.

# Will be used if no SIMBAD match for the given RA and DEC
user_source_name = "MySourceName"
user_class_name = "MyClassName"

PORB = None
PORB_ERR = None

PSPIN = None
PSPIN_ERR = None

# ------------------------------------------------------------
# Copy over files
# ------------------------------------------------------------

# ------------------------------------------------------------
# Prepare output database
# ------------------------------------------------------------

if not MASTER_DB.exists():
    raise FileNotFoundError(f"Input database not found: {MASTER_DB}")

OUTPUT_DB.parent.mkdir(parents=True, exist_ok=True)

if OUTPUT_DB.exists():
    OUTPUT_DB.unlink()

db_path = OUTPUT_DB

print("Created new updates-only database:")
print(f"  {OUTPUT_DB}")


# ------------------------------------------------------------
# Source Table
# ------------------------------------------------------------

primary_source_id = maxid() + 1

df_primary = pd.DataFrame({
    "PRIMARY_FLAG": [True],
    "RA": [ra],
    "DEC": [dec],
})

df_vlass_neighbors = vlass_neighbors(ra, dec)
df_vlass_neighbors["PRIMARY_FLAG"] = False

df_racs_neighbors = racs_neighbors(ra, dec)
df_racs_neighbors["PRIMARY_FLAG"] = False

df_source = pd.concat(
    [df_primary, df_vlass_neighbors, df_racs_neighbors],
    ignore_index=True,
)

df_source = add_simbad_fields(df_source)

df_source = add_gaia_dist(df_source)

df_source = pd.concat([
    df_source.iloc[[0]],
    df_source.iloc[1:].dropna(subset=["SOURCE_NAME"])
]).reset_index(drop=True)

if pd.isna(df_source.at[0, "SOURCE_NAME"]):
    no_simbad = True
else:
    no_simbad = False
    
if no_simbad == True:
    df_source.at[0, "SOURCE_NAME"] = user_source_name    

mask_duplicate_neighbor = (
    (df_source["PRIMARY_FLAG"] == False)
    & df_source.duplicated(subset="SOURCE_NAME", keep="first")
)

df_source = df_source.loc[~mask_duplicate_neighbor].reset_index(drop=True)


# # ------------------------------------------------------------
# # Remove neighbors already present in input database
# # ------------------------------------------------------------

conn = sqlite3.connect(MASTER_DB)

try:
    df_current_names = pd.read_sql_query(
        """
        SELECT s.PRIMARY_FLAG, n.NAME
        FROM source_table AS s
        JOIN name_table AS n
            ON s.SOURCE_ID = n.SOURCE_ID
        WHERE s.PRIMARY_FLAG = 0
          AND n.DEFAULT_NAME = 1
        """,
        conn,
    )

finally:
    conn.close()

current_names = df_current_names["NAME"].dropna()

mask_existing_neighbor = (
    (df_source["PRIMARY_FLAG"] == False)
    & df_source["SOURCE_NAME"].isin(current_names)
)

df_source = df_source.loc[~mask_existing_neighbor].reset_index(drop=True)

df_source["SOURCE_ID"] = np.arange(
    primary_source_id,
    primary_source_id + len(df_source),
)

df_source_for_name = df_source.copy()

df_source = df_source.drop(columns=["SOURCE_NAME"])

print(f"New source rows: {len(df_source)}")


# # ------------------------------------------------------------
# # Name Table
# # ------------------------------------------------------------

if no_simbad==True: 
    
    df_name = pd.DataFrame({
        "SOURCE_ID": [primary_source_id],
        "NAME": [user_source_name],
        "DEFAULT_NAME": [True]
    })

else: 
    df_name = name(df_source_for_name)

print(f"New name rows: {len(df_name)}")


# ------------------------------------------------------------
# Measurement Table
# ------------------------------------------------------------

df_vlass_1 = pd.DataFrame()
df_vlass_23 = pd.DataFrame()
df_racs = pd.DataFrame()

for _, row in df_source.iterrows():

    source_id = row["SOURCE_ID"]
    ra_source = row["RA"]
    dec_source = row["DEC"]

    df_temp = vlass(1, source_id, ra_source, dec_source)
    df_vlass_1 = pd.concat([df_vlass_1, df_temp], ignore_index=True)

    df_temp = vlass(23, source_id, ra_source, dec_source)
    df_vlass_23 = pd.concat([df_vlass_23, df_temp], ignore_index=True)

    df_temp = racs(source_id, ra_source, dec_source)
    df_racs = pd.concat([df_racs, df_temp], ignore_index=True)

df_vlass = pd.concat([df_vlass_1, df_vlass_23], ignore_index=True)

df_primary_new = df_source[df_source["PRIMARY_FLAG"] == True].copy()


primary_row = df_primary_new.iloc[0]

primary_id = primary_row["SOURCE_ID"]
primary_ra = primary_row["RA"]
primary_dec = primary_row["DEC"]

df_tmass = tmass(primary_id, primary_ra, primary_dec)
df_gaia = gaia(primary_id, primary_ra, primary_dec)
df_galex = galex(primary_id, primary_ra, primary_dec)
df_chandra = chandra(primary_id, primary_ra, primary_dec)
df_rosat = rosat(primary_id, primary_ra, primary_dec)
df_xmm = xmm(primary_id, primary_ra, primary_dec)


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

df_measurement = df_measurement.rename(columns={"CAT": "EPOCH"})

print(f"New measurement rows: {len(df_measurement)}")


# ------------------------------------------------------------
# Class Table
# ------------------------------------------------------------
df_class = pd.DataFrame()

for _, row in df_source.iterrows():

    source_id = row["SOURCE_ID"]
    ra_source = row["RA"]
    dec_source = row["DEC"]

    df_temp = source_class(source_id, ra_source, dec_source)
    df_class = pd.concat([df_class, df_temp], ignore_index=True)
    
    
if no_simbad:

    df_class = df_class[
        df_class["SOURCE_ID"] != primary_source_id
    ].reset_index(drop=True)

    df_primary_class = pd.DataFrame({
        "SOURCE_ID": [primary_source_id],
        "CLASS": [user_class_name],
        "DEFAULT_CLASS": [True],
    })

    df_class = pd.concat(
        [df_primary_class, df_class],
        ignore_index=True
    )    


print(f"New class rows: {len(df_class)}")



# ------------------------------------------------------------
# Period Table
# ------------------------------------------------------------


df_period = pd.DataFrame({
    "SOURCE_ID": [primary_id, primary_id],
    "TYPE": ["Orbital", "Spin"],
    "UNITS": ["Hours", "Seconds"],
    "PERIOD": [PORB, PSPIN],
    "PERIOD_ERROR": [PORB_ERR, PSPIN_ERR],
})


print(f"New period rows: {len(df_period)}")

##############################  df_summary ########################################


# ------------------------------------------------------------
# Default name
# ------------------------------------------------------------

df_default_name = (
    df_name.loc[
        df_name["DEFAULT_NAME"].isin([1, True]),
        ["SOURCE_ID", "NAME"]
    ]
    .drop_duplicates(subset="SOURCE_ID", keep="first")
)

# ------------------------------------------------------------
# Default class
# ------------------------------------------------------------

df_default_class = (
    df_class.loc[
        df_class["DEFAULT_CLASS"].isin([1, True]),
        ["SOURCE_ID", "CLASS"]
    ]
    .drop_duplicates(subset="SOURCE_ID", keep="first")
)

# ------------------------------------------------------------
# Orbital period
# ------------------------------------------------------------

df_porb = (
    df_period.loc[
        df_period["TYPE"].astype(str).str.casefold() == "orbit",
        ["SOURCE_ID", "PERIOD"]
    ]
    .drop_duplicates(subset="SOURCE_ID", keep="first")
    .rename(columns={"PERIOD": "PORB"})
)

# ------------------------------------------------------------
# Ensure measurement values are numeric
#
# Non-numeric values are converted to NaN. The original
# df_measurement is not modified.
# ------------------------------------------------------------

df_measurement_numeric = df_measurement.copy()

df_measurement_numeric["MEASUREMENT_NUMERIC"] = pd.to_numeric(
    df_measurement_numeric["MEASUREMENT"],
    errors="coerce"
)

# ------------------------------------------------------------
# Gaia G magnitude
# ------------------------------------------------------------

gaia_mask = (
    df_measurement_numeric["OBSERVATORY"].astype(str).str.casefold().eq("gaia")
    & df_measurement_numeric["BAND"].astype(str).str.casefold().eq("g")
)

df_gaia = (
    df_measurement_numeric.loc[
        gaia_mask,
        ["SOURCE_ID", "MEASUREMENT_NUMERIC"]
    ]
    .drop_duplicates(subset="SOURCE_ID", keep="first")
    .rename(columns={"MEASUREMENT_NUMERIC": "GAIA_G_MAG"})
)

# ------------------------------------------------------------
# RACS flux
# ------------------------------------------------------------

racs_mask = (
    df_measurement_numeric["OBSERVATORY"]
    .astype(str)
    .str.casefold()
    .eq("askap")
)

df_racs = (
    df_measurement_numeric.loc[
        racs_mask,
        ["SOURCE_ID", "MEASUREMENT_NUMERIC"]
    ]
    .drop_duplicates(subset="SOURCE_ID", keep="first")
    .rename(columns={"MEASUREMENT_NUMERIC": "RACS_FLUX"})
)

# ------------------------------------------------------------
# VLASS maximum flux and number of numerical detections
# ------------------------------------------------------------

vlass_mask = (
    df_measurement_numeric["OBSERVATORY"]
    .astype(str)
    .str.casefold()
    .eq("vla")
)

df_vlass = (
    df_measurement_numeric.loc[vlass_mask]
    .groupby("SOURCE_ID", as_index=False)
    .agg(
        VLASS_MAX_FLUX=("MEASUREMENT_NUMERIC", "max"),
        VLASS_NUM_DETECTIONS=("MEASUREMENT_NUMERIC", "count")
    )
)

# ------------------------------------------------------------
# Construct summary dataframe
# ------------------------------------------------------------

df_summary = (
    df_source[["SOURCE_ID", "RA", "DEC"]]
    .merge(df_default_name, on="SOURCE_ID", how="left")
    .merge(df_default_class, on="SOURCE_ID", how="left")
    .merge(df_porb, on="SOURCE_ID", how="left")
    .merge(df_gaia, on="SOURCE_ID", how="left")
    .merge(df_racs, on="SOURCE_ID", how="left")
    .merge(df_vlass, on="SOURCE_ID", how="left")
)

# Sources with no VLA measurement rows should have zero detections,
# rather than NaN.
df_summary["VLASS_NUM_DETECTIONS"] = (
    df_summary["VLASS_NUM_DETECTIONS"]
    .fillna(0)
    .astype(int)
)

# Explicit final column order
df_summary = df_summary[
    [
        "SOURCE_ID",
        "RA",
        "DEC",
        "NAME",
        "CLASS",
        "PORB",
        "GAIA_G_MAG",
        "RACS_FLUX",
        "VLASS_MAX_FLUX",
        "VLASS_NUM_DETECTIONS",
    ]
]

df_summary = df_summary.where(pd.notna(df_summary), None)

df_summary.to_csv(
    OUTPUT_SUMMARY,
    index=False
)

print(df_summary)



######################################################################################

# ------------------------------------------------------------
# Add only new rows to output database
# ------------------------------------------------------------

tables = {
    "source_table": df_source,
    "name_table": df_name,
    "measurement_table": df_measurement,
    "class_table": df_class,
    "period_table": df_period,
}

conn = sqlite3.connect(db_path)

try:
    for table_name, df in tables.items():

        df = df.replace({np.nan: None})

        if len(df) == 0:
            print(f"Skipping {table_name}: no rows")
            continue

        df.to_sql(
            table_name,
            conn,
            if_exists="append",
            index=False,
        )

        print(f"Appended {len(df)} rows to {table_name}")

    conn.commit()

finally:
    conn.close()


# ------------------------------------------------------------
# Final cleanup
# ------------------------------------------------------------

gc.collect()

print("Update complete.")
print(f"Updated database written to: {OUTPUT_DB}")

sys.exit()