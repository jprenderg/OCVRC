from pathlib import Path

# ------------------------------------------------------------
# Project root
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent

# ------------------------------------------------------------
# CASDA Login
# ------------------------------------------------------------

OPAL_USERNAME = "jprenderg@outlook.com"

# ------------------------------------------------------------
# Main directories
# ------------------------------------------------------------

DATA_DIR = PROJECT_ROOT / "data"
LOOKUP_DIR = DATA_DIR / "lookup"
CACHE_DIR = DATA_DIR / "cache" / "casda"

MASTER_DIR = PROJECT_ROOT / "master"
OUTPUT_DIR = PROJECT_ROOT / "output"
DOCS_DIR = PROJECT_ROOT / "docs"
SUBMISSIONS_DIR = PROJECT_ROOT / "submissions"

# ------------------------------------------------------------
# Lookup files
# ------------------------------------------------------------

# RA_DEC_FITS = LOOKUP_DIR / "ra_dec_fits.pkl"

VLASS1_LIST = LOOKUP_DIR / "VLASS1_img_list_10Aug25.txt"

VLASS23_LIST = LOOKUP_DIR / "VLASS_2_3_img_list_10Aug25.txt"

# ------------------------------------------------------------
# Databases and Summary
# ------------------------------------------------------------

MASTER_DB = MASTER_DIR / "OCVRC_Master.db"
MASTER_SUMMARY = MASTER_DIR / "Summary_Table_Master.csv"

OUTPUT_DB = OUTPUT_DIR / "OCVRC_Build.db"
OUTPUT_SUMMARY = OUTPUT_DIR / "Summary_Table_Build.csv"

# ------------------------------------------------------------
# Table directories
# ------------------------------------------------------------

SOURCE_DIR = PROJECT_ROOT / "source_table"
NAME_DIR = PROJECT_ROOT / "name_table"
MEASUREMENT_DIR = PROJECT_ROOT / "measurement_table"
CLASS_DIR = PROJECT_ROOT / "class_table"
PERIOD_DIR = PROJECT_ROOT / "period_table"

# ------------------------------------------------------------
# Python directories
# ------------------------------------------------------------

SOURCE_PY_DIR = SOURCE_DIR / "python"
NAME_PY_DIR = NAME_DIR / "python"
MEASUREMENT_PY_DIR = MEASUREMENT_DIR / "python"
CLASS_PY_DIR = CLASS_DIR / "python"
PERIOD_PY_DIR = PERIOD_DIR / "python"

# ------------------------------------------------------------
# Pickle outputs
# ------------------------------------------------------------

SOURCE_PKL = SOURCE_DIR / "pickle" / "Source.pkl"
NAME_PKL = NAME_DIR / "pickle" / "Name.pkl"
MEASUREMENT_PKL = MEASUREMENT_DIR / "pickle" / "Measurement.pkl"
CLASS_PKL = CLASS_DIR / "pickle" / "Class.pkl"
PERIOD_PKL = PERIOD_DIR / "pickle" / "Period.pkl"

# ------------------------------------------------------------
# CSV outputs
# ------------------------------------------------------------

SOURCE_CSV = SOURCE_DIR / "csv" / "Source.csv"
NAME_CSV = NAME_DIR / "csv" / "Name.csv"
MEASUREMENT_CSV = MEASUREMENT_DIR / "csv" / "Measurement.csv"
CLASS_CSV = CLASS_DIR / "csv" / "Class.csv"
PERIOD_CSV = PERIOD_DIR / "csv" / "Period.csv"

# ------------------------------------------------------------
# Create writable directories if missing
# ------------------------------------------------------------

for directory in [
    CACHE_DIR,
    MASTER_DIR,
    OUTPUT_DIR,
    SOURCE_DIR / "pickle",
    SOURCE_DIR / "csv",
    NAME_DIR / "pickle",
    NAME_DIR / "csv",
    MEASUREMENT_DIR / "pickle",
    MEASUREMENT_DIR / "csv",
    CLASS_DIR / "pickle",
    CLASS_DIR / "csv",
    PERIOD_DIR / "pickle",
    PERIOD_DIR / "csv",
]:
    directory.mkdir(parents=True, exist_ok=True)