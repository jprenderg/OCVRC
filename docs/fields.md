# Catalog Fields

This is the authoritative **final SQLite schema** implemented by `database/Build_SQLite_Database.py` and `database/Build_Final_Database.py` in this release. Table names and field names are exact.

The five tables are linked by `SOURCE_ID`. SQLite `INTEGER` flags encode True as `1`, False as `0`, and unknown/not-applicable values as `NULL` where supported.

## `source_table`

| Field | SQLite type | Units | SQL NULL allowed | Description |
|---|---|---|---|---|
| `SOURCE_ID` | INTEGER | — | No | Primary key, identifies the source |
| `RA` | REAL | degrees | No | Adopted right ascension (0 ≤ RA < 360) |
| `DEC` | REAL | degrees | No | Adopted declination (−90 ≤ DEC ≤ 90) |
| `SIMBAD` | INTEGER | 0/1 | No | Whether a SIMBAD positional match was found |
| `MAJOR_AXIS` | REAL | mas | Yes | SIMBAD coordinate-uncertainty ellipse axis value, where available |
| `MINOR_AXIS` | REAL | mas | Yes | SIMBAD coordinate-uncertainty ellipse axis value, where available |
| `POSITION_ANGLE` | REAL | degrees | Yes | SIMBAD coordinate-uncertainty position angle |
| `DISTANCE` | REAL | kpc | Yes | Distance estimated using Gaia parallax |
| `DISTANCE_ERROR` | REAL | kpc | Yes | Uncertainty in Gaia-derived distance |
| `PM_RA` | REAL | mas/yr | Yes | Gaia proper motion in RA, including cos(DEC) convention |
| `PM_DEC` | REAL | mas/yr | Yes | Gaia proper motion in declination |

## `name_table`

| Field | SQLite type | Units | SQL NULL allowed | Description |
|---|---|---|---|---|
| `SOURCE_ID` | INTEGER | — | No | Foreign key to source_table |
| `NAME` | TEXT | — | No | Source name or alias |
| `DEFAULT_NAME` | INTEGER | 0/1 | No | Exactly one preferred name (1) per source |

## `class_table`

| Field | SQLite type | Units | SQL NULL allowed | Description |
|---|---|---|---|---|
| `SOURCE_ID` | INTEGER | — | No | Foreign key to source_table |
| `CLASS` | TEXT | — | No | Source classification |
| `DEFAULT_CLASS` | INTEGER | 0/1 | No | Exactly one preferred class for each classified source |

## `period_table`

| Field | SQLite type | Units | SQL NULL allowed | Description |
|---|---|---|---|---|
| `SOURCE_ID` | INTEGER | — | No | Foreign key to source_table |
| `TYPE` | TEXT | — | No | Period type (Orbital or Spin) |
| `UNITS` | TEXT | hours or seconds | No | Hours for orbital, seconds for spin |
| `PERIOD` | REAL | as UNITS | No | User-provided period |
| `PERIOD_ERROR` | REAL | as UNITS | Yes | User-provided period uncertainty |
| `CONFIRMED` | INTEGER | 0/1/NULL | Yes | User-provided confirmation flag where available |

## `measurement_table`

| Field | SQLite type | Units | SQL NULL allowed | Description |
|---|---|---|---|---|
| `SOURCE_ID` | INTEGER | — | No | Foreign key to source_table |
| `OBSERVATORY` | TEXT | — | Yes¹ | Survey / observatory identifier |
| `EPOCH` | TEXT | — | Yes¹ | Observation epoch or catalog release |
| `BAND` | TEXT | varies | Yes¹ | Spectral, photometric or energy band |
| `UNITS` | TEXT | varies | Yes¹ | Units of measured value and associated uncertainty |
| `MEASUREMENT` | REAL | as UNITS | Yes¹ | Measured flux or catalog measurement; negative radio fits retained |
| `ERROR` | REAL | as UNITS | Yes¹ | Symmetric uncertainty when supplied |
| `POS_ERROR` | REAL | as UNITS | Yes¹ | Positive asymmetric uncertainty when supplied |
| `NEG_ERROR` | REAL | as UNITS | Yes¹ | Negative asymmetric uncertainty when supplied |
| `NEIGHBOR` | INTEGER | 0/1/NULL | Yes | Radio neighbor-contamination flag: 0 not flagged; 1 flagged; NULL not applicable |

¹ The SQLite schema permits NULL in measurement columns, but the current validation suite imposes stricter completeness checks on processed measurement rows. An observation can use a symmetric `ERROR` **or** a pair of asymmetric `POS_ERROR` and `NEG_ERROR` values.

## Measurement interpretation

The `OBSERVATORY`, `EPOCH`, `BAND` and `UNITS` fields must be interpreted together; unlike fixed radio bands, ancillary surveys have survey-specific measurement units and bands. VLASS (`VLA`) radio measurements use `2-4 GHz` and `mJy`; RACS (`ASKAP`) radio measurements use `887.5 MHz` and `mJy`. `NEIGHBOR` applies to those two radio surveys, and is NULL for all other surveys.

Repeated VLASS records for the same source and epoch may arise when image tiles overlap. They are retained as distinct fits rather than silently deduplicated.

## Constraints

* `SOURCE_ID` is unique and non-NULL in `source_table`.
* Each source has exactly one preferred `NAME`.
* Each classified source has exactly one preferred `CLASS`; unclassified sources may have no class row.
* The period table allows one row per `(SOURCE_ID, TYPE)`.
* The SQLite database includes foreign-key constraints and indexes for names, classes and measurements.
