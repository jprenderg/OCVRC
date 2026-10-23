# Input catalog

The default input is `input/sources.csv`, specified by `[PATHS] INPUT_CSV` in `config.txt`. The ingestion code also supports Excel files when `[PATHS] INPUT_FILE` points to one.

| Column | Required | Meaning |
|---|---|---|
| `RA` | Yes | Right ascension in decimal degrees, `[0,360)` |
| `DEC` | Yes | Declination in decimal degrees, `[-90,90]` |
| `SOURCE_ID` | No | Unique input source identifier; generated sequentially if absent |
| `NAME` | No | User-preferred source name |
| `CLASS` | No | User-preferred source classification |
| `PORB` | No | Orbital period in hours |
| `PORB_ERROR` | No | Orbital-period error, in hours |
| `PORB_CONFIRMED` | No | Optional confirmed flag |
| `PSPIN` | No | Spin period in seconds |
| `PSPIN_ERROR` | No | Spin-period error, in seconds |
| `PSPIN_CONFIRMED` | No | Optional confirmed flag |

A user-provided name or class takes precedence over SIMBAD as the default; SIMBAD alternatives are retained as additional name/class rows. Sources without a name get an automatically generated `User_name_{SOURCE_ID}` when necessary. Periods come from user input, not SIMBAD.

An unavailable SIMBAD service is **not** treated as a negative SIMBAD match. Successfully completed results may be reused according to the ingestion cache and verified-output checks.
