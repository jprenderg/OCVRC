# Configuration

All editable settings are in the root `config.txt` file, with uppercase keys. Relative paths resolve relative to the project root. Key sections:

| Section | Purpose |
|---|---|
| `[PATHS]` | Input file, output directory, VLASS image URL lists, RACS cache directory |
| `[SIMBAD]` | Matching radius, timeout and retry controls |
| `[PROCESSING]` | `SKIP_SIMBAD` (local testing only) |
| `[IMAGE_RETRIEVAL]` | Source row range and optional plotting |
| `[VLASS]`, `[RACS]` | Radio cutout parameters and RACS cache access |
| Survey sections | Coordinate matching radius, query row ranges and retry settings |
| `[QUERY_MANAGEMENT]` | `RESUME`, `SAVE_EVERY`, `CONSECUTIVE_FAILURE_LIMIT` |
| `[PIPELINE]` | Stage switches, including `VALIDATION = TRUE` |

Do not commit passwords or authentication tokens. Set `CASDA_USERNAME` in your local `config.txt` or use the `CASDA_USERNAME` environment variable as supported by the retrieval script. Before a full run, ensure `START_ROW = 0` and `END_ROW = -1` for stages where the options apply.

The repository ships both VLASS URL lists and precomputed VLASS/RACS GLS correlation files. These files must remain in their configured locations.

### Diagnostic plots

Production suppresses optional image and neighbor-detection diagnostic figures by default. To enable them for investigation, set `[DIAGNOSTICS] SHOW_PLOTS = TRUE` and select an image using `[IMAGE_RETRIEVAL] DISPLAY_INDEX` in `config.txt`. Plot settings do not change the underlying measurements.
