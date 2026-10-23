# Building a catalog

Run `python Run_OCVRC.py` from the repository root. The master runner calls: source ingestion; VLASS/RACS image retrieval; neighbor detection; GLS radio PSF fitting; radio table assembly; Gaia, 2MASS, GALEX, LoTSS, ROSAT, XMM, Chandra and eROSITA integrations; final five-table assembly; SQLite creation; and optional validation.

`[PIPELINE] VALIDATION = TRUE` runs the same validation suite after database creation. Source queries reuse completed checkpoints with `[QUERY_MANAGEMENT] RESUME = TRUE`. Failed service requests and genuine no-matches are tracked separately.

Outputs include `output/final/OCVRC.db` and five final pandas DataFrames. The SQLite builder can also run independently as `python database/Build_SQLite_Database.py` using existing final DataFrames.

The SQLite database includes **all** retained radio fit results, including negative fluxes and measurements flagged as affected by neighbors. Where VLASS tiles overlap, more than one measurement per source/epoch may appear.
