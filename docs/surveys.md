# Additional surveys

The pipeline queries Gaia DR3, 2MASS, GALEX, LoTSS DR2, ROSAT, XMM-Newton, Chandra CSC 2.1 and eROSITA. Each integration populates the shared `measurement_table` format using survey-specific `OBSERVATORY`, `EPOCH`, `BAND` and `UNITS`. Gaia also supplies adopted distances and proper motions for `source_table`.

Survey matching radii and retry settings are defined independently in `config.txt`. The query manager stores status and checkpoint files under `output/`, permitting the pipeline to distinguish true `NO_MATCH` results from interrupted or unavailable-service queries.

A `NO_MATCH` does not prove an intrinsic non-detection or absence from the sky: footprints, catalog selection criteria, sensitivity and release coverage matter.
