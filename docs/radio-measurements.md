# Radio measurements

VLASS (2–4 GHz) and ASKAP RACS (887.5 MHz) survey images are sampled at catalog positions. The pipeline applies neighbor-detection logic and generalized least-squares (GLS) PSF fitting using precomputed spatial-correlation calibrations in `calibration/`.

The final `measurement_table` retains measurements even when fitted flux is negative or the neighbor flag is `1`. The `NEIGHBOR` flag indicates a possible nearby-source contribution; it is **not** a detection significance indicator. Interpret flux using its uncertainty and image context.

Repeated source/epoch VLASS entries can reflect overlapping tiles. Use the actual measurement rows, not a presumed single row per epoch. See [Catalog Fields](fields.md).
