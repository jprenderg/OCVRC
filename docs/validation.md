# Validation

The master pipeline can run the combined suite automatically with `[PIPELINE] VALIDATION = TRUE` or you can execute:

```bash
python validation/Run_Validation.py
```

Four independent checks test (1) SQLite/database structure and integrity, (2) astrometry and ingestion consistency, (3) measurement schemas and uncertainties and (4) query statuses and cached survey output consistency. The validator is read-only and performs no network queries.

The October 23, 2026 six-source integration test reported **292 PASS, 2 WARNING, 0 FAIL**. Both warnings describe the same four groups of overlapping VLASS source/epoch measurements. These are example-test results, **not** a guarantee that arbitrary input catalogs will pass.
