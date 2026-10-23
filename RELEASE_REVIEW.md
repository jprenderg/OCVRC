# Release review — 2026-10-23

This is a **release candidate for review**, not a published repository. The provided project archive contained 12 root README variants with conflicting historical instructions, generated output files (including SQLite DB), local caches and compiled Python files. This candidate consolidates current-use documentation into `README.md` and `docs/`; it preserves the original scientific Python code and both GLS calibration files and both VLASS URL inventories.

## Confirmed code-authoritative schema

Source: 11 fields; name: 3; class: 3; period: 6; measurement: 10. See `docs/fields.md` for exact names and SQL NULL permissions drawn directly from `database/Build_SQLite_Database.py`.

## Explicit release considerations

- The supplied `config.txt` contained a personal CASDA email address. It has been blanked in the candidate. Configure your own credentials locally; **never publish passwords**.
- No license was included. Pick a license before public release; this candidate intentionally does not invent one.
- The old public website sources were not included in the archive. The new `docs/fields.md` is a code-grounded replacement; any unique historical commentary or scientific references should be migrated separately.
- `requirements.txt` is an inferred dependency list from imports, not a tested pinned environment. Check CASDA installation on a clean system.
- Outputs, image caches and SQLite test databases are excluded from this package; `output/` is created by running scripts.
- The existing first-stage `README.md` is outdated (states the complete pipeline is not implemented), and other historical root READMEs say validation is not integrated; neither statement is correct for the tested pipeline. They were deliberately not carried forward.
- Scientific Python scripts have not been changed or executed in this candidate.
- Fresh SIMBAD empty-cache behavior, including outage recovery and correct source provenance, requires a dedicated test.
