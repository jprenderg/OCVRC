# Getting started

Install Python and the dependencies from the repository root:

```bash
python -m pip install -r requirements.txt
```

Inspect `config.txt`, especially `[PATHS]`, `[SIMBAD]`, `[RACS]` and `[PIPELINE]`. Replace the demonstration `input/sources.csv` with your own source coordinates. The default output location is `output/`.

Run the complete workflow from the project root:

```bash
python Run_OCVRC.py
```

On completion, the SQLite database is at `output/final/OCVRC.db`; five final pandas tables are also saved to `output/final/*.pkl`. Run `python validation/Run_Validation.py` to validate existing outputs without contacting astronomical services.

**External dependencies:** SIMBAD, NRAO VLASS quicklook services, CSIRO CASDA, survey catalog services and their authentication/access rules can change or be unavailable. The sample data are not a scientifically representative sample.
