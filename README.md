# OCVRC

Source-list-to-SQLite astronomical catalog construction pipeline, including VLASS/RACS radio PSF fitting, SIMBAD identification, ancillary survey crossmatching and internal validation.

**[Documentation](https://jprenderg.github.io/OCVRC/)** · **[Catalog Fields](https://jprenderg.github.io/OCVRC/fields/)**

## Quick start

```bash
python -m pip install -r requirements.txt
python Run_OCVRC.py
```

Update `config.txt` and supply your own `input/sources.csv` first. The SQLite output is `output/final/OCVRC.db`. To validate existing outputs without querying services: `python validation/Run_Validation.py`.

Full installation, fields, validation and scientific caveats are in the bundled MkDocs sources in `docs/`. Build locally with `python -m pip install -r requirements-docs.txt` and `mkdocs serve`.

**Note:** This candidate has not been smoke-tested from a clean installation. Source code is copied unmodified from the supplied working project, aside from making the example CASDA username empty in `config.txt`.
