# Open Cataclysmisic Variable Radio Catalog (OCVRC)

The **Open Cataclysmic Variable Radio Vatalog (OCVRC)** is a database of radio observations of known and candidate cataclysmic variables (CVs). It combines measurements from the **Very Large Array Sky Survey (VLASS)** and the **Rapid ASKAP Continuum Survey (RACS)** together with supplementary measurements from major optical, infrared, ultraviolet, X-ray, and astrometric surveys.

This repository contains:

- The OCVRC master database.
- A CSV version of a Source Summary Table.
- Python software for building a user-defined database.
- Documentation describing the catalog and its database schema.

---

## Documentation

Complete documentation is provided in the docs/ directory and is published through the GitHub Pages site:

**https://jprenderg.github.io/OCVRC/**

The documentation includes:

- Catalog overview
- Database field definitions
- Download instructions
- Instructions for creating new source entries

---

## Downloads

The latest catalog data products include:

- SQLite OCVRC database
- Source Summary Table (CSV)

See the **Downloads** page in the documentation for the latest files.

---

## Updating the Catalog

The repository includes a Python utility for building a new database with user-defined sources.

The master database itself is maintained separately and cannot be updated directly by users.

Users should contact the authors if they wish to contribute to the database or have questions about the database.

Complete instructions are provided in the **Updating the Catalog** section of the documentation.

---

## Repository Structure

```text
OcvcVlassRacs/
│
├── docs/                  For documentation website
│
├── downloads/             Full database file and summary table
│
├── master/                Master catalog database and source summary table
│
├── output/                User created database and corresponding CSV source summary table.
│
├── data/
│   ├── lookup/            VLASS lookup tables
│   └── cache/             CASDA download cache
│
├── source_table/          Source table generation
├── measurement_table/     Measurement table generation
├── name_table/            Name table generation
├── class_table/           Class table generation
├── period_table/          Period table generation
│
├── update_database.py     Main update program
├── config.py              User configuration
├── requirements.txt       Python package requirements
├── mkdocs.yml             MkDocs configuration
└── README.md
```

---

## Citation

If you use the OCVRC Catalog in your research, please cite the accompanying publication.