# Downloads

This directory contains the primary data products for the **OCVRC**.

## Files

### `OcvcVlassRacs.sqlite`

The complete SQLite database containing the following tables:

- `Source_Table`
- `Name_Table`
- `Measurement_Table`
- `Class_Table`
- `Period_Table`


The tables are linked through the `SOURCE_ID` field.

This file is recommended for users who wish to perform SQL queries or make use of the complete relational database.

---

### `Summary_Table.csv`

A CSV version of the `Source Summary Table`.

This file contains one row for each primary catalog source and is intended for users who require only the primary source catalog without the complete relational database.

---

## Documentation

Complete documentation, including descriptions of all database fields, download instructions, and update procedures, is available through the project's GitHub Pages documentation site.

---

## Version

The files in this directory correspond to the current catalog release.

Please ensure that the downloaded files and the documentation are from the same release.

---

## Citation

If you use the **OCVRC** in your research, please cite the accompanying publication.
