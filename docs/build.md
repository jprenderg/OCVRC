# Building a User Defined Catalog

This page describes how to use the Open Cataclysmic Variable Radio Catalog (OCVRC) Build Utility to build a new user catalog by adding sources and their associated measurements.

The build utility **never modifies the master database**. Instead, all new records are written to an output database. Users should contact the authors if they wish to contribute to the master database or have questions about the database.

---

## Requirements

Before running the update utility, ensure that you have:

- Python 3.11 or later
- An OPAL account with access to CASDA
- An internet connection
- The Python packages listed in `requirements.txt`

Install the required packages:

```bash
pip install -r requirements.txt
```

## OPAL / CASDA Login

Some measurements require downloading data from the CASDA archive.

Before running the program, edit `config.py` by entering your OPAL username:

```python
OPAL_USERNAME = "your_email@example.com"
```

When the program runs, you will be prompted to enter your OPAL password in the terminal or console.

Your password is **never stored** by this project.

---

## Running the Program

Open `build_database.py`.

Near the beginning of the file is a section labeled **User Inputs**.

Enter the following information:

- Right Ascension (RA)
- Declination (Dec)
- Orbital period (optional)
- Orbital period uncertainty (optional)
- Spin period (optional)
- Spin period uncertainty (optional)
- A user provided source and class name for use in the event that a SIMBAD match is not available

Run the program:

```bash
python build_database.py
```

The program automatically:

- Assigns a new `SOURCE_ID`
- Creates records in `Source_Table`
- Creates records in `Name_Table`
- Creates records in `Measurement_Table`
- Creates records in `Class_Table`
- Creates records in `Period_Table`

---

## Output

The output database contains only user created records and is written to:

```text
output/OCVRC_Build.db
```
The master database is **never modified**.

A Source Summary Table is located at:

```text
output/Summary_Table_Build.csv
```
This file contains data regarding primary sources only (i.e., radio-loud neighbors are excluded)

The user can control whether the current version of `OCVRC_Build.db` and `summary_table_build.csv` are deleted before adding a new source or new source information is appended to those files. This is accomplished by setting `strt_fresh = True or False`

---

## CASDA Cache

Downloaded CASDA files are stored in:

```text
data/cache/
```

Previously downloaded files are reused whenever possible, reducing download time on future runs.

The cache may be deleted at any time. Files will be downloaded again automatically when needed.

---


