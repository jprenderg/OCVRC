# Updating the Catalog

This page describes how to use the Open Cataclysmic Variable Radio Catalog (OCVRC) Update Utility to create a new source and its associated measurements to the catalog.

The update utility **never modifies the master database**. Instead, all new records are written to an output database. Users should contact the authors if they wish to contribute to the master database or have questions about the database.

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

---

## OPAL / CASDA Login

Some measurements require downloading data from the CASDA archive.

Before running the program, edit `config.py` and enter your OPAL username:

```python
OPAL_USERNAME = "your_email@example.com"
```

When the program runs, you will be prompted to enter your OPAL password in the terminal or console.

Your password is **never stored** by this project.

---

## Running the Program

Open `update_database.py`.

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
python update_database.py
```

The program automatically:

- Assigns a new `SOURCE_ID`
- Creates records in `Source_Table`
- Creates records in `Name_Table`
- Creates records in `Measurement_Table`
- Creates records in `Class_Table`
- Creates records in `Period_Table`
- Creates records in `Summary_Table`

The output database contains only the newly created record and is written to:

```text
output/VLASS_RACS_Updated.db
```

The master database is **never modified**.

---

## Output

After the program finishes, the ouptut database is located at:

```text
output/VLASS_RACS_Updated.db
```
This database contains only the newly added records and their associated measurements.

The update Source Summary Table CSV file is located at:

```text
output/summary_table.csv
```
---

## CASDA Cache

Downloaded CASDA files are stored in:

```text
data/cache/
```

Previously downloaded files are reused whenever possible, reducing download time on future runs.

The cache may be deleted at any time. Files will be downloaded again automatically when needed.

---

## Submitting Updates

After verifying the new entries, upload `VLASS_RACS_Update.db` to a file-sharing service such as OneDrive or Google Drive.

Email the download link to:

**jprenderg@outlook.com**
