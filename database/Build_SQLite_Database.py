"""Build and validate output/final/OCVRC.db from the five final pandas tables.

Run standalone in Spyder or as the last Run_OCVRC.py stage.
Requires Python stdlib sqlite3, pandas and numpy.
"""
from pathlib import Path
import configparser
import os
import sqlite3
import tempfile
from contextlib import closing
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
config = configparser.ConfigParser(interpolation=None)
if not config.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')
out_setting = Path(config.get('PATHS', 'OUTPUT_DIR', fallback='output')).expanduser()
OUTPUT = out_setting if out_setting.is_absolute() else ROOT / out_setting
FINAL = OUTPUT / 'final'
DB_PATH = FINAL / 'OCVRC.db'

SCHEMAS = {
    'source': ['SOURCE_ID', 'RA', 'DEC', 'SIMBAD', 'MAJOR_AXIS',
               'MINOR_AXIS', 'POSITION_ANGLE', 'DISTANCE', 'DISTANCE_ERROR',
               'PM_RA', 'PM_DEC'],
    'name': ['SOURCE_ID', 'NAME', 'DEFAULT_NAME'],
    'class': ['SOURCE_ID', 'CLASS', 'DEFAULT_CLASS'],
    'period': ['SOURCE_ID', 'TYPE', 'UNITS', 'PERIOD', 'PERIOD_ERROR', 'CONFIRMED'],
    'measurement': ['SOURCE_ID', 'OBSERVATORY', 'EPOCH', 'BAND', 'UNITS',
                    'MEASUREMENT', 'ERROR', 'POS_ERROR', 'NEG_ERROR', 'NEIGHBOR'],
}


def load_table(name):
    path = FINAL / f'{name}.pkl'
    if not path.is_file():
        raise FileNotFoundError(f'Missing final table: {path}')
    frame = pd.read_pickle(path)
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f'{path}: expected pandas DataFrame')
    expected = SCHEMAS[name]
    if list(frame.columns) != expected:
        raise ValueError(f'{path.name}: wrong columns\nExpected: {expected}\nActual: {list(frame.columns)}')
    return frame.copy()


def sql_value(value):
    if value is None or (not isinstance(value, (tuple, list, dict, np.ndarray)) and pd.isna(value)):
        return None
    if isinstance(value, (np.bool_, bool)):
        return int(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError('Nonfinite value encountered in input')
    if isinstance(value, (str, int, float, bytes)):
        return value
    raise TypeError(f'Unsupported SQLite value type: {type(value).__name__}')


def checked_boolean(value, column, nullable=False):
    value = sql_value(value)
    if value is None and nullable:
        return None
    if isinstance(value, str):
        raise ValueError(f'{column}: string boolean {value!r} is not permitted')
    if value not in (0, 1):
        raise ValueError(f'{column}: expected 0 or 1, got {value!r}')
    return int(value)


def insert_table(conn, name, frame):
    columns = SCHEMAS[name]
    table = f'{name}_table'
    placeholders = ', '.join(['?'] * len(columns))
    sql = f'INSERT INTO {table} ({", ".join(columns)}) VALUES ({placeholders})'
    rows = []
    flags = {'source': ('SIMBAD',), 'name': ('DEFAULT_NAME',),
             'class': ('DEFAULT_CLASS',), 'period': ('CONFIRMED',),
             'measurement': ('NEIGHBOR',)}[name]
    for row in frame.itertuples(index=False, name=None):
        values = [sql_value(x) for x in row]
        for col in flags:
            pos = columns.index(col)
            values[pos] = checked_boolean(values[pos], col,
                                           nullable=(col in ('CONFIRMED', 'NEIGHBOR')))
        rows.append(tuple(values))
    conn.executemany(sql, rows)
    return len(rows)


frames = {name: load_table(name) for name in SCHEMAS}
FINAL.mkdir(parents=True, exist_ok=True)
fd, temp_name = tempfile.mkstemp(prefix='.OCVRC_build_', suffix='.db', dir=FINAL)
os.close(fd)
temp_db = Path(temp_name)

try:
    # sqlite3 connection context managers commit/rollback but DO NOT close.
    # Explicit closing is required on Windows before os.replace().
    with closing(sqlite3.connect(temp_db)) as conn, conn:
        conn.execute('PRAGMA foreign_keys = ON')
        conn.executescript('''
        CREATE TABLE source_table (
            SOURCE_ID INTEGER PRIMARY KEY,
            RA REAL NOT NULL CHECK (RA >= 0 AND RA < 360),
            DEC REAL NOT NULL CHECK (DEC >= -90 AND DEC <= 90),
            SIMBAD INTEGER NOT NULL CHECK (SIMBAD IN (0, 1)),
            MAJOR_AXIS REAL, MINOR_AXIS REAL, POSITION_ANGLE REAL,
            DISTANCE REAL, DISTANCE_ERROR REAL, PM_RA REAL, PM_DEC REAL
        );
        CREATE TABLE name_table (
            SOURCE_ID INTEGER NOT NULL REFERENCES source_table(SOURCE_ID),
            NAME TEXT NOT NULL CHECK (length(trim(NAME)) > 0),
            DEFAULT_NAME INTEGER NOT NULL CHECK (DEFAULT_NAME IN (0, 1)),
            UNIQUE (SOURCE_ID, NAME)
        );
        CREATE TABLE class_table (
            SOURCE_ID INTEGER NOT NULL REFERENCES source_table(SOURCE_ID),
            CLASS TEXT NOT NULL CHECK (length(trim(CLASS)) > 0),
            DEFAULT_CLASS INTEGER NOT NULL CHECK (DEFAULT_CLASS IN (0, 1)),
            UNIQUE (SOURCE_ID, CLASS)
        );
        CREATE TABLE period_table (
            SOURCE_ID INTEGER NOT NULL REFERENCES source_table(SOURCE_ID),
            TYPE TEXT NOT NULL,
            UNITS TEXT NOT NULL,
            PERIOD REAL NOT NULL,
            PERIOD_ERROR REAL,
            CONFIRMED INTEGER CHECK (CONFIRMED IN (0, 1) OR CONFIRMED IS NULL),
            UNIQUE (SOURCE_ID, TYPE)
        );
        CREATE TABLE measurement_table (
            SOURCE_ID INTEGER NOT NULL REFERENCES source_table(SOURCE_ID),
            OBSERVATORY TEXT, EPOCH TEXT, BAND TEXT, UNITS TEXT,
            MEASUREMENT REAL, ERROR REAL, POS_ERROR REAL, NEG_ERROR REAL,
            NEIGHBOR INTEGER CHECK (NEIGHBOR IN (0, 1) OR NEIGHBOR IS NULL)
        );
        CREATE UNIQUE INDEX idx_name_one_default
            ON name_table(SOURCE_ID) WHERE DEFAULT_NAME = 1;
        CREATE UNIQUE INDEX idx_class_one_default
            ON class_table(SOURCE_ID) WHERE DEFAULT_CLASS = 1;
        CREATE INDEX idx_name_name ON name_table(NAME);
        CREATE INDEX idx_class_class ON class_table(CLASS);
        CREATE INDEX idx_period_source ON period_table(SOURCE_ID);
        CREATE INDEX idx_measurement_source ON measurement_table(SOURCE_ID);
        CREATE INDEX idx_measurement_survey_epoch
            ON measurement_table(OBSERVATORY, EPOCH);
        ''')
        counts = {name: insert_table(conn, name, frame)
                  for name, frame in frames.items()}
        missing_names = conn.execute('''
            SELECT SOURCE_ID FROM source_table s WHERE
            (SELECT COUNT(*) FROM name_table n
             WHERE n.SOURCE_ID=s.SOURCE_ID AND n.DEFAULT_NAME=1) <> 1
            LIMIT 10''').fetchall()
        if missing_names:
            raise ValueError(f'Sources without exactly one default name: {missing_names}')
        bad_classes = conn.execute('''
            SELECT SOURCE_ID FROM class_table GROUP BY SOURCE_ID
            HAVING SUM(DEFAULT_CLASS) <> 1 LIMIT 10''').fetchall()
        if bad_classes:
            raise ValueError(f'Classified sources without exactly one default class: {bad_classes}')
        fk_errors = conn.execute('PRAGMA foreign_key_check').fetchall()
        integrity = conn.execute('PRAGMA integrity_check').fetchone()[0]
        if fk_errors or integrity != 'ok':
            raise ValueError(f'SQLite validation failed: FK={fk_errors[:10]}, integrity={integrity}')
        for name, expected in counts.items():
            count = conn.execute(f'SELECT COUNT(*) FROM {name}_table').fetchone()[0]
            if count != expected:
                raise ValueError(f'{name} row count mismatch: {count} != {expected}')
        conn.commit()
    # Replace the previous database only when construction and checks succeeded.
    try:
        os.replace(temp_db, DB_PATH)
    except PermissionError as exc:
        raise PermissionError(
            f"Cannot replace {DB_PATH}. Close any SQLite browser or other "
            "application holding the existing database open, then rerun. "
            "The existing database has not been overwritten."
        ) from exc
except BaseException:
    try:
        temp_db.unlink(missing_ok=True)
    except PermissionError:
        print(f"WARNING: Could not remove temporary database: {temp_db}")
    raise

print('\nSQLITE DATABASE CREATED')
print('=' * 55)
for name, count in counts.items():
    print(f'{name + "_table":22s} {count:8,d} rows')
print('Foreign keys: OK')
print('Integrity check: OK')
print('Database:', DB_PATH)
