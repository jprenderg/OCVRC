"""OCVRC_Github: build source, name, class and period tables from CSV/XLSX and SIMBAD.

Required columns: RA, DEC (decimal degrees).
Optional: SOURCE_ID, NAME, CLASS, PORB, PORB_ERROR, PORB_CONFIRMED,
PSPIN, PSPIN_ERROR, PSPIN_CONFIRMED. Missing SOURCE_ID is assigned sequentially.
Install: pip install pandas numpy astropy astroquery
Run in Spyder using the current console namespace (Run > Configuration per file).
"""
from __future__ import annotations

import configparser
import time
import hashlib
import json
import os
import uuid
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CFG = configparser.ConfigParser(interpolation=None)
if not CFG.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')


def setting_path(key):
    p = Path(CFG.get('PATHS', key)).expanduser()
    return p if p.is_absolute() else ROOT / p


INPUT = setting_path('INPUT_FILE') if CFG.has_option('PATHS', 'INPUT_FILE') else setting_path('INPUT_CSV')
OUT = setting_path('OUTPUT_DIR')
RADIUS = CFG.getfloat('SIMBAD', 'MATCH_RADIUS_ARCSEC')
TIMEOUT = CFG.getint('SIMBAD', 'TIMEOUT_SECONDS')
RETRIES = CFG.getint('SIMBAD', 'MAX_RETRIES')
DELAY = CFG.getfloat('SIMBAD', 'RETRY_DELAY_SECONDS')
SKIP = CFG.getboolean('PROCESSING', 'SKIP_SIMBAD')
RESUME = CFG.getboolean('QUERY_MANAGEMENT', 'RESUME', fallback=True)
CACHE = OUT / 'simbad_cache'
CACHE_VERSION = 1



class SimbadUnavailableError(RuntimeError):
    """SIMBAD cannot be reached and no verified offline results are available."""


def _file_sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as file:
        for part in iter(lambda: file.read(1024 * 1024), b''):
            digest.update(part)
    return digest.hexdigest()


def _input_fingerprint():
    # Hash input bytes AND SIMBAD-related configuration to avoid stale reuse.
    settings = {'RADIUS': RADIUS, 'SKIP_SIMBAD': SKIP, 'CACHE_VERSION': CACHE_VERSION}
    return hashlib.sha256((_file_sha256(INPUT) + json.dumps(settings, sort_keys=True)).encode()).hexdigest()


def _output_paths():
    return {key: OUT / (key + '.pkl') for key in ('source', 'name', 'class', 'period')}


def _legacy_outputs_match_input(df, source, name, klass, period):
    """Conservative offline check of pre-manifest, already completed outputs.

    Does not purport to reverify SIMBAD remotely or recover raw SIMBAD query data.
    """
    expected_ids = df['SOURCE_ID'].astype('int64').tolist()
    if len(source) != len(df) or source['SOURCE_ID'].tolist() != expected_ids:
        return False
    if not np.array_equal(source['RA'].to_numpy(dtype=float), df['RA'].to_numpy(dtype=float)):
        return False
    if not np.array_equal(source['DEC'].to_numpy(dtype=float), df['DEC'].to_numpy(dtype=float)):
        return False
    required = {'SIMBAD_FLAG', 'SIMBAD_MAIN_ID', 'MAJOR_AXIS', 'MINOR_AXIS', 'POSITION_ANGLE'}
    if not required.issubset(source.columns):
        return False
    if source['SIMBAD_FLAG'].isna().any():
        return False
    names_by_id = name.groupby('SOURCE_ID', sort=False)
    classes_by_id = klass.groupby('SOURCE_ID', sort=False)
    if set(name['SOURCE_ID']) != set(expected_ids) or not set(klass['SOURCE_ID']).issubset(set(expected_ids)):
        return False
    if name.duplicated(['SOURCE_ID', 'NAME']).any() or klass.duplicated(['SOURCE_ID', 'CLASS']).any():
        return False
    if not names_by_id['DEFAULT_NAME'].sum().eq(1).all():
        return False
    if len(klass) and not classes_by_id['DEFAULT_CLASS'].sum().eq(1).all():
        return False
    for row in df.itertuples(index=False):
        sid = int(row.SOURCE_ID)
        srow = source.loc[source.SOURCE_ID == sid].iloc[0]
        matched = bool(srow.SIMBAD_FLAG)
        if matched != (clean(srow.SIMBAD_MAIN_ID) is not None):
            return False
        if not matched and any(pd.notna(srow[x]) for x in ('MAJOR_AXIS', 'MINOR_AXIS', 'POSITION_ANGLE')):
            return False
        default_name = name.loc[(name.SOURCE_ID == sid) & name.DEFAULT_NAME.isin([1, True]), 'NAME']
        expected_name = clean(row.NAME) or (f'User_name_{sid}' if not matched else None)
        if expected_name is not None and (len(default_name) != 1 or clean(default_name.iloc[0]) != expected_name):
            return False
        supplied_class = clean(row.CLASS)
        if supplied_class:
            default_class = klass.loc[(klass.SOURCE_ID == sid) & klass.DEFAULT_CLASS.isin([1, True]), 'CLASS']
            if len(default_class) != 1 or clean(default_class.iloc[0]) != supplied_class:
                return False
    expected_period = make_period_table(df)
    if set(period.columns) != set(expected_period.columns) or len(period) != len(expected_period):
        return False
    try:
        pd.testing.assert_frame_equal(
            period[expected_period.columns].reset_index(drop=True),
            expected_period.reset_index(drop=True), check_dtype=False, check_exact=False,
            rtol=1e-12, atol=1e-12)
    except AssertionError:
        return False
    return True


def load_reusable_outputs(df):
    """Reuse validated prior outputs without querying SIMBAD when RESUME=TRUE.

    A stored manifest binds output hashes to exact input bytes and query settings.
    Older, pre-manifest runs undergo a conservative input/outputs consistency audit.
    """
    if not RESUME or SKIP:
        return None
    paths = _output_paths()
    if not all(path.is_file() for path in paths.values()):
        return None
    manifest_path = OUT / 'source_ingestion_manifest.json'
    identity = _input_fingerprint()
    if manifest_path.exists():
        try:
            record = json.loads(manifest_path.read_text(encoding='utf-8'))
            if record.get('input_fingerprint') != identity:
                print('[SOURCE_INGESTION] Input/configuration changed: previous outputs cannot be reused.')
                return None
            if record.get('files') != {k: _file_sha256(v) for k, v in paths.items()}:
                print('[SOURCE_INGESTION] Output hashes changed: cache manifest invalid.')
                return None
        except (OSError, ValueError) as exc:
            print(f'[SOURCE_INGESTION] Invalid saved manifest: {exc}')
            return None
    try:
        tables = [pd.read_pickle(paths[key]) for key in ('source', 'name', 'class', 'period')]
        # Compatibility: previously verified ingestion tables use SIMBAD's old column names.
        legacy = {'COO_ERR_MAJ': 'MAJOR_AXIS', 'COO_ERR_MIN': 'MINOR_AXIS',
                  'COO_ERR_ANGLE': 'POSITION_ANGLE'}
        tables[0] = tables[0].rename(columns=legacy)
        if not _legacy_outputs_match_input(df, *tables):
            print('[SOURCE_INGESTION] Existing source tables do not match current input.')
            return None
        validate(*tables)
    except Exception as exc:
        print(f'[SOURCE_INGESTION] Previous tables failed validation: {type(exc).__name__}: {exc}')
        return None
    if not manifest_path.exists():
        # Manifest bootstrapping is limited to audited, previously complete outputs.
        print('[SOURCE_INGESTION] WARNING: audited legacy outputs; original SIMBAD responses')
        print('  cannot be independently verified offline. No SIMBAD cache is fabricated.')
        save_manifest()
    print(f'[SOURCE_INGESTION] Reusing validated outputs for {len(df)} sources; SIMBAD not contacted.')
    return tables


def save_manifest():
    OUT.mkdir(parents=True, exist_ok=True)
    record = {'input_fingerprint': _input_fingerprint(),
              'files': {k: _file_sha256(v) for k, v in _output_paths().items()}}
    path = OUT / 'source_ingestion_manifest.json'
    temp = OUT / f'.source_ingestion_manifest_{uuid.uuid4().hex}.tmp'
    try:
        temp.write_text(json.dumps(record, indent=2), encoding='utf-8')
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def clean(v):
    if v is None or pd.isna(v):
        return None
    s = str(v).strip()
    return None if not s or s.lower() in ('nan', 'none', 'null') else s


def unique(values):
    return list(dict.fromkeys(v for v in (clean(x) for x in values) if v))


def read_input():
    if INPUT.suffix.lower() == '.csv':
        df = pd.read_csv(INPUT, dtype={'NAME': 'string', 'CLASS': 'string'})
    elif INPUT.suffix.lower() == '.xlsx':
        df = pd.read_excel(INPUT, dtype={'NAME': 'string', 'CLASS': 'string'})
    else:
        raise ValueError('Input file must be .csv or .xlsx')
    df.columns = df.columns.str.strip().str.upper()
    for col in ('RA', 'DEC'):
        if col not in df:
            raise ValueError(f'Missing required column {col}')
        df[col] = pd.to_numeric(df[col], errors='coerce')
    bad = (~np.isfinite(df['RA']) | ~np.isfinite(df['DEC']) |
           ~df['RA'].between(0, 360, inclusive='left') | ~df['DEC'].between(-90, 90))
    if bad.any():
        raise ValueError(f'Invalid RA/DEC in CSV row numbers: {(df.index[bad] + 2).tolist()}')
    if 'SOURCE_ID' not in df:
        df.insert(0, 'SOURCE_ID', np.arange(1, len(df) + 1, dtype=int))
    else:
        ids = pd.to_numeric(df['SOURCE_ID'], errors='coerce')
        if ids.isna().any() or (ids % 1 != 0).any() or (ids <= 0).any() or ids.duplicated().any():
            raise ValueError('SOURCE_ID must contain unique positive integers.')
        df['SOURCE_ID'] = ids.astype('int64')
    for col in ('NAME', 'CLASS'):
        if col not in df:
            df[col] = None
        df[col] = df[col].map(clean)
    # Periods must be supplied or explicitly represented by blank values.
    # Each period is in a fixed unit: PORB hours, PSPIN seconds.
    for prefix in ('PORB', 'PSPIN'):
        for col in (prefix, f'{prefix}_ERROR'):
            if col not in df:
                df[col] = np.nan
            original = df[col]
            df[col] = pd.to_numeric(original, errors='coerce')
            invalid = original.notna() & original.astype(str).str.strip().ne('') & df[col].isna()
            if invalid.any():
                raise ValueError(f'{col}: invalid numeric values on CSV/Excel data rows '
                                 f'{(df.index[invalid] + 2).tolist()}')
        if (df[prefix].notna() & (df[prefix] <= 0)).any():
            raise ValueError(f'{prefix} must be positive where supplied.')
        if (df[f'{prefix}_ERROR'].notna() & (df[f'{prefix}_ERROR'] < 0)).any():
            raise ValueError(f'{prefix}_ERROR must be nonnegative.')
        bad = df[f'{prefix}_ERROR'].notna() & df[prefix].isna()
        if bad.any():
            raise ValueError(f'{prefix}_ERROR supplied without {prefix} on rows '
                             f'{(df.index[bad] + 2).tolist()}')
        confirmed = f'{prefix}_CONFIRMED'
        if confirmed not in df:
            df[confirmed] = False
        def parse_confirmation(value):
            if pd.isna(value) or str(value).strip() == '':
                return False
            if isinstance(value, (bool, np.bool_)):
                return bool(value)
            val = str(value).strip().lower()
            if val in ('true', 'yes', '1', 'y'):
                return True
            if val in ('false', 'no', '0', 'n'):
                return False
            raise ValueError(f'Invalid {confirmed}: {value!r}; use TRUE/FALSE.')
        df[confirmed] = df[confirmed].map(parse_confirmation)
        if (df[confirmed] & df[prefix].isna()).any():
            raise ValueError(f'{confirmed}=TRUE without {prefix}.')
    return df


def make_period_table(df):
    records = []
    for row in df.itertuples(index=False):
        for prefix, period_type, units in [('PORB', 'Orbital', 'Hours'),
                                           ('PSPIN', 'Spin', 'Seconds')]:
            val = getattr(row, prefix)
            if pd.notna(val):
                error = getattr(row, f'{prefix}_ERROR')
                records.append({'SOURCE_ID': int(row.SOURCE_ID), 'TYPE': period_type,
                                'UNITS': units, 'PERIOD': float(val),
                                'PERIOD_ERROR': float(error) if pd.notna(error) else np.nan,
                                'CONFIRMED': bool(getattr(row, f'{prefix}_CONFIRMED'))})
    return pd.DataFrame(records, columns=['SOURCE_ID', 'TYPE', 'UNITS', 'PERIOD',
                                          'PERIOD_ERROR', 'CONFIRMED'])


def make_client():
    from astroquery.simbad import Simbad
    s = Simbad()
    s.TIMEOUT = TIMEOUT
    s.add_votable_fields('otype', 'coo_err_maj', 'coo_err_min', 'coo_err_angle')
    return s


def query_simbad(simbad, ra, dec):
    """Return identifiers, classes, separation, and SIMBAD error ellipse.

    Query errors are fatal after retries; they are NOT interpreted as no match.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    position = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    def positional_error(result, row, key, unit):
        # Astroquery's SIMBAD fields are normally in mas, mas, and degrees.
        if key not in result.colnames:
            raise ValueError(f'SIMBAD query did not return {key}')
        col = result[key]
        value = col[row]
        if np.ma.is_masked(value) or value is None:
            return np.nan
        try:
            if getattr(col, 'unit', None) is not None:
                return float((float(value) * col.unit).to_value(unit))
            return float(value)
        except (TypeError, ValueError):
            return np.nan
    for attempt in range(1, RETRIES + 1):
        try:
            result = simbad.query_region(position, radius=RADIUS * u.arcsec)
            if result is None or len(result) == 0:
                return None, [], None, [], np.nan, np.nan, np.nan, np.nan
            coords = SkyCoord(ra=np.array(result['ra'], dtype=float) * u.deg,
                              dec=np.array(result['dec'], dtype=float) * u.deg)
            separations = position.separation(coords).arcsec
            j = int(np.argmin(separations))
            if separations[j] > RADIUS:
                return None, [], None, [], np.nan, np.nan, np.nan, np.nan
            main_id = clean(result['main_id'][j])
            primary_class = clean(result['otype'][j]) if 'otype' in result.colnames else None
            ids = simbad.query_objectids(main_id)
            aliases = unique([] if ids is None else list(ids['id']))
            # Retrieve every additional otype for the matched SIMBAD object.
            # Escape single quotes for ADQL literals.
            escaped = main_id.replace("'", "''")
            tap_query = ("SELECT DISTINCT otypes.otype AS alternative_class "
                         "FROM basic JOIN otypes ON basic.oid = otypes.oidref "
                         f"WHERE basic.main_id = '{escaped}'")
            types = simbad.query_tap(tap_query)
            alt_classes = unique([] if types is None else list(types['alternative_class']))
            err_maj = positional_error(result, j, 'coo_err_maj', u.mas)
            err_min = positional_error(result, j, 'coo_err_min', u.mas)
            err_angle = positional_error(result, j, 'coo_err_angle', u.deg)
            return (main_id, aliases, primary_class, alt_classes,
                    float(separations[j]), err_maj, err_min, err_angle)
        except Exception as exc:
            message = str(exc)
            message = 'HTTP 503: Service Unavailable' if re.search(r'\b503\b', message) else (message.splitlines()[0][:180] if message else type(exc).__name__)
            print(f'  SIMBAD query failed (attempt {attempt}/{RETRIES}): {message}')
            if attempt == RETRIES:
                raise RuntimeError(f'SIMBAD query failed for RA={ra}, DEC={dec}. '
                                   'No output was written; retry after service recovery.') from exc
            time.sleep(DELAY)


def cached_resolver():
    """Reuse verified per-source SIMBAD results; never mistake outages for NO_MATCH.

    Cache identity includes SOURCE_ID, coordinates, search radius, and format version.
    Completed rows are persisted atomically, even if a later lookup fails.
    """
    client = None
    CACHE.mkdir(parents=True, exist_ok=True)

    def resolve(sid, ra, dec):
        nonlocal client
        identity = {'version': CACHE_VERSION, 'SOURCE_ID': int(sid),
                    'RA': float(ra), 'DEC': float(dec), 'RADIUS': float(RADIUS)}
        digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        path = CACHE / f'SIMBAD_{sid}_{digest[:20]}.pkl'
        if RESUME and path.exists():
            payload = pd.read_pickle(path)
            if payload.get('identity') == identity and len(payload.get('result', ())) == 8:
                print(f'SIMBAD SOURCE_ID={sid}: CACHED')
                return tuple(payload['result'])
            raise ValueError(f'Invalid SIMBAD cache record: {path}')
        if client is None:
            try:
                client = make_client()  # Only initialize for genuinely uncached sources.
            except Exception as exc:
                message = str(exc).splitlines()[0][:160]
                raise SimbadUnavailableError(
                    'SIMBAD initialization failed: ' + message +
                    '. Existing outputs were preserved; retry when service recovers.') from exc
        try:
            result = query_simbad(client, ra, dec)
        except Exception as exc:
            raise SimbadUnavailableError(
                f'SIMBAD query unavailable for SOURCE_ID={sid}. Existing outputs preserved.') from exc
        payload = {'identity': identity, 'result': tuple(result)}
        temp = path.with_name(f'.{path.name}.{uuid.uuid4().hex}.tmp')
        try:
            pd.to_pickle(payload, temp)
            os.replace(temp, path)
        finally:
            temp.unlink(missing_ok=True)
        return result

    return resolve


def make_tables(df, resolver):
    sources, names, classes, generated = [], [], [], []
    for i, row in df.iterrows():
        sid = int(row['SOURCE_ID'])
        user_name, user_class = clean(row['NAME']), clean(row['CLASS'])
        (primary, aliases, sim_class, other_classes, sep,
         err_maj, err_min, err_angle) = resolver(sid, float(row['RA']), float(row['DEC']))
        all_names = unique([user_name, primary, *aliases])
        if not all_names:
            auto = f'User_name_{sid}'
            all_names = [auto]
            generated.append((sid, auto))
        all_classes = unique([user_class, sim_class, *other_classes])
        names.extend({'SOURCE_ID': sid, 'NAME': name, 'DEFAULT_NAME': int(j == 0)}
                     for j, name in enumerate(all_names))
        classes.extend({'SOURCE_ID': sid, 'CLASS': klass, 'DEFAULT_CLASS': int(j == 0)}
                       for j, klass in enumerate(all_classes))
        sources.append({'SOURCE_ID': sid, 'RA': row['RA'], 'DEC': row['DEC'],
                        'SIMBAD_FLAG': primary is not None, 'SIMBAD_MAIN_ID': primary,
                        'SIMBAD_SEPARATION_ARCSEC': sep,
                        'MAJOR_AXIS': err_maj,
                        'MINOR_AXIS': err_min,
                        'POSITION_ANGLE': err_angle})
        print(f'Processed {i+1}/{len(df)} SOURCE_ID={sid}', end='\r')
    print()
    source = pd.DataFrame(sources, columns=['SOURCE_ID', 'RA', 'DEC', 'SIMBAD_FLAG',
                                             'SIMBAD_MAIN_ID', 'SIMBAD_SEPARATION_ARCSEC',
                                             'MAJOR_AXIS', 'MINOR_AXIS', 'POSITION_ANGLE'])
    name = pd.DataFrame(names, columns=['SOURCE_ID', 'NAME', 'DEFAULT_NAME'])
    klass = pd.DataFrame(classes, columns=['SOURCE_ID', 'CLASS', 'DEFAULT_CLASS'])
    return source, name, klass, generated


def validate(source, name, klass, period):
    ids = set(source.SOURCE_ID)
    assert source.SOURCE_ID.is_unique
    assert set(name.SOURCE_ID) == ids
    assert name.groupby('SOURCE_ID').DEFAULT_NAME.sum().eq(1).all()
    assert not name.duplicated(['SOURCE_ID', 'NAME']).any()
    assert set(klass.SOURCE_ID).issubset(ids)
    assert klass.groupby('SOURCE_ID').DEFAULT_CLASS.sum().eq(1).all()
    assert not klass.duplicated(['SOURCE_ID', 'CLASS']).any()
    assert set(period.SOURCE_ID).issubset(ids)
    assert not period.duplicated(['SOURCE_ID', 'TYPE']).any()


# Top-level variables are deliberately assigned outside functions so Spyder's
# Variable Explorer can show them when run in the current IPython console.
df_input = read_input()
prior = load_reusable_outputs(df_input)
if prior is not None:
    df_source, df_name, df_class, df_period = prior
    generated_names = []
else:
    if SKIP:
        print('WARNING: SKIP_SIMBAD=TRUE: SIMBAD matching disabled for testing.')
        resolver = lambda sid, ra, dec: (None, [], None, [], np.nan, np.nan, np.nan, np.nan)
    else:
        resolver = cached_resolver()
    df_source, df_name, df_class, generated_names = make_tables(df_input, resolver)
    df_period = make_period_table(df_input)
    validate(df_source, df_name, df_class, df_period)
    OUT.mkdir(parents=True, exist_ok=True)
    for filename, dataframe in [('source.pkl', df_source), ('name.pkl', df_name),
                                ('class.pkl', df_class), ('period.pkl', df_period)]:
        dataframe.to_pickle(OUT / filename)
    # Never create a reusable manifest for a deliberately SIMBAD-skipped run.
    if not SKIP:
        save_manifest()

print('\n' + '=' * 50)
print('AUTOMATIC SOURCE NAME ASSIGNMENTS')
print('=' * 50)
if prior is None:
    for sid, value in generated_names:
        print(f'SOURCE_ID {sid}: Assigned name {value}')
    print(f'Total automatically assigned names: {len(generated_names)}')
else:
    print('Reused existing name assignments without modification.')
print('=' * 50)
print(f'Sources: {len(df_source)}; names: {len(df_name)}; '
      f'classes: {len(df_class)}; periods: {len(df_period)}')
print(f'Pickle output directory: {OUT}')
