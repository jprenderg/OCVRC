# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:54:00 2026

@author: Joe
"""

# GET_VLASS_Images.py

import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy import wcs
from astropy.wcs.utils import proj_plane_pixel_scales


# ============================================================
# SETTINGS AND PATHS (config.txt)
# ============================================================
# Paths are relative to the project root unless absolute.
from pathlib import Path
import configparser

ROOT = Path(__file__).resolve().parents[1]
CFG = configparser.ConfigParser(interpolation=None)
if not CFG.read(ROOT / "config.txt", encoding="utf-8"):
    raise FileNotFoundError(ROOT / "config.txt")

def cfg_path(section, key):
    value = CFG.get(section, key, fallback="").strip()
    if not value:
        raise ValueError(f"Set [{section}] {key} in config.txt before running.")
    p = Path(value).expanduser()
    return str(p if p.is_absolute() else ROOT / p)

START_ROW = CFG.getint("IMAGE_RETRIEVAL", "START_ROW", fallback=0)
END_ROW_TEXT = CFG.get("IMAGE_RETRIEVAL", "END_ROW", fallback="NONE").strip().upper()
END_ROW = None if END_ROW_TEXT in ("NONE", "") else int(END_ROW_TEXT)
INDEX = CFG.getint("IMAGE_RETRIEVAL", "DISPLAY_INDEX", fallback=-1)
path_catalog = cfg_path("PATHS", "SOURCE_PKL")
output_dir = cfg_path("PATHS", "IMAGE_OUTPUT_DIR")

BOX_SIZE = CFG.getint("VLASS", "BOX_SIZE", fallback=11)
RA_LIMIT = CFG.getfloat("VLASS", "RA_LIMIT", fallback=0.5)
DEC_LIMIT = CFG.getfloat("VLASS", "DEC_LIMIT", fallback=0.5)
path_vlass1 = cfg_path("PATHS", "VLASS1_LINKS")
path_vlass23 = cfg_path("PATHS", "VLASS23_LINKS")
if BOX_SIZE % 2 == 0 or BOX_SIZE < 3:
    raise ValueError("BOX_SIZE must be odd and >=3")
half_box = BOX_SIZE // 2

# ============================================================
# READ CATALOG
# ============================================================

df_catalog = pd.read_pickle(path_catalog).copy()
df_catalog = df_catalog.reset_index(drop=True)

required = ["SOURCE_ID", "RA", "DEC"]

missing_columns = [
    col for col in required
    if col not in df_catalog.columns
]

if missing_columns:
    raise ValueError(
        f"Missing catalog columns: {missing_columns}"
    )

n_catalog = len(df_catalog)

# END_ROW = None or -1 means process all sources
if END_ROW is None or END_ROW <= 0:
    END_ROW = n_catalog

# Prevent END_ROW from exceeding catalog size
END_ROW = min(END_ROW, n_catalog)

if not (0 <= START_ROW < END_ROW):
    raise ValueError(
        f"Invalid loop bounds: START_ROW={START_ROW}, "
        f"END_ROW={END_ROW}, catalog size={n_catalog}"
    )

df_processed = (
    df_catalog.iloc[START_ROW:END_ROW]
    .copy()
    .reset_index(drop=True)
)

partial_run = (
    START_ROW != 0 or END_ROW != n_catalog
)

suffix = "_test" if partial_run else ""

path_out = os.path.join(
    output_dir,
    f"VLASS_images{suffix}.pkl"
)

path_bad = os.path.join(
    output_dir,
    f"VLASS_bad_images{suffix}.pkl"
)

path_missing = os.path.join(
    output_dir,
    f"VLASS_missing_sources{suffix}.pkl"
)

print()
print("CATALOG")
print("========================================")
print("Total sources:", n_catalog)
print("Start row:", START_ROW)
print("End row:", END_ROW)
print("Sources to process:", len(df_processed))


# ============================================================
# READ VLASS IMAGE LISTS
# ============================================================

def read_vlass_list(path, epoch_code):

    with open(path) as f:
        img_list = [
            line.rstrip("\r\n")
            for line in f
            if line.strip()
        ]

    if epoch_code == 1:
        ra_slice = slice(83, 89)
        dec_slice = slice(89, 96)

    elif epoch_code == 23:
        ra_slice = slice(81, 87)
        dec_slice = slice(87, 94)

    else:
        raise ValueError("Unexpected epoch_code")

    rows = []
    rejected = 0

    for link in img_list:

        try:

            # --------------------------------------------
            # KEEP INTENSITY IMAGES
            # --------------------------------------------

            image_type = link[-14:-11]

            if image_type != "tt0":
                continue

            epoch = link[45:53]

            # --------------------------------------------
            # RA: HHMMSS
            # --------------------------------------------

            ra_str = link[ra_slice]

            ra = (
                int(ra_str[:2]) * 15
                + int(ra_str[2:4]) / 60 * 15
                + float(ra_str[4:6]) / 3600 * 15
            ) % 360

            # --------------------------------------------
            # DEC: +/-DDMMSS OR +/-DDMM
            # --------------------------------------------

            dec_str = link[dec_slice].strip()

            sign = -1 if dec_str.startswith("-") else 1
            dec_abs = dec_str.lstrip("+-")

            if len(dec_abs) not in (4, 6):
                raise ValueError(
                    f"Unexpected DEC format: {dec_str}"
                )

            degrees = int(dec_abs[:2])
            minutes = int(dec_abs[2:4])

            seconds = (
                float(dec_abs[4:6])
                if len(dec_abs) == 6
                else 0.0
            )

            if minutes >= 60 or seconds >= 60:
                raise ValueError(
                    f"Invalid DEC components: {dec_str}"
                )

            dec = sign * (
                degrees
                + minutes / 60
                + seconds / 3600
            )

            if not -90 <= dec <= 90:
                raise ValueError(
                    f"Invalid DEC: {dec_str} -> {dec}"
                )

            rows.append({
                "img_link": link,
                "VL_radeg": ra,
                "VL_decdeg": dec,
                "epoch": epoch
            })

        except Exception as error:

            rejected += 1

            print(
                "Invalid image-list entry:",
                error
            )

    print()
    print("Image list:", os.path.basename(path))
    print("Valid images:", len(rows))
    print("Invalid entries:", rejected)

    return pd.DataFrame(
        rows,
        columns=[
            "img_link",
            "VL_radeg",
            "VL_decdeg",
            "epoch"
        ]
    )


# ============================================================
# BUILD VLASS IMAGE TABLE
# ============================================================

vlass1 = read_vlass_list(
    path_vlass1,
    1
)

vlass23 = read_vlass_list(
    path_vlass23,
    23
)

vlass_data = pd.concat(
    [vlass1, vlass23],
    ignore_index=True
)

vlass_data = (
    vlass_data
    .drop_duplicates(subset=["img_link"])
    .reset_index(drop=True)
)

if vlass_data.empty:
    raise ValueError("No VLASS images found.")

print()
print("Total VLASS images:", len(vlass_data))


# ============================================================
# OUTPUT STORAGE
# ============================================================

results = []
bad_images = []
no_match_ids = []

start_time = time.time()


# ============================================================
# PROCESS SOURCES
# ============================================================

for j, source in df_processed.iterrows():

    source_id = source["SOURCE_ID"]

    radeg = float(source["RA"])
    decdeg = float(source["DEC"])

    # --------------------------------------------------------
    # PROGRESS
    # --------------------------------------------------------

    if j % 100 == 0:

        elapsed = (time.time() - start_time) / 60

        print(
            f"Row {START_ROW + j} / {END_ROW} "
            f"SOURCE_ID {source_id} "
            f"Elapsed: {elapsed:.1f} min"
        )

    # --------------------------------------------------------
    # CORRECTED VLASS TILE MATCHING
    #
    # Corrects:
    #   - Declination-dependent RA distances
    #   - RA wraparound at 0/360 degrees
    #
    # Retains original 0.5-degree search limits.
    # --------------------------------------------------------

    dra = (
        (vlass_data["VL_radeg"] - radeg + 180) % 360
    ) - 180

    ddec = vlass_data["VL_decdeg"] - decdeg

    cosdec = np.cos(np.radians(decdeg))

    match = vlass_data[
        (np.abs(dra * cosdec) <= RA_LIMIT) &
        (np.abs(ddec) <= DEC_LIMIT)
    ]

    # --------------------------------------------------------
    # NO MATCHING TILES
    # --------------------------------------------------------

    if match.empty:

        no_match_ids.append(source_id)

        bad_images.append({
            "SOURCE_ID": source_id,
            "EPOCH": None,
            "REASON": "No matching tiles",
            "NAN_PIXELS": np.nan,
            "INF_PIXELS": np.nan,
            "FILENAME": None
        })

        continue

    # --------------------------------------------------------
    # PROCESS ALL MATCHING IMAGES
    #
    # Each image has its own exception handler.
    # --------------------------------------------------------

    for _, row in match.iterrows():

        link = row["img_link"]
        epoch = row["epoch"]

        try:

            with fits.open(
                link,
                memmap=True
            ) as hdul:

                header = hdul[0].header
                data = hdul[0].data

                # ----------------------------------------
                # IMAGE
                # ----------------------------------------

                if data.ndim == 4:

                    image = data[0, 0, :, :]

                elif data.ndim == 2:

                    image = data

                else:

                    raise ValueError(
                        f"Unexpected dimensions: {data.ndim}"
                    )

                ny, nx = image.shape

                # ----------------------------------------
                # WCS
                # ----------------------------------------

                mywcs = (
                    wcs.WCS(header)
                    .celestial
                )

                # ----------------------------------------
                # EXACT SOURCE COORDINATES
                # ----------------------------------------

                x_wcs, y_wcs = (
                    mywcs.world_to_pixel_values(
                        radeg,
                        decdeg
                    )
                )

                if not (
                    np.isfinite(x_wcs)
                    and np.isfinite(y_wcs)
                ):

                    raise ValueError(
                        "Non-finite WCS coordinates"
                    )

                # ----------------------------------------
                # NEAREST PIXEL
                # ----------------------------------------

                xpix = int(np.round(x_wcs))
                ypix = int(np.round(y_wcs))

                # ----------------------------------------
                # VERIFY COMPLETE CUTOUT FITS
                # ----------------------------------------

                valid = (
                    half_box <= xpix
                    and xpix + half_box < nx
                    and half_box <= ypix
                    and ypix + half_box < ny
                )

                if not valid:

                    bad_images.append({
                        "SOURCE_ID": source_id,
                        "EPOCH": epoch,
                        "REASON": "Outside image boundaries",
                        "NAN_PIXELS": np.nan,
                        "INF_PIXELS": np.nan,
                        "FILENAME": link
                    })

                    continue

                # ----------------------------------------
                # SOURCE POSITION WITHIN CUTOUT
                # ----------------------------------------

                source_x = (
                    half_box + x_wcs - xpix
                )

                source_y = (
                    half_box + y_wcs - ypix
                )

                # ----------------------------------------
                # EXTRACT 11x11 CUTOUT
                # ----------------------------------------

                reduced_image = image[
                    ypix-half_box:ypix+half_box+1,
                    xpix-half_box:xpix+half_box+1
                ].copy()

                # ----------------------------------------
                # CHECK DIMENSIONS
                # ----------------------------------------

                if reduced_image.shape != (
                    BOX_SIZE,
                    BOX_SIZE
                ):

                    raise ValueError(
                        "Incomplete cutout"
                    )

                # ----------------------------------------
                # CHECK FINITE PIXELS
                # ----------------------------------------

                n_nan = int(
                    np.isnan(reduced_image).sum()
                )

                n_inf = int(
                    np.isinf(reduced_image).sum()
                )

                if n_nan or n_inf:

                    bad_images.append({
                        "SOURCE_ID": source_id,
                        "EPOCH": epoch,
                        "REASON": "Non-finite pixels",
                        "NAN_PIXELS": n_nan,
                        "INF_PIXELS": n_inf,
                        "FILENAME": link
                    })

                    continue

                # ----------------------------------------
                # BEAM PARAMETERS
                # ----------------------------------------

                bmaj = float(header["BMAJ"])
                bmin = float(header["BMIN"])
                bpa = float(header["BPA"])

                # ----------------------------------------
                # PIXEL SCALES
                # ----------------------------------------

                pixscale_x, pixscale_y = (
                    proj_plane_pixel_scales(mywcs)
                )

                pixscale_x = float(abs(pixscale_x))
                pixscale_y = float(abs(pixscale_y))

                # ----------------------------------------
                # STORE RESULT
                # ----------------------------------------

                results.append({
                    "SOURCE_ID": source_id,
                    "RA_DEG": radeg,
                    "DEC_DEG": decdeg,
                    "EPOCH": epoch,
                    "FILENAME": link,
                    "SOURCE_X": source_x,
                    "SOURCE_Y": source_y,
                    "BMAJ": bmaj,
                    "BMIN": bmin,
                    "BPA": bpa,
                    "PIXSCALE_X": pixscale_x,
                    "PIXSCALE_Y": pixscale_y,
                    "IMAGE": reduced_image
                })

        # ----------------------------------------------------
        # IMAGE EXCEPTION
        #
        # Continue with remaining images/epochs.
        # ----------------------------------------------------

        except Exception as error:

            bad_images.append({
                "SOURCE_ID": source_id,
                "EPOCH": epoch,
                "REASON": str(error),
                "NAN_PIXELS": np.nan,
                "INF_PIXELS": np.nan,
                "FILENAME": link
            })

            print(
                "Image exception:",
                source_id,
                epoch,
                error
            )


# ============================================================
# FINAL DATAFRAME
# ============================================================

columns = [
    "SOURCE_ID",
    "RA_DEG",
    "DEC_DEG",
    "EPOCH",
    "FILENAME",
    "SOURCE_X",
    "SOURCE_Y",
    "BMAJ",
    "BMIN",
    "BPA",
    "PIXSCALE_X",
    "PIXSCALE_Y",
    "IMAGE"
]

df_final = pd.DataFrame(
    results,
    columns=columns
)

df_final = (
    df_final
    .sort_values(
        ["SOURCE_ID", "EPOCH"]
    )
    .reset_index(drop=True)
)


# ============================================================
# BAD IMAGE DATAFRAME
# ============================================================

df_bad = pd.DataFrame(
    bad_images,
    columns=[
        "SOURCE_ID",
        "EPOCH",
        "REASON",
        "NAN_PIXELS",
        "INF_PIXELS",
        "FILENAME"
    ]
)


# ============================================================
# MISSING SOURCES
#
# Only sources within the selected loop bounds.
# ============================================================

df_missing = df_processed[
    ~df_processed["SOURCE_ID"].isin(
        df_final["SOURCE_ID"]
    )
].copy()

df_missing = (
    df_missing
    .sort_values("SOURCE_ID")
    .reset_index(drop=True)
)


# ============================================================
# VALIDATE STORED IMAGES
# ============================================================

bad_final = df_final["IMAGE"].apply(
    lambda image: (
        image.shape != (BOX_SIZE, BOX_SIZE)
        or not np.all(np.isfinite(image))
    )
)


# ============================================================
# SUMMARY
# ============================================================

elapsed = (time.time() - start_time) / 60

print()
print("VLASS IMAGE SUMMARY")
print("========================================")

print(
    "Catalog sources:",
    n_catalog
)

print(
    "Processed sources:",
    len(df_processed)
)

print(
    "Sources without matching tiles:",
    len(no_match_ids)
)

print(
    "Stored cutouts:",
    len(df_final)
)

print(
    "Unique sources with images:",
    df_final["SOURCE_ID"].nunique()
)

print(
    "Sources without images:",
    len(df_missing)
)

print(
    "Rejected images:",
    len(df_bad)
)

print(
    "Invalid stored images:",
    int(bad_final.sum())
)

print(
    "Elapsed time:",
    f"{elapsed:.1f} min"
)


# ============================================================
# REJECTION REASONS
# ============================================================

print()
print("REJECTION REASONS")
print("----------------------------------------")

if not df_bad.empty:

    print(
        df_bad["REASON"].value_counts()
    )


# ============================================================
# MISSING SOURCE DECLINATIONS
# ============================================================

print()
print("MISSING SOURCE DECLINATIONS")
print("----------------------------------------")

if not df_missing.empty:

    print(
        "DEC >= -40:",
        int((df_missing["DEC"] >= -40).sum())
    )

    print(
        "DEC < -40:",
        int((df_missing["DEC"] < -40).sum())
    )


# ============================================================
# SAVE
# ============================================================

os.makedirs(
    output_dir,
    exist_ok=True
)

df_final.to_pickle(path_out)
df_bad.to_pickle(path_bad)
df_missing.to_pickle(path_missing)

print()
print("SAVED FILES")
print("----------------------------------------")

print(path_out)
print(path_bad)
print(path_missing)


# ============================================================
# DISPLAY EXAMPLE IMAGE
# ============================================================

if 0 <= INDEX < len(df_final):

    row = df_final.iloc[INDEX]

    image = row["IMAGE"]

    source_x = row["SOURCE_X"]
    source_y = row["SOURCE_Y"]

    plt.figure(figsize=(6, 5))

    plt.imshow(
        image * 1000,
        origin="lower",
        cmap="hot"
    )

    plt.scatter(
        source_x,
        source_y,
        marker="+",
        s=150,
        linewidths=2
    )

    plt.colorbar(
        label="Flux (mJy/beam)"
    )

    plt.xlabel("X pixel")
    plt.ylabel("Y pixel")

    plt.title(
        f"SOURCE_ID {row['SOURCE_ID']}\n"
        f"Epoch {row['EPOCH']}   "
        f"Position = "
        f"({source_x:.3f}, {source_y:.3f})"
    )

    plt.tight_layout()
    plt.show()

    # --------------------------------------------------------
    # IMAGE INFORMATION
    # --------------------------------------------------------

    print()
    print("PLOTTED IMAGE")
    print("----------------------------------------")

    print("Index       =", INDEX)
    print("SOURCE_ID   =", row["SOURCE_ID"])
    print("Epoch       =", row["EPOCH"])
    print("Filename    =", row["FILENAME"])

    print(
        "Position    =",
        f"({source_x:.3f}, {source_y:.3f})"
    )

    print(
        "Beam        =",
        f"{row['BMAJ']*3600:.3f} x "
        f"{row['BMIN']*3600:.3f} arcsec FWHM"
    )

    print(
        "Pixel scale =",
        f"{row['PIXSCALE_X']*3600:.3f} x "
        f"{row['PIXSCALE_Y']*3600:.3f} arcsec/pixel"
    )

    print(
        "Image shape =",
        image.shape
    )

else:

    print()
    print(
        f"No image available at INDEX = {INDEX}. "
        f"df_final contains {len(df_final)} rows."
    )
# Friendly aliases for Spyder Variable Explorer
df_vlass_images = df_final
df_vlass_bad_images = df_bad
df_vlass_missing = df_missing
