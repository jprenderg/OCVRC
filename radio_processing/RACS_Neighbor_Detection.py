# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 09:20:17 2026

@author: Joe
"""

# RACS Bright neighbor and streak detection.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skew


# ============================================================
# SETTINGS
# ============================================================
from pathlib import Path
import configparser
ROOT = Path(__file__).resolve().parents[1]
CFG = configparser.ConfigParser(interpolation=None)
if not CFG.read(ROOT / "config.txt", encoding="utf-8"):
    raise FileNotFoundError(ROOT / "config.txt")
image_dir = Path(CFG.get("PATHS", "IMAGE_OUTPUT_DIR", fallback="output/images"))
image_dir = image_dir if image_dir.is_absolute() else ROOT / image_dir
image_dir.mkdir(parents=True, exist_ok=True)
path_in = image_dir / "RACS_images.pkl"
path_out = image_dir / "RACS_neighbor_detection.pkl"
if not path_in.exists():
    raise FileNotFoundError(f"Missing {path_in}; run Get_RACS_Images.py first")

S = 0.6260
FIT_RADIUS = 13.0
INDEX = 0
SHOW_PLOTS = CFG.getboolean("DIAGNOSTICS", "SHOW_PLOTS", fallback=False)

# ============================================================
# READ DATA
# ============================================================

df = pd.read_pickle(path_in).copy()

df["SKEWNESS"] = np.nan
df["NEIGHBOR_DETECTION"] = np.nan


# ============================================================
# ELLIPTICAL RACS PSF
#
# SOURCE_X and SOURCE_Y are always fixed at the catalog/WCS
# coordinates.
#
# BMAJ/BMIN     = degrees, FWHM
# PIXSCALE_X/Y  = degrees/pixel
# ============================================================

def gaussian(
    x, y, x0, y0,
    bmaj, bmin, bpa,
    pixscale_x, pixscale_y
):

    sigma_maj = bmaj / (2.35482 * pixscale_x)
    sigma_min = bmin / (2.35482 * pixscale_y)

    theta = np.deg2rad(bpa)

    dx = x - x0
    dy = y - y0

    major = dx*np.sin(theta) + dy*np.cos(theta)
    minor = dx*np.cos(theta) - dy*np.sin(theta)

    return np.exp(
        -0.5 * (
            (major / sigma_maj)**2 +
            (minor / sigma_min)**2
        )
    )


# ============================================================
# FIT AND SUBTRACT PRIMARY
#
# Position is fixed at SOURCE_X, SOURCE_Y.
#
# Fit:
#
#     image = a*G + b
#
# a = peak flux
# b = background
#
# Only a*G is subtracted.
# Background remains in residual image.
# ============================================================

def fit_primary(
    image,
    source_x,
    source_y,
    bmaj,
    bmin,
    bpa,
    pixscale_x,
    pixscale_y
):

    ny, nx = image.shape
    y, x = np.indices((ny, nx))

    g = gaussian(
        x, y,
        source_x,
        source_y,
        bmaj,
        bmin,
        bpa,
        pixscale_x,
        pixscale_y
    )

    r = np.hypot(
        x - source_x,
        y - source_y
    )

    fit_mask = (
        (r <= FIT_RADIUS) &
        np.isfinite(image) &
        np.isfinite(g)
    )

    if fit_mask.sum() < 3:
        raise ValueError(
            f"Too few finite fit pixels: {fit_mask.sum()}"
        )

    X = np.column_stack((
        g[fit_mask],
        np.ones(fit_mask.sum())
    ))

    z = image[fit_mask]

    peak_flux, background = np.linalg.lstsq(
        X,
        z,
        rcond=None
    )[0]

    fitted_primary = peak_flux * g

    subtracted_image = (
        image -
        fitted_primary
    )

    return (
        subtracted_image,
        peak_flux,
        background
    )


# ============================================================
# PROCESS IMAGES
# ============================================================

START_ROW = 0
STOP_ROW = len(df)

# STOP_ROW = 300

for i in range(
    START_ROW,
    STOP_ROW
):

    try:

        image = np.asarray(
            df.at[i, "IMAGE"],
            dtype=float
        )

        source_x = float(
            df.at[i, "SOURCE_X"]
        )

        source_y = float(
            df.at[i, "SOURCE_Y"]
        )

        bmaj = float(
            df.at[i, "BMAJ"]
        )

        bmin = float(
            df.at[i, "BMIN"]
        )

        bpa = float(
            df.at[i, "BPA"]
        )

        pixscale_x = float(
            df.at[i, "PIXSCALE_X"]
        )

        pixscale_y = float(
            df.at[i, "PIXSCALE_Y"]
        )

        (
            subtracted_image,
            peak_flux,
            background
        ) = fit_primary(
            image,
            source_x,
            source_y,
            bmaj,
            bmin,
            bpa,
            pixscale_x,
            pixscale_y
        )

        hist_skew = skew(
            subtracted_image.ravel(),
            nan_policy="omit"
        )

        df.at[
            i,
            "SKEWNESS"
        ] = hist_skew

        df.at[
            i,
            "NEIGHBOR_DETECTION"
        ] = int(
            hist_skew > S
        )

        print(
            i,
            "SOURCE_ID =", df.at[i, "SOURCE_ID"],
            "POSITION =",
            f"({source_x:.2f}, {source_y:.2f})",
            "PEAK =",
            f"{peak_flux*1000:.3f}",
            "BEAM =",
            f"{bmaj*3600:.2f} x {bmin*3600:.2f}",
            "BPA =",
            f"{bpa:.1f}",
            "SKEW =",
            round(hist_skew, 3),
            "NEIGHBOR =",
            int(hist_skew > S)
        )

    except Exception as error:

        print(
            "Exception:",
            i,
            error
        )


# ============================================================
# SAVE
# ============================================================

df.to_pickle(
    path_out
)

print()

print(
    "Processed rows:",
    START_ROW,
    "to",
    STOP_ROW - 1
)

print(
    "Neighbor detections:",
    int(
        df.loc[
            START_ROW:STOP_ROW-1,
            "NEIGHBOR_DETECTION"
        ]
        .fillna(0)
        .sum()
    )
)

print(
    "Saved:",
    path_out
)


# Optional plots and selected-image diagnostics. Disabled by default.
if SHOW_PLOTS and len(df) > 0:
    # ============================================================
    # DISPLAY SELECTED IMAGE
    # ============================================================

    if INDEX == -1:
        print("Image display disabled.")
    else:
        row = df.iloc[INDEX]
    
    image = np.asarray(
        row["IMAGE"],
        dtype=float
    )

    source_x = float(
        row["SOURCE_X"]
    )

    source_y = float(
        row["SOURCE_Y"]
    )

    bmaj = float(
        row["BMAJ"]
    )

    bmin = float(
        row["BMIN"]
    )

    bpa = float(
        row["BPA"]
    )

    pixscale_x = float(
        row["PIXSCALE_X"]
    )

    pixscale_y = float(
        row["PIXSCALE_Y"]
    )

    (
        subtracted_image,
        peak_flux,
        background
    ) = fit_primary(
        image,
        source_x,
        source_y,
        bmaj,
        bmin,
        bpa,
        pixscale_x,
        pixscale_y
    )

    hist_skew = skew(
        subtracted_image.ravel(),
        nan_policy="omit"
    )


    # ============================================================
    # ORIGINAL IMAGE
    # ============================================================

    plt.figure(
        figsize=(6, 5)
    )

    plt.imshow(
        image * 1000,
        origin="lower",
        cmap="hot"
    )

    plt.colorbar(
        label="Flux (mJy/beam)"
    )

    plt.scatter(
        source_x,
        source_y,
        marker="+",
        s=150,
        linewidths=2,
        label="Catalog position"
    )

    plt.xlabel(
        "X pixel"
    )

    plt.ylabel(
        "Y pixel"
    )

    plt.title(
        f"RACS Original Image — Index {INDEX}"
    )

    plt.legend()
    plt.tight_layout()
    plt.show()


    # ============================================================
    # PRIMARY-SUBTRACTED IMAGE
    # ============================================================

    plt.figure(
        figsize=(6, 5)
    )

    plt.imshow(
        subtracted_image * 1000,
        origin="lower",
        cmap="hot"
    )

    plt.colorbar(
        label="Flux (mJy/beam)"
    )

    plt.xlabel(
        "X pixel"
    )

    plt.ylabel(
        "Y pixel"
    )

    plt.title(
        f"RACS Primary Removed — Index {INDEX}"
    )

    plt.tight_layout()
    plt.show()


    # ============================================================
    # HISTOGRAM OF PRIMARY-SUBTRACTED PIXELS
    # ============================================================

    plt.figure(
        figsize=(7, 5)
    )

    plt.hist(
        subtracted_image.ravel() * 1000,
        bins=40,
        edgecolor="black"
    )

    plt.axvline(
        0,
        linestyle="--",
        linewidth=1
    )

    plt.xlabel(
        "Residual Flux (mJy/beam)"
    )

    plt.ylabel(
        "Number of Pixels"
    )

    plt.title(
        f"RACS Residual Pixel Distribution — Index {INDEX}\n"
        f"Skewness = {hist_skew:.3f}"
    )

    plt.tight_layout()
    plt.show()


    # ============================================================
    # DIAGNOSTIC INFORMATION
    # ============================================================

    sigma_maj = (
        bmaj /
        (
            2.35482 *
            pixscale_x
        )
    )

    sigma_min = (
        bmin /
        (
            2.35482 *
            pixscale_y
        )
    )

    print()
    print("DIAGNOSTIC IMAGE")
    print("--------------------------------")

    print(
        f"Index            = {INDEX}"
    )

    print(
        f"Catalog position = "
        f"({source_x:.3f}, {source_y:.3f})"
    )

    print(
        f"Beam             = "
        f"{bmaj*3600:.3f} x "
        f"{bmin*3600:.3f} arcsec FWHM"
    )

    print(
        f"BPA              = "
        f"{bpa:.3f} deg"
    )

    print(
        f"Sigma major/minor= "
        f"{sigma_maj:.3f}, "
        f"{sigma_min:.3f} pixels"
    )

    print(
        f"Fit radius       = "
        f"{FIT_RADIUS:.3f} pixels"
    )

    print(
        f"Peak flux        = "
        f"{peak_flux*1000:.3f} mJy/beam"
    )

    print(
        f"Background       = "
        f"{background*1000:.3f} mJy/beam"
    )

    print(
        f"Skewness         = "
        f"{hist_skew:.3f}"
    )

    print(
        f"Neighbor         = "
        f"{int(hist_skew > S)}"
    )

# ============================================================
# FLAGGED NEIGHBORS
# ============================================================

print()
print("FLAGGED NEIGHBORS")
print("--------------------------------")

flagged = df.loc[
    df["NEIGHBOR_DETECTION"] == 1,
    ["SOURCE_ID"]
]

# for index, row in flagged.iterrows():

    # print(
    #     "INDEX =",
    #     index,
    #     "SOURCE_ID =",
    #     row["SOURCE_ID"]
    # )