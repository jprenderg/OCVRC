# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 09:27:21 2026

@author: Joe
"""

# RACS PSF Fit.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# SETTINGS
# ============================================================
from pathlib import Path
import configparser
ROOT = Path(__file__).resolve().parents[1]
CFG = configparser.ConfigParser(interpolation=None)
if not CFG.read(ROOT / 'config.txt', encoding='utf-8'):
    raise FileNotFoundError(ROOT / 'config.txt')
image_dir = Path(CFG.get('PATHS','IMAGE_OUTPUT_DIR',fallback='output/images'))
image_dir = image_dir if image_dir.is_absolute() else ROOT / image_dir
image_dir.mkdir(parents=True, exist_ok=True)
path_in = image_dir / 'RACS_neighbor_detection.pkl'
path_out = image_dir / 'RACS_flux_measurements.pkl'
correlation_path = ROOT / 'calibration' / 'RACS_GLS_Correlation.npz'
FIT_RADIUS = 13.0
NOISE_MASK_RADIUS = 15.0
MAX_CORR_LAG = 15
EIGENVALUE_FLOOR = 0.01
INDEX_TO_PLOT = -1
PEAK_MIN, PEAK_MAX = -1.0, 3.0
RMS_MIN, RMS_MAX = 0.0, 0.6
SN_MIN, SN_MAX = -5.0, 15.0
N_BINS = 50
if not correlation_path.exists():
    raise FileNotFoundError(f'Missing precomputed GLS calibration: {correlation_path}. '
                            'Copy calibration from full OCVRC run; do not estimate from a few user targets.')


# ============================================================
# READ DATA
# ============================================================

df = pd.read_pickle(path_in).copy()

# Keep all images, including NEIGHBOR_DETECTION == 1.
# The flag remains in the dataframe for later filtering.

print("Number of RACS images:", len(df))


# ============================================================
# ELLIPTICAL RACS PSF
#
# BMAJ/BMIN     = degrees, FWHM
# PIXSCALE_X/Y  = degrees/pixel
#
# model = a*G + b
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
# VALIDATE BEAM PARAMETERS
# ============================================================

def validate_beam(row):

    values = np.array([
        row["BMAJ"],
        row["BMIN"],
        row["BPA"],
        row["PIXSCALE_X"],
        row["PIXSCALE_Y"]
    ], dtype=float)

    if not np.all(np.isfinite(values)):
        raise ValueError(
            "Non-finite beam or pixel-scale parameter"
        )

    if row["BMAJ"] <= 0:
        raise ValueError("Invalid BMAJ")

    if row["BMIN"] <= 0:
        raise ValueError("Invalid BMIN")

    if row["PIXSCALE_X"] <= 0:
        raise ValueError("Invalid PIXSCALE_X")

    if row["PIXSCALE_Y"] <= 0:
        raise ValueError("Invalid PIXSCALE_Y")


# ============================================================
# ORDINARY FIT
#
# Used only to estimate global spatial correlation.
#
# Position fixed at catalog/WCS coordinates.
#
# model = a*G + b
# ============================================================

def ordinary_fit(row):

    validate_beam(row)

    image = np.asarray(
        row["IMAGE"],
        dtype=float
    )

    x0 = float(row["SOURCE_X"])
    y0 = float(row["SOURCE_Y"])

    if (
        not np.isfinite(x0)
        or
        not np.isfinite(y0)
    ):
        raise ValueError(
            "Non-finite source position"
        )

    y, x = np.indices(
        image.shape
    )

    g = gaussian(
        x, y,
        x0, y0,
        row["BMAJ"],
        row["BMIN"],
        row["BPA"],
        row["PIXSCALE_X"],
        row["PIXSCALE_Y"]
    )

    r = np.hypot(
        x - x0,
        y - y0
    )

    fit_mask = (
        (r <= FIT_RADIUS) &
        np.isfinite(image) &
        np.isfinite(g)
    )

    if fit_mask.sum() < 3:
        raise ValueError(
            f"Too few finite fit pixels: "
            f"{fit_mask.sum()}"
        )

    X = np.column_stack((
        g[fit_mask],
        np.ones(
            fit_mask.sum()
        )
    ))

    z = image[
        fit_mask
    ]

    beta = np.linalg.lstsq(
        X,
        z,
        rcond=None
    )[0]

    residual = image - (
        beta[0] * g +
        beta[1]
    )

    return residual, x, y


# ============================================================
# DIAGNOSTIC FUNCTION
# ============================================================

def print_row_diagnostic(
    i,
    row,
    error
):

    image = np.asarray(
        row["IMAGE"],
        dtype=float
    )

    print()
    print("DIAGNOSTIC")
    print("--------------------------------")
    print("Index          =", i)
    print("SOURCE_ID      =", row.get("SOURCE_ID"))
    print("Error          =", error)
    print("Image shape    =", image.shape)
    print("NaN pixels     =", np.isnan(image).sum())
    print("Inf pixels     =", np.isinf(image).sum())
    print("Finite pixels  =", np.isfinite(image).sum())

    print(
        "SOURCE_X/Y     =",
        row.get("SOURCE_X"),
        row.get("SOURCE_Y")
    )

    print(
        "BMAJ/BMIN      =",
        row.get("BMAJ"),
        row.get("BMIN")
    )

    print(
        "BPA            =",
        row.get("BPA")
    )

    print(
        "PIXSCALE_X/Y   =",
        row.get("PIXSCALE_X"),
        row.get("PIXSCALE_Y")
    )

    print()


# ============================================================
# ESTIMATE GLOBAL RACS SPATIAL CORRELATION
#
# Match VLASS: estimate global noise correlation from all images,
# including those flagged NEIGHBOR_DETECTION == 1.
# The neighbor flag remains available for downstream filtering.
# Correlation calculation is vectorized for the large RACS cutouts.
# ============================================================

with np.load(correlation_path, allow_pickle=False) as cal:
    correlation = np.asarray(cal['correlation'], dtype=float)
    lag = int(cal['max_corr_lag'])
if lag != MAX_CORR_LAG or correlation.shape != (2*lag+1, 2*lag+1):
    raise ValueError(f'RACS calibration lag/shape mismatch: {lag}, {correlation.shape}')
if not np.isfinite(correlation).all():
    raise ValueError('RACS calibration contains nonfinite correlations')
df_corr = df
correlation_failures = []
print('Loaded precomputed RACS correlation:', correlation_path)


# ============================================================
# CORRELATION MATRIX CACHE
#
# For RACS the fit contains roughly 530 pixels. Repeatedly
# reconstructing and inverting a ~530x530 matrix is expensive.
#
# The matrix depends only on which integer pixels are included
# in the radius-13 fit. Cache matrices with identical geometry.
# ============================================================

correlation_cache = {}

# Center index of the spatial correlation matrix
center = MAX_CORR_LAG


def make_correlation_matrix(
    coords
):

    key = tuple(
        map(
            tuple,
            coords.astype(int)
        )
    )

    if key in correlation_cache:

        return correlation_cache[
            key
        ]


    dx = (
        coords[:, None, 0] -
        coords[None, :, 0]
    ).astype(int)

    dy = (
        coords[:, None, 1] -
        coords[None, :, 1]
    ).astype(int)


    R = np.zeros(
        dx.shape,
        dtype=float
    )

    valid = (
        (np.abs(dx) <= MAX_CORR_LAG) &
        (np.abs(dy) <= MAX_CORR_LAG)
    )

    R[
        valid
    ] = correlation[
        dy[valid] + center,
        dx[valid] + center
    ]


    # Force exact numerical symmetry
    R = 0.5 * (
        R +
        R.T
    )


    # --------------------------------------------------------
    # EIGENVALUE REGULARIZATION
    # --------------------------------------------------------

    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            R
        )
    )

    eigenvalue_binding = int(
        np.any(
            eigenvalues <
            EIGENVALUE_FLOOR
        )
    )

    eigenvalues = np.maximum(
        eigenvalues,
        EIGENVALUE_FLOOR
    )

    R = (
        eigenvectors @
        np.diag(
            eigenvalues
        ) @
        eigenvectors.T
    )


    # Renormalize diagonal to 1
    d = np.sqrt(
        np.diag(
            R
        )
    )

    R = (
        R /
        np.outer(
            d,
            d
        )
    )


    # Inverse cached because geometry is reused
    R_inv = np.linalg.inv(
        R
    )

    result = (
        R_inv,
        eigenvalue_binding
    )

    correlation_cache[
        key
    ] = result

    return result


# ============================================================
# GLS FIT
#
# model = a*G + b
#
# Position fixed at SOURCE_X, SOURCE_Y.
#
# PEAK_FLUX  = a
# RMS        = standard error of a
# BACKGROUND = b
# SN         = a / RMS
# ============================================================

for col in [
    "PEAK_FLUX",
    "RMS",
    "BACKGROUND"
]:

    df[col] = np.nan


df[
    "EIGENVALUE_BINDING"
] = np.nan


print()
print(
    "Performing RACS GLS fits..."
)
print()


gls_failures = []


for i, row in df.iterrows():

    try:

        validate_beam(
            row
        )

        image = np.asarray(
            row["IMAGE"],
            dtype=float
        )

        x0 = float(
            row["SOURCE_X"]
        )

        y0 = float(
            row["SOURCE_Y"]
        )

        if (
            not np.isfinite(x0)
            or
            not np.isfinite(y0)
        ):

            raise ValueError(
                "Non-finite source position"
            )


        y, x = np.indices(
            image.shape
        )


        g = gaussian(
            x, y,
            x0, y0,
            row["BMAJ"],
            row["BMIN"],
            row["BPA"],
            row["PIXSCALE_X"],
            row["PIXSCALE_Y"]
        )


        r = np.hypot(
            x - x0,
            y - y0
        )


        fit_mask = (
            (r <= FIT_RADIUS) &
            np.isfinite(image) &
            np.isfinite(g)
        )


        noise_mask = (
            (r > NOISE_MASK_RADIUS) &
            np.isfinite(image)
        )


        if fit_mask.sum() < 3:

            raise ValueError(
                f"Too few finite fit pixels: "
                f"{fit_mask.sum()}"
            )


        if noise_mask.sum() < 3:

            raise ValueError(
                f"Too few finite noise pixels: "
                f"{noise_mask.sum()}"
            )


        z = image[
            fit_mask
        ].ravel()


        X = np.column_stack((
            g[
                fit_mask
            ].ravel(),

            np.ones(
                fit_mask.sum()
            )
        ))


        coords = np.column_stack((
            x[
                fit_mask
            ].ravel(),

            y[
                fit_mask
            ].ravel()
        ))


        if not np.all(
            np.isfinite(X)
        ):

            raise ValueError(
                "Non-finite design matrix"
            )


        if not np.all(
            np.isfinite(z)
        ):

            raise ValueError(
                "Non-finite fit data"
            )


        # ----------------------------------------------------
        # LOCAL IMAGE NOISE
        #
        # Used only to scale the covariance matrix.
        # ----------------------------------------------------

        noise_std = np.std(
            image[
                noise_mask
            ],
            ddof=1
        )

        if (
            not np.isfinite(noise_std)
            or
            noise_std <= 0
        ):

            raise ValueError(
                "Invalid local noise standard deviation"
            )


        # ----------------------------------------------------
        # SPATIAL CORRELATION MATRIX
        # ----------------------------------------------------

        (
            R_inv,
            eigenvalue_binding
        ) = make_correlation_matrix(
            coords
        )


        df.at[
            i,
            "EIGENVALUE_BINDING"
        ] = eigenvalue_binding


        # ----------------------------------------------------
        # GLS SOLUTION
        # ----------------------------------------------------

        XtRiX = (
            X.T @
            R_inv @
            X
        )


        beta = np.linalg.solve(
            XtRiX,
            X.T @
            R_inv @
            z
        )


        peak_flux, background = (
            beta
        )


        # ----------------------------------------------------
        # COVARIANCE OF FITTED PARAMETERS
        #
        # RMS = standard error of a
        # ----------------------------------------------------

        covariance = (
            noise_std**2 *
            np.linalg.inv(
                XtRiX
            )
        )


        rms = np.sqrt(
            max(
                covariance[
                    0,
                    0
                ],
                0
            )
        )


        # ----------------------------------------------------
        # Jy/beam -> mJy/beam
        # ----------------------------------------------------

        (
            peak_flux,
            rms,
            background
        ) = (
            1000 *
            np.array([
                peak_flux,
                rms,
                background
            ])
        )


        # ----------------------------------------------------
        # STORE
        # ----------------------------------------------------

        df.loc[
            i,
            [
                "PEAK_FLUX",
                "RMS",
                "BACKGROUND"
            ]
        ] = (
            peak_flux,
            rms,
            background
        )


    except Exception as error:

        gls_failures.append(
            (
                i,
                row.get(
                    "SOURCE_ID"
                ),
                str(error)
            )
        )

        print()
        print(
            "GLS exception:",
            i,
            error
        )

        print_row_diagnostic(
            i,
            row,
            error
        )


# ============================================================
# SIGNAL-TO-NOISE
# ============================================================

df["SN"] = (
    df["PEAK_FLUX"] /
    df["RMS"]
)


# ============================================================
# EIGENVALUE SUMMARY
# ============================================================

n_binding = int(
    (
        df[
            "EIGENVALUE_BINDING"
        ] == 1
    ).sum()
)

n_not_binding = int(
    (
        df[
            "EIGENVALUE_BINDING"
        ] == 0
    ).sum()
)

n_tested = int(
    df[
        "EIGENVALUE_BINDING"
    ]
    .notna()
    .sum()
)


# ============================================================
# MEASUREMENT TABLE OUTPUT
# ============================================================

df_out = df[
    [
        "SOURCE_ID",
        "PEAK_FLUX",
        "RMS",
        "NEIGHBOR_DETECTION"
    ]
].copy()


df_out[
    "OBSERVATORY"
] = "ASKAP"

df_out[
    "EPOCH"
] = "RACS"

df_out[
    "BAND"
] = "887.5 MHz"

df_out[
    "UNITS"
] = "mJy"


df_out = df_out.rename(
    columns={
        "PEAK_FLUX":
            "MEASUREMENT",
        "RMS":
            "ERROR"
    }
)


df_out = df_out[
    [
        "SOURCE_ID",
        "OBSERVATORY",
        "EPOCH",
        "BAND",
        "UNITS",
        "MEASUREMENT",
        "ERROR",
        "NEIGHBOR_DETECTION"
    ]
]


df_out = df_out.sort_values("SOURCE_ID", ascending=True).reset_index(drop=True)

df_out.to_pickle(
    path_out
)


print(f"Saved {len(df_out)} fitted observations to {path_out}")
print(f"Neighbor-flagged observations: {(df_out['NEIGHBOR_DETECTION'] == 1).sum()}")
