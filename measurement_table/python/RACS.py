# -*- coding: utf-8 -*-
"""
Created on Wed May 27 12:01:55 2026

@author: Joe
"""

import os
import time
import numpy as np
import pandas as pd

from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import wcs

from astroquery.casda import Casda

from config import CACHE_DIR, OPAL_USERNAME


# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

cache_dir = CACHE_DIR
os.makedirs(cache_dir, exist_ok=True)

MAX_CUTOUT_TRIES = 5
MAX_DOWNLOAD_TRIES = 3

CUTOUT_RETRY_DELAY = 10
DOWNLOAD_RETRY_DELAY = 5


# ------------------------------------------------------------
# CASDA login
# ------------------------------------------------------------

casda = Casda()
casda.login(username=OPAL_USERNAME)


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def make_output(source_id, measurement=None, error=None):

    return pd.DataFrame([{
        "SOURCE_ID": source_id,
        "OBSERVATORY": "ASKAP",
        "CAT": "RACS-DR1",
        "BAND": "887.5 MHz",
        "UNITS": "mJy",
        "MEASUREMENT": measurement,
        "ERROR": error
    }])


def is_valid_fits(path):

    try:

        with open(path, "rb") as f:
            return f.read(8).startswith(b"SIMPLE")

    except Exception:
        return False


# ------------------------------------------------------------
# Download FITS file with retries
# ------------------------------------------------------------

def download_valid_fits(
    casda,
    urls,
    savedir,
    label,
    max_tries=MAX_DOWNLOAD_TRIES
):

    for attempt in range(1, max_tries + 1):

        print(
            f"RACS {label} download attempt "
            f"{attempt}/{max_tries}"
        )

        try:

            files = casda.download_files(
                urls,
                savedir=savedir
            )

            fits_files = [
                str(f)
                for f in files
                if str(f).lower().endswith(".fits")
            ]

            for path in fits_files:

                if is_valid_fits(path):

                    print(f"RACS {label} downloaded:")
                    print(f"  {path}")

                    return path

                # Invalid FITS file -- remove it
                try:
                    os.remove(path)
                except Exception:
                    pass

        except Exception as e:

            print(
                f"RACS {label} download attempt "
                f"{attempt} failed: {e}"
            )

        if attempt < max_tries:
            time.sleep(DOWNLOAD_RETRY_DELAY)

    return None


# ------------------------------------------------------------
# Request CASDA cutout with retries
# ------------------------------------------------------------

def get_cutout_with_retries(
    selected,
    coord,
    radius,
    label,
    max_tries=MAX_CUTOUT_TRIES
):

    for attempt in range(1, max_tries + 1):

        print(
            f"RACS {label} cutout attempt "
            f"{attempt}/{max_tries}"
        )

        try:

            urls = casda.cutout(
                selected,
                coordinates=coord,
                radius=radius
            )

            if urls is not None and len(urls) > 0:
                return urls

            print(
                f"RACS {label} cutout returned no URLs."
            )

        except Exception as e:

            print(
                f"RACS {label} cutout attempt "
                f"{attempt} failed: {e}"
            )

        if attempt < max_tries:
            time.sleep(CUTOUT_RETRY_DELAY)

    return None


# ------------------------------------------------------------
# Get a CASDA cutout and download it
# ------------------------------------------------------------

def get_racs_fits(
    selected,
    coord,
    radius,
    label
):

    # --------------------------------------------------------
    # Stage/request cutout
    # --------------------------------------------------------

    urls = get_cutout_with_retries(
        selected=selected,
        coord=coord,
        radius=radius,
        label=label
    )

    if urls is None:

        print(
            f"RACS {label} cutout failed after "
            f"{MAX_CUTOUT_TRIES} attempts."
        )

        return None

    # --------------------------------------------------------
    # Download cutout
    # --------------------------------------------------------

    path = download_valid_fits(
        casda=casda,
        urls=urls,
        savedir=cache_dir,
        label=label
    )

    if path is None:

        print(
            f"RACS {label} download failed after "
            f"{MAX_DOWNLOAD_TRIES} attempts."
        )

        return None

    return path


# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------

def racs(source_id, ra_deg, dec_deg):

    print("\nStart racs\n")

    search_radius = 30 * u.arcmin
    cutout_radius = 1 * u.arcmin

    dim = 5

    measurement = None
    error = None

    coord = SkyCoord(
        ra=ra_deg * u.deg,
        dec=dec_deg * u.deg,
        frame="icrs"
    )

    try:

        # ----------------------------------------------------
        # Query CASDA near coordinate
        # ----------------------------------------------------

        result = casda.query_region(
            coord,
            radius=search_radius
        )

        public_data = casda.filter_out_unreleased(
            result
        )

        filenames = public_data[
            "filename"
        ].astype(str)

        # ----------------------------------------------------
        # Select RACS DR1 flux image
        #
        # Example:
        # RACS-DR1_0202-37A.fits
        # ----------------------------------------------------

        flux_subset = public_data[
            (
                public_data["obs_collection"]
                == "The Rapid ASKAP Continuum Survey"
            )
            & np.char.startswith(
                filenames,
                "RACS-DR1_"
            )
            & np.char.endswith(
                filenames,
                "A.fits"
            )
            & ~np.char.endswith(
                filenames,
                "_RMS.fits"
            )
        ]

        if len(flux_subset) == 0:

            print(
                "No RACS DR1 flux FITS image found."
            )

            return make_output(source_id)

        flux_selected = flux_subset[:1]

        flux_filename = str(
            flux_selected["filename"][0]
        )

        # ----------------------------------------------------
        # Find corresponding RMS image
        #
        # Example:
        #
        # RACS-DR1_0202-37A.fits
        # RACS-DR1_0202-37A_RMS.fits
        # ----------------------------------------------------

        rms_filename = flux_filename.replace(
            ".fits",
            "_RMS.fits"
        )

        rms_subset = public_data[
            public_data[
                "filename"
            ].astype(str) == rms_filename
        ]

        if len(rms_subset) == 0:

            print(
                f"No matching RMS FITS image found: "
                f"{rms_filename}"
            )

            return make_output(source_id)

        rms_selected = rms_subset[:1]

        print("Flux file:", flux_filename)
        print("RMS file: ", rms_filename)

        # ----------------------------------------------------
        # Get flux cutout
        # ----------------------------------------------------

        flux_path = get_racs_fits(
            selected=flux_selected,
            coord=coord,
            radius=cutout_radius,
            label="flux"
        )

        if flux_path is None:

            print(
                "No valid flux FITS file downloaded."
            )

            return make_output(source_id)

        # ----------------------------------------------------
        # Get RMS cutout
        # ----------------------------------------------------

        rms_path = get_racs_fits(
            selected=rms_selected,
            coord=coord,
            radius=cutout_radius,
            label="RMS"
        )

        if rms_path is None:

            print(
                "No valid RMS FITS file downloaded."
            )

            return make_output(source_id)

        # ----------------------------------------------------
        # Measure flux from flux image
        # ----------------------------------------------------

        with fits.open(
            flux_path,
            memmap=False
        ) as hdul:

            flux_data = hdul[0].data
            flux_header = hdul[0].header.copy()

            flux_image = np.squeeze(
                flux_data
            )

            flux_wcs = wcs.WCS(
                flux_header
            ).celestial

            pix = flux_wcs.wcs_world2pix(
                [(ra_deg, dec_deg)],
                0
            )[0]

            xpix = int(round(pix[0]))
            ypix = int(round(pix[1]))

            nrows, ncols = flux_image.shape

            y0 = ypix - dim
            y1 = ypix + dim
            x0 = xpix - dim
            x1 = xpix + dim

            if not (
                0 <= xpix < ncols
                and 0 <= ypix < nrows
                and x0 >= 0
                and y0 >= 0
                and x1 <= ncols
                and y1 <= nrows
            ):

                print(
                    "Flux cutout lies outside image."
                )

                return make_output(source_id)

            flux_cutout = flux_image[
                y0:y1,
                x0:x1
            ]

            if np.all(
                np.isnan(flux_cutout)
            ):

                print(
                    "Flux cutout contains only NaN values."
                )

                return make_output(source_id)

            max_y_cutout, max_x_cutout = (
                np.unravel_index(
                    np.nanargmax(flux_cutout),
                    flux_cutout.shape
                )
            )

            y_flux = (
                y0 + max_y_cutout
            )

            x_flux = (
                x0 + max_x_cutout
            )

            measurement = float(
                flux_image[
                    y_flux,
                    x_flux
                ]
                * 1000.0
            )

            flux_peak_sky = (
                flux_wcs.wcs_pix2world(
                    [(x_flux, y_flux)],
                    0
                )[0]
            )

            peak_ra = flux_peak_sky[0]
            peak_dec = flux_peak_sky[1]

        # ----------------------------------------------------
        # Measure RMS at same sky position as flux peak
        # ----------------------------------------------------

        with fits.open(
            rms_path,
            memmap=False
        ) as hdul:

            rms_data = hdul[0].data
            rms_header = hdul[0].header.copy()

            rms_image = np.squeeze(
                rms_data
            )

            rms_wcs = wcs.WCS(
                rms_header
            ).celestial

            rms_pix = rms_wcs.wcs_world2pix(
                [(peak_ra, peak_dec)],
                0
            )[0]

            x_rms = int(
                round(rms_pix[0])
            )

            y_rms = int(
                round(rms_pix[1])
            )

            rms_rows, rms_cols = (
                rms_image.shape
            )

            if (
                0 <= x_rms < rms_cols
                and
                0 <= y_rms < rms_rows
            ):

                error = float(
                    rms_image[
                        y_rms,
                        x_rms
                    ]
                    * 1000.0
                )

            else:

                print(
                    "RMS coordinate lies outside "
                    "RMS image."
                )

                return make_output(source_id)

    except Exception as e:

        print(
            f"ASKAP exception for SOURCE_ID "
            f"{source_id}: {e}"
        )

        return make_output(source_id)

    return make_output(
        source_id,
        measurement=measurement,
        error=error
    )


# ------------------------------------------------------------
# Test
# ------------------------------------------------------------

# unique_id = 1
# ra_deg = 0.0370542
# dec_deg = -77.3388
#
# temp = racs(
#     unique_id,
#     ra_deg,
#     dec_deg
# )
#
# print(temp)