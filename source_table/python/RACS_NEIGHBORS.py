# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 11:55:16 2026

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

CACHE_DIR.mkdir(
    parents=True,
    exist_ok=True
)

cache_dir = str(CACHE_DIR)


# ------------------------------------------------------------
# CASDA settings
# ------------------------------------------------------------

search_radius = 30 * u.arcmin
cutout_radius = 2 * u.arcmin

max_cutout_tries = 5
cutout_retry_wait = 10

max_download_tries = 3
download_retry_wait = 5


# ------------------------------------------------------------
# RACS neighbor-search settings
# ------------------------------------------------------------

search_radius_pix = 120

peak_box_size = 5
half_box = peak_box_size // 2

threshold_mjy = 1.0
threshold_jy = threshold_mjy / 1000.0

# Minimum angular separation from primary source
min_neighbor_sep_arcsec = 2.0


# ------------------------------------------------------------
# CASDA login
# ------------------------------------------------------------

casda = Casda()

casda.login(
    username=OPAL_USERNAME
)


# ------------------------------------------------------------
# Helper: convert FITS data to 2-D image
# ------------------------------------------------------------

def get_2d_image(data):

    image = np.squeeze(data)

    if image.ndim == 2:
        return image

    return None


# ------------------------------------------------------------
# Helper: verify downloaded file is a FITS file
# ------------------------------------------------------------

def is_valid_fits(path):

    try:

        with open(path, "rb") as f:
            return f.read(8).startswith(b"SIMPLE")

    except Exception:
        return False


# ------------------------------------------------------------
# Helper: retry CASDA cutout/staging operation
# ------------------------------------------------------------

def casda_cutout_with_retry(
    casda,
    selected,
    coord,
    radius,
    label="CASDA file",
    max_tries=max_cutout_tries,
    wait_seconds=cutout_retry_wait,
):

    for attempt in range(
        1,
        max_tries + 1
    ):

        try:

            print(
                f"{label} cutout attempt "
                f"{attempt}/{max_tries}"
            )

            urls = casda.cutout(
                selected,
                coordinates=coord,
                radius=radius
            )

            return urls

        except Exception as e:

            print(
                f"{label} cutout attempt "
                f"{attempt}/{max_tries} failed:"
            )

            print(
                f"  {e}"
            )

            if attempt < max_tries:

                print(
                    f"Retrying in "
                    f"{wait_seconds} seconds..."
                )

                time.sleep(
                    wait_seconds
                )

    print()

    print(
        f"WARNING: {label} cutout failed "
        f"after {max_tries} attempts."
    )

    print(
        "This is a CASDA staging/download failure, "
        "not evidence that no RACS source exists."
    )

    print()

    return None


# ------------------------------------------------------------
# Helper: download and validate FITS file
# ------------------------------------------------------------

def download_valid_fits(
    casda,
    urls,
    savedir,
    label="CASDA file",
    max_tries=max_download_tries,
    wait_seconds=download_retry_wait,
):

    if urls is None:
        return None

    for attempt in range(
        1,
        max_tries + 1
    ):

        try:

            print(
                f"{label} download attempt "
                f"{attempt}/{max_tries}"
            )

            files = casda.download_files(
                urls,
                savedir=savedir
            )

        except Exception as e:

            print(
                f"{label} download attempt "
                f"{attempt}/{max_tries} failed:"
            )

            print(
                f"  {e}"
            )

            if attempt < max_tries:

                print(
                    f"Retrying in "
                    f"{wait_seconds} seconds..."
                )

                time.sleep(
                    wait_seconds
                )

            continue

        fits_files = [
            f
            for f in files
            if str(f).lower().endswith(
                ".fits"
            )
        ]

        for path in fits_files:

            if is_valid_fits(path):

                print(
                    f"{label} downloaded:"
                )

                print(
                    f"  {path}"
                )

                return path

            print(
                f"Downloaded {label} is not "
                f"a valid FITS file:"
            )

            print(
                f"  {path}"
            )

            try:
                os.remove(path)

            except Exception:
                pass

        if attempt < max_tries:

            print(
                f"Retrying download in "
                f"{wait_seconds} seconds..."
            )

            time.sleep(
                wait_seconds
            )

    print()

    print(
        f"WARNING: Unable to download "
        f"a valid {label}."
    )

    print()

    return None


# ------------------------------------------------------------
# Get RACS flux and RMS cutouts
# ------------------------------------------------------------

def get_racs_cutouts(
    ra_deg,
    dec_deg
):
    """
    Query CASDA for RACS-DR1 data.

    Finds matching flux and RMS images,
    stages cutouts through CASDA,
    downloads them,
    verifies that they are valid FITS files,
    and returns their local paths.

    Returns
    -------
    flux_path
        Local path to RACS flux FITS cutout.

    rms_path
        Local path to RACS RMS FITS cutout.

    Returns (None, None) if the operation fails.
    """

    coord = SkyCoord(
        ra=ra_deg * u.deg,
        dec=dec_deg * u.deg,
        frame="icrs"
    )


    # --------------------------------------------------------
    # Query CASDA
    # --------------------------------------------------------

    try:

        result = casda.query_region(
            coord,
            radius=search_radius
        )

    except Exception as e:

        print()

        print(
            "WARNING: CASDA query failed:"
        )

        print(
            f"  {e}"
        )

        print()

        return None, None


    # --------------------------------------------------------
    # Keep only public data
    # --------------------------------------------------------

    try:

        public_data = (
            casda.filter_out_unreleased(
                result
            )
        )

    except Exception as e:

        print()

        print(
            "WARNING: Could not filter "
            "CASDA search results:"
        )

        print(
            f"  {e}"
        )

        print()

        return None, None


    # --------------------------------------------------------
    # Find RACS-DR1 flux image
    # --------------------------------------------------------

    filenames = (
        public_data["filename"]
        .astype(str)
    )

    flux_subset = public_data[
        (
            public_data["obs_collection"]
            ==
            "The Rapid ASKAP Continuum Survey"
        )
        &
        np.char.startswith(
            filenames,
            "RACS-DR1_"
        )
        &
        np.char.endswith(
            filenames,
            "A.fits"
        )
        &
        ~np.char.endswith(
            filenames,
            "_RMS.fits"
        )
    ]


    if len(flux_subset) == 0:

        print(
            "No RACS-DR1 flux image found."
        )

        return None, None


    flux_selected = flux_subset[:1]

    flux_filename = str(
        flux_selected["filename"][0]
    )


    # --------------------------------------------------------
    # Find matching RMS image
    # --------------------------------------------------------

    rms_filename = (
        flux_filename.replace(
            ".fits",
            "_RMS.fits"
        )
    )

    rms_subset = public_data[
        public_data["filename"]
        .astype(str)
        ==
        rms_filename
    ]


    if len(rms_subset) == 0:

        print(
            "Matching RACS RMS image "
            "was not found."
        )

        print(
            f"Expected RMS file: "
            f"{rms_filename}"
        )

        return None, None


    rms_selected = rms_subset[:1]


    print()

    print(
        "Flux file:",
        flux_filename
    )

    print(
        "RMS file: ",
        rms_filename
    )


    # --------------------------------------------------------
    # Stage flux cutout
    # --------------------------------------------------------

    flux_urls = casda_cutout_with_retry(
        casda,
        flux_selected,
        coord,
        cutout_radius,
        label="RACS flux",
    )


    if flux_urls is None:

        print(
            "RACS flux cutout could not "
            "be staged."
        )

        return None, None


    # --------------------------------------------------------
    # Download flux cutout
    # --------------------------------------------------------

    flux_path = download_valid_fits(
        casda,
        flux_urls,
        cache_dir,
        label="RACS flux",
    )


    if flux_path is None:

        print(
            "RACS flux cutout could not "
            "be downloaded."
        )

        return None, None


    # --------------------------------------------------------
    # Stage RMS cutout
    # --------------------------------------------------------

    rms_urls = casda_cutout_with_retry(
        casda,
        rms_selected,
        coord,
        cutout_radius,
        label="RACS RMS",
    )


    if rms_urls is None:

        print(
            "RACS RMS cutout could not "
            "be staged."
        )

        return None, None


    # --------------------------------------------------------
    # Download RMS cutout
    # --------------------------------------------------------

    rms_path = download_valid_fits(
        casda,
        rms_urls,
        cache_dir,
        label="RACS RMS",
    )


    if rms_path is None:

        print(
            "RACS RMS cutout could not "
            "be downloaded."
        )

        return None, None


    return (
        flux_path,
        rms_path
    )


# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------

def racs_neighbors(
    ra_deg,
    dec_deg
):
    """
    Find RACS local-peak neighbors around input RA/DEC.

    A detected source is only considered a neighbor if
    its angular separation from the primary source is
    greater than 2 arcsec.

    Downloads RACS-DR1 flux and RMS cutouts from CASDA.

    Searches for local flux-density peaks above the
    configured threshold within search_radius_pix of
    the input position.

    Parameters
    ----------
    ra_deg : float
        RA in degrees.

    dec_deg : float
        DEC in degrees.

    Returns
    -------
    df_neighbors : pandas.DataFrame
        Columns:
            RA
            DEC
    """

    print("\nStart racs_neighbors\n")

    empty = pd.DataFrame(
        columns=[
            "RA",
            "DEC"
        ]
    )


    # --------------------------------------------------------
    # Validate input coordinates
    # --------------------------------------------------------

    if (
        pd.isna(ra_deg)
        or
        pd.isna(dec_deg)
    ):
        return empty


    if (
        not np.isfinite(ra_deg)
        or
        not np.isfinite(dec_deg)
    ):
        return empty


    # --------------------------------------------------------
    # Define primary-source coordinate
    # --------------------------------------------------------

    primary_coord = SkyCoord(
        ra=ra_deg * u.deg,
        dec=dec_deg * u.deg,
        frame="icrs"
    )


    # --------------------------------------------------------
    # Download RACS cutouts
    # --------------------------------------------------------

    try:

        (
            fits_path_flux,
            fits_path_rms
        ) = get_racs_cutouts(
            ra_deg,
            dec_deg
        )

    except Exception as e:

        print()

        print(
            "WARNING: Unexpected CASDA "
            "exception:"
        )

        print(
            f"  {e}"
        )

        print(
            "RACS neighbor search could "
            "not be completed."
        )

        print()

        return empty


    if (
        fits_path_flux is None
        or
        fits_path_rms is None
    ):

        print(
            "RACS neighbor search skipped "
            "because required FITS data "
            "could not be obtained."
        )

        return empty


    # --------------------------------------------------------
    # Open flux FITS
    # --------------------------------------------------------

    try:

        with fits.open(
            fits_path_flux,
            memmap=False
        ) as hdul:

            data = (
                hdul[0].data
            )

            header = (
                hdul[0].header.copy()
            )

            image = get_2d_image(
                data
            )


            if image is None:

                print(
                    "RACS flux FITS does not "
                    "contain a usable 2-D image."
                )

                return empty


    except Exception as e:

        print(
            "Could not open flux FITS:"
        )

        print(
            f"  {e}"
        )

        return empty


    # --------------------------------------------------------
    # Open RMS FITS
    # --------------------------------------------------------

    try:

        with fits.open(
            fits_path_rms,
            memmap=False
        ) as hdul_rms:

            data_rms = (
                hdul_rms[0].data
            )

            image_rms = get_2d_image(
                data_rms
            )


            if image_rms is None:

                print(
                    "RACS RMS FITS does not "
                    "contain a usable 2-D image."
                )

                return empty


    except Exception as e:

        print(
            "Could not open RMS FITS:"
        )

        print(
            f"  {e}"
        )

        return empty


    # --------------------------------------------------------
    # Flux and RMS images must have same dimensions
    # --------------------------------------------------------

    if (
        image.shape
        !=
        image_rms.shape
    ):

        print(
            "Flux and RMS FITS image "
            "dimensions do not match."
        )

        return empty


    image_rows, image_columns = (
        image.shape
    )


    # --------------------------------------------------------
    # Build celestial WCS
    # --------------------------------------------------------

    try:

        mywcs = (
            wcs.WCS(header)
            .celestial
        )

    except Exception as e:

        print(
            "WCS exception:"
        )

        print(
            f"  {e}"
        )

        return empty


    # --------------------------------------------------------
    # Convert input RA/DEC to pixel
    # --------------------------------------------------------

    try:

        pix = mywcs.wcs_world2pix(
            [
                (
                    ra_deg,
                    dec_deg
                )
            ],
            0
        )[0]

    except Exception as e:

        print(
            "Could not convert RA/DEC "
            "to image pixel:"
        )

        print(
            f"  {e}"
        )

        return empty


    if (
        not np.isfinite(pix[0])
        or
        not np.isfinite(pix[1])
    ):

        return empty


    xpix = int(
        round(pix[0])
    )

    ypix = int(
        round(pix[1])
    )


    if (
        xpix < 0
        or
        xpix >= image_columns
        or
        ypix < 0
        or
        ypix >= image_rows
    ):

        print(
            "Requested coordinate falls "
            "outside the downloaded "
            "RACS cutout."
        )

        return empty


    # --------------------------------------------------------
    # Define search region
    # --------------------------------------------------------

    y0 = max(
        half_box,
        ypix - search_radius_pix
    )

    y1 = min(
        image_rows - half_box,
        ypix + search_radius_pix + 1
    )

    x0 = max(
        half_box,
        xpix - search_radius_pix
    )

    x1 = min(
        image_columns - half_box,
        xpix + search_radius_pix + 1
    )


    if (
        y1 <= y0
        or
        x1 <= x0
    ):
        return empty


    center = image[
        y0:y1,
        x0:x1
    ]


    if center.size == 0:
        return empty


    if not np.isfinite(
        center
    ).any():
        return empty


    yy, xx = np.mgrid[
        y0:y1,
        x0:x1
    ]


    inside_radius = (
        (xx - xpix) ** 2
        +
        (yy - ypix) ** 2
    ) <= (
        search_radius_pix ** 2
    )


    # --------------------------------------------------------
    # Find local peaks above threshold
    # --------------------------------------------------------

    peak_mask = (
        np.isfinite(center)
        &
        (center > threshold_jy)
        &
        inside_radius
    )


    for dy in range(
        -half_box,
        half_box + 1
    ):

        for dx in range(
            -half_box,
            half_box + 1
        ):

            if (
                dy == 0
                and
                dx == 0
            ):
                continue


            neighbor_y0 = (
                y0 + dy
            )

            neighbor_y1 = (
                y1 + dy
            )

            neighbor_x0 = (
                x0 + dx
            )

            neighbor_x1 = (
                x1 + dx
            )


            if (
                neighbor_y0 < 0
                or
                neighbor_x0 < 0
                or
                neighbor_y1 > image_rows
                or
                neighbor_x1 > image_columns
            ):

                peak_mask[:, :] = False

                continue


            neighbor = image[
                neighbor_y0:neighbor_y1,
                neighbor_x0:neighbor_x1
            ]


            if (
                neighbor.shape
                !=
                center.shape
            ):

                peak_mask[:, :] = False

                continue


            peak_mask &= (
                center > neighbor
            )


    peak_y = yy[
        peak_mask
    ]

    peak_x = xx[
        peak_mask
    ]


    # --------------------------------------------------------
    # No peaks found
    # --------------------------------------------------------

    if len(peak_x) == 0:

        print(
            "No RACS peaks above "
            f"{threshold_mjy:.1f} mJy "
            "were found in the search region."
        )

        return empty


    # --------------------------------------------------------
    # Bounds and edge checks
    # --------------------------------------------------------

    rms_rows, rms_columns = (
        image_rms.shape
    )


    good_peak = (
        np.isfinite(peak_x)
        &
        np.isfinite(peak_y)
        &
        (peak_x >= 0)
        &
        (peak_x < image_columns)
        &
        (peak_y >= 0)
        &
        (peak_y < image_rows)
        &
        (peak_x < rms_columns)
        &
        (peak_y < rms_rows)
        &
        (
            peak_x - half_box
            >= 0
        )
        &
        (
            peak_x + half_box
            < image_columns
        )
        &
        (
            peak_y - half_box
            >= 0
        )
        &
        (
            peak_y + half_box
            < image_rows
        )
        &
        (
            peak_x + half_box
            < rms_columns
        )
        &
        (
            peak_y + half_box
            < rms_rows
        )
    )


    peak_x = (
        peak_x[
            good_peak
        ]
        .astype(int)
    )

    peak_y = (
        peak_y[
            good_peak
        ]
        .astype(int)
    )


    if len(peak_x) == 0:
        return empty


    # --------------------------------------------------------
    # Require finite flux and RMS at peak pixels
    # --------------------------------------------------------

    peak_flux_mjy = (
        image[
            peak_y,
            peak_x
        ]
        *
        1000.0
    )

    peak_rms_mjy = (
        image_rms[
            peak_y,
            peak_x
        ]
        *
        1000.0
    )


    finite_flux_rms = (
        np.isfinite(
            peak_flux_mjy
        )
        &
        np.isfinite(
            peak_rms_mjy
        )
    )


    peak_x = peak_x[
        finite_flux_rms
    ]

    peak_y = peak_y[
        finite_flux_rms
    ]


    if len(peak_x) == 0:
        return empty


    # --------------------------------------------------------
    # Convert peak pixels back to RA/DEC
    # --------------------------------------------------------

    try:

        world = (
            mywcs.wcs_pix2world(
                np.column_stack(
                    [
                        peak_x,
                        peak_y
                    ]
                ),
                0
            )
        )

    except Exception as e:

        print(
            "Could not convert peak pixels "
            "to RA/DEC:"
        )

        print(
            f"  {e}"
        )

        return empty


    peak_ra = world[:, 0]
    peak_dec = world[:, 1]


    finite_world = (
        np.isfinite(peak_ra)
        &
        np.isfinite(peak_dec)
    )


    peak_ra = peak_ra[
        finite_world
    ]

    peak_dec = peak_dec[
        finite_world
    ]


    if len(peak_ra) == 0:
        return empty


    # --------------------------------------------------------
    # Calculate angular separation from primary
    # --------------------------------------------------------

    peak_coords = SkyCoord(
        ra=peak_ra * u.deg,
        dec=peak_dec * u.deg,
        frame="icrs"
    )

    separation_arcsec = (
        primary_coord
        .separation(peak_coords)
        .arcsec
    )


    # --------------------------------------------------------
    # Keep only true neighbors
    #
    # Neighbor must be MORE than 2 arcsec from primary.
    # --------------------------------------------------------

    neighbor_mask = (
        separation_arcsec
        >
        min_neighbor_sep_arcsec
    )


    peak_ra = peak_ra[
        neighbor_mask
    ]

    peak_dec = peak_dec[
        neighbor_mask
    ]


    # --------------------------------------------------------
    # No neighbors remain after 2-arcsec exclusion
    # --------------------------------------------------------

    if len(peak_ra) == 0:

        print(
            "RACS neighbors found: 0"
        )

        return empty


    # --------------------------------------------------------
    # Construct output dataframe
    # --------------------------------------------------------

    df_neighbors = pd.DataFrame(
        {
            "RA": peak_ra,
            "DEC": peak_dec,
        }
    )


    # --------------------------------------------------------
    # Remove exact duplicate coordinates
    # --------------------------------------------------------

    df_neighbors = (
        df_neighbors
        .drop_duplicates(
            subset=[
                "RA",
                "DEC"
            ]
        )
        .reset_index(
            drop=True
        )
    )


    print(
        f"RACS neighbors found: "
        f"{len(df_neighbors)}"
    )


    return df_neighbors