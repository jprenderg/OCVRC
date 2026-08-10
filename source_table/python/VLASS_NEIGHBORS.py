# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 17:23:59 2026

@author: Joe
"""

import re
import numpy as np
import pandas as pd

import astropy.io.fits as fits
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u

from config import VLASS1_LIST


# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

path_vlass1_links = VLASS1_LIST

search_radius_pix = 300

peak_box_size = 5
half_box = peak_box_size // 2

threshold_mjy = 1.0
threshold_jy = threshold_mjy / 1000.0

# ------------------------------------------------------------
# Neighbor separation limits
# ------------------------------------------------------------

min_neighbor_sep_arcsec = 2.0
max_neighbor_sep_arcsec = 300.0


# ------------------------------------------------------------
# Read and parse VLASS image-link list once
# ------------------------------------------------------------

def make_vlass_table(path_vlass_links):

    with open(path_vlass_links) as f:
        img_list = f.read().splitlines()

    ra_slice = slice(83, 89)
    dec_slice = slice(89, 96)

    good_rows = []

    for line in img_list:

        try:

            tile_match = re.findall(r"T\w+", line)

            if len(tile_match) == 0:
                continue

            drop_flag = line[54]

            if drop_flag == "Q":
                continue

            vl_ra = line[ra_slice]
            vl_dec = line[dec_slice]
            vl_tile = tile_match[0]
            img_type = line[-14:-11]

            vl_radeg = (
                int(vl_ra[:2]) * 15.0
                + int(vl_ra[2:4]) / 60.0 * 15.0
                + float(vl_ra[4:]) / 3600.0 * 15.0
            )

            if vl_dec[0] == "-":

                vl_decdeg = (
                    int(vl_dec[:3])
                    - int(vl_dec[3:5]) / 60.0
                )

            else:

                vl_decdeg = (
                    int(vl_dec[:3])
                    + int(vl_dec[3:5]) / 60.0
                )

            good_rows.append({
                "IMG_LINK": line,
                "VL_RA": vl_ra,
                "VL_DEC": vl_dec,
                "VL_TILE": vl_tile,
                "TYPE": img_type,
                "VL_RADEG": vl_radeg,
                "VL_DECDEG": vl_decdeg,
                "VL_RADEG_H": vl_radeg + 0.5,
                "VL_RADEG_L": vl_radeg - 0.5,
                "VL_DECDEG_H": vl_decdeg + 0.5,
                "VL_DECDEG_L": vl_decdeg - 0.5
            })

        except Exception:
            continue

    vlass_data = pd.DataFrame(good_rows)

    if len(vlass_data) == 0:
        raise RuntimeError(
            "No usable VLASS image-link rows were parsed."
        )

    return vlass_data


vlass_data = make_vlass_table(
    path_vlass1_links
)


# ------------------------------------------------------------
# Function
# ------------------------------------------------------------

def vlass_neighbors(ra, dec):
    """
    Find VLASS local-peak neighbors around the primary source.

    A source is considered a neighbor only if its angular
    separation from the primary satisfies:

        separation > 2 arcsec
        separation < 300 arcsec

    Parameters
    ----------
    ra : float
        Primary-source RA in degrees.

    dec : float
        Primary-source DEC in degrees.

    Returns
    -------
    df_neighbors : pandas.DataFrame
        Columns:
            RA
            DEC
    """

    print("\nStart vlass_neighbors\n")

    empty = pd.DataFrame(
        columns=["RA", "DEC"]
    )

    # --------------------------------------------------------
    # Validate input coordinate
    # --------------------------------------------------------

    if pd.isna(ra) or pd.isna(dec):
        return empty

    if not np.isfinite(ra) or not np.isfinite(dec):
        return empty

    # --------------------------------------------------------
    # Primary coordinate
    # --------------------------------------------------------

    primary_coord = SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        frame="icrs"
    )

    # --------------------------------------------------------
    # Find VLASS tiles containing primary coordinate
    # --------------------------------------------------------

    match = vlass_data.loc[
        (
            (vlass_data["VL_RADEG_L"] <= ra)
            &
            (ra <= vlass_data["VL_RADEG_H"])
        )
        &
        (
            (vlass_data["VL_DECDEG_L"] <= dec)
            &
            (dec <= vlass_data["VL_DECDEG_H"])
        )
    ].reset_index(drop=True)

    if len(match) == 0:
        return empty

    all_peaks = []

    # --------------------------------------------------------
    # Process science/RMS image pairs
    # --------------------------------------------------------

    for k in range(
        2,
        len(match) + 1,
        2
    ):

        try:

            match2 = match.iloc[
                k - 2:k
            ].copy()

            if not any(
                match2["TYPE"] == "tt0"
            ):
                continue

            if not any(
                match2["TYPE"] == "rms"
            ):
                continue

            link = match2.loc[
                match2["TYPE"] == "tt0",
                "IMG_LINK"
            ].iloc[0]

            link_rms = match2.loc[
                match2["TYPE"] == "rms",
                "IMG_LINK"
            ].iloc[0]

            # ------------------------------------------------
            # Read science image
            # ------------------------------------------------

            with fits.open(
                link,
                memmap=True
            ) as hdul:

                header = (
                    hdul[0]
                    .header
                    .copy()
                )

                data = hdul[0].data

                if data.ndim == 4:

                    image = data[
                        0,
                        0,
                        :,
                        :
                    ]

                elif data.ndim == 2:

                    image = data

                else:
                    continue

            # ------------------------------------------------
            # Read RMS image
            # ------------------------------------------------

            with fits.open(
                link_rms,
                memmap=True
            ) as hdul_rms:

                data_rms = (
                    hdul_rms[0].data
                )

                if data_rms.ndim == 4:

                    image_rms = data_rms[
                        0,
                        0,
                        :,
                        :
                    ]

                elif data_rms.ndim == 2:

                    image_rms = data_rms

                else:
                    continue

            # ------------------------------------------------
            # Validate image dimensions
            # ------------------------------------------------

            ny, nx = image.shape

            if image_rms.shape != image.shape:
                continue

            if ny <= 0 or nx <= 0:
                continue

            # ------------------------------------------------
            # Build 2D WCS
            # ------------------------------------------------

            h2 = header.copy()

            keys_to_remove = [
                "CRPIX4",
                "CRVAL4",
                "CDELT4",
                "CUNIT4",
                "CTYPE4",

                "CRPIX3",
                "CRVAL3",
                "CDELT3",
                "CUNIT3",
                "CTYPE3",

                "NAXIS3",
                "NAXIS4",

                "PC1_3",
                "PC2_3",
                "PC3_3",
                "PC4_3",

                "PC1_4",
                "PC2_4",
                "PC3_4",
                "PC4_4",

                "PC3_1",
                "PC4_1",

                "PC3_2",
                "PC4_2"
            ]

            for key in keys_to_remove:

                if key in h2:
                    h2.remove(key)

            h2["NAXIS"] = 2

            mywcs = wcs.WCS(h2)

            # ------------------------------------------------
            # Convert primary RA/DEC to image pixel
            # ------------------------------------------------

            pix = mywcs.wcs_world2pix(
                [(ra, dec)],
                0
            )[0]

            if (
                not np.isfinite(pix[0])
                or
                not np.isfinite(pix[1])
            ):
                continue

            xpix = int(
                round(pix[0])
            )

            ypix = int(
                round(pix[1])
            )

            if (
                xpix < 0
                or
                xpix >= nx
                or
                ypix < 0
                or
                ypix >= ny
            ):
                continue

            # ------------------------------------------------
            # Define candidate search region
            # ------------------------------------------------

            y0 = max(
                half_box,
                ypix - search_radius_pix
            )

            y1 = min(
                ny - half_box,
                ypix + search_radius_pix + 1
            )

            x0 = max(
                half_box,
                xpix - search_radius_pix
            )

            x1 = min(
                nx - half_box,
                xpix + search_radius_pix + 1
            )

            if (
                y1 <= y0
                or
                x1 <= x0
            ):
                continue

            center = image[
                y0:y1,
                x0:x1
            ]

            if center.size == 0:
                continue

            if not np.isfinite(
                center
            ).any():
                continue

            yy, xx = np.mgrid[
                y0:y1,
                x0:x1
            ]

            # ------------------------------------------------
            # Initial pixel-radius search
            #
            # This only limits the image area searched.
            # The exact 2-300 arcsec criterion is applied
            # later using SkyCoord.
            # ------------------------------------------------

            inside_radius = (
                (xx - xpix) ** 2
                +
                (yy - ypix) ** 2
            ) <= search_radius_pix ** 2

            # ------------------------------------------------
            # Initial peak mask
            # ------------------------------------------------

            peak_mask = (
                np.isfinite(center)
                &
                (center > threshold_jy)
                &
                inside_radius
            )

            # ------------------------------------------------
            # Require each candidate to be the maximum
            # within the 5x5 box
            # ------------------------------------------------

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
                        neighbor_y1 > ny
                        or
                        neighbor_x1 > nx
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

            # ------------------------------------------------
            # Candidate peak pixels
            # ------------------------------------------------

            peak_y = yy[
                peak_mask
            ]

            peak_x = xx[
                peak_mask
            ]

            if len(peak_x) == 0:
                continue

            # ------------------------------------------------
            # Validate peak pixels
            # ------------------------------------------------

            good_peak = (
                np.isfinite(peak_x)
                &
                np.isfinite(peak_y)
                &
                (peak_x >= 0)
                &
                (peak_x < nx)
                &
                (peak_y >= 0)
                &
                (peak_y < ny)
                &
                (
                    peak_x - half_box
                    >= 0
                )
                &
                (
                    peak_x + half_box
                    < nx
                )
                &
                (
                    peak_y - half_box
                    >= 0
                )
                &
                (
                    peak_y + half_box
                    < ny
                )
            )

            peak_x = peak_x[
                good_peak
            ].astype(int)

            peak_y = peak_y[
                good_peak
            ].astype(int)

            if len(peak_x) == 0:
                continue

            # ------------------------------------------------
            # Keep only finite science/RMS values
            # ------------------------------------------------

            peak_flux = image[
                peak_y,
                peak_x
            ]

            peak_rms = image_rms[
                peak_y,
                peak_x
            ]

            finite_value = (
                np.isfinite(peak_flux)
                &
                np.isfinite(peak_rms)
            )

            peak_x = peak_x[
                finite_value
            ]

            peak_y = peak_y[
                finite_value
            ]

            if len(peak_x) == 0:
                continue

            # ------------------------------------------------
            # Convert peak pixels to RA/DEC
            # ------------------------------------------------

            world = mywcs.wcs_pix2world(
                np.column_stack([
                    peak_x,
                    peak_y
                ]),
                0
            )

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
                continue

            # ------------------------------------------------
            # Calculate exact angular separation from primary
            # ------------------------------------------------

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

            # ------------------------------------------------
            # NEIGHBOR DEFINITION
            #
            # Must be:
            #
            #     > 2 arcsec
            #     < 300 arcsec
            #
            # from primary source.
            # ------------------------------------------------

            neighbor_mask = (
                (
                    separation_arcsec
                    >
                    min_neighbor_sep_arcsec
                )
                &
                (
                    separation_arcsec
                    <
                    max_neighbor_sep_arcsec
                )
            )

            peak_ra = peak_ra[
                neighbor_mask
            ]

            peak_dec = peak_dec[
                neighbor_mask
            ]

            if len(peak_ra) == 0:
                continue

            # ------------------------------------------------
            # Add retained neighbors
            # ------------------------------------------------

            df_peaks = pd.DataFrame({
                "RA": peak_ra,
                "DEC": peak_dec
            })

            all_peaks.append(
                df_peaks
            )

        except Exception as e:

            print(
                "VLASS image pair skipped:",
                repr(e)
            )

            continue

    # --------------------------------------------------------
    # No neighbors found
    # --------------------------------------------------------

    if len(all_peaks) == 0:

        print(
            "VLASS neighbors found: 0"
        )

        return empty

    # --------------------------------------------------------
    # Combine peaks from all image pairs
    # --------------------------------------------------------

    df_neighbors = pd.concat(
        all_peaks,
        ignore_index=True
    )

    # --------------------------------------------------------
    # Remove exact duplicate coordinates
    # --------------------------------------------------------

    df_neighbors = (
        df_neighbors
        .drop_duplicates(
            subset=["RA", "DEC"]
        )
        .reset_index(drop=True)
    )

    # --------------------------------------------------------
    # Final safety check using exact sky separation
    #
    # This ensures that even after combining image pairs,
    # every returned row satisfies:
    #
    #     2 < separation < 300 arcsec
    # --------------------------------------------------------

    neighbor_coords = SkyCoord(
        ra=df_neighbors["RA"].to_numpy() * u.deg,
        dec=df_neighbors["DEC"].to_numpy() * u.deg,
        frame="icrs"
    )

    separation_arcsec = (
        primary_coord
        .separation(neighbor_coords)
        .arcsec
    )

    keep = (
        (
            separation_arcsec
            >
            min_neighbor_sep_arcsec
        )
        &
        (
            separation_arcsec
            <
            max_neighbor_sep_arcsec
        )
    )

    df_neighbors = (
        df_neighbors.loc[
            keep,
            ["RA", "DEC"]
        ]
        .reset_index(drop=True)
    )

    print(
       f"VLASS neighbors found: "
       f"{len(df_neighbors)}"
    )

    return df_neighbors


# ------------------------------------------------------------
# Test
# ------------------------------------------------------------

# df_neighbors = vlass_neighbors(
#     ra=272.73,
#     dec=7.678333
# )
#
# print(df_neighbors)