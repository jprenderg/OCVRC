# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 07:26:48 2026

@author: Joe
"""

import pandas as pd
import numpy as np
import re

from astropy.io import fits
from astropy import wcs

from config import VLASS1_LIST
from config import VLASS23_LIST


def vlass(epoch_code, source_id, ra_deg, dec_deg):
    """
    Measure VLASS flux/RMS for one source.

    Returns dataframe with columns:
        SOURCE_ID
        OBSERVATORY
        CAT
        BAND
        UNITS
        MEASUREMENT
        ERROR
    """
    
    print("\nStart vlass\n")

    # ------------------------------------------------------------
    # Controls
    # ------------------------------------------------------------

    dim = 5

    if epoch_code == 1:
        img_list_path = VLASS1_LIST
        
        ra_slice = slice(83, 89)
        dec_slice = slice(89, 96)

    elif epoch_code == 23:
        img_list_path = VLASS23_LIST
        ra_slice = slice(81, 87)
        dec_slice = slice(87, 94)

    else:
        raise ValueError("epoch_code must be 1 or 23")

    empty_out = pd.DataFrame(
        columns=[
            "SOURCE_ID",
            "OBSERVATORY",
            "CAT",
            "BAND",
            "UNITS",
            "MEASUREMENT",
            "ERROR"
        ]
    )

    # ------------------------------------------------------------
    # Read image list
    # ------------------------------------------------------------

    with open(img_list_path) as f:
        img_list = [line.strip() for line in f if line.strip()]

    # ------------------------------------------------------------
    # Parse VLASS image list
    # ------------------------------------------------------------

    rows_img = []

    for link in img_list:

        try:
            tile = re.findall(r"T\w+", link)[0]

            drop_flag = link[54]
            if drop_flag == "Q":
                continue

            img_type = link[-14:-11]
            epoch = link[45:53]

            ra_str = link[ra_slice]
            dec_str = link[dec_slice]

            ra_center = (
                int(ra_str[:2]) * 15.0
                + int(ra_str[2:4]) / 60.0 * 15.0
                + float(ra_str[4:]) / 3600.0 * 15.0
            )

            if dec_str[0] == "-":
                dec_center = int(dec_str[:3]) - int(dec_str[3:5]) / 60.0
            else:
                dec_center = int(dec_str[:3]) + int(dec_str[3:5]) / 60.0

            rows_img.append({
                "img_link": link,
                "VL_tile": tile,
                "type": img_type,
                "epoch": epoch,
                "VL_radeg": ra_center,
                "VL_decdeg": dec_center,
                "VL_radeg_l": ra_center - 0.5,
                "VL_radeg_h": ra_center + 0.5,
                "VL_decdeg_l": dec_center - 0.5,
                "VL_decdeg_h": dec_center + 0.5
            })

        except Exception as e:
            print(f"Skipping image-list row because of parse error: {e}")
            continue

    vlass_data = pd.DataFrame(rows_img)

    if len(vlass_data) == 0:
        return empty_out

    # ------------------------------------------------------------
    # Match RA/DEC to possible VLASS tiles
    # ------------------------------------------------------------

    match = vlass_data.loc[
        (vlass_data["VL_radeg_l"] <= ra_deg)
        & (ra_deg <= vlass_data["VL_radeg_h"])
        & (vlass_data["VL_decdeg_l"] <= dec_deg)
        & (dec_deg <= vlass_data["VL_decdeg_h"])
    ].reset_index(drop=True)

    if len(match) == 0:
        return empty_out

    # ------------------------------------------------------------
    # Measure all matching epoch/tile tt0/rms pairs
    # ------------------------------------------------------------

    rows = []

    grouped = match.groupby(["epoch", "VL_tile"])

    for (epoch, tile), match_group in grouped:

        try:
            tt0_rows = match_group[match_group["type"] == "tt0"]
            rms_rows = match_group[match_group["type"] == "rms"]

            if len(tt0_rows) == 0 or len(rms_rows) == 0:
                print(f"Skipping {epoch} {tile}: missing tt0 or rms image")
                continue

            link = tt0_rows["img_link"].iloc[0]
            link_rms = rms_rows["img_link"].iloc[0]

            # ----------------------------------------------------
            # Open flux image
            # ----------------------------------------------------

            with fits.open(link, memmap=True) as hdul:
                header = hdul[0].header.copy()
                data = hdul[0].data

                if data.ndim == 4:
                    image = data[0, 0, :, :]
                elif data.ndim == 2:
                    image = data
                else:
                    print(f"Skipping {epoch} {tile}: unexpected data shape {data.shape}")
                    continue

            ny, nx = image.shape

            # ----------------------------------------------------
            # Reduce header to 2D WCS
            # ----------------------------------------------------

            header_2d = header.copy()

            remove_keys = [
                "CRPIX4", "CRVAL4", "CDELT4", "CUNIT4", "CTYPE4",
                "CRPIX3", "CRVAL3", "CDELT3", "CUNIT3", "CTYPE3",
                "NAXIS3", "NAXIS4",
                "PC1_3", "PC2_3", "PC3_3", "PC4_3",
                "PC1_4", "PC2_4", "PC3_4", "PC4_4",
                "PC3_1", "PC4_1",
                "PC3_2", "PC4_2"
            ]

            for key in remove_keys:
                if key in header_2d:
                    header_2d.remove(key)

            header_2d["NAXIS"] = 2
            mywcs = wcs.WCS(header_2d)

            # ----------------------------------------------------
            # Convert RA/DEC to pixel
            # ----------------------------------------------------

            pix = mywcs.wcs_world2pix([(ra_deg, dec_deg)], 0)[0]

            if not np.isfinite(pix[0]) or not np.isfinite(pix[1]):
                print(f"Skipping {epoch} {tile}: non-finite WCS pixel")
                continue

            xpix = int(round(pix[0]))
            ypix = int(round(pix[1]))

            # ----------------------------------------------------
            # Explicit pixel bounds check
            # ------------------------------------------------------------

            if (
                xpix - dim < 0
                or xpix + dim >= nx
                or ypix - dim < 0
                or ypix + dim >= ny
            ):
                print(
                    f"Skipping {epoch} {tile}: pixel out of bounds "
                    f"x={xpix}, y={ypix}, image shape={image.shape}"
                )
                continue

            # ----------------------------------------------------
            # 11 x 11 cutout centered on source pixel
            # ------------------------------------------------------------

            cutout = image[
                ypix - dim:ypix + dim + 1,
                xpix - dim:xpix + dim + 1
            ]

            if cutout.size == 0 or np.all(np.isnan(cutout)):
                print(f"Skipping {epoch} {tile}: empty or all-NaN cutout")
                continue

            max_y_cutout, max_x_cutout = np.unravel_index(
                np.nanargmax(cutout),
                cutout.shape
            )

            i_flux = ypix - dim + max_y_cutout
            j_flux = xpix - dim + max_x_cutout

            flux = image[i_flux, j_flux] * 1000.0

            # ----------------------------------------------------
            # Open RMS image
            # ------------------------------------------------------------

            with fits.open(link_rms, memmap=True) as hdul_rms:
                data_rms = hdul_rms[0].data

                if data_rms.ndim == 4:
                    image_rms = data_rms[0, 0, :, :]
                elif data_rms.ndim == 2:
                    image_rms = data_rms
                else:
                    print(f"Skipping {epoch} {tile}: unexpected RMS shape {data_rms.shape}")
                    continue

            rms_ny, rms_nx = image_rms.shape

            if (
                i_flux < 0
                or i_flux >= rms_ny
                or j_flux < 0
                or j_flux >= rms_nx
            ):
                print(
                    f"Skipping {epoch} {tile}: RMS pixel out of bounds "
                    f"x={j_flux}, y={i_flux}, RMS shape={image_rms.shape}"
                )
                continue

            rms = image_rms[i_flux, j_flux] * 1000.0

            if not np.isfinite(flux) or not np.isfinite(rms):
                print(f"Skipping {epoch} {tile}: non-finite flux or RMS")
                continue

            rows.append({
                "SOURCE_ID": source_id,
                "EPOCH": epoch,
                "FLUX": flux,
                "RMS": rms
            })

        except Exception as e:
            print(f"Exception for {epoch} {tile}: {e}")
            continue

    # ------------------------------------------------------------
    # Build intermediate dataframe
    # ------------------------------------------------------------

    df_interm = pd.DataFrame(
        rows,
        columns=[
            "SOURCE_ID",
            "EPOCH",
            "FLUX",
            "RMS"
        ]
    )

    if len(df_interm) == 0:
        return empty_out

    df_interm = df_interm.sort_values(
        ["SOURCE_ID", "EPOCH"]
    ).reset_index(drop=True)

    # ------------------------------------------------------------
    # Final measurement dataframe
    # ------------------------------------------------------------

    df_out = pd.DataFrame({
        "SOURCE_ID": df_interm["SOURCE_ID"],
        "OBSERVATORY": "VLA",
        "CAT": df_interm["EPOCH"],
        "BAND": "2-4 GHz",
        "UNITS": "mJy",
        "MEASUREMENT": df_interm["FLUX"],
        "ERROR": df_interm["RMS"]
    })

    return df_out


# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------

unique_id = 1
ra_deg = 268.794667
dec_deg = -30.175139
epoch_code = 23

# temp = vlass(epoch_code, unique_id, ra_deg, dec_deg)
