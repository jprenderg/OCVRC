# -*- coding: utf-8 -*-
"""
Created on Tue May 26 08:43:55 2026

@author: Joe
"""

from astroquery.heasarc import Heasarc
from astroquery.exceptions import NoResultsWarning
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import warnings
import pandas as pd


# ------------------------------------------------------------
# Input coordinates
# ------------------------------------------------------------

def rosat(source_id, ra_deg, dec_deg): 
    
    print("\nStart rosat\n")

    # ------------------------------------------------------------
    # Defaults
    # ------------------------------------------------------------
    
    ROSAT_COUNT_RATE = None
    ROSAT_COUNT_RATE_ERR = None
    
    # ------------------------------------------------------------
    # HEASARC setup
    # ------------------------------------------------------------
    
    heasarc = Heasarc()
    radius = 30 * u.arcsec
    
    warnings.simplefilter("ignore", NoResultsWarning)
    
    # ------------------------------------------------------------
    # Helper
    # ------------------------------------------------------------
    
    def get_col(table, candidates):
        cols_lower = {c.lower(): c for c in table.colnames}
    
        for c in candidates:
            if c.lower() in cols_lower:
                return cols_lower[c.lower()]
    
        return None
    
    # ------------------------------------------------------------
    # Query
    # ------------------------------------------------------------
    
    try:
        coord = SkyCoord(
            ra=ra_deg * u.deg,
            dec=dec_deg * u.deg,
            frame="icrs"
        )
    
        res = heasarc.query_region(
            coord,
            catalog="rass2rxs",
            radius=radius
        )
    
        if len(res) > 0:
    
            ra_col = get_col(res, ["RA"])
            dec_col = get_col(res, ["DEC"])
    
            cr_col = get_col(res, ["Count_Rate"])
            crerr_col = get_col(res, ["Count_Rate_Error"])
    
            if ra_col is not None and dec_col is not None:
    
                src_coords = SkyCoord(
                    np.array(res[ra_col], dtype=float) * u.deg,
                    np.array(res[dec_col], dtype=float) * u.deg,
                    frame="icrs"
                )
    
                sep = coord.separation(src_coords)
                idx = np.argmin(sep.arcsec)
    
                if cr_col is not None:
                    ROSAT_COUNT_RATE = float(res[cr_col][idx])
    
                if crerr_col is not None:
                    ROSAT_COUNT_RATE_ERR = float(res[crerr_col][idx])
    
    except Exception as e:
        print(f"ROSAT query failed for UNIQUE_ID {unique_id}: {e}")
    
    # ------------------------------------------------------------
    # Output dataframe
    # ------------------------------------------------------------
    
    df_out = pd.DataFrame([{
        "SOURCE_ID": source_id,
        "OBSERVATORY": "ROSAT",
        "CAT": "rass2rxs",
        "BAND": "0.1-2 keV",
        "UNITS": "counts s-1",
        "MEASUREMENT": ROSAT_COUNT_RATE,
        "ERROR": ROSAT_COUNT_RATE_ERR
    }])
    
    return df_out
    
# unique_id = 1
# ra_deg = 0.47958
# dec_deg = -67.12869

# temp = rosat(unique_id, ra_deg, dec_deg)