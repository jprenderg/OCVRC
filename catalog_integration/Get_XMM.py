"""
XMM XSA EPIC source catalogue via ESA TAP.

HTTP 503 errors are normalized so the shared query manager
can recognize service outages without printing HTML.
"""

from _xray_common import *
from astroquery.esa.xmm_newton import XMMNewton
import re


# ============================================================
# XMM QUERY
# ============================================================

def query(ra, dec, radius, timeout):

    radius_deg = radius / 3600

    q = f"""
        SELECT ra, dec,
               ep_1_flux, ep_1_flux_err,
               ep_2_flux, ep_2_flux_err,
               ep_3_flux, ep_3_flux_err,
               ep_4_flux, ep_4_flux_err,
               ep_5_flux, ep_5_flux_err
        FROM xsa.v_epic_source_cat
        WHERE 1 = CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra}, {dec}, {radius_deg})
        )
    """

    try:
        return XMMNewton.query_xsa_tap(
            q,
            output_format='votable'
        )

    except Exception as exc:

        error_text = str(exc)

        # ESA service unavailable: suppress returned HTML.
        if re.search(r'\b503\b', error_text):
            raise RuntimeError(
                'HTTP 503: Service Unavailable'
            ) from None

        # Other errors: suppress multiline HTML and long messages.
        message = error_text.splitlines()[0][:150]

        raise RuntimeError(
            f'{type(exc).__name__}: {message}'
        ) from None


# ============================================================
# XMM CATALOG FIELDS
# ============================================================

fields = {
    f'ep_{i}_{s}': [f'ep_{i}_{s}']
    for i in range(1, 6)
    for s in ('flux', 'flux_err')
}


# ============================================================
# ENERGY BANDS
# ============================================================

bands = [
    (
        energy,
        f'ep_{i}_flux',
        f'ep_{i}_flux_err',
        None,
        None
    )
    for i, energy in enumerate(
        [
            '0.2–0.5 keV',
            '0.5–1.0 keV',
            '1.0–2.0 keV',
            '2.0–4.5 keV',
            '4.5–12.0 keV'
        ],
        1
    )
]


# ============================================================
# EXECUTE QUERY
# ============================================================

df_xmm, df_xmm_measurement, df_xmm_failures = run(
    'XMM',
    query,
    fields,
    bands,
    'XMM',
    'xsa.v_epic_source_cat',
    'erg cm-2 s-1',
    5
)


# ============================================================
# EXPOSE STATUS IN SPYDER
# ============================================================

df_query_status = pd.read_pickle(
    OUTPUT / 'XMM_query_status.pkl'
)