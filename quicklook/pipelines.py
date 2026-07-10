"""Single source of truth for the TESS light-curve pipelines we support.

Two orthogonal classifications live here:

* ``FULL_FRAME_TESS_PIPELINES`` — pipelines that derive light curves from
  Full Frame Images (FFI). These need separate handling in the analysis
  path (cadence detection, aperture choice, etc.).
* ``HLSP_PIPELINES`` — pipelines whose saved output filenames follow the
  HLSP naming convention ``<target>_s<NN>_<pipeline>_<cadence>``. This
  set drives the filename → pipeline reverse-lookup used by the GUI.

These overlap but are not identical: TASOC is not full-frame derived but
its output filenames use the HLSP convention, so it appears in
``HLSP_PIPELINES`` only.

Adding a new pipeline:
  1. Append the lowercase name to ``ALL_TESS_PIPELINES`` (preserves the
     order shown in the GUI form).
  2. If it derives from FFI, add to ``FULL_FRAME_TESS_PIPELINES``.
  3. If its saved-file stems embed the pipeline name as the lctype
     token, add to ``HLSP_PIPELINES``.
"""

from __future__ import annotations

# Order is meaningful — the GUI form renders pipelines in this order.
FULL_FRAME_TESS_PIPELINES: list[str] = [
    "tess-spoc",
    "qlp",
    "tglc",
    "cdips",
    "pathos",
    "eleanor",
    "t16",
    "gsfc-eleanor-lite",
    "tequila",
    "tica",
    "diamante",
]

ALL_TESS_PIPELINES: list[str] = ["spoc", "tasoc"] + FULL_FRAME_TESS_PIPELINES

# Pipelines whose saved-file stems put the pipeline name as the lctype
# token (e.g. "TIC123_s12_qlp_lc_tls.h5"). Used by app.py's filename
# parser when recovering pipeline/cadence from a stem.
HLSP_PIPELINES: frozenset[str] = frozenset(
    {
        "qlp",
        "tglc",
        "tasoc",
        "cdips",
        "pathos",
        "tess-spoc",
        "t16",
    }
)
