"""Benchmark the local TGLC pipeline against the MAST TGLC HLSP.

The local light curve comes from :func:`quicklook.tglc.get_tglc_lc` (TESScut
cutout + ePSF fit run here).  The reference comes from the TGLC High Level
Science Product hosted at MAST (https://archive.stsci.edu/hlsp/tglc),
downloaded through ``lightkurve``.

The two products are aligned cadence-by-cadence on BTJD time and each
photometry column (``cal_aper_flux``, ``cal_psf_flux``, ``aperture_flux``,
``psf_flux``) is compared with a handful of robust metrics.  Raw fluxes are
median-normalised before comparison because the MAST product is built from a
150x150 px FFI cut while the local pipeline uses a <=99 px TESScut cutout, so
absolute background levels are not expected to agree exactly.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field

import numpy as np

logger = logging.getLogger(__name__)

COLUMNS = ("cal_aper_flux", "cal_psf_flux", "aperture_flux", "psf_flux")
CALIBRATED_COLUMNS = ("cal_aper_flux", "cal_psf_flux")

# Half of the shortest FFI cadence (200 s) in days: two cadences closer than
# this are the same exposure.
DEFAULT_TIME_TOL_DAYS = 100.0 / 86400.0
MAD_TO_STD = 1.4826

# Pass criteria applied to the calibrated (normalised, detrended) columns.
DEFAULT_MIN_PEARSON = 0.9
DEFAULT_MAX_RESID_TO_SCATTER = 0.5
DEFAULT_MIN_OVERLAP_FRAC = 0.8


@dataclass(frozen=True)
class ColumnMetrics:
    """Agreement metrics for one photometry column."""

    column: str
    n_compared: int
    pearson_r: float
    median_ratio: float
    resid_std: float
    local_scatter: float
    mast_scatter: float
    resid_to_scatter: float
    scatter_ratio: float


@dataclass(frozen=True)
class BenchmarkResult:
    """Outcome of a local-vs-MAST TGLC comparison."""

    target: str
    sector: int
    n_local: int
    n_mast: int
    n_matched: int
    overlap_frac: float
    columns: dict = field(default_factory=dict)
    passed: bool = False
    failures: tuple = ()

    def to_dict(self) -> dict:
        out = asdict(self)
        out["failures"] = list(self.failures)
        return out


def robust_std(x: np.ndarray) -> float:
    """MAD-based standard deviation, ignoring non-finite values."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(MAD_TO_STD * np.median(np.abs(x - np.median(x))))


def point_to_point_scatter(x: np.ndarray) -> float:
    """Robust point-to-point scatter, insensitive to slow trends."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return float("nan")
    return robust_std(np.diff(x)) / np.sqrt(2.0)


def match_cadences(
    t_local: np.ndarray, t_mast: np.ndarray, tol: float = DEFAULT_TIME_TOL_DAYS
) -> tuple[np.ndarray, np.ndarray]:
    """Return index arrays pairing each MAST cadence with its nearest local one.

    Pairs further apart than ``tol`` days are dropped, and each local cadence
    is used at most once.
    """
    t_local = np.asarray(t_local, dtype=float)
    t_mast = np.asarray(t_mast, dtype=float)
    if t_local.size == 0 or t_mast.size == 0:
        empty = np.array([], dtype=int)
        return empty, empty

    order = np.argsort(t_local)
    t_sorted = t_local[order]
    if t_sorted.size == 1:
        nearest = np.zeros(t_mast.size, dtype=int)
    else:
        pos = np.clip(np.searchsorted(t_sorted, t_mast), 1, t_sorted.size - 1)
        left, right = t_sorted[pos - 1], t_sorted[pos]
        nearest = np.where(np.abs(t_mast - left) <= np.abs(t_mast - right), pos - 1, pos)
    dt = np.abs(t_sorted[nearest] - t_mast)
    ok = dt <= tol

    idx_local = order[nearest[ok]]
    idx_mast = np.flatnonzero(ok)
    _, first = np.unique(idx_local, return_index=True)
    return idx_local[first], idx_mast[first]


def compare_series(local: np.ndarray, mast: np.ndarray, column: str = "") -> ColumnMetrics:
    """Compare two already-aligned flux series.

    Both series are divided by their own median so that a constant scale
    offset (e.g. different background estimates) does not count as
    disagreement; ``median_ratio`` reports that offset separately.
    """
    local = np.asarray(local, dtype=float)
    mast = np.asarray(mast, dtype=float)
    if local.shape != mast.shape:
        raise ValueError(f"shape mismatch: local {local.shape} vs mast {mast.shape}")

    good = np.isfinite(local) & np.isfinite(mast)
    n = int(good.sum())
    nan = float("nan")
    if n < 3:
        return ColumnMetrics(column, n, nan, nan, nan, nan, nan, nan, nan)

    a, b = local[good], mast[good]
    med_a, med_b = np.median(a), np.median(b)
    median_ratio = float(med_a / med_b) if med_b != 0 else nan
    a_n = a / med_a if med_a != 0 else a
    b_n = b / med_b if med_b != 0 else b

    if np.std(a_n) == 0 or np.std(b_n) == 0:
        pearson = nan
    else:
        pearson = float(np.corrcoef(a_n, b_n)[0, 1])

    resid_std = robust_std(a_n - b_n)
    local_scatter = point_to_point_scatter(a_n)
    mast_scatter = point_to_point_scatter(b_n)
    resid_to_scatter = resid_std / mast_scatter if mast_scatter > 0 else nan
    scatter_ratio = local_scatter / mast_scatter if mast_scatter > 0 else nan

    return ColumnMetrics(
        column=column,
        n_compared=n,
        pearson_r=pearson,
        median_ratio=median_ratio,
        resid_std=resid_std,
        local_scatter=local_scatter,
        mast_scatter=mast_scatter,
        resid_to_scatter=resid_to_scatter,
        scatter_ratio=scatter_ratio,
    )


def _evaluate(
    metrics: dict,
    overlap_frac: float,
    min_pearson: float,
    max_resid_to_scatter: float,
    min_overlap_frac: float,
) -> tuple[str, ...]:
    """Return human-readable failure reasons (empty tuple means pass)."""
    failures = []
    if not overlap_frac >= min_overlap_frac:
        failures.append(f"cadence overlap {overlap_frac:.2f} < {min_overlap_frac}")
    # A column with too few finite cadences (e.g. MAST PSF fluxes are all NaN
    # for some HLSP files) is unavailable rather than wrong; it is skipped, but
    # at least one calibrated column must be comparable.
    compared = [c for c in CALIBRATED_COLUMNS if c in metrics and metrics[c].n_compared >= 3]
    if not compared:
        failures.append("no calibrated column has enough finite cadences to compare")
    for col in compared:
        m = metrics[col]
        if not m.pearson_r >= min_pearson:
            failures.append(f"{col}: pearson_r {m.pearson_r:.3f} < {min_pearson}")
        if not m.resid_to_scatter <= max_resid_to_scatter:
            failures.append(
                f"{col}: resid/scatter {m.resid_to_scatter:.3f} > {max_resid_to_scatter}"
            )
    return tuple(failures)


def mast_good_mask(mast_lc) -> np.ndarray:
    """Cadences with neither TESS nor TGLC quality flags set."""
    good = np.ones(len(mast_lc), dtype=bool)
    for col in ("tess_flags", "tglc_flags"):
        if col in mast_lc.colnames:
            good &= np.asarray(mast_lc[col]).astype(int) == 0
    return good


def _values(lc, col: str) -> np.ndarray:
    data = lc[col]
    return np.asarray(getattr(data, "value", data), dtype=float)


def compare_light_curves(
    local_lc,
    mast_lc,
    target: str = "",
    sector: int = -1,
    columns: tuple = COLUMNS,
    time_tol: float = DEFAULT_TIME_TOL_DAYS,
    min_pearson: float = DEFAULT_MIN_PEARSON,
    max_resid_to_scatter: float = DEFAULT_MAX_RESID_TO_SCATTER,
    min_overlap_frac: float = DEFAULT_MIN_OVERLAP_FRAC,
) -> BenchmarkResult:
    """Align a local TGLC light curve with the MAST one and score agreement."""
    good = mast_good_mask(mast_lc)
    mast = mast_lc[good]
    idx_local, idx_mast = match_cadences(local_lc.time.value, mast.time.value, tol=time_tol)
    n_matched = int(idx_local.size)
    overlap_frac = n_matched / len(mast) if len(mast) else 0.0

    metrics = {}
    for col in columns:
        if col not in local_lc.colnames or col not in mast.colnames:
            logger.warning(f"Column {col!r} missing from one light curve; skipped.")
            continue
        metrics[col] = compare_series(
            _values(local_lc, col)[idx_local], _values(mast, col)[idx_mast], column=col
        )

    failures = _evaluate(metrics, overlap_frac, min_pearson, max_resid_to_scatter, min_overlap_frac)
    return BenchmarkResult(
        target=target,
        sector=int(sector),
        n_local=len(local_lc),
        n_mast=len(mast),
        n_matched=n_matched,
        overlap_frac=float(overlap_frac),
        columns={k: asdict(v) for k, v in metrics.items()},
        passed=not failures,
        failures=failures,
    )


def download_mast_tglc(target: str, sector: int, download_dir: str | None = None):
    """Download the MAST TGLC HLSP light curve for ``target`` in ``sector``.

    ``quality_bitmask='none'`` is used on purpose: ``'hardest'`` overflows the
    int16 ``tess_flags`` column under NumPy 2 (see
    ``LIGHTKURVE_TGLC_BUG_REPORT.txt``); flags are applied later in
    :func:`mast_good_mask` instead.
    """
    import lightkurve as lk

    search = lk.search_lightcurve(target, author="TGLC", sector=int(sector))
    if len(search) == 0:
        raise LookupError(f"No MAST TGLC light curve for {target} sector {sector}")
    return search[0].download(quality_bitmask="none", download_dir=download_dir)


def plot_comparison(local_lc, mast_lc, result: BenchmarkResult, outfile: str):
    """Overlay local and MAST light curves for each compared column."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mast = mast_lc[mast_good_mask(mast_lc)]
    cols = [c for c, m in result.columns.items() if m["n_compared"] >= 3]
    if not cols:
        logger.warning("Nothing to plot: no column has comparable cadences.")
        return None
    fig, axes = plt.subplots(len(cols), 1, figsize=(10, 2.6 * len(cols)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, col in zip(axes, cols):
        for lc, label, kw in (
            (mast, "MAST", {"color": "k", "alpha": 0.6}),
            (local_lc, "local", {"color": "C1", "alpha": 0.6}),
        ):
            y = _values(lc, col)
            med = np.nanmedian(y)
            ax.plot(lc.time.value, y / med if med else y, ".", ms=2, label=label, **kw)
        m = result.columns[col]
        ax.set_ylabel(col)
        ax.set_title(
            f"r={m['pearson_r']:.3f}  resid/scatter={m['resid_to_scatter']:.2f}  "
            f"scatter ratio={m['scatter_ratio']:.2f}",
            fontsize=9,
        )
        ax.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("BTJD")
    status = "PASS" if result.passed else "FAIL"
    fig.suptitle(f"{result.target} s{result.sector}: local TGLC vs MAST HLSP [{status}]")
    fig.tight_layout()
    fig.savefig(outfile, dpi=120)
    plt.close(fig)
    return outfile


def run_benchmark(
    target: str,
    sector: int,
    size: int = 99,
    limit_mag: float = 16,
    cache_dir: str | None = None,
    use_cache: bool = True,
    verbose: bool = False,
    **compare_kwargs,
):
    """Run the local pipeline, download MAST, and compare.

    Returns ``(result, local_lc, mast_lc)``.
    """
    from quicklook.tglc import get_tglc_lc

    mast_lc = download_mast_tglc(target, sector)
    local_lc = get_tglc_lc(
        target,
        sector=sector,
        size=size,
        limit_mag=limit_mag,
        cache_dir=cache_dir,
        use_cache=use_cache,
        verbose=verbose,
    )
    result = compare_light_curves(local_lc, mast_lc, target=target, sector=sector, **compare_kwargs)
    return result, local_lc, mast_lc
