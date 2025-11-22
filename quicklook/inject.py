#!/usr/bin/env python
"""
Experimental transit injection and recovery
used to determine the best window length
to flatten the baseline model to avoid
overfitting the transit.
A grid search over a list of window lengths
and transit parameters.
This is used when --window_length is set to 0 in ql script.
"""
from dataclasses import dataclass
from typing import List
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from wotan import flatten
from transitleastsquares import transitleastsquares


# --------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------
@dataclass
class InjectionParams:
    t0: float
    period: float
    duration: float
    depth: float


@dataclass
class InjectionResult:
    window: float
    depth_ratio: float
    duration_ratio: float
    asymmetry: float
    recovered_sde: float


# --------------------------------------------------------------------
# Transit injector
# --------------------------------------------------------------------
def inject_transit(time: np.ndarray, flux: np.ndarray, params: InjectionParams) -> np.ndarray:

    phase = ((time - params.t0) / params.period) % 1
    in_transit = phase < (params.duration / params.period)

    out = flux.copy()
    out[in_transit] *= 1 - params.depth
    return out


# --------------------------------------------------------------------
# Transit symmetry metric
# --------------------------------------------------------------------
def phase_asymmetry(time, flux, period, t0):
    phase = ((time - t0) / period) % 1
    order = np.argsort(phase)
    folded = flux[order]
    n = len(folded)
    left = folded[: n // 2]
    right = folded[n // 2 :]
    return abs(np.median(left) - np.median(right))


# --------------------------------------------------------------------
# One test evaluation
# --------------------------------------------------------------------
def evaluate_window(time, flux, window, inj: InjectionParams, method="biweight"):

    injected = inject_transit(time, flux, inj)

    flat, _ = flatten(
        time,
        injected,
        method=method,
        window_length=window,
        return_trend=True,
    )

    tls = transitleastsquares(time, flat)
    res = tls.power(show_progress_bar=False, verbose=False)

    depth_ratio = res.depth / inj.depth if inj.depth > 0 else np.nan
    duration_ratio = res.duration / inj.duration if inj.duration > 0 else np.nan
    asym = phase_asymmetry(time, flat, inj.period, inj.t0)

    return InjectionResult(
        window=window,
        depth_ratio=depth_ratio,
        duration_ratio=duration_ratio,
        asymmetry=asym,
        recovered_sde=res.SDE,
    )


# --------------------------------------------------------------------
# Sequential grid search
# --------------------------------------------------------------------
def run_grid(
    time: np.ndarray,
    flux: np.ndarray,
    windows: List[float],
    injections: List[InjectionParams],
    method="biweight",
):

    results = []

    for w in tqdm(windows, desc="window"):
        for inj in tqdm(injections, desc="transit", leave=False):
            results.append(evaluate_window(time, flux, w, inj, method=method))

    return pd.DataFrame([r.__dict__ for r in results])


# --------------------------------------------------------------------
# Select optimal window
# --------------------------------------------------------------------
def select_best_window(df, depth_min=0.9, duration_tol=0.2, asym_max=0.002):

    good = df[
        (df.depth_ratio >= depth_min)
        & (df.duration_ratio.between(1 - duration_tol, 1 + duration_tol))
        & (df.asymmetry <= asym_max)
    ]

    if good.empty:
        return None, df

    good["score"] = good.depth_ratio - 2 * abs(good.duration_ratio - 1) - 1.5 * good.asymmetry

    best = good.sort_values("score", ascending=False).iloc[0]
    return best.window, good


if __name__ == "__main__":
    import lightkurve as lk

    res = lk.search_lightcurve("TOI-837", author="tglc")
    lc = res[0].download()
    time, flux = lc.time.value, lc.flux.value

    np.random.default_rng(1)
    baseline = time[-1] - time[0]
    half_baseline = baseline / 2
    Nmodels = 10
    t0s = np.random.uniform(low=0, high=half_baseline, size=Nmodels)
    periods = np.random.uniform(low=0.1, high=half_baseline, size=Nmodels)
    durations = np.random.uniform(low=0.05, high=0.5, size=Nmodels)
    depths = np.random.uniform(low=0.005, high=0.05, size=Nmodels)

    windows = np.linspace(0.1, 0.5, 5)

    # Several injections
    injections = []
    for t0, p, dur, depth in zip(t0s, periods, durations, depths):
        injections.append(InjectionParams(t0=time[0] + t0, period=p, duration=dur, depth=depth))

    # Run grid
    method = "biweight"
    df = run_grid(time, flux, windows, injections, method=method)

    # Select optimal window
    best_window, valid = select_best_window(df)
    print("Selected window length:", best_window)

    # Run grid with different flatten method
    # method2 = 'rspline'
    # df2 = run_grid(time, flux, windows, injections, method=method2)

    # best_window2, valid = select_best_window(df2)
    # print(f"Best window length {best_window2} using method={method2}")

    flat, trend = flatten(
        time,
        flux,
        method="biweight",
        window_length=best_window,
        return_trend=True,
    )

    # flat2, trend2 = flatten(
    #         time,
    #         flux,
    #         method="biweight",
    #         window_length=best_window2,
    #         return_trend=True,
    #     )

    ax = lc.scatter()
    ax.plot(time, trend, c="C3")
    # ax.plot(time, trend2, c='C0')

    # Run TLS
    tls = transitleastsquares(time, flat)
    results = tls.power(verbose=False)

    # tls2 = transitleastsquares(time, flat2)
    # results2 = tls2.power()
