"""Duotransits: pair two single events and list the periods they allow.

Two transits separated by ``dT`` allow periods ``P = dT / n`` for whole numbers
``n``. An alias is ruled out when any of its other predicted transits lands in
covered data (Cooke et al. 2020; Hawthorn et al. 2024). The survivors are
ranked by a geometric/period prior and by how well the observed duration fits
a transit at that period around the host star (Osborn et al. 2016, 2022).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np

from quicklook.single_event import (
    SingleEvent,
    _cadence,
    _prepare,
    central_duration,
    period_is_allowed,
)

DEPTH_SIGMA_MAX = 3.0  # depths must agree within this many combined sigma
DEPTH_FRACTION_FLOOR = 0.05  # systematic depth error per event (box shape/fit)
DURATION_RATIO_MAX = 1.5
PERIOD_PRIOR_INDEX = -5 / 3  # transit probability ~P^-2/3 x log-uniform prior ~P^-1
_DURATION_TAIL_WIDTH = 0.2  # eccentric-orbit tail beyond the central duration
_X_CLIP = 0.99  # avoid the integrable singularity at D = T_central
_MIN_PERIOD_DURATIONS = 2  # shortest alias considered, in transit durations
MAX_ALIASES = 5000


@dataclass(frozen=True)
class DuoAlias:
    period: float
    n: int
    weight: float  # normalised over the allowed aliases


@dataclass(frozen=True)
class DuoPair:
    first: SingleEvent
    second: SingleEvent
    delta_t: float
    aliases: tuple[DuoAlias, ...]
    depth_ratio: float  # second / first
    duration_ratio: float  # longer / shorter

    @property
    def best(self) -> DuoAlias:
        return max(self.aliases, key=lambda a: a.weight)

    @property
    def p_min(self) -> float:
        return min(a.period for a in self.aliases)


def events_consistent(a: SingleEvent, b: SingleEvent) -> bool:
    """Same depth within ``DEPTH_SIGMA_MAX`` sigma and similar duration.

    Each event's depth uncertainty is ``depth / snr`` plus a
    ``DEPTH_FRACTION_FLOOR`` systematic term: at high SNR the box fit's bias
    (duration and centre on a grid) exceeds the statistical error. Matching
    depth and shape between the two events is the main eclipsing-binary check.
    """

    def sigma_of(e: SingleEvent) -> float:
        return float(np.hypot(e.depth / e.snr, DEPTH_FRACTION_FLOOR * e.depth))

    sigma = np.hypot(sigma_of(a), sigma_of(b))
    if abs(a.depth - b.depth) > DEPTH_SIGMA_MAX * sigma:
        return False
    ratio = max(a.duration, b.duration) / min(a.duration, b.duration)
    return ratio <= DURATION_RATIO_MAX


def duration_likelihood(duration: float, t_central: float) -> float:
    """Relative likelihood of an observed duration given the central duration.

    For a circular orbit with uniform impact parameter, p(D) = D / (Tc sqrt(Tc^2 - D^2))
    for D < Tc. Longer durations need an eccentric orbit and get a Gaussian tail.
    """
    x = duration / t_central
    if x < 1:
        xc = min(x, _X_CLIP)
        return float(xc / np.sqrt(1 - xc**2) / t_central)
    edge = _X_CLIP / np.sqrt(1 - _X_CLIP**2) / t_central
    return float(edge * np.exp(-(((x - 1) / _DURATION_TAIL_WIDTH) ** 2)))


def duo_aliases(time, a: SingleEvent, b: SingleEvent, rho_star: float = 1.0) -> tuple:
    """Allowed periods ``dT/n`` for a pair of events, with normalised weights."""
    time, _ = _prepare(time, np.zeros_like(np.asarray(time, dtype=float)))
    first, second = sorted((a, b), key=lambda e: e.t0)
    delta_t = second.t0 - first.t0
    duration = (first.duration + second.duration) / 2
    cadence = _cadence(time)
    n_max = min(int(delta_t / (_MIN_PERIOD_DURATIONS * duration)), MAX_ALIASES)
    kept = []
    for n in range(1, n_max + 1):
        period = delta_t / n
        if period_is_allowed(time, first.t0, duration, period, cadence, observed=(0, n)):
            t_central = central_duration(period, rho_star)
            weight = period**PERIOD_PRIOR_INDEX * duration_likelihood(duration, t_central)
            kept.append((period, n, weight))
    total = sum(w for *_, w in kept)
    if not kept or not total > 0:
        return tuple(DuoAlias(p, n, 1 / len(kept)) for p, n, _ in kept)
    return tuple(DuoAlias(p, n, w / total) for p, n, w in kept)


def find_duo(time, events, rho_star: float = 1.0) -> DuoPair | None:
    """The strongest pair of consistent events with at least one allowed alias."""
    pairs = sorted(combinations(events, 2), key=lambda ab: -(ab[0].snr + ab[1].snr))
    for a, b in pairs:
        if not events_consistent(a, b):
            continue
        aliases = duo_aliases(time, a, b, rho_star)
        if not aliases:
            continue
        first, second = sorted((a, b), key=lambda e: e.t0)
        return DuoPair(
            first=first,
            second=second,
            delta_t=second.t0 - first.t0,
            aliases=aliases,
            depth_ratio=second.depth / first.depth,
            duration_ratio=max(a.duration, b.duration) / min(a.duration, b.duration),
        )
    return None
