# Quicklook Roadmap: Improving Transit Detection Efficiency for Young Variable Stars

Based on literature from `projects/transit-search/` (*Fernandes et al. 2022*, *Hippke & David 2019*, *Rizzuto et al. 2017*, *Garcia et al. 2024*, *Kunimoto et al. 2021*, *Heller et al. 2019a/2020*) and analysis of the `quicklook` package.

---

## 1. Algorithmic Insights & Best Practices

### A. Stellar Variability Detrending
- **Avoid Savitzky-Golay (SG) Filters:** SG overfits transit profiles, underestimating transit depths and reducing recovery rates.
- **Penalized Splines / Adaptive Knots:** Use `pspline` or `rspline` in `Wotan` with adaptive knot density based on $P_{\text{rot}}$ derived from GLS.
  - $P_{\text{rot}} \le 2\text{ d} \implies \text{max\_splines} = 200$
  - $2 < P_{\text{rot}} \le 10\text{ d} \implies \text{max\_splines} = 200 / P_{\text{rot}}$
  - $P_{\text{rot}} \ge 10\text{ d} \implies \text{max\_splines} = 25$
- **Adaptive Sigma Clipping:** Use $3\sigma$ clipping instead of $2\sigma$ for high-amplitude rotators ($>2\%$ variability) to prevent clipping rotational peaks.
- **Edge Trimming:** Increase `edge_cutoff` to $0.5\text{ days}$ to trim edge-detrending systematics.
- **Sector-Independent Detrending:** Flatten sectors/segments individually prior to stitching into a single light curve.

### B. Search & Vetting Enhancements
- **TLS over BLS:** Default to TLS/GTLS for transit shape and limb-darkening matching.
- **Advanced Vetting Metrics:**
  - Individual transit depth variance ratio ($\sigma_{\text{depth}} / \bar{\text{depth}} < 1.5$)
  - Duration consistency ($0.75 < D_{\text{obs}} / D_{\text{exp}} < 1.33$)
  - Odd-even depth mismatch significance ($\sigma$)
  - Secondary eclipse search over phases 0.1 to 0.9 (flag if secondary SDE > 3 or depth > 10% of primary)
  - Rotation harmonic flags ($P_{\text{orb}} \approx P_{\text{rot}}$ within 20% tolerance)
- **Iterative Search:** Mask primary candidate transits and re-run TLS with `minimum_transits = 1` for outer/smaller planets.

---

## 2. Implementation Tasks

- [ ] **Task 1: Evaluate `detrend_and_search_tls` vs. Existing Pipelines**
  - Compare integrated `detrend_and_search_tls` helper against current `TessQuickLook` workflow in `tql.py`.
- [ ] **Task 2: Implement Rotation-Aware Window Tuning**
  - Integrate GLS periodogram $P_{\text{rot}}$ into default window selection when `window_length` is default/unspecified in `flatten_raw_lc()`.
  - Dynamically set `max_splines` for spline methods and adaptive window length ($0.5\text{ d}$ for $P_{\text{rot}} \le 2\text{ d}$, $1.0\text{ d}$ for $P_{\text{rot}} > 2\text{ d}$).
- [ ] **Task 3: Sector-Independent Detrending Workflow**
  - Update multi-sector light curve handling to detrend individual sectors before stitching.
- [ ] **Task 4: Implement Advanced Vetting Metrics & Flags**
  - Add depth variance ratio, duration consistency, odd-even mismatch, and secondary eclipse metrics to `measure.py` / `tql.py`.
- [ ] **Task 5: Add Iterative Transit Search**
  - Support automatic transit masking and secondary planet search execution.
