"""Descriptive variance-ratio efficiency scores and cross-industry dispersion.

Run .venv/bin/python analyze_industry_efficiency.py. Uses frozen local data.
E = 1 - mean_q(abs(VR(q)-1)), q in (2,5,10,20). This is a declared
descriptive score, not a universal efficiency index or probability.
"""
from pathlib import Path
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
from arch.bootstrap import StationaryBootstrap
from arch.unitroot import VarianceRatio
from analyze_variance_ratio import direct_vr, holm

HORIZONS = (2, 5, 10, 20)
NAMES = {"NoDur": "Consumer nondurables", "Durbl": "Consumer durables",
         "Manuf": "Manufacturing", "Enrgy": "Energy", "Chems": "Chemicals",
         "BusEq": "Computers, software & electronics", "Telcm": "Telecommunications",
         "Utils": "Utilities", "Shops": "Wholesale, retail & selected services",
         "Hlth": "Healthcare", "Money": "Finance", "Other": "Other industries"}
OUT = ROOT / "results/industry_efficiency"


def variance_ratios(r):
    """All columns and horizons; log-return input, overlapping and debiased."""
    n, k = r.shape
    assert n > max(HORIZONS) and np.isfinite(r).all()
    mu = r.mean(axis=0)
    variance1 = np.sum((r-mu)**2, axis=0)/(n-1)
    assert np.all(variance1 > 0)
    y = np.vstack([np.zeros(k), np.cumsum(r, axis=0)])
    values = []
    for q in HORIZONS:
        varianceq = np.sum((y[q:]-y[:-q]-q*mu)**2, axis=0)/(q*(n-q+1)*(1-q/n))
        values.append(varianceq/variance1)
    return np.array(values).T


def efficiency(vr):
    return 1-np.mean(np.abs(vr-1), axis=-1)


def dispersion_bounds(lower, upper):
    """Exact box maximum for 12 industries; convex minimization for the minimum.

    SD is convex, so its maximum over a box occurs at a vertex. A common
    point in all score intervals implies the minimum possible SD is zero.
    """
    lower, upper = np.asarray(lower), np.asarray(upper)
    k = len(lower)
    assert k <= 12 and np.all(lower <= upper)
    if max(lower) <= min(upper):
        minimum = 0.0
    else:
        objective = lambda mean: np.mean((np.clip(mean, lower, upper)-mean)**2)
        result = minimize_scalar(objective, bounds=(float(min(lower)), float(max(upper))), method="bounded",
                                 options={"xatol": 1e-12})
        assert result.success
        minimum = float(np.sqrt(result.fun))
    mask = ((np.arange(2**k)[:, None] >> np.arange(k)) & 1).astype(bool)
    corners = np.where(mask, upper, lower)
    return minimum, float(np.max(np.std(corners, axis=1, ddof=0)))


def bootstrap_bands(r, point, block_length, reps, seed):
    """Resample whole daily vectors to retain cross-industry dependence.

    Joint bands are constructed for signed VR parameters before propagating
    them through absolute values and SD. This avoids using a naive percentile
    SD interval to test zero at a nonnegative boundary.
    """
    bs = StationaryBootstrap(block_length, r, seed=seed)
    draws = np.array([variance_ratios(sample[0][0]) for sample in bs.bootstrap(reps)])
    se = draws.std(axis=0, ddof=1)
    assert np.all(se > 0)
    maximum_errors = np.max(np.abs((draws-point)/se), axis=(1, 2))
    # 97.5% simultaneous coverage within each period, Bonferroni over 2 periods.
    critical = float(np.quantile(maximum_errors, .975, method="higher"))
    vr_lower, vr_upper = np.maximum(0, point-critical*se), point+critical*se
    absolute_lower = np.maximum(0, np.maximum(vr_lower-1, 1-vr_upper))
    absolute_upper = np.maximum(np.abs(vr_lower-1), np.abs(vr_upper-1))
    score_lower, score_upper = 1-absolute_upper.mean(axis=1), 1-absolute_lower.mean(axis=1)
    sd_bounds = dispersion_bounds(score_lower, score_upper)
    common_lower, common_upper = float(score_lower.max()), float(score_upper.min())
    return {"mean_block_length": block_length, "replicates": reps, "seed": seed,
            "within_period_joint_coverage_target": .975,
            "two_period_joint_coverage_target": .95,
            "critical_max_standardized_error": critical,
            "vr_lower": vr_lower.tolist(), "vr_upper": vr_upper.tolist(),
            "score_lower": score_lower.tolist(), "score_upper": score_upper.tolist(),
            "population_sd_bounds": list(sd_bounds),
            "common_score_interval": [common_lower, common_upper] if common_lower <= common_upper else None,
            "equal_scores_excluded_by_joint_bands": common_lower > common_upper}


def validate_helpers():
    assert np.all(efficiency(np.ones((3, 4))) == 1)
    np.testing.assert_allclose(efficiency(np.full((3, 4), .8)), .8)
    np.testing.assert_allclose(efficiency(np.full((3, 4), 1.2)), .8)
    np.testing.assert_allclose(dispersion_bounds([0, 0], [1, 1]), [0, .5])
    np.testing.assert_allclose(dispersion_bounds([0, .8], [.2, 1]), [.3, .5], atol=1e-8)
    assert np.std(np.full(12, .8), ddof=0) < 1e-14


def analyze(number):
    source_path = ROOT / f"data/industry_portfolios/{number}_industry_value_weighted_study_windows.json"
    data = json.loads(source_path.read_text())
    columns = data["metadata"]["columns"]
    result = {"industry_count": number, "source_metadata": data["metadata"],
              "prepared_data_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(), "periods": {}}
    full_tests = []
    for period_index, (period, rows) in enumerate(data["data"].items()):
        raw = np.array([[row["returns_pct"][c] for c in columns] for row in rows])
        assert np.isfinite(raw).all() and not np.isin(raw, [-99.99, -999]).any()
        assert (raw > -100).all()
        r = np.log1p(raw/100)
        vr = variance_ratios(r)
        scores = efficiency(vr)
        tests = []
        for i, c in enumerate(columns):
            y = np.r_[0, np.cumsum(r[:, i])]
            for j, q in enumerate(HORIZONS):
                lm = VarianceRatio(y, lags=q, trend="c", robust=True, overlap=True, debiased=True)
                np.testing.assert_allclose(vr[i, j], lm.vr, rtol=1e-10, atol=1e-12)
                np.testing.assert_allclose([lm.vr, lm.stat, lm.pvalue], direct_vr(y, q), rtol=1e-10, atol=1e-12)
                tests.append({"industry": c, "period": period, "q": q, "vr": float(lm.vr),
                              "robust_z": float(lm.stat), "robust_p": float(lm.pvalue)})
        for row, p in zip(tests, holm([t["robust_p"] for t in tests])):
            row["holm_p_all_industries_and_horizons_in_period"] = p
        full_tests.extend(tests)
        records = [{"code": c, "name": NAMES.get(c, c) if number == 12 else c,
                    "efficiency_score": float(scores[i]), "inefficiency_distance": float(1-scores[i]),
                    "vr_by_horizon": dict(zip(map(str, HORIZONS), map(float, vr[i])))}
                   for i, c in enumerate(columns)]
        summary = {"daily_returns_per_industry": len(rows), "first_return": rows[0]["date"],
                   "last_return": rows[-1]["date"], "mean_efficiency_score": float(scores.mean()),
                   "population_sd": float(scores.std(ddof=0)), "sample_sd": float(scores.std(ddof=1)),
                   "population_sd_is_exactly_zero": bool(scores.std(ddof=0) == 0),
                   "min_efficiency_score": float(scores.min()), "max_efficiency_score": float(scores.max()),
                   "single_horizon_population_sds": dict(zip(map(str, HORIZONS), map(float, np.std(1-np.abs(vr-1), axis=0))))}
        result["periods"][period] = {"summary": summary, "industries": records, "individual_tests": tests}
        print(f"{number} industries, {period}: SD={summary['population_sd']:.6f}; mean={scores.mean():.6f}", flush=True)
        if number == 12:
            bands = []
            for length in (40, 20, 60):
                reps = 4999 if length == 40 else 2499
                bands.append(bootstrap_bands(r, vr, length, reps, 20260917+1000*period_index+length))
                print(f"  bootstrap block={length}: SD bounds={bands[-1]['population_sd_bounds']}; equality excluded={bands[-1]['equal_scores_excluded_by_joint_bands']}", flush=True)
            result["periods"][period]["uncertainty_main"] = bands[0]
            result["periods"][period]["uncertainty_block_length_sensitivity"] = bands[1:]
    for row, p in zip(full_tests, holm([t["robust_p"] for t in full_tests])):
        row["holm_p_all_tests_both_periods_this_dataset"] = p
    return result


def report(results):
    main, detail = results["datasets"]
    periods = list(main["periods"])
    lines = ["# US Market Standard Deviation by Sectors", "",
             "## Definition and interpretation", "",
             "This analysis uses all 12 Kenneth French broad industry portfolios as its primary panel and all 49 narrower industries as a descriptive supplement. They are SIC-based portfolios of US stocks, not the 11 GICS sectors. The two taxonomies are analyzed separately, never pooled.", "",
             "Define D_i = (1/4) * sum over q in {2,5,10,20} of |VR_i(q) − 1|, and the efficiency score E_i = 1 − D_i. Higher E means closer to the random-walk variance restriction. E = 1 is the benchmark. This transparent, study-specific score gives equal weight to each horizon; it is not a standardized market-efficiency index, probability, or percentage. It can be negative; no clipping is applied. Opposite-signed VR deviations of equal size receive the same score. Averages can hide differences across horizons.", "",
             "Cross-industry population SD = sqrt(sum_i (E_i − mean(E))² / K), where K=12 (or 49 in the supplement). All industries in each chosen classification are included. The alternative sample SD with denominator K−1 is also reported. Both SDs are calculated from full precision scores before rounding. SD(E) equals SD(D).", "",
             "SD = 0 means identical estimated scores. Equal scores below 1 mean equal estimated average distances from the VR benchmark, not proof of equal inefficiency. Nonzero estimated SD may arise from estimation error even when population scores are identical. Neither result proves or disproves the full efficient-market hypothesis.", "",
             "## Primary results: 12 industries", "",
             "| Industry | E: 2008–2010 | E: 2020–2021 |", "|---|---:|---:|"]
    for left, right in zip(main["periods"][periods[0]]["industries"], main["periods"][periods[1]]["industries"]):
        lines.append(f'| {left["name"]} ({left["code"]}) | {left["efficiency_score"]:.6f} | {right["efficiency_score"]:.6f} |')
    for key, label in [("mean_efficiency_score", "Mean score"), ("population_sd", "Population SD (K denominator)"), ("sample_sd", "Sample SD (K−1 denominator)")]:
        values = [main["periods"][p]["summary"][key] for p in periods]
        lines.append(f"| **{label}** | **{values[0]:.6f}** | **{values[1]:.6f}** |")
    lines += ["", "The population SD is nonzero in both periods: the estimated industry scores are not identical. Comparing these two SD point estimates is descriptive; no formal test of a change in dispersion between periods is claimed.", "", "![Industry efficiency scores](efficiency_scores.png)", "",
              "## Sampling uncertainty and equality", "",
              "To avoid treating nonzero estimated SD as a significant difference, stationary block bootstrap draws resample the whole daily industry-return vector. This retains contemporaneous cross-industry dependence and dependence within sampled blocks. The primary mean block length is 40 trading days with 4,999 replicates; 20 and 60 days with 2,499 replicates each are sensitivity checks. All random seeds are saved.", "",
              "For each period, the maximum absolute standardized bootstrap VR estimation error across 12 industries × 4 horizons determines joint 97.5% bands for the signed VR parameters. Bootstrap standard deviations supply the fixed scales; the bootstrap errors are centered at the original estimates. Bonferroni across the two periods targets approximate 95% simultaneous coverage. Bounds are intersected with the nonnegative VR parameter space. The bands are then propagated through the absolute-distance score. This handles the nondifferentiability at VR=1 without directly bootstrapping the absolute-value score for inference.", "",
              "The minimum possible population SD within the score bounds is zero if the score intervals share a common point; otherwise it is computed by convex minimization. The maximum is obtained over all 2^12 interval corners. These propagated bounds are conservative. Equality is excluded only if there is no score common to all industry intervals. A naive percentile interval for a nonnegative SD is not used to test zero.", "",
              "| Period | Estimated population SD | Propagated SD bounds | Equal scores excluded? |", "|---|---:|---|---|"]
    for p in periods:
        obj = main["periods"][p]
        band = obj["uncertainty_main"]
        lo, hi = band["population_sd_bounds"]
        lines.append(f'| {p} | {obj["summary"]["population_sd"]:.6f} | [{lo:.6f}, {hi:.6f}] | {"Yes" if band["equal_scores_excluded_by_joint_bands"] else "No"} |')
    lines += ["", "These are approximate bootstrap bounds under weak-dependence and stationarity assumptions, not exact finite-sample confidence guarantees. Crisis windows mix volatile regimes; small samples, block choice, and the conservativeness of band propagation limit inference. Failure to exclude equality is not evidence that equality is true. Individual robust VR tests and their multiple-testing adjustments are recorded separately in JSON; their p-values are not efficiency scores.", "", "### Block-length sensitivity", "",
              "| Period | Mean block length | Propagated SD bounds | Equal scores excluded? |", "|---|---:|---|---|"]
    for p in periods:
        obj = main["periods"][p]
        for band in [obj["uncertainty_main"]]+obj["uncertainty_block_length_sensitivity"]:
            lo, hi = band["population_sd_bounds"]
            lines.append(f'| {p} | {band["mean_block_length"]} | [{lo:.6f}, {hi:.6f}] | {"Yes" if band["equal_scores_excluded_by_joint_bands"] else "No"} |')
    lines += ["", "## Supplement: all 49 industries", "", "The following is descriptive; the bootstrap equality analysis above applies to the 12-industry primary panel only. Codes follow the [49-industry definitions](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_49_ind_port.html).", "",
              "| Industry code | E: 2008–2010 | E: 2020–2021 |", "|---|---:|---:|"]
    for left, right in zip(detail["periods"][periods[0]]["industries"], detail["periods"][periods[1]]["industries"]):
        lines.append(f'| {left["code"]} | {left["efficiency_score"]:.6f} | {right["efficiency_score"]:.6f} |')
    for key, label in [("mean_efficiency_score", "Mean score"), ("population_sd", "Population SD"), ("sample_sd", "Sample SD")]:
        values = [detail["periods"][p]["summary"][key] for p in periods]
        lines.append(f"| **{label}** | **{values[0]:.6f}** | **{values[1]:.6f}** |")
    lines += ["", "## Data, method, and validation", "",
              "Value-weighted daily industry returns come from the frozen July 2026 CRSP vintage in the Kenneth French library. Each industry has 757 returns for 2008–2010 and 505 for 2020–2021. Raw percentage simple returns are transformed with log(1+R/100); log wealth starts at zero and accumulates these returns. The industry analysis uses total returns and is not directly interchangeable with the earlier S&P 500 price-index analysis.", "",
              "Lo–MacKinlay variance ratios use nonzero drift, overlapping returns, finite-sample debiasing, and heteroskedasticity-robust normal inference. Robust p-values are Holm-adjusted within each period across all industries/horizons, and additionally across both periods; the 12- and 49-industry analyses have separate testing families (96 and 392 tests).", "",
              "All 488 observed VR estimates match arch.unitroot.VarianceRatio, and every robust Z and p-value is checked against the independent scalar implementation. Validation also covers symmetric distances around VR=1, zero dispersion for identical scores, and known minimum/maximum dispersion bounds. Full precision estimates, source hashes, individual tests, and bootstrap settings are in industry_efficiency_results.json.", "",
              "Reproduce with `.venv/bin/python analyze_industry_efficiency.py` from the project root after installing the saved requirements. No network requests are made.", "",
              "Sources: [12-industry construction](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_12_ind_port.html); [49-industry construction](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_49_ind_port.html); [variance-ratio implementation](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.VarianceRatio.html); [stationary bootstrap implementation](https://arch.readthedocs.io/en/latest/bootstrap/generated/arch.bootstrap.StationaryBootstrap.html); [Politis and Romano (1994)](https://doi.org/10.1080/01621459.1994.10476870). The particular four-horizon score and propagation of joint bands are choices documented for this study."]
    (OUT / "us_market_standard_deviation_by_sectors.md").write_text("\n".join(lines)+"\n")


def plot(main):
    fig, ax = plt.subplots(figsize=(10.5, 7))
    labels = [r["name"] for r in next(iter(main["periods"].values()))["industries"]]
    y = np.arange(len(labels))
    for i, (period, result) in enumerate(main["periods"].items()):
        scores = [r["efficiency_score"] for r in result["industries"]]
        ax.barh(y+(i-.5)*.34, scores, height=.3, color=["#237a95", "#c55d31"][i], label=period)
    ax.set(yticks=y, yticklabels=labels, xlabel="Descriptive efficiency score (higher = closer to VR benchmark)", xlim=(0, 1.03))
    ax.invert_yaxis()
    ax.axvline(1, color="#566273", linestyle="--", linewidth=1)
    ax.set_title("US industry efficiency scores", fontsize=18, loc="left", fontweight="bold", pad=16)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=.15)
    ax.set_axisbelow(True)
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc="lower center", bbox_to_anchor=(.65, .09), ncol=2)
    fig.subplots_adjust(left=.34, right=.96, top=.90, bottom=.24)
    fig.text(.06, .025, "E = 1 − mean |VR(q) − 1|, q = 2, 5, 10, 20 trading days. Scores are not probabilities.\nKenneth French 12 value-weighted industry portfolios · See report for sampling uncertainty.", fontsize=9, color="#4a5565")
    fig.savefig(OUT / "efficiency_scores.png", dpi=180)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    validate_helpers()
    results = {"score_formula": "E = 1 - mean_q(abs(VR(q)-1))", "horizons": HORIZONS,
               "primary_industry_count": 12, "bootstrap_main_mean_block_length": 40,
               "bootstrap_main_replicates": 4999, "datasets": [analyze(12), analyze(49)]}
    (OUT / "industry_efficiency_results.json").write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")
    plot(results["datasets"][0])
    report(results)
    print("All calculations and validation completed.", flush=True)


if __name__ == "__main__":
    main()
