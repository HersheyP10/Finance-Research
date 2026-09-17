"""Extend the original S&P 500 variance-ratio analysis back to 1950.

Reproduce: .venv/bin/python analyze_historical_variance_ratio.py
All input data are frozen locally; no network requests are made.
"""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib
import json
import math
import os

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from arch.bootstrap import StationaryBootstrap
from arch.unitroot import VarianceRatio
from analyze_variance_ratio import direct_vr, holm

ROOT = Path(__file__).resolve().parent
HORIZONS = (2, 5, 10, 20)
CUTOFF = "2026-09-15"  # Same completed-session cutoff as the original run.
PERIODS = (
    ("1950s", "1950-01-01", "1959-12-31", "historical_test"),
    ("End of Cold War", "1988-01-01", "1990-12-31", "historical_test"),
    ("Dot-com run-up and bust", "1995-01-01", "2002-12-31", "historical_test"),
    ("Dot-com crash only", "2000-01-01", "2002-12-31", "dotcom_sensitivity"),
    ("2008–2010 baseline", "2008-01-01", "2010-12-31", "original_baseline"),
    ("2020–2021 baseline", "2020-01-01", "2021-12-31", "original_baseline"),
)
RAW_PATH = ROOT / "data/yahoo_gspc_1950_2026_raw.json"
OLD_PATH = ROOT / "data/yahoo_gspc_raw.json"
OUT = ROOT / "results/historical_variance_ratio"


def load_prices(path):
    content = path.read_bytes()
    chart = json.loads(content)["chart"]
    assert chart["error"] is None
    series = chart["result"][0]
    assert series["meta"]["symbol"] == "^GSPC"
    timezone = ZoneInfo(series["meta"]["exchangeTimezoneName"])
    dates = [datetime.fromtimestamp(t, timezone).date().isoformat() for t in series["timestamp"]]
    closes = series["indicators"]["adjclose"][0]["adjclose"]
    assert dates == sorted(set(dates)) and len(dates) == len(closes)
    prices = {day: float(price) for day, price in zip(dates, closes)
              if day <= CUTOFF and price is not None}
    assert len(prices) and min(prices.values()) > 0 and np.isfinite(list(prices.values())).all()
    return content, prices, dates


def run_test(log_prices, q):
    test = VarianceRatio(log_prices, lags=q, trend="c", robust=True, overlap=True, debiased=True)
    vr, z, p = direct_vr(log_prices, q)
    np.testing.assert_allclose([test.vr, test.stat, test.pvalue], [vr, z, p], rtol=1e-10, atol=1e-12)
    return {"horizon_trading_days": q, "variance_ratio": float(test.vr),
            "robust_z": float(test.stat), "robust_p": float(test.pvalue),
            "absolute_distance_from_random_walk": float(abs(test.vr-1))}


def mean_absolute_distance(r):
    y = np.r_[0.0, np.cumsum(r)]
    return float(np.mean([abs(VarianceRatio(y, lags=q, trend="c", robust=True,
                                            overlap=True, debiased=True).vr-1)
                          for q in HORIZONS]))


def bootstrap_distances(r, block_length, reps, seed):
    sampler = StationaryBootstrap(block_length, np.asarray(r).reshape(-1, 1), seed=seed)
    return np.fromiter((mean_absolute_distance(sample[0][0][:, 0])
                        for sample in sampler.bootstrap(reps)), dtype=float, count=reps)


def compare_historical_periods(series):
    """Compare each requested old window with the mean of both original-run windows."""
    historical = ["1950s", "End of Cold War", "Dot-com run-up and bust", "Dot-com crash only"]
    baselines = ["2008–2010 baseline", "2020–2021 baseline"]
    point = {name: mean_absolute_distance(r) for name, r in series.items()}
    output = {length: [] for length in (20, 40, 60)}
    for length in (40, 20, 60):
        reps, commonseed = (4999 if length == 40 else 2499), 20300917+length
        draws = {name: bootstrap_distances(series[name], length, reps, commonseed+i*100)
                 for i, name in enumerate(historical+baselines)}
        for name in historical:
            baseline_draws = .5*(draws[baselines[0]]+draws[baselines[1]])
            delta = point[name]-.5*(point[baselines[0]]+point[baselines[1]])
            errors = draws[name]-baseline_draws-delta
            se = float(np.std(errors, ddof=1))
            assert se > 0
            z = delta/se
            p = math.erfc(abs(z)/math.sqrt(2))
            # 95% Bonferroni simultaneous Wald bands for four comparisons.
            critical = float(norm.ppf(1-.05/(2*len(historical))))
            output[length].append({"historical_period":name,
                "historical_mean_abs_vr_distance":point[name],
                "baseline_periods":baselines,
                "baseline_mean_abs_vr_distance":float(.5*(point[baselines[0]]+point[baselines[1]])),
                "difference_historical_minus_baseline":float(delta),
                "stationary_bootstrap_standard_error":se,
                "wald_z":float(z),"two_sided_p":p,
                "holm_p_within_four_period_comparisons":None,
                "familywise_95pct_simultaneous_interval":[float(delta-critical*se),float(delta+critical*se)],
                "earlier_period_more_distant_point_estimate":bool(delta>0),
                "reject_equal_mean_distances_holm_5pct":False,
                "mean_block_length":length,"replicates_per_series":reps,
                "critical_value_familywise_95pct":critical,"seed":commonseed})
        for row, adj in zip(output[length], holm([r["two_sided_p"] for r in output[length]])):
            row["holm_p_within_four_period_comparisons"]=adj
            row["reject_equal_mean_distances_holm_5pct"]=adj<.05
    return {"quantity_compared":"mean absolute distance from VR=1 across q=2,5,10,20",
            "difference_sign":"historical minus equally weighted mean of 2008–2010 and 2020–2021 baseline distances; positive means the earlier-window point estimate is more distant.",
            "method":"Independent stationary block bootstrap of each non-overlapping period return series. Each daily block contains the whole return path for that period; the two baseline distance estimates are averaged equally.",
            "inference":"Bootstrap standard error and two-sided Wald p-value; Holm adjustment across four comparisons within each block-length specification. 95% Bonferroni simultaneous normal intervals cover the four comparisons.",
            "assumptions":"Approximate inference assumes weak dependence within each period and independent non-overlapping period samples. Structural breaks, small samples, and block choice can affect validity.",
            "period_comparisons_by_mean_block_length":{str(k):v for k,v in output.items()}}


def main():
    source_bytes, all_prices, all_dates = load_prices(RAW_PATH)
    _, original_prices, _ = load_prices(OLD_PATH)
    common = sorted(set(original_prices).intersection(all_prices))
    assert set(original_prices) <= set(all_prices)
    max_abs = max(abs(original_prices[d]-all_prices[d]) for d in common)
    assert max_abs == 0.0, f"Yahoo's completed-session history revised in the overlapping sample (max diff {max_abs})"
    assert "1950-01-03" in all_prices and CUTOFF in all_prices

    results = {"symbol": "^GSPC", "source": "Yahoo Finance historical S&P 500 (^GSPC)",
               "history_source_url": "https://finance.yahoo.com/quote/%5EGSPC/history/",
    "chart_endpoint": "https://query2.finance.yahoo.com/v8/finance/chart/%5EGSPC?period1=-632620800&period2=1789603200&interval=1d",
               "retrieved_utc_date": "2026-09-16", "last_complete_session": CUTOFF,
               "intraday_session_excluded": "2026-09-16",
               "historical_data_sha256": hashlib.sha256(source_bytes).hexdigest(),
               "original_run_data_sha256": hashlib.sha256(OLD_PATH.read_bytes()).hexdigest(),
               "overlapping_completed_sessions_exactly_match": True,
               "overlap_sessions": len(common), "maximum_absolute_price_difference": max_abs,
               "price_series": "Yahoo adjusted close (equal to close throughout the completed-session overlap)",
               "return_type": "Daily close-to-close log price-index returns; not dividend-reinvested sector returns",
               "method": "Same as original run: Lo–MacKinlay VR, drift, overlapping returns, debiased variances, heteroskedasticity-robust two-sided normal inference",
               "horizons_trading_days": HORIZONS,
               "multiple_testing": "Holm adjustment across four horizons within each period; also across all 24 estimates.",
               "periods": []}
    return_series = {}
    for label, start, end, role in PERIODS:
        in_period = sorted(d for d in all_prices if start <= d <= end)
        base_dates = [d for d in all_prices if d < start]
        assert len(in_period) > max(HORIZONS) and base_dates, f"Insufficient data: {label}"
        dates = [max(base_dates), *in_period]
        logp = np.log(np.array([all_prices[d] for d in dates], dtype=float))
        return_series[label] = np.diff(logp)
        tests = [run_test(logp, q) for q in HORIZONS]
        per_period = holm([r["robust_p"] for r in tests])
        for row, adjusted in zip(tests, per_period):
            row["holm_p_within_period"] = adjusted
            row["reject_5pct_unadjusted"] = row["robust_p"] < .05
            row["reject_5pct_holm_within_period"] = adjusted < .05
        stats = {"period": label, "requested_start": start, "requested_end": end,
                 "first_return_date": dates[1], "last_return_date": dates[-1],
                 "daily_returns": len(in_period),
                 "mean_absolute_vr_distance": float(np.mean([r["absolute_distance_from_random_walk"] for r in tests])),
                 "holm_rejections_within_period": [r["horizon_trading_days"] for r in tests if r["reject_5pct_holm_within_period"]],
                 "unadjusted_rejections": [r["horizon_trading_days"] for r in tests if r["reject_5pct_unadjusted"]]}
        results["periods"].append({"summary": stats, "tests": tests, "role": role})
    all_tests = [row for period in results["periods"] for row in period["tests"]]
    for row, adjusted in zip(all_tests, holm([r["robust_p"] for r in all_tests])):
        row["holm_p_all_24_tests"] = adjusted
        row["reject_5pct_holm_all_periods"] = adjusted < .05
    results["conclusion_limit"] = "Historical-versus-baseline contrasts use an exploratory stationary block bootstrap. The point estimates are lower in all historical windows, but evidence is sensitive to block length, and inference assumes weak dependence and approximate stationarity within each period."
    results["historical_vs_baseline_comparisons"] = compare_historical_periods(return_series)
    results["validation"] = "All 24 library VR, robust Z, and p-value estimates match an independent scalar implementation (rtol=1e-10, atol=1e-12); extended and original Yahoo adjusted closes match exactly on all overlapping completed sessions."
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "historical_variance_ratio_results.json").write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")
    report(results)
    chart(results)
    for period in results["periods"]:
        s=period["summary"]
        print(f'{s["period"]}: n={s["daily_returns"]}, mean|VR-1|={s["mean_absolute_vr_distance"]:.6f}, Holm rejects={s["holm_rejections_within_period"]}', flush=True)
    print(results["validation"], flush=True)


def report(results):
    periods=results["periods"]
    lines=["# Earlier US market variance ratio tests", "",
           "This extends the original ^GSPC analysis using Yahoo Finance historical adjusted-close data. The expanded series begins January 3, 1950. The original and expanded series match to the cent on every overlapping completed trading session (4,706 observations through September 15, 2026). September 16, 2026 is omitted because its quote was still in progress in the extended download.", "",
           "## Method", "",
           "This uses the original run's Lo–MacKinlay variance ratio test: daily log price returns, a drift, overlapping q-day returns for q=2, 5, 10, 20 trading days, finite-sample debiasing, heteroskedasticity-robust inference, and two-sided p-values. The null is VR(q)=1. VR below 1 indicates negative aggregate serial correlation over that horizon; VR above 1 indicates positive serial correlation. The robust test uses large-sample normal inference, as in the [arch VarianceRatio documentation](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.VarianceRatio.html).", "",
           "Holm p-values correct for four horizons within each period. A separate adjustment covers all 24 estimates across the six windows. Period windows are inclusive calendar dates; returns begin on the first trading session and end on the last session in each window. The period-specific mean absolute distance mean(|VR(q)-1|) is a descriptive summary; direct between-period comparisons appear below.", "",
           "## Results", "",
           "| Period | Daily returns | Mean abs(VR−1) | q=2 VR (Holm p) | q=5 VR (Holm p) | q=10 VR (Holm p) | q=20 VR (Holm p) | Reject any horizon after within-window Holm? |", "|---|---:|---:|---:|---:|---:|---:|---|"]
    for item in periods:
        s=item["summary"]; tests=item["tests"]
        vals=[f'{r["variance_ratio"]:.4f} ({r["holm_p_within_period"]:.4f})' for r in tests]
        lines.append(f'| {s["period"]} | {s["daily_returns"]:,} | {s["mean_absolute_vr_distance"]:.4f} | '+" | ".join(vals)+f' | {"Yes: "+", ".join(map(str,s["holm_rejections_within_period"]))+" days" if s["holm_rejections_within_period"] else "No"} |')
    lines += ["", "Within-window adjusted p-values below 0.05 reject that horizon's random-walk restriction. For the 1950s, q=2 also survives Holm correction across all 24 tests (global adjusted p=0.0194); it is the only global rejection.", "",
              "![S&P 500 variance ratios by period](variance_ratios.png)", "",
              "## How to read this as a test of historical inefficiency", "",
              "A rejection says this series departs from the particular random-walk variance restriction at the chosen horizon, under this test and sample. It does not establish all forms of market inefficiency, nor does failure to reject prove efficiency. The data are a broad S&P 500 price index, not the full US equity market or a dividend-reinvested total-return portfolio.", "",
              "1950–1959, 1988–1990, the full 1995–2002 dot-com run-up and bust, and the nested 2000–2002 crash window show the requested historical tests. The 2008–2010 and 2020–2021 periods repeat the original analysis as comparison baselines. The run-up-and-bust window spans more than twice as many years as the 2000–2002 window; their estimates do not isolate a causal bubble effect.", "",
              "### Direct comparison with the 2008–2021 baseline windows", "",
              "For each historical window, the table compares its mean absolute VR distance with the equally weighted average for 2008–2010 and 2020–2021 (0.2421). Negative differences mean the historical point estimate is closer to VR=1. The primary inference uses stationary bootstrap blocks averaging 40 trading days, 4,999 resamples per period, Holm correction over these four historical comparisons, and family-wise 95% intervals.", "",
              "| Earlier window | Historical distance | Difference vs baseline average | Holm p | Family-wise 95% interval for difference |", "|---|---:|---:|---:|---:|"]
    primary=results["historical_vs_baseline_comparisons"]["period_comparisons_by_mean_block_length"]["40"]
    for row in primary:
        ci=row["familywise_95pct_simultaneous_interval"]
        lines.append(f'| {row["historical_period"]} | {row["historical_mean_abs_vr_distance"]:.4f} | {row["difference_historical_minus_baseline"]:.4f} | {row["holm_p_within_four_period_comparisons"]:.4f} | [{ci[0]:.4f}, {ci[1]:.4f}] |')
    lines += ["", "No historical-versus-baseline contrast is significant after Holm correction at the primary 40-day block length. Sensitivity checks using 20-day blocks also find no significant contrast; with 60-day blocks, only the 1950s comparison crosses 5% (Holm p=0.0484). Since that result depends on block length, the bootstrap evidence does not establish that any older window was less inefficient than the recent baseline. The point estimates are lower in all four windows, but they are not conclusive evidence of a difference.", "",
              "This comparison is exploratory: it assumes weak dependence within each period and independent non-overlapping period samples. Structural breaks, market regime changes, and the block-length choice can affect the bootstrap inference.", "",
              "## Validation and provenance", "",
              results["validation"], "The raw extended response, exact source hash, cutoff, requested bounds, all full-precision test values, period roles, Holm adjustments, and bootstrap comparisons are recorded in `historical_variance_ratio_results.json`.", "",
              "Sources: [Yahoo Finance S&P 500 (^GSPC) history](https://finance.yahoo.com/quote/%5EGSPC/history/); [arch VarianceRatio documentation](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.VarianceRatio.html). Reproduce with `.venv/bin/python analyze_historical_variance_ratio.py`; no network access is needed."]
    (OUT/"report.md").write_text("\n".join(lines)+"\n")


def chart(results):
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    colors=["#237a95","#c55d31","#7a5cab","#599861","#9b653d","#7988a9"]
    qs=list(HORIZONS)
    for i,item in enumerate(results["periods"]):
        label=item["summary"]["period"]
        vals=[r["variance_ratio"] for r in item["tests"]]
        axes[0].plot(qs,vals,marker="o",linewidth=2,markersize=5,label=label,color=colors[i])
        dist=[r["absolute_distance_from_random_walk"] for r in item["tests"]]
        axes[1].plot(qs,dist,marker="o",linewidth=2,markersize=5,label=label,color=colors[i])
    axes[0].axhline(1,color="#555f6e",linestyle="--",linewidth=1)
    axes[0].set_ylabel("Variance ratio")
    axes[1].set_ylabel("|VR − 1|")
    axes[1].set_xlabel("Holding horizon (trading days)")
    axes[1].set_xticks(qs)
    axes[0].set_title("S&P 500 variance ratios across historical periods",loc="left",fontweight="bold",fontsize=15,pad=12)
    for ax in axes:
        ax.grid(axis="y",alpha=.18)
        ax.spines[["top","right"]].set_visible(False)
    fig.legend(*axes[0].get_legend_handles_labels(),loc="lower center",bbox_to_anchor=(.52,.015),ncol=3,frameon=False,fontsize=9)
    fig.subplots_adjust(left=.11,right=.98,top=.91,bottom=.19,hspace=.16)
    fig.text(.06,.002,"Source: Yahoo Finance ^GSPC daily adjusted close; completed sessions through Sep 15, 2026. See report for test details.",fontsize=8,color="#4a5565")
    fig.savefig(OUT/"variance_ratios.png",dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
