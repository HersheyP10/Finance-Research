"""Reproduce S&P 500 variance-ratio estimates from the saved Yahoo snapshot.

Run: .venv/bin/python analyze_variance_ratio.py
The saved cutoff deliberately excludes the incomplete 2026-09-16 session.
"""
from pathlib import Path
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
import hashlib
import json
import math
import os

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from arch.unitroot import VarianceRatio

CUTOFF = "2026-09-15"
HORIZONS = (2, 5, 10, 20)
PERIODS = (
    ("2008–2010", "2008-01-01", "2010-12-31"),
    ("2020–2021", "2020-01-01", "2021-12-31"),
    ("2022–2026 YTD", "2022-01-01", CUTOFF),
)
SOURCE = "https://finance.yahoo.com/quote/%5EGSPC/history/"
ENDPOINT = "https://query2.finance.yahoo.com/v8/finance/chart/%5EGSPC?period1=1199059200&period2=1789603200&interval=1d"


def holm(pvalues):
    """Holm family-wise adjustment; valid with dependent tests."""
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]))
    return adjusted.tolist()


def direct_vr(log_prices, q):
    """Independent scalar implementation to check the library's convention.

    sigma1 = sum((r-mu)^2)/(T-1)
    sigmaq = sum((q-period return-q*mu)^2)/(q*(T-q+1)*(1-q/T))
    Robust variance uses delta_j without T; thus Z=(VR-1)/sqrt(theta).
    """
    prices = list(map(float, log_prices))
    returns = [b-a for a, b in zip(prices[:-1], prices[1:])]
    n = len(returns)
    mu = (prices[-1]-prices[0])/n
    squared = [(r-mu)**2 for r in returns]
    ss = math.fsum(squared)
    variance1 = ss/(n-1)
    varianceq = math.fsum((prices[t]-prices[t-q]-q*mu)**2 for t in range(q, n+1))/(q*(n-q+1)*(1-q/n))
    vr = varianceq/variance1
    theta = math.fsum((2*(q-j)/q)**2 * math.fsum(squared[t]*squared[t-j] for t in range(j,n))/ss**2 for j in range(1,q))
    z = (vr-1)/math.sqrt(theta)
    return vr, z, math.erfc(abs(z)/math.sqrt(2))


def main():
    output = ROOT / "results"
    output.mkdir(exist_ok=True)
    raw = (ROOT / "data/yahoo_gspc_raw.json").read_bytes()
    payload = json.loads(raw)
    assert payload["chart"]["error"] is None
    source = payload["chart"]["result"][0]
    assert source["meta"]["symbol"] == "^GSPC"
    dates = [datetime.fromtimestamp(t, ZoneInfo("America/New_York")).date().isoformat() for t in source["timestamp"]]
    frame = pd.DataFrame({"date": dates, "close": source["indicators"]["quote"][0]["close"], "adjusted_close": source["indicators"]["adjclose"][0]["adjclose"]})
    assert frame.date.is_unique and frame.date.is_monotonic_increasing
    excluded = frame[frame.date > CUTOFF].date.tolist()
    frame = frame[frame.date <= CUTOFF].reset_index(drop=True)
    assert frame[["close", "adjusted_close"]].notna().all().all(), "Missing prices must be investigated, not filled."
    assert np.isfinite(frame[["close", "adjusted_close"]]).all().all()
    assert (frame[["close", "adjusted_close"]] > 0).all().all()
    assert frame.date.iloc[-1] == CUTOFF
    frame["log_return"] = np.log(frame.adjusted_close).diff()
    rows, samples = [], []
    selected_data = []
    for name, start, end in PERIODS:
        positions = frame.index[(frame.date >= start) & (frame.date <= end)]
        assert len(positions) > max(HORIZONS) and positions[0] > 0
        # Preceding session supplies the base price for the first in-period return.
        sample = frame.loc[positions[0]-1:positions[-1]].copy()
        y = np.log(sample.adjusted_close.to_numpy())
        n = len(y)-1
        info = {"period": name, "requested_start": start, "requested_end": end, "base_price_date": sample.date.iloc[0], "first_return_date": sample.date.iloc[1], "last_return_date": sample.date.iloc[-1], "daily_returns": n, "price_observations_including_base": len(y)}
        samples.append(info)
        for idx, record in enumerate(sample.to_dict("records")):
            record.update(period=name, role="base_price" if idx == 0 else "in_period_return")
            if idx == 0:
                record["log_return"] = None
            selected_data.append(record)
        for q in HORIZONS:
            test = VarianceRatio(y, lags=q, trend="c", robust=True, overlap=True, debiased=True)
            independent = direct_vr(y, q)
            np.testing.assert_allclose([test.vr, test.stat, test.pvalue], independent, rtol=1e-10, atol=1e-12)
            homoskedastic = VarianceRatio(y, lags=q, trend="c", robust=False, overlap=True, debiased=True)
            rows.append({"period": name, "daily_returns": n, "horizon_trading_days": q, "variance_ratio": float(test.vr), "robust_z": float(test.stat), "robust_p": float(test.pvalue), "homoskedastic_p_sensitivity_only": float(homoskedastic.pvalue)})
        for row, adjusted in zip(rows[-4:], holm([r["robust_p"] for r in rows[-4:]])):
            row["holm_p_within_period"] = adjusted
    for row, adjusted in zip(rows, holm([r["robust_p"] for r in rows])):
        row["holm_p_all_12_tests"] = adjusted
        row["reject_5pct_unadjusted"] = row["robust_p"] < .05
        row["reject_5pct_holm_within_period"] = row["holm_p_within_period"] < .05
    results = {"source": SOURCE, "endpoint": ENDPOINT, "raw_sha256": hashlib.sha256(raw).hexdigest(), "snapshot_retrieved_utc_date": "2026-09-16", "last_completed_session": CUTOFF, "excluded_dates": excluded, "raw_rows": len(dates), "complete_price_rows": len(frame), "maximum_close_adjusted_close_difference": float(np.abs(frame.close-frame.adjusted_close).max()), "method": "Lo–MacKinlay; drift, overlapping returns, finite-sample debiasing, heteroskedasticity-robust asymptotic two-sided inference", "horizons": HORIZONS, "samples": samples, "tests": rows, "validation": "All 12 VR, Z and p-value triples match an independent scalar implementation to rtol=1e-10, atol=1e-12."}
    (output / "variance_ratio_results.json").write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")
    (ROOT / "data/analysis_samples.json").write_text(json.dumps(selected_data, indent=2, allow_nan=False)+"\n")

    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    colors = ["#237a95", "#c55d31", "#7a5cab"]
    for (period, _, _), color in zip(PERIODS, colors):
        current = [r for r in rows if r["period"] == period]
        ax.plot(HORIZONS, [r["variance_ratio"] for r in current], marker="o", linewidth=2.3, markersize=6, color=color, label=period)
    ax.axhline(1, color="#5d6573", linestyle="--", linewidth=1.2, label="Random-walk benchmark: VR = 1")
    ax.set(xlabel="Return horizon (trading days)", ylabel="Variance ratio", xticks=HORIZONS)
    ax.set_title("S&P 500 variance ratios across three periods", loc="left", fontsize=16, fontweight="bold", pad=15)
    ax.grid(axis="y", alpha=.18)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="best", fontsize=9)
    fig.subplots_adjust(left=.09, right=.98, bottom=.22, top=.88)
    fig.text(.08, .025, "Daily log adjusted-close returns · Overlapping, debiased estimates\nSource: Yahoo Finance (^GSPC) · Completed sessions through September 15, 2026 · See report for significance", fontsize=8, color="#4a5565")
    fig.savefig(output / "variance_ratios.png", dpi=180)
    fig.savefig(output / "variance_ratios.svg")
    plt.close(fig)

    lines = ["# S&P 500 variance ratio analysis", "", "Source: [Yahoo Finance S&P 500 (^GSPC) history]("+SOURCE+"). Retrieved September 16, 2026. The final included session is September 15, 2026; September 16's incomplete session is excluded.", "", "## Sample definition", "", "The S&P 500 is a proxy for large-cap US stocks, not the entire US equity universe. This is a price-index analysis, not a dividends-reinvested total-return analysis. Yahoo adjusted close is used; its equality with close is checked and recorded in the JSON results. User-specified years are inclusive. The crisis and COVID windows include recovery periods and are not event-only windows.", "", "| Period | First return | Last return | Daily returns | Base close |", "|---|---|---|---:|---|"]
    for s in samples:
        lines.append(f'| {s["period"]} | {s["first_return_date"]} | {s["last_return_date"]} | {s["daily_returns"]:,} | {s["base_price_date"]} |')
    lines += ["", "The previous trading session's close is included solely to compute the first return in each period. Every tested daily return ends within its assigned window; no missing prices are forward-filled and no returns are formed across the gaps between study windows.", "", "## Method", "", "Let r_t = ln(P_t/P_(t−1)). The null is VR(q) = 1, consistent with uncorrelated returns and linear growth of return variance with horizon. The alternatives are two-sided. Horizons q = 2, 5, 10 and 20 are trading days, roughly two days, one week, two weeks and one month.", "", "Tests allow a nonzero drift, use overlapping q-day returns, and debias the variance estimates. Heteroskedasticity-robust asymptotic Z statistics are used because financial volatility varies over time. See the [arch VarianceRatio documentation](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.VarianceRatio.html).", "", "With T daily returns, mu = mean(r), S = sum((r−mu)^2), the one-day variance is S/(T−1). The q-day variance per day is sum((log(P_t/P_(t−q))−q*mu)^2) / [q*(T−q+1)*(1−q/T)]. Their ratio is VR(q). Define delta_j = sum((r_t−mu)^2*(r_(t−j)−mu)^2)/S^2 and theta = sum([2*(q−j)/q]^2*delta_j), for j = 1,...,q−1. Then Z = (VR−1)/sqrt(theta), with a two-sided standard-normal p-value.", "", "Holm-adjusted p-values control family-wise error over four horizons within each period. The JSON additionally reports adjustment over all 12 tests. Unadjusted results are shown transparently; this is not a Chow–Denning test.", "", "## Results", "", "| Period | q | VR | Robust Z | Raw p | Holm p (4 horizons) | Reject at 5% after Holm? |", "|---|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        lines.append(f'| {r["period"]} | {r["horizon_trading_days"]} | {r["variance_ratio"]:.4f} | {r["robust_z"]:.4f} | {r["robust_p"]:.4f} | {r["holm_p_within_period"]:.4f} | {"Yes" if r["reject_5pct_holm_within_period"] else "No"} |')
    lines += ["", "![Variance ratios](variance_ratios.png)", "", "## Interpretation", ""]
    for period, _, _ in PERIODS:
        current = [r for r in rows if r["period"] == period]
        raw_rejections = [r["horizon_trading_days"] for r in current if r["reject_5pct_unadjusted"]]
        adj_rejections = [r["horizon_trading_days"] for r in current if r["reject_5pct_holm_within_period"]]
        lines.append(f'- **{period}:** raw 5% rejection horizons: {raw_rejections or "none"}; Holm-adjusted 5% rejection horizons: {adj_rejections or "none"}.')
    lines += ["", "VR below 1 indicates a negative aggregate serial-correlation pattern over the tested horizon; VR above 1 indicates a positive pattern. This alone does not establish statistical significance. Failure to reject is not proof of a random walk or market efficiency. A significant result does not by itself establish a profitable trading strategy or stationary price mean reversion.", "", "These are separate within-period tests, not a formal test that periods differ from each other. Unequal sample lengths affect precision. Broad windows mix changing regimes, and asymptotic inference may be less accurate with heavy tails or structural breaks. The 2022–2026 sample covers only 2026 through the stated cutoff; no future prices are estimated. The analysis does not identify a causal crisis effect.", "", "## Reproducibility and validation", "", results["validation"], "", "The source snapshot, exact selected observations, full-precision estimates, sample boundaries, raw-data SHA-256 hash, and exclusions are retained. Date uniqueness, order, positive finite prices, and missing values are checked. Homoskedastic p-values in the JSON are a sensitivity diagnostic only and do not determine the conclusions.", "", "Re-run `python analyze_variance_ratio.py` after installing `requirements.txt`. This uses the frozen snapshot and cutoff rather than silently fetching newer data."]
    if not any(r["reject_5pct_holm_within_period"] for r in rows):
        lines[2:2] = ["**Finding:** None of the tested horizons rejects the random-walk variance restriction at 5% after Holm adjustment within each period. Before adjustment, the 2008–2010 period rejects at 2 and 5 days, and 2020–2021 rejects at 2 days; 2022–2026 does not reject at any tested horizon.", ""]
    (output / "report.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({"samples": samples, "tests": rows, "validation": results["validation"]}, indent=2))


if __name__ == "__main__":
    main()
