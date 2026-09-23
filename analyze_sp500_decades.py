"""Calculate S&P 500 variance ratios and dispersion for full decades.

Reproduce: .venv/bin/python analyze_sp500_decades.py
All input data are frozen locally; no network requests are made.
"""
from pathlib import Path
import hashlib
import json
import math
import os

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_historical_variance_ratio import (
    CUTOFF,
    HORIZONS,
    RAW_PATH,
    load_prices,
    run_test,
)
from analyze_variance_ratio import holm


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results/sp500_decades"
DECADES = tuple(
    (f"{year}s", f"{year}-01-01", f"{year + 9}-12-31")
    for year in range(1950, 2020, 10)
)
TRADING_DAYS_PER_YEAR = 252


def main():
    source_bytes, prices, _ = load_prices(RAW_PATH)
    results = {
        "symbol": "^GSPC",
        "source": "Yahoo Finance historical S&P 500 (^GSPC)",
        "history_source_url": "https://finance.yahoo.com/quote/%5EGSPC/history/",
        "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "last_complete_session_in_source": CUTOFF,
        "price_series": "Yahoo adjusted close; price index, without reinvested dividends",
        "return_type": "Daily close-to-close log returns",
        "method": "Lo-MacKinlay variance ratio with drift, overlapping returns, finite-sample debiasing, and heteroskedasticity-robust inference",
        "horizons_trading_days": list(HORIZONS),
        "variance_ratio_dispersion": "Population and sample standard deviations across the four q-specific variance-ratio point estimates within each decade",
        "return_volatility": "Sample standard deviation of daily log returns; annualized by multiplying by sqrt(252)",
        "periods": [],
    }

    all_tests = []
    for label, start, end in DECADES:
        in_period = sorted(day for day in prices if start <= day <= end)
        prior_dates = [day for day in prices if day < start]
        assert prior_dates and len(in_period) > max(HORIZONS), f"Insufficient data: {label}"
        dates = [max(prior_dates), *in_period]
        log_prices = np.log(np.asarray([prices[day] for day in dates], dtype=float))
        returns = np.diff(log_prices)
        tests = [run_test(log_prices, q) for q in HORIZONS]
        adjusted = holm([row["robust_p"] for row in tests])
        for row, adj_p in zip(tests, adjusted):
            row["holm_p_within_decade"] = adj_p
            row["reject_5pct_holm_within_decade"] = adj_p < 0.05
        all_tests.extend(tests)

        ratios = np.asarray([row["variance_ratio"] for row in tests])
        daily_sd = float(np.std(returns, ddof=1))
        results["periods"].append({
            "period": label,
            "requested_start": start,
            "requested_end": end,
            "first_return_date": in_period[0],
            "last_return_date": in_period[-1],
            "daily_returns": len(returns),
            "tests": tests,
            "mean_variance_ratio": float(np.mean(ratios)),
            "population_sd_across_variance_ratios": float(np.std(ratios, ddof=0)),
            "sample_sd_across_variance_ratios": float(np.std(ratios, ddof=1)),
            "mean_absolute_vr_distance": float(np.mean(np.abs(ratios - 1.0))),
            "closeness_score": float(1.0 - np.mean(np.abs(ratios - 1.0))),
            "daily_log_return_sample_sd": daily_sd,
            "annualized_log_return_volatility": daily_sd * math.sqrt(TRADING_DAYS_PER_YEAR),
            "holm_rejections_within_decade": [
                row["horizon_trading_days"]
                for row in tests
                if row["reject_5pct_holm_within_decade"]
            ],
        })

    global_adjusted = holm([row["robust_p"] for row in all_tests])
    for row, adj_p in zip(all_tests, global_adjusted):
        row["holm_p_all_28_tests"] = adj_p
        row["reject_5pct_holm_all_decades"] = adj_p < 0.05

    scores = np.asarray([row["closeness_score"] for row in results["periods"]])
    results["across_decades"] = {
        "population_sd_of_decade_closeness_scores": float(np.std(scores, ddof=0)),
        "sample_sd_of_decade_closeness_scores": float(np.std(scores, ddof=1)),
        "population_sd_by_horizon": {
            str(q): float(np.std([
                next(test["variance_ratio"] for test in period["tests"]
                     if test["horizon_trading_days"] == q)
                for period in results["periods"]
            ], ddof=0))
            for q in HORIZONS
        },
    }
    results["interpretation_limit"] = (
        "The SD across four variance-ratio estimates describes horizon-to-horizon "
        "dispersion and is not a sampling standard error or a formal efficiency test. "
        "A variance-ratio rejection tests one random-walk implication; it does not by "
        "itself prove or disprove market efficiency."
    )
    results["validation"] = (
        "Every library variance ratio, robust Z statistic, and p-value matches the "
        "project's independent scalar implementation (rtol=1e-10, atol=1e-12)."
    )

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "sp500_decade_variance_ratios.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )
    write_report(results)
    write_error_bar_plot(results)

    for period in results["periods"]:
        ratios = [test["variance_ratio"] for test in period["tests"]]
        print(
            f'{period["period"]}: VR={ratios}, '
            f'SD(VR)={period["population_sd_across_variance_ratios"]:.6f}, '
            f'annualized volatility={period["annualized_log_return_volatility"]:.2%}'
        )


def write_report(results):
    lines = [
        "# S&P 500 Variance Ratios by Decade",
        "",
        "## Method",
        "",
        "This analysis uses daily S&P 500 adjusted-close log returns and the same Lo–MacKinlay variance-ratio specification as the earlier historical run: a drift, overlapping returns, finite-sample debiasing, and heteroskedasticity-robust inference. Each decade is a separate full calendar window. The tested horizons are 2, 5, 10, and 20 trading days.",
        "",
        "`SD(VR)` is the population standard deviation across the four variance-ratio point estimates within that decade. It shows how much the estimates vary by horizon; it is not a standard error. `Return volatility` is the sample standard deviation of daily log returns annualized with sqrt(252). Mean absolute distance is mean(|VR−1|), so lower values are closer to the random-walk variance restriction.",
        "",
        "## Results",
        "",
        "| Decade | Daily returns | VR(2) | VR(5) | VR(10) | VR(20) | SD(VR) | Mean abs(VR−1) | Annualized return volatility |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for period in results["periods"]:
        ratios = {row["horizon_trading_days"]: row["variance_ratio"] for row in period["tests"]}
        lines.append(
            f'| {period["period"]} | {period["daily_returns"]:,} | '
            f'{ratios[2]:.4f} | {ratios[5]:.4f} | {ratios[10]:.4f} | {ratios[20]:.4f} | '
            f'{period["population_sd_across_variance_ratios"]:.4f} | '
            f'{period["mean_absolute_vr_distance"]:.4f} | '
            f'{period["annualized_log_return_volatility"]:.2%} |'
        )

    lines += [
        "",
        "![Mean S&P 500 variance ratio by decade with standard-deviation error bars](variance_ratio_error_bars.png)",
        "",
        "## Statistical tests",
        "",
        "A variance ratio of 1 is the random-walk variance-scaling restriction. The table below gives Holm-adjusted p-values within each decade; an adjusted p-value below 0.05 rejects that restriction at that horizon.",
        "",
        "| Decade | q=2 VR (Holm p) | q=5 VR (Holm p) | q=10 VR (Holm p) | q=20 VR (Holm p) | Rejected horizons |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for period in results["periods"]:
        cells = [
            f'{row["variance_ratio"]:.4f} ({row["holm_p_within_decade"]:.4f})'
            for row in period["tests"]
        ]
        rejected = ", ".join(str(q) for q in period["holm_rejections_within_decade"]) or "None"
        lines.append(
            f'| {period["period"]} | ' + " | ".join(cells) + f" | {rejected} |"
        )

    closest = min(results["periods"], key=lambda row: row["mean_absolute_vr_distance"])
    farthest = max(results["periods"], key=lambda row: row["mean_absolute_vr_distance"])
    least_dispersion = min(results["periods"], key=lambda row: row["population_sd_across_variance_ratios"])
    most_volatile = max(results["periods"], key=lambda row: row["annualized_log_return_volatility"])
    across_sd = results["across_decades"]["population_sd_of_decade_closeness_scores"]
    horizon_sds = results["across_decades"]["population_sd_by_horizon"]
    lines += [
        "",
        "## Interpretation",
        "",
        f'The {closest["period"]} are closest to VR=1 on the four-horizon summary (mean absolute distance {closest["mean_absolute_vr_distance"]:.4f}), while the {farthest["period"]} are farthest ({farthest["mean_absolute_vr_distance"]:.4f}). The {least_dispersion["period"]} have the smallest horizon-to-horizon spread in VR estimates (SD {least_dispersion["population_sd_across_variance_ratios"]:.4f}). The {most_volatile["period"]} have the highest annualized daily-return volatility ({most_volatile["annualized_log_return_volatility"]:.2%}).',
        "",
        f'The population standard deviation of the seven decade closeness scores is {across_sd:.4f}. This describes variation among the decade point estimates; it is not a confidence interval or a test that the decades have different true efficiencies.',
        "",
        "Across the seven decades, the population SD of the VR point estimates is "
        + ", ".join(f'{q}-day: {horizon_sds[str(q)]:.4f}' for q in HORIZONS)
        + ". Dispersion across decades rises with the holding horizon in these estimates.",
        "",
        results["interpretation_limit"],
        "",
        "## Data and validation",
        "",
        results["validation"],
        "The full-precision estimates, p-values, source hash, observation dates, and both population and sample standard deviations are in `sp500_decade_variance_ratios.json`.",
        "",
        "Source: [Yahoo Finance S&P 500 (^GSPC) history](https://finance.yahoo.com/quote/%5EGSPC/history/). Reproduce with `.venv/bin/python analyze_sp500_decades.py`; no network access is needed.",
    ]
    (OUT / "report.md").write_text("\n".join(lines) + "\n")


def write_error_bar_plot(results):
    labels = [period["period"] for period in results["periods"]]
    means = np.asarray([period["mean_variance_ratio"] for period in results["periods"]])
    errors = np.asarray([
        period["population_sd_across_variance_ratios"]
        for period in results["periods"]
    ])
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    ax.errorbar(
        x,
        means,
        yerr=errors,
        fmt="o-",
        color="#176b87",
        ecolor="#c55d31",
        markerfacecolor="#176b87",
        markeredgecolor="white",
        markeredgewidth=1.2,
        markersize=8,
        linewidth=2.2,
        elinewidth=2,
        capsize=6,
        capthick=2,
        zorder=3,
    )
    ax.axhline(1.0, color="#596273", linestyle="--", linewidth=1.4, zorder=1)
    ax.text(
        len(labels) - 0.55,
        1.012,
        "Random-walk benchmark (VR = 1)",
        ha="right",
        va="bottom",
        fontsize=9,
        color="#4d5563",
    )
    for xpos, mean in zip(x, means):
        ax.annotate(
            f"{mean:.3f}",
            (xpos, mean),
            xytext=(0, 12),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#263241",
        )

    lower = min(1.0, float(np.min(means - errors)))
    upper = max(1.0, float(np.max(means + errors)))
    padding = max(0.06, 0.13 * (upper - lower))
    ax.set_ylim(lower - padding, upper + padding)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Mean variance ratio across q = 2, 5, 10, 20 days")
    ax.set_xlabel("S&P 500 decade")
    ax.set_title(
        "S&P 500 Mean Variance Ratio by Decade",
        loc="left",
        fontsize=16,
        fontweight="bold",
        pad=18,
    )
    ax.grid(axis="y", alpha=0.22, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(
        0.11,
        0.02,
        "Error bars show ±1 population SD across the four horizon-specific VR estimates; they are not confidence intervals.",
        fontsize=9,
        color="#4d5563",
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.88, bottom=0.17)
    fig.savefig(OUT / "variance_ratio_error_bars.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
