"""Compare small- and large-cap variance ratios by decade.

Reproduce: .venv/bin/python analyze_size_decades.py
The official Kenneth French source archive is frozen in data/.
"""
from pathlib import Path
from zipfile import ZipFile
import csv
import hashlib
import io
import json
import math
import os

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import numpy as np
from arch.bootstrap import StationaryBootstrap
from arch.unitroot import VarianceRatio
from scipy.stats import norm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_industry_efficiency import efficiency, variance_ratios
from analyze_variance_ratio import direct_vr, holm


SOURCE = ROOT / "data/size_portfolios_me_daily.csv.zip"
OUT = ROOT / "results/size_decades"
HORIZONS = (2, 5, 10, 20)
PORTFOLIOS = (("small", "Lo 30"), ("large", "Hi 30"))
DECADES = tuple(
    (f"{year}s", f"{year}0101", f"{year + 9}1231")
    for year in range(1960, 2020, 10)
)


def read_value_weighted_returns():
    """Read the first (value-weighted) daily return table from the source zip."""
    with ZipFile(SOURCE) as archive:
        names = archive.namelist()
        assert names == ["Portfolios_Formed_on_ME_daily.csv"]
        text = archive.read(names[0]).decode("utf-8")
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines)
                 if line.strip() == "Average Value Weighted Returns -- Daily")
    reader = csv.reader(io.StringIO("\n".join(lines[start + 1:])))
    header = [cell.strip() for cell in next(reader)]
    assert header[:5] == ["", "<= 0", "Lo 30", "Med 40", "Hi 30"]
    indices = {label: header.index(label) for _, label in PORTFOLIOS}
    rows = []
    for fields in reader:
        if not fields or not fields[0].strip():
            break
        date = fields[0].strip()
        assert len(date) == 8 and date.isdigit()
        values = [float(fields[indices[label]]) for _, label in PORTFOLIOS]
        assert all(value not in (-99.99, -999.0) and value > -100 for value in values)
        rows.append((date, *values))
    assert rows and rows == sorted(rows) and len(rows) == len({row[0] for row in rows})
    assert rows[0][0] == "19260701" and rows[-1][0] >= "20260731"
    return rows


def run_individual_tests(log_returns, labels):
    vr = variance_ratios(log_returns)
    tests = []
    for i, label in enumerate(labels):
        log_wealth = np.r_[0.0, np.cumsum(log_returns[:, i])]
        for j, q in enumerate(HORIZONS):
            test = VarianceRatio(
                log_wealth, lags=q, trend="c", robust=True,
                overlap=True, debiased=True,
            )
            np.testing.assert_allclose(vr[i, j], test.vr, rtol=1e-10, atol=1e-12)
            np.testing.assert_allclose(
                [test.vr, test.stat, test.pvalue],
                direct_vr(log_wealth, q), rtol=1e-10, atol=1e-12,
            )
            tests.append({
                "portfolio": label,
                "horizon_trading_days": q,
                "variance_ratio": float(test.vr),
                "robust_z": float(test.stat),
                "robust_p": float(test.pvalue),
            })
    for row, adjusted in zip(tests, holm([row["robust_p"] for row in tests])):
        row["holm_p_within_decade_eight_tests"] = adjusted
        row["reject_5pct_holm_within_decade"] = adjusted < 0.05
    return vr, tests


def score_difference(log_returns):
    scores = efficiency(variance_ratios(log_returns))
    return float(scores[1] - scores[0])


def bootstrap_differences(log_returns, block_length, reps, seed):
    sampler = StationaryBootstrap(block_length, log_returns, seed=seed)
    return np.fromiter(
        (score_difference(sample[0][0]) for sample in sampler.bootstrap(reps)),
        dtype=float,
        count=reps,
    )


def add_comparison_inference(periods):
    """Paired block-bootstrap inference for large score minus small score."""
    comparisons = {length: [] for length in (40, 20, 60)}
    for length in (40, 20, 60):
        reps = 4999 if length == 40 else 2499
        for i, period in enumerate(periods):
            returns = np.asarray(period.pop(f"_returns_{length}"))
            point = period["comparison"]["large_minus_small_score"]
            draws = bootstrap_differences(
                returns, length, reps, 20260922 + 1000 * i + length,
            )
            errors = draws - point
            se = float(np.std(errors, ddof=1))
            assert se > 0
            z = point / se
            p = math.erfc(abs(z) / math.sqrt(2))
            comparisons[length].append({
                "decade": period["decade"],
                "large_minus_small_score": point,
                "stationary_bootstrap_standard_error": se,
                "wald_z": float(z),
                "two_sided_p": p,
                "mean_block_length": length,
                "replicates": reps,
                "seed": 20260922 + 1000 * i + length,
            })
        adjusted = holm([row["two_sided_p"] for row in comparisons[length]])
        critical = float(norm.ppf(1 - 0.05 / (2 * len(periods))))
        for row, adj_p in zip(comparisons[length], adjusted):
            delta = row["large_minus_small_score"]
            se = row["stationary_bootstrap_standard_error"]
            row["holm_p_across_six_decades"] = adj_p
            row["familywise_95pct_simultaneous_interval"] = [
                float(delta - critical * se),
                float(delta + critical * se),
            ]
            row["reject_equal_scores_holm_5pct"] = adj_p < 0.05
    for period, row in zip(periods, comparisons[40]):
        period["comparison"]["inference_main"] = row
    return {
        "quantity": "large-cap closeness score minus small-cap closeness score",
        "sign": "Positive means the large-cap portfolio is closer to VR=1 across the four horizons.",
        "method": "Paired stationary block bootstrap of the two daily return series, retaining contemporaneous dependence; Wald inference with Holm p-values and Bonferroni simultaneous intervals across six decades.",
        "main_mean_block_length": 40,
        "comparisons_by_mean_block_length": {
            str(length): rows for length, rows in comparisons.items()
        },
    }


def main():
    raw_rows = read_value_weighted_returns()
    labels = [name for name, _ in PORTFOLIOS]
    periods = []
    all_tests = []

    for decade, start, end in DECADES:
        selected = [row for row in raw_rows if start <= row[0] <= end]
        assert len(selected) > max(HORIZONS)
        simple_pct = np.asarray([row[1:] for row in selected], dtype=float)
        log_returns = np.log1p(simple_pct / 100.0)
        assert np.isfinite(log_returns).all()
        vr, tests = run_individual_tests(log_returns, labels)
        all_tests.extend(tests)
        scores = efficiency(vr)
        portfolio_results = []
        for i, label in enumerate(labels):
            portfolio_tests = [row for row in tests if row["portfolio"] == label]
            annualized_volatility = float(np.std(log_returns[:, i], ddof=1) * math.sqrt(252))
            portfolio_results.append({
                "portfolio": label,
                "source_column": PORTFOLIOS[i][1],
                "variance_ratios": {
                    str(q): float(vr[i, j]) for j, q in enumerate(HORIZONS)
                },
                "mean_variance_ratio": float(np.mean(vr[i])),
                "population_sd_across_variance_ratios": float(np.std(vr[i], ddof=0)),
                "sample_sd_across_variance_ratios": float(np.std(vr[i], ddof=1)),
                "mean_absolute_vr_distance": float(1.0 - scores[i]),
                "closeness_score": float(scores[i]),
                "annualized_log_return_volatility": annualized_volatility,
                "tests": portfolio_tests,
            })
        difference = float(scores[1] - scores[0])
        period = {
            "decade": decade,
            "requested_start": f"{start[:4]}-01-01",
            "requested_end": f"{end[:4]}-12-31",
            "first_return_date": f"{selected[0][0][:4]}-{selected[0][0][4:6]}-{selected[0][0][6:]}",
            "last_return_date": f"{selected[-1][0][:4]}-{selected[-1][0][4:6]}-{selected[-1][0][6:]}",
            "daily_returns": len(selected),
            "portfolios": portfolio_results,
            "comparison": {
                "large_minus_small_score": difference,
                "closer_point_estimate": "large" if difference > 0 else "small" if difference < 0 else "equal",
            },
        }
        for length in (40, 20, 60):
            period[f"_returns_{length}"] = log_returns.tolist()
        periods.append(period)

    for row, adjusted in zip(all_tests, holm([row["robust_p"] for row in all_tests])):
        row["holm_p_all_48_tests"] = adjusted
        row["reject_5pct_holm_all_tests"] = adjusted < 0.05

    inference = add_comparison_inference(periods)
    results = {
        "source": "Kenneth R. French Data Library, Portfolios Formed on Size",
        "source_url": "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_port_form_sz.html",
        "download_url": "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/Portfolios_Formed_on_ME_Daily_CSV.zip",
        "source_archive": str(SOURCE.relative_to(ROOT)),
        "source_archive_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "source_vintage": "202607 CRSP database",
        "retrieved_date": "2026-09-22",
        "portfolio_definition": {
            "small": "Bottom 30% by market equity using NYSE breakpoints",
            "large": "Top 30% by market equity using NYSE breakpoints",
            "construction": "Value-weighted portfolios formed at the end of each June; includes eligible NYSE, AMEX, and NASDAQ stocks",
        },
        "return_type": "Daily value-weighted total returns converted from percentage simple returns to log returns",
        "method": "Lo-MacKinlay variance ratio with drift, overlapping returns, finite-sample debiasing, and heteroskedasticity-robust inference",
        "horizons_trading_days": list(HORIZONS),
        "score": "1 - mean(abs(VR(q)-1)) across q=2,5,10,20; higher means closer to the random-walk variance restriction",
        "periods": periods,
        "large_vs_small_inference": inference,
        "validation": "All 48 variance ratios, robust Z statistics, and p-values match an independent scalar implementation (rtol=1e-10, atol=1e-12).",
        "interpretation_limit": "Closeness to VR=1 is evidence about one random-walk implication, not a complete measure of market efficiency or investable predictability.",
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "size_decade_variance_ratios.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )
    write_report(results)
    write_plot(results)
    for period in periods:
        comparison = period["comparison"]
        print(
            f'{period["decade"]}: large-small score={comparison["large_minus_small_score"]:.6f}, '
            f'Holm p={comparison["inference_main"]["holm_p_across_six_decades"]:.6f}'
        )


def write_report(results):
    periods = results["periods"]
    lines = [
        "# Large-Cap versus Small-Cap Variance Ratios by Decade",
        "",
        "## Definitions and method",
        "",
        "Small-cap is the value-weighted bottom 30% portfolio and large-cap is the value-weighted top 30% portfolio from the Kenneth French Portfolios Formed on Size dataset. The size cutoffs use NYSE market-equity breakpoints; eligible NYSE, AMEX, and Nasdaq stocks are assigned to portfolios formed at the end of each June. These are research portfolios, not fixed constituent indexes.",
        "",
        "Daily percentage total returns are converted to log returns. Each decade uses the same Lo–MacKinlay specification as the earlier S&P 500 analysis: q = 2, 5, 10, and 20 trading days, a drift, overlapping returns, finite-sample debiasing, and heteroskedasticity-robust inference.",
        "",
        "The closeness score is E = 1 − mean(|VR(q)−1|). Higher values mean the four estimates are closer to the random-walk variance-scaling restriction. The score is not a percentage or a complete test of market efficiency.",
        "",
        "## Direct comparison",
        "",
        "The difference is large-cap score minus small-cap score. Positive values favor large-cap. The interval is a family-wise 95% simultaneous interval across the six decades from the paired 40-day stationary-block bootstrap; Holm p-values also adjust across the six comparisons.",
        "",
        "| Decade | Small score | Large score | Large − small | Closer point estimate | Family-wise 95% interval | Holm p |",
        "|---|---:|---:|---:|---|---:|---:|",
    ]
    for period in periods:
        portfolios = {row["portfolio"]: row for row in period["portfolios"]}
        comparison = period["comparison"]
        inference = comparison["inference_main"]
        interval = inference["familywise_95pct_simultaneous_interval"]
        lines.append(
            f'| {period["decade"]} | {portfolios["small"]["closeness_score"]:.4f} | '
            f'{portfolios["large"]["closeness_score"]:.4f} | '
            f'{comparison["large_minus_small_score"]:+.4f} | '
            f'{comparison["closer_point_estimate"].title()} | '
            f'[{interval[0]:.4f}, {interval[1]:.4f}] | '
            f'{format_p(inference["holm_p_across_six_decades"])} |'
        )
    lines += [
        "",
        "![Large-cap and small-cap mean variance ratios by decade](size_variance_ratio_error_bars.png)",
        "",
        "## Variance ratios and dispersion",
        "",
        "SD(VR) is the population standard deviation across the four horizon-specific ratios. Annualized volatility is the daily log-return sample standard deviation multiplied by sqrt(252).",
        "",
        "| Decade | Portfolio | VR(2) | VR(5) | VR(10) | VR(20) | SD(VR) | Mean abs(VR−1) | Annualized volatility |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for period in periods:
        for portfolio in period["portfolios"]:
            vr = portfolio["variance_ratios"]
            lines.append(
                f'| {period["decade"]} | {portfolio["portfolio"].title()} | '
                f'{vr["2"]:.4f} | {vr["5"]:.4f} | {vr["10"]:.4f} | {vr["20"]:.4f} | '
                f'{portfolio["population_sd_across_variance_ratios"]:.4f} | '
                f'{portfolio["mean_absolute_vr_distance"]:.4f} | '
                f'{portfolio["annualized_log_return_volatility"]:.2%} |'
            )

    significant = [
        period for period in periods
        if period["comparison"]["inference_main"]["reject_equal_scores_holm_5pct"]
    ]
    large_wins = [
        period["decade"] for period in periods
        if period["comparison"]["closer_point_estimate"] == "large"
    ]
    small_wins = [
        period["decade"] for period in periods
        if period["comparison"]["closer_point_estimate"] == "small"
    ]
    lines += [
        "",
        "## Answer to the large-versus-small question",
        "",
        f'Large-cap is closer to VR=1 by point estimate in {", ".join(large_wins) if large_wins else "no decade"}. Small-cap is closer in {", ".join(small_wins) if small_wins else "no decade"}.',
        "",
        (
            "After the paired bootstrap and Holm adjustment, statistically detectable size differences occur in "
            + ", ".join(period["decade"] for period in significant)
            + "."
            if significant else
            "After the paired bootstrap and Holm adjustment, no decade has a statistically detectable large-versus-small score difference at 5%."
        ),
        "",
        "The ranking can change because VR values above and below 1 both count as departures; a lower raw VR is not automatically more or less random. Compare the absolute distance from 1 or the closeness score.",
        "",
        "The very high older small-cap ratios indicate strong positive aggregate serial correlation in the portfolio returns. Nonsynchronous trading and stale prices can contribute to this pattern in less-liquid small stocks, so it should not automatically be interpreted as an exploitable trading opportunity.",
        "",
        "The 20-day and 60-day block-length sensitivity results are stored in the JSON. A result that changes under those alternatives should be treated cautiously.",
        "",
        "## Data and validation",
        "",
        results["validation"],
        "The full-precision results, individual tests, multiple-testing adjustments, bootstrap settings, sensitivity checks, and source hash are in `size_decade_variance_ratios.json`.",
        "",
        "Source: [Kenneth French, Portfolios Formed on Size](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_port_form_sz.html). Reproduce with `.venv/bin/python analyze_size_decades.py`; the saved source archive means no network request is needed.",
    ]
    (OUT / "report.md").write_text("\n".join(lines) + "\n")


def format_p(value):
    return "<0.0001" if value < 0.0001 else f"{value:.4f}"


def write_plot(results):
    periods = results["periods"]
    labels = [period["decade"] for period in periods]
    x = np.arange(len(labels), dtype=float)
    fig, ax = plt.subplots(figsize=(11.0, 6.6))
    styles = {
        "small": {"color": "#c55d31", "marker": "s", "offset": -0.08},
        "large": {"color": "#176b87", "marker": "o", "offset": 0.08},
    }
    for label in ("small", "large"):
        rows = [
            next(row for row in period["portfolios"] if row["portfolio"] == label)
            for period in periods
        ]
        means = np.asarray([row["mean_variance_ratio"] for row in rows])
        errors = np.asarray([row["population_sd_across_variance_ratios"] for row in rows])
        style = styles[label]
        ax.errorbar(
            x + style["offset"], means, yerr=errors,
            fmt=style["marker"] + "-", label=label.title(),
            color=style["color"], ecolor=style["color"],
            markerfacecolor=style["color"], markeredgecolor="white",
            markeredgewidth=1.1, markersize=7.5, linewidth=2,
            elinewidth=1.8, capsize=5, capthick=1.8, zorder=3,
        )
    ax.axhline(1.0, color="#596273", linestyle="--", linewidth=1.4, zorder=1)
    ax.text(
        len(labels) - 0.55, 1.012, "Random-walk benchmark (VR = 1)",
        ha="right", va="bottom", fontsize=9, color="#4d5563",
    )
    ax.set_xticks(x, labels)
    ax.set_xlabel("Decade")
    ax.set_ylabel("Mean variance ratio across q = 2, 5, 10, 20 days")
    ax.set_title(
        "Large-Cap and Small-Cap Variance Ratios by Decade",
        loc="left", fontsize=16, fontweight="bold", pad=18,
    )
    ax.grid(axis="y", alpha=0.22, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, loc="upper right")
    fig.text(
        0.10, 0.02,
        "Points are four-horizon means. Error bars show ±1 population SD across horizons; they are not confidence intervals.\nKenneth French value-weighted bottom 30% and top 30% size portfolios.",
        fontsize=9, color="#4d5563",
    )
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.19)
    fig.savefig(OUT / "size_variance_ratio_error_bars.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
