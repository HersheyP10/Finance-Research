"""Test whether large-cap returns lead small-cap returns by decade.

Reproduce: .venv/bin/python analyze_large_leads_small.py
The official Kenneth French source archive is frozen in data/.
"""
from pathlib import Path
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
from scipy.stats import norm
import statsmodels.api as sm

from analyze_size_decades import DECADES, SOURCE, read_value_weighted_returns
from analyze_variance_ratio import holm


OUT = ROOT / "results/large_leads_small"
HAC_LAGS = 5
ALPHA = 0.05


def _independent_hac_covariance(x, residuals, max_lags):
    """Newey-West covariance used to validate statsmodels' HAC result."""
    nobs, nparams = x.shape
    scores = x * residuals[:, None]
    middle = scores.T @ scores
    for lag in range(1, max_lags + 1):
        weight = 1.0 - lag / (max_lags + 1.0)
        gamma = scores[lag:].T @ scores[:-lag]
        middle += weight * (gamma + gamma.T)
    bread = np.linalg.inv(x.T @ x)
    covariance = bread @ middle @ bread
    return covariance * nobs / (nobs - nparams)


def fit_predictive_regression(target, target_lag, predictor_lag):
    """Estimate target(t+1) on its own and the other portfolio's t return."""
    x = sm.add_constant(np.column_stack([target_lag, predictor_lag]))
    fit = sm.OLS(target, x).fit(
        cov_type="HAC",
        cov_kwds={"maxlags": HAC_LAGS, "use_correction": True},
        use_t=False,
    )
    restricted = sm.OLS(target, sm.add_constant(target_lag)).fit()

    beta_direct = np.linalg.solve(x.T @ x, x.T @ target)
    np.testing.assert_allclose(beta_direct, fit.params, rtol=1e-11, atol=1e-13)
    cov_direct = _independent_hac_covariance(x, fit.resid, HAC_LAGS)
    np.testing.assert_allclose(cov_direct, fit.cov_params(), rtol=1e-10, atol=1e-13)

    beta = float(fit.params[2])
    se = float(fit.bse[2])
    z_value = beta / se
    p_value = float(2.0 * norm.sf(abs(z_value)))
    np.testing.assert_allclose(p_value, fit.pvalues[2], rtol=1e-12, atol=1e-15)

    return {
        "predictor_coefficient": beta,
        "hac_standard_error": se,
        "z_statistic": z_value,
        "two_sided_p_value": p_value,
        "target_own_lag_coefficient": float(fit.params[1]),
        "intercept": float(fit.params[0]),
        "full_model_r_squared": float(fit.rsquared),
        "restricted_model_r_squared": float(restricted.rsquared),
        "incremental_r_squared": float(fit.rsquared - restricted.rsquared),
        "standardized_predictor_coefficient": float(
            beta * np.std(predictor_lag, ddof=1) / np.std(target, ddof=1)
        ),
    }


def calculate_results():
    raw_rows = read_value_weighted_returns()
    periods = []
    for decade, start, end in DECADES:
        selected = [row for row in raw_rows if start <= row[0] <= end]
        returns = np.log1p(np.asarray([row[1:] for row in selected], dtype=float) / 100.0)
        small = returns[:, 0]
        large = returns[:, 1]

        large_to_small = fit_predictive_regression(
            target=small[1:], target_lag=small[:-1], predictor_lag=large[:-1]
        )
        small_to_large = fit_predictive_regression(
            target=large[1:], target_lag=large[:-1], predictor_lag=small[:-1]
        )
        periods.append({
            "decade": decade,
            "first_return_date": (
                f"{selected[0][0][:4]}-{selected[0][0][4:6]}-{selected[0][0][6:]}"
            ),
            "last_return_date": (
                f"{selected[-1][0][:4]}-{selected[-1][0][4:6]}-{selected[-1][0][6:]}"
            ),
            "regression_observations": int(len(small) - 1),
            "large_t_to_small_t_plus_1": large_to_small,
            "small_t_to_large_t_plus_1_diagnostic": small_to_large,
        })

    main_p_values = [
        row["large_t_to_small_t_plus_1"]["two_sided_p_value"] for row in periods
    ]
    reverse_p_values = [
        row["small_t_to_large_t_plus_1_diagnostic"]["two_sided_p_value"]
        for row in periods
    ]
    familywise_critical = float(norm.ppf(1.0 - ALPHA / (2.0 * len(periods))))
    for period, adjusted in zip(periods, holm(main_p_values)):
        row = period["large_t_to_small_t_plus_1"]
        row["holm_p_value_across_decades"] = adjusted
        row["reject_zero_holm_5pct"] = adjusted < ALPHA
        row["familywise_95pct_bonferroni_interval"] = [
            row["predictor_coefficient"] - familywise_critical * row["hac_standard_error"],
            row["predictor_coefficient"] + familywise_critical * row["hac_standard_error"],
        ]
    for period, adjusted in zip(periods, holm(reverse_p_values)):
        row = period["small_t_to_large_t_plus_1_diagnostic"]
        row["holm_p_value_across_decades"] = adjusted
        row["reject_zero_holm_5pct"] = adjusted < ALPHA

    return {
        "question": "Does large-cap return at t help predict small-cap return at t+1?",
        "source": "Kenneth R. French Data Library, Portfolios Formed on Size",
        "source_url": (
            "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/"
            "Data_Library/det_port_form_sz.html"
        ),
        "source_archive": str(SOURCE.relative_to(ROOT)),
        "source_archive_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "portfolio_definition": {
            "small": "Value-weighted bottom 30% by market equity using NYSE breakpoints",
            "large": "Value-weighted top 30% by market equity using NYSE breakpoints",
        },
        "return_type": "Daily total returns converted from percentage simple returns to log returns",
        "main_model": (
            "small(t+1) = intercept + phi*small(t) + beta*large(t) + error(t+1)"
        ),
        "null_hypothesis": "beta = 0",
        "inference": {
            "covariance": "Newey-West HAC with Bartlett kernel",
            "hac_max_lags": HAC_LAGS,
            "individual_test": "Two-sided asymptotic z test",
            "multiple_testing": "Holm adjustment across the six decade tests",
            "familywise_interval": (
                "Bonferroni simultaneous 95% interval across six decades"
            ),
            "familywise_normal_critical_value": familywise_critical,
        },
        "periods": periods,
        "validation": (
            "Every OLS coefficient and Newey-West covariance matrix matches an "
            "independent matrix implementation."
        ),
    }


def format_p(value):
    return "<0.0001" if value < 0.0001 else f"{value:.4f}"


def write_report(results):
    lines = [
        "# Do Large-Cap Stocks Lead Small-Cap Stocks?",
        "",
        "## Method",
        "",
        "For each decade, the test estimates the daily predictive regression "
        "`small(t+1) = intercept + phi*small(t) + beta*large(t) + error(t+1)`. "
        "Controlling for `small(t)` asks whether the large-cap return adds "
        "information beyond the small-cap portfolio's own one-day persistence.",
        "",
        "The coefficient uses log total returns. A beta of 0.10 means that a 1% "
        "large-cap return today predicts about an additional 0.10% small-cap "
        "return on the next trading day, conditional on today's small-cap return. "
        "Inference uses Newey-West HAC standard errors with five lags, two-sided "
        "tests, and Holm adjustment across the six decades.",
        "",
        "## Results",
        "",
        "| Decade | Beta on large(t) | HAC SE | Family-wise 95% interval | Holm p | Incremental R² | Evidence large leads small? |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for period in results["periods"]:
        row = period["large_t_to_small_t_plus_1"]
        low, high = row["familywise_95pct_bonferroni_interval"]
        answer = "Yes" if row["reject_zero_holm_5pct"] and row["predictor_coefficient"] > 0 else "No"
        lines.append(
            f'| {period["decade"]} | {row["predictor_coefficient"]:+.4f} | '
            f'{row["hac_standard_error"]:.4f} | [{low:+.4f}, {high:+.4f}] | '
            f'{format_p(row["holm_p_value_across_decades"])} | '
            f'{row["incremental_r_squared"]:.4%} | {answer} |'
        )
    lines.extend([
        "",
        "![Large-cap return coefficient by decade with family-wise intervals](large_to_small_coefficients.png)",
        "",
        "## Interpretation",
        "",
        "After controlling for the small-cap portfolio's own lag and correcting "
        "for six decade tests, large-cap returns have positive next-day predictive "
        "content for small-cap returns in the 1970s and 1980s. The other decades "
        "do not reject a zero incremental effect at the 5% family-wise level.",
        "",
        "The 1970s coefficient is smaller but precisely estimated. The 1980s effect "
        "is both larger and more economically meaningful by incremental R². This "
        "matches the earlier variance-ratio evidence that older small-cap returns "
        "contained substantial positive serial dependence.",
        "",
        "As a direction check, the reverse regressions test whether small-cap(t) "
        "predicts large-cap(t+1), controlling for large-cap(t). No reverse test is "
        "significant after Holm correction. That asymmetry supports a large-to-small "
        "lead in the 1970s and 1980s, although it does not establish economic causality.",
        "",
        "Statistical predictability is not automatically a profitable strategy: the "
        "test does not subtract trading costs, and older small-stock prices were more "
        "exposed to infrequent trading and stale-price effects.",
        "",
        "## Reverse-direction diagnostic",
        "",
        "| Decade | Beta on small(t) | Holm p | Significant after correction? |",
        "|---|---:|---:|---|",
    ])
    for period in results["periods"]:
        row = period["small_t_to_large_t_plus_1_diagnostic"]
        lines.append(
            f'| {period["decade"]} | {row["predictor_coefficient"]:+.4f} | '
            f'{format_p(row["holm_p_value_across_decades"])} | '
            f'{"Yes" if row["reject_zero_holm_5pct"] else "No"} |'
        )
    lines.extend([
        "",
        "## Data and validation",
        "",
        "The portfolios are the Kenneth French value-weighted bottom 30% and top "
        "30% size portfolios. They use NYSE size breakpoints and include eligible "
        "NYSE, AMEX, and Nasdaq stocks.",
        "",
        "Every regression coefficient and HAC covariance matrix matches an independent "
        "matrix implementation. Full-precision results, dates, model settings, and the "
        "source hash are stored in `large_to_small_lead_lag.json`.",
        "",
        "Reproduce with `.venv/bin/python analyze_large_leads_small.py`; no network "
        "request is needed.",
    ])
    (OUT / "report.md").write_text("\n".join(lines) + "\n")


def write_plot(results):
    labels = [row["decade"] for row in results["periods"]]
    estimates = np.asarray([
        row["large_t_to_small_t_plus_1"]["predictor_coefficient"]
        for row in results["periods"]
    ])
    intervals = np.asarray([
        row["large_t_to_small_t_plus_1"]["familywise_95pct_bonferroni_interval"]
        for row in results["periods"]
    ])
    significant = np.asarray([
        row["large_t_to_small_t_plus_1"]["reject_zero_holm_5pct"]
        for row in results["periods"]
    ])
    x = np.arange(len(labels))
    yerr = np.vstack([estimates - intervals[:, 0], intervals[:, 1] - estimates])

    fig, ax = plt.subplots(figsize=(11, 6.5))
    colors = np.where(significant, "#c75b2a", "#24748d")
    for i in range(len(labels)):
        ax.errorbar(
            x[i], estimates[i], yerr=yerr[:, i:i+1], fmt="o", markersize=8,
            color=colors[i], ecolor=colors[i], elinewidth=2, capsize=6,
        )
    ax.axhline(0.0, color="#64748b", linestyle="--", linewidth=1.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Coefficient on large-cap return at t")
    ax.set_xlabel("Decade")
    ax.set_title("Does the Large-Cap Return Predict the Next-Day Small-Cap Return?")
    ax.grid(axis="y", alpha=0.25)
    ax.text(
        0.0, -0.19,
        "Model controls for small-cap return at t. Bars are family-wise 95% Bonferroni intervals; orange points reject after Holm correction.",
        transform=ax.transAxes, fontsize=9, color="#475569",
    )
    fig.subplots_adjust(bottom=0.22, left=0.11, right=0.98, top=0.88)
    fig.savefig(OUT / "large_to_small_coefficients.png", dpi=200)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = calculate_results()
    (OUT / "large_to_small_lead_lag.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    write_report(results)
    write_plot(results)
    for period in results["periods"]:
        row = period["large_t_to_small_t_plus_1"]
        print(
            f'{period["decade"]}: beta={row["predictor_coefficient"]:+.6f}, '
            f'Holm p={row["holm_p_value_across_decades"]:.6f}, '
            f'incremental R2={row["incremental_r_squared"]:.6f}'
        )


if __name__ == "__main__":
    main()
