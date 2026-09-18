"""Cross-country GICS-sector variance-ratio comparison using frozen JKP files.

Run: .venv/bin/python analyze_global_market_sectors.py

The JKP GICS industry files are monthly and are licensed CC BY-NC 4.0. This
script keeps only classified developed, emerging, and frontier markets; it
requires a complete monthly return series inside each window. No network calls
are made.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import date
import csv
import hashlib
import json
import math
from pathlib import Path
import zipfile
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from arch.unitroot import VarianceRatio

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data/global_industry"
RAW = DATA / "raw"
OUT = ROOT / "results/global_market_efficiency"
HORIZONS = (2, 3, 6)  # months; compatible with the shortest 24-month window
PERIODS = {
    "2008-2010": (2008, 2010),
    "2020-2021": (2020, 2021),
    "2022-2025": (2022, 2025),
}
GROUP_LABELS = {
    "developed": "Developed",
    "emerging": "Developing (MSCI Emerging proxy)",
    "frontier": "Frontier",
}
SECTORS = {
    "10": "Energy", "15": "Materials", "20": "Industrials",
    "25": "Consumer Discretionary", "30": "Consumer Staples",
    "35": "Health Care", "40": "Financials",
    "45": "Information Technology", "50": "Communication Services",
    "55": "Utilities", "60": "Real Estate",
}

# Familiar broad local-market index counterparts. These examples are context;
# the analysis below uses the matched JKP country-sector portfolio returns.
BENCHMARKS = [
    ("developed", "USA", "United States", "S&P 500", "https://www.spglobal.com/spdji/en/indices/equity/sp-500/"),
    ("developed", "CAN", "Canada", "S&P/TSX Composite", "https://www.tsx.com/en/listings/tsx-and-tsxv-issuer-resources/tsx-issuer-resources/sp-tsx-index-eligibility"),
    ("developed", "JPN", "Japan", "TOPIX", "https://www.jpx.co.jp/english/markets/indices/topix/"),
    ("emerging", "IND", "India", "NIFTY 50", "https://www.nseindia.com/static/products-services/indices-nifty50-index"),
    ("emerging", "BRA", "Brazil", "Ibovespa", "https://www.b3.com.br/en_us/market-data-and-indices/indices/broad-indices/ibovespa.htm"),
    ("emerging", "ZAF", "South Africa", "FTSE/JSE All Share", "https://www.jse.co.za/headline"),
    ("frontier", "VNM", "Vietnam", "VN-Index", "https://ssc.gov.vn/webcenter/portal/ssc/pages_r/l/chitit?dDocName=APPSSCGOVVN1620162951"),
    ("frontier", "PAK", "Pakistan", "KSE-100", "https://www.psx.com.pk/psx/product-and-services/indices"),
    ("frontier", "MAR", "Morocco", "MASI", "https://www.casablanca-bourse.com/en/live-market/indices"),
]


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_classifications(path: Path) -> dict[str, str]:
    """Read the small JKP XLSX classification file with the stdlib XML tools."""
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = ["".join(t.text or "" for t in item.findall(".//x:t", ns))
                   for item in shared_root.findall("x:si", ns)]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    mapping = {}
    for row in sheet.findall(".//x:sheetData/x:row", ns)[1:]:
        values = {}
        for cell in row.findall("x:c", ns):
            col = "".join(c for c in cell.attrib["r"] if c.isalpha())
            val = cell.find("x:v", ns)
            if val is None:
                continue
            values[col] = strings[int(val.text)] if cell.attrib.get("t") == "s" else val.text
        group = values.get("B", "").lower()
        code = values.get("A", "").upper()
        if group in GROUP_LABELS and code:
            mapping[code] = group
    return mapping


def month_ends(start_year: int, end_year: int) -> list[str]:
    out = []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            last_day = (date(year + (month == 12), month % 12 + 1, 1) - date.resolution).day
            out.append(f"{year:04d}-{month:02d}-{last_day:02d}")
    return out


def read_country_file(path: Path) -> dict[str, dict[str, float]]:
    with zipfile.ZipFile(path) as archive:
        csv_name = next(name for name in archive.namelist() if name.endswith(".csv"))
        records = csv.DictReader(archive.read(csv_name).decode("utf-8-sig").splitlines())
        result: dict[str, dict[str, float]] = defaultdict(dict)
        for row in records:
            sector = str(row["gics"])
            month = row["date"]
            if sector not in SECTORS:
                continue
            try:
                ret = float(row["ret"])
                n_stocks = int(row["n"])
            except (TypeError, ValueError):
                continue
            if n_stocks > 0 and math.isfinite(ret):
                result[sector][month] = ret
        return result


def variance_ratios(log_returns: np.ndarray) -> list[float]:
    """Lo–MacKinlay overlapping, debiased estimates, matching arch convention."""
    n = len(log_returns)
    mu = float(log_returns.mean())
    centered = log_returns - mu
    variance_1 = float(centered @ centered) / (n - 1)
    if not math.isfinite(variance_1) or variance_1 <= 0:
        raise ValueError("Zero or non-finite one-month variance")
    log_prices = np.r_[0.0, np.cumsum(log_returns)]
    ratios = []
    for q in HORIZONS:
        q_returns = log_prices[q:] - log_prices[:-q]
        variance_q = float(np.sum((q_returns - q * mu) ** 2)) / (
            q * (n - q + 1) * (1 - q / n)
        )
        ratios.append(variance_q / variance_1)
    return ratios


def validate_vr(series: np.ndarray, values: list[float]) -> None:
    y = np.r_[0.0, np.cumsum(series)]
    for q, value in zip(HORIZONS, values):
        test = VarianceRatio(y, lags=q, trend="c", robust=True,
                             overlap=True, debiased=True)
        np.testing.assert_allclose(value, test.vr, rtol=1e-10, atol=1e-12)


def analyze() -> dict:
    class_path = DATA / "jkp_country_classification.xlsx"
    classes = read_classifications(class_path)
    source_files = sorted(RAW.glob("*_gics_monthly_vw.zip"))
    countries = {p.name[:3] for p in source_files}
    classified = {c: classes[c] for c in countries if c in classes}
    assert set(classified.values()) == set(GROUP_LABELS)
    grouped_codes = {g: sorted(c for c, value in classified.items() if value == g)
                     for g in GROUP_LABELS}
    # The extracted files are the intersection of GICS data availability and
    # the classification snapshot (not every market in the source's 93-country set).
    assert {g: len(v) for g, v in grouped_codes.items()} == {
        "developed": 23, "emerging": 23, "frontier": 18
    }

    country_rows = []
    coverage = {g: {} for g in GROUP_LABELS}
    validation_count = 0
    hashes = {"classification_xlsx_sha256": source_hash(class_path),
              "country_gics_zip_sha256": {}}
    for path in source_files:
        iso3 = path.name[:3]
        group = classified.get(iso3)
        if not group:
            continue
        hashes["country_gics_zip_sha256"][iso3] = source_hash(path)
        panel = read_country_file(path)
        for period, (start_year, end_year) in PERIODS.items():
            expected = month_ends(start_year, end_year)
            coverage[group].setdefault(period, {})
            for code, observations in panel.items():
                if not all(day in observations for day in expected):
                    continue
                values = np.asarray([observations[day] for day in expected], dtype=float)
                if not np.isfinite(values).all() or np.any(values <= -1):
                    continue
                log_returns = np.log1p(values)
                ratios = variance_ratios(log_returns)
                # Check every usable series against the published arch implementation.
                validate_vr(log_returns, ratios)
                validation_count += len(HORIZONS)
                score = 1.0 - float(np.mean(np.abs(np.asarray(ratios) - 1.0)))
                row = {
                    "market_class": group,
                    "market_class_label": GROUP_LABELS[group],
                    "country_iso3": iso3,
                    "gics_code": code,
                    "sector": SECTORS[code],
                    "period": period,
                    "monthly_returns": len(values),
                    "first_return_month": expected[0],
                    "last_return_month": expected[-1],
                    "efficiency_score": score,
                    "inefficiency_distance": 1 - score,
                }
                row.update({f"vr_{q}m": value for q, value in zip(HORIZONS, ratios)})
                country_rows.append(row)
                coverage[group][period].setdefault(code, set()).add(iso3)

    # First average country-sector scores equally inside each class and sector;
    # then average sector means equally for the class-level aggregate.
    group_sector = []
    for group in GROUP_LABELS:
        for period in PERIODS:
            codes = sorted({r["gics_code"] for r in country_rows
                            if r["market_class"] == group and r["period"] == period})
            for code in codes:
                rows = [r for r in country_rows if r["market_class"] == group
                        and r["period"] == period and r["gics_code"] == code]
                record = {
                    "market_class": group,
                    "market_class_label": GROUP_LABELS[group],
                    "period": period,
                    "gics_code": code,
                    "sector": SECTORS[code],
                    "country_count": len(rows),
                    "countries_iso3": ",".join(sorted(r["country_iso3"] for r in rows)),
                    "mean_efficiency_score": float(np.mean([r["efficiency_score"] for r in rows])),
                    "population_sd_country_scores": float(np.std([r["efficiency_score"] for r in rows], ddof=0)),
                }
                for q in HORIZONS:
                    record[f"mean_vr_{q}m"] = float(np.mean([r[f"vr_{q}m"] for r in rows]))
                group_sector.append(record)

    summary = {"group": {}, "balanced_common_sector_set": {}}
    for group in GROUP_LABELS:
        summary["group"][group] = {}
        for period in PERIODS:
            rows = [r for r in group_sector if r["market_class"] == group and r["period"] == period]
            sector_scores = [r["mean_efficiency_score"] for r in rows]
            sector_codes = [r["gics_code"] for r in rows]
            country_codes = sorted({r["country_iso3"] for r in country_rows
                                    if r["market_class"] == group and r["period"] == period})
            summary["group"][group][period] = {
                "countries_in_classification_snapshot": len(grouped_codes[group]),
                "countries_with_at_least_one_complete_sector": len(country_codes),
                "countries_with_at_least_one_complete_sector_iso3": country_codes,
                "complete_country_sector_cells": sum(r["country_count"] for r in rows),
                "sectors_represented": len(rows),
                "sector_codes": sector_codes,
                "mean_of_sector_means": float(np.mean(sector_scores)) if sector_scores else None,
                "population_sd_across_sector_means": float(np.std(sector_scores, ddof=0)) if sector_scores else None,
            }
    common_codes = set(SECTORS)
    for group in GROUP_LABELS:
        for period in PERIODS:
            common_codes &= {r["gics_code"] for r in group_sector
                             if r["market_class"] == group and r["period"] == period}
    common_codes = sorted(common_codes)
    summary["balanced_common_sector_set"]["gics_codes"] = common_codes
    summary["balanced_common_sector_set"]["sectors"] = [SECTORS[c] for c in common_codes]
    summary["balanced_common_sector_set"]["periods"] = {}
    for period in PERIODS:
        summary["balanced_common_sector_set"]["periods"][period] = {}
        for group in GROUP_LABELS:
            rows = [r for r in group_sector if r["market_class"] == group
                    and r["period"] == period and r["gics_code"] in common_codes]
            values = [r["mean_efficiency_score"] for r in rows]
            summary["balanced_common_sector_set"]["periods"][period][group] = {
                "mean_of_common_sector_means": float(np.mean(values)) if values else None,
                "population_sd_across_common_sector_means": float(np.std(values, ddof=0)) if values else None,
                "sector_means": {r["sector"]: r["mean_efficiency_score"] for r in rows},
                "country_counts_by_sector": {r["sector"]: r["country_count"] for r in rows},
            }

    for benchmark in BENCHMARKS:
        group, iso3, *_ = benchmark
        assert classified.get(iso3) == group

    return {
        "title": "Global Market Efficiency by Sector",
        "source": "Jensen, Kelly, and Pedersen (2023) Global Factor Data, GICS monthly value-weighted country-industry returns",
        "source_url": "https://www.jkpfactors.com/data",
        "data_license": "Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC 4.0)",
        "data_updated_through": "2025-12",
        "requested_latest_window": "2022-2026",
        "available_latest_window": "2022-2025",
        "market_classification_source": "Downloaded JKP country-classification XLSX; msci_development field, used as a fixed classification snapshot for all windows",
        "groups": GROUP_LABELS,
        "country_counts_in_gics_data_classification_intersection": {g: len(v) for g, v in grouped_codes.items()},
        "periods_monthly": {p: {"first": month_ends(*years)[0], "last": month_ends(*years)[-1], "n": len(month_ends(*years))} for p, years in PERIODS.items()},
        "horizons_months": list(HORIZONS),
        "efficiency_score": "1 - mean(abs(VR(q)-1)) over q=2,3,6 months",
        "aggregation": "Equal-weight country efficiency scores within each market-class x sector x period; then equal-weight sector means within each class-period. Country portfolios are value-weighted in the source data.",
        "complete_monthly_coverage_rule": "All calendar month-end observations in a window must exist with n>0 and finite return; no gaps are imputed.",
        "country_sector_cells_validated_against_arch_unitroot": validation_count,
        "source_hashes": hashes,
        "coverage": {g: {p: {c: len(ids) for c, ids in cov.items()}
                          for p, cov in per.items()} for g, per in coverage.items()},
        "summary": summary,
        "representative_benchmark_indices": [
            {"classification_group": g, "country_iso3": code, "country": country,
             "broad_market_index": index, "official_source": url}
            for g, code, country, index, url in BENCHMARKS
        ],
        "country_sector_scores": country_rows,
        "group_sector_means": group_sector,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_plot(result: dict) -> None:
    group_sector = result["group_sector_means"]
    periods = list(PERIODS)
    fig, axes = plt.subplots(1, 3, figsize=(16, 7), sharey=True)
    colors = ["#237a95", "#c55d31", "#7a5cab"]
    for ax, group in zip(axes, GROUP_LABELS):
        # Keep the same 11-sector axis in all three panels. A missing sector
        # remains visibly blank rather than shifting labels between groups.
        sector_codes = sorted(SECTORS)
        labels = [SECTORS[c] for c in sector_codes]
        y = np.arange(len(labels))
        for i, period in enumerate(periods):
            look = {(r["gics_code"]): r["mean_efficiency_score"] for r in group_sector
                    if r["market_class"] == group and r["period"] == period}
            vals = [look.get(c, np.nan) for c in sector_codes]
            ax.barh(y + (i - 1) * .24, vals, height=.21, color=colors[i], label=period)
        ax.set_title(GROUP_LABELS[group], loc="left", fontsize=13, fontweight="bold")
        ax.set_xlabel("Mean sector efficiency score")
        ax.set_yticks(y, labels)
        ax.axvline(1, color="#566273", linestyle="--", linewidth=1)
        ax.grid(axis="x", alpha=.15)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_axisbelow(True)
        ax.invert_yaxis()
    axes[0].set_ylabel("")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper center", bbox_to_anchor=(.5, .99), ncol=3)
    fig.suptitle("Country-averaged GICS-sector variance-ratio scores", x=.06, y=1.05,
                 ha="left", fontsize=17, fontweight="bold")
    fig.text(.06, -.03, "Higher scores mean VR estimates are closer to 1. Monthly horizons: 2, 3, and 6 months.\nCountry scores are averaged equally within each sector; incomplete country-sector histories are omitted.", fontsize=9, color="#4a5565")
    fig.tight_layout(rect=(0, .04, 1, .92), w_pad=2.2)
    fig.savefig(OUT / "group_sector_scores.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def render_report(result: dict) -> None:
    lines = [
        "# Global Market Efficiency by Sector", "",
        "## What was compared", "",
        "The analysis groups country GICS industry portfolios into Developed, Developing (the source's MSCI Emerging category), and Frontier. It uses every market in the intersection of the JKP classification file and available GICS files: 23 developed, 23 emerging, and 18 frontier markets. A country-sector-period observation is included only when every monthly return in that period exists and the portfolio has at least one stock. No missing returns are filled.", "",
        "For well-known national broad-market benchmarks, this report includes representative equivalents: S&P 500 (United States), S&P/TSX Composite (Canada), TOPIX (Japan), NIFTY 50 (India), Ibovespa (Brazil), FTSE/JSE All Share (South Africa), VN-Index (Vietnam), KSE-100 (Pakistan), and MASI (Morocco). These index names help orient the country comparison; the calculated sector scores use the JKP country-sector return portfolios, not those named headline index price series.", "",
        "The JKP source provides GICS industry returns for all countries, but these industry files are monthly. They are value-weighted sector portfolios, not daily headline benchmark index closes. Because the latest data end in December 2025, the requested 2022–2026 period is represented by 2022–2025 only; no 2026 returns are present.", "",
        "## Method", "",
        "For each country, sector, and window, monthly simple portfolio returns are converted to log returns. We estimate overlapping, finite-sample-debiased Lo–MacKinlay variance ratios at q = 2, 3, and 6 months. With random-walk variance scaling, VR(q) is 1. The descriptive score is E = 1 − mean(|VR(q) − 1|); higher scores are closer to that restriction. Scores may be below zero and are not probabilities or formal tests of market efficiency.", "",
        "First, country scores are averaged equally within each market group and GICS sector. Then sector means are averaged equally to obtain a group-period mean and the population standard deviation across sector means. Value weighting happens inside the source's country portfolios; countries receive equal weight when averaging their sector scores. The standard deviation is descriptive and does not itself test whether sector efficiencies differ statistically.", "",
        "Monthly windows contain only 36 observations (2008–2010), 24 (2020–2021), and 48 (2022–2025). In particular, the two-year window gives noisy variance-ratio estimates. The horizon set is adapted to monthly data and is not numerically comparable with the earlier U.S. daily q = 2, 5, 10, 20-day scores.", "",
        "## Group averages and across-sector dispersion", "",
        "| Market group | Period | Countries with ≥1 complete sector | Country-sector cells | Sectors represented | Mean of sector means | SD across sector means |", "|---|---|---:|---:|---:|---:|---:|"]
    for group in GROUP_LABELS:
        for period in PERIODS:
            row = result["summary"]["group"][group][period]
            lines.append(f"| {GROUP_LABELS[group]} | {period} | {row['countries_with_at_least_one_complete_sector']} | {row['complete_country_sector_cells']} | {row['sectors_represented']} | {row['mean_of_sector_means']:.6f} | {row['population_sd_across_sector_means']:.6f} |")
    lines += ["", "These scores average the country results inside each sector before aggregating across sectors. Frontier estimates use fewer countries and sectors, especially in 2008–2010, so the group summaries do not have identical underlying coverage.", "", "## Sector means by class", "", "Each value is the equal-country mean score for that class, sector, and period. The `N` column is the number of complete country-sector return histories contributing to that cell. Blank cells have no complete history and are not zero scores.", ""]
    for group in GROUP_LABELS:
        lines += [f"### {GROUP_LABELS[group]}", "", "| GICS sector | N 2008–2010 | Score 2008–2010 | N 2020–2021 | Score 2020–2021 | N 2022–2025 | Score 2022–2025 |", "|---|---:|---:|---:|---:|---:|---:|"]
        # Keep the same 11-sector axis in all three panels. A missing sector
        # remains visibly blank rather than shifting labels between groups.
        sector_codes = sorted(SECTORS)
        for code in sector_codes:
            cells = {(r["period"]): r for r in result["group_sector_means"]
                     if r["market_class"] == group and r["gics_code"] == code}
            items = []
            for period in PERIODS:
                row = cells.get(period)
                items.extend([str(row["country_count"]) if row else "—",
                              f"{row['mean_efficiency_score']:.6f}" if row else "—"])
            lines.append(f"| {SECTORS[code]} | " + " | ".join(items) + " |")
        lines.append("")
    balanced = result["summary"]["balanced_common_sector_set"]
    common = balanced["gics_codes"]
    lines += ["## Like-for-like sector subset", "", f"To hold sector composition fixed, the script finds GICS sectors with data in every group and window. The common set is {', '.join(SECTORS[c] for c in common) if common else 'empty'}. This balanced subset is still sparse in Frontier markets; inspect country counts below.", "", "| Period | Group | Mean of common-sector means | SD across common-sector means | Country N by sector |", "|---|---|---:|---:|---|"]
    for period in PERIODS:
        for group in GROUP_LABELS:
            row = balanced["periods"][period][group]
            counts = ", ".join(f"{name}: {n}" for name, n in row["country_counts_by_sector"].items()) or "—"
            mean_value = row["mean_of_common_sector_means"]
            sd_value = row["population_sd_across_common_sector_means"]
            lines.append(f"| {period} | {GROUP_LABELS[group]} | {mean_value:.6f} | {sd_value:.6f} | {counts} |")
    lines += ["", "## Results", "", "![Country-averaged sector scores by market class](group_sector_scores.png)", "", "The country-level scores and variance ratios are in `country_sector_scores.csv`; equal-country sector means and contributing country counts are in `group_sector_means.csv`. The JSON includes source file hashes, country coverage counts, complete source metadata, and the dynamically determined balanced sector set.", "", "## Data and limitations", "", "The data are Jensen, Kelly, and Pedersen (2023) Global Factor Data, specifically monthly value-weighted GICS industry portfolio returns. GICS sectors are available across countries, whereas the source's Fama–French 49-industry portfolios are U.S.-only. Data are updated through December 2025. The released source dataset is CC BY-NC 4.0; attribution and non-commercial terms apply.", "", "The classification is the `msci_development` field in the downloaded JKP country-classification workbook and is held fixed across historical periods. The groups are therefore a consistent analysis bucket, not a reconstruction of each country's historical MSCI classification. Sector definitions and membership can also change over time; this analysis uses GICS labels as delivered in the source data. Frontier coverage is particularly thin in the first window; country and sector counts are shown rather than imputing missing markets. Broad-market benchmark index names are examples, not return inputs for this calculation.", "", "Reproduce from the project root with `.venv/bin/python analyze_global_market_sectors.py`. The script reads only the frozen files in `data/global_industry/`, makes no network calls, and validates every usable VR estimate against `arch.unitroot.VarianceRatio`.", "", "Sources: [JKP Global Factor Data download and coverage](https://www.jkpfactors.com/data); [JKP release update through December 2025](https://github.com/bkelly-lab/jkp-data/discussions/89); [Lo and MacKinlay (1988)](https://web.mit.edu/~alo/www/Papers/lo-mackinlay-88.html); [S&P/TSX Composite description](https://www.tsx.com/en/listings/tsx-and-tsxv-issuer-resources/tsx-issuer-resources/sp-tsx-index-eligibility); [TOPIX description](https://www.jpx.co.jp/english/markets/indices/topix/); [NIFTY 50](https://www.nseindia.com/static/products-services/indices-nifty50-index); [Ibovespa](https://www.b3.com.br/en_us/market-data-and-indices/indices/broad-indices/ibovespa.htm); [FTSE/JSE All Share](https://www.jse.co.za/headline); [VN-Index data](https://ssc.gov.vn/webcenter/portal/ssc/pages_r/l/chitit?dDocName=APPSSCGOVVN1620162951); [KSE-100](https://www.psx.com.pk/psx/product-and-services/indices); [MASI](https://www.casablanca-bourse.com/en/live-market/indices)."]
    (OUT / "global_market_efficiency_by_sector.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    result = analyze()
    write_csv(OUT / "country_sector_scores.csv", result["country_sector_scores"])
    write_csv(OUT / "group_sector_means.csv", result["group_sector_means"])
    benchmark_rows = [{"market_class": g, "country_iso3": code, "country": country,
                       "broad_market_index": index, "official_source": url}
                      for g, code, country, index, url in BENCHMARKS]
    write_csv(OUT / "benchmark_index_examples.csv", benchmark_rows)
    metadata = {k: v for k, v in result.items()
                if k not in ("country_sector_scores", "group_sector_means")}
    (OUT / "global_market_efficiency_results.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    make_plot(result)
    render_report(result)
    print(json.dumps({
        "market_counts": result["country_counts_in_gics_data_classification_intersection"],
        "country_sector_score_rows": len(result["country_sector_scores"]),
        "group_sector_mean_rows": len(result["group_sector_means"]),
        "validated_vr_estimates": result["country_sector_cells_validated_against_arch_unitroot"],
        "balanced_common_sectors": result["summary"]["balanced_common_sector_set"]["sectors"],
        "summary": result["summary"]["group"],
    }, indent=2))


if __name__ == "__main__":
    main()
