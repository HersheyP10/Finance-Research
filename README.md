# Finance Research

Reproducible S&P 500 variance ratio analysis using Yahoo Finance daily history.

Study windows: 2008–2010, 2020–2021, and 2022–September 15, 2026.

Read [the analysis report](results/report.md) for methods, sample sizes, results, and limitations. Full precision results are in [JSON](results/variance_ratio_results.json).

## Reproduce

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python analyze_variance_ratio.py
```

The script runs against the frozen Yahoo source snapshot in `data/yahoo_gspc_raw.json`; it makes no network requests. Dependencies used in the original run are recorded in `requirements-lock.txt`.

Yahoo endpoint and source provenance are included in the results JSON. September 16, 2026 was an incomplete market session at download time and is excluded. The adjusted-close price index is not a dividends-reinvested total-return index.

## Industry data for follow-up research

[Industry data sources and validation](results/sector_data_sources.md) documents downloaded Kenneth French 12- and 49-industry daily datasets. Both have complete observations for 2008–2010 and 2020–2021. Run `python3 prepare_industry_data.py` to regenerate the selected study windows. These are SIC-based industry portfolios.

## Industry efficiency scores

[US Market Standard Deviation by Sectors](results/industry_efficiency/us_market_standard_deviation_by_sectors.md) contains scores for every industry, cross-industry standard deviations, and a conservative bootstrap uncertainty analysis for the 12-industry panel. The declared score is `E = 1 - mean(abs(VR(q)-1))` across 2-, 5-, 10-, and 20-day horizons. Estimated standard deviations are nonzero, but the uncertainty analysis does not exclude equal underlying scores.

Reproduce with `.venv/bin/python analyze_industry_efficiency.py`. The full-precision JSON includes all individual variance ratio tests and a separate 49-industry supplement. Scores describe proximity to the random-walk variance restriction; they are not probabilities or proof of market efficiency.
