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
