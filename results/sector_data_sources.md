# Data for comparing industry-level return predictability

## Recommended public source

The Kenneth R. French Data Library provides daily US industry-portfolio returns covering both requested windows. Start with the 12-industry value-weighted series for a broad comparison. Use the 49-industry series when separating banking, insurance, real estate, software, hardware, and other narrower industries is important.

- [12-industry methodology](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_12_ind_port.html)
- [12-industry daily TXT archive](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/12_Industry_Portfolios_daily_TXT.zip)
- [49-industry methodology](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_49_ind_port.html)
- [49-industry daily TXT archive](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/49_Industry_Portfolios_daily_TXT.zip)
- [12-industry SIC definitions](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/Siccodes12.zip)
- [Data library and revision notes](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html)

Stocks from NYSE, AMEX, and NASDAQ are assigned using their SIC codes at the annual June formation date. These are industry research portfolios, not the current 11 GICS sectors or portfolios restricted to S&P 500 members. Constituents change over time; using one classification scheme does not imply fixed holdings. The downloaded vintage uses the July 2026 CRSP database and extends through July 31, 2026. Both value-weighted and equal-weighted return tables are retained in the original archives; prepared files select value-weighted returns only.

## Verified study-window coverage

| Dataset | 2008-01-02 to 2010-12-31 | 2020-01-02 to 2021-12-31 | Missing return cells in these windows |
|---|---:|---:|---:|
| 12 industries | 757 daily returns per industry | 505 daily returns per industry | 0 |
| 49 industries | 757 daily returns per industry | 505 daily returns per industry | 0 |

Dates are unique, ordered, and match the corresponding dates in the project's saved S&P 500 sample. All study-window returns are finite, avoid the source's missing-data codes, and exceed -100%. No missing observations were filled. These checks concern the requested windows, not completeness of all industries back to 1926.

## Broad industry groups

| Code | Group |
|---|---|
| NoDur | Consumer nondurables |
| Durbl | Consumer durables |
| Manuf | Manufacturing |
| Enrgy | Energy |
| Chems | Chemicals |
| BusEq | Computers, software, and electronic equipment |
| Telcm | Telephone and television transmission |
| Utils | Utilities |
| Shops | Wholesale, retail, and selected services |
| Hlth | Healthcare, medical equipment, and pharmaceuticals |
| Money | Finance |
| Other | Remaining industries |

These labels summarize the source definitions; use the SIC definition archive for exact membership rules. Do not relabel these groups as equivalent GICS sectors.

## Prepared local files

- [12-industry study samples](../data/industry_portfolios/12_industry_value_weighted_study_windows.json)
- [49-industry study samples](../data/industry_portfolios/49_industry_value_weighted_study_windows.json)
- [Coverage, source URLs, and SHA-256 hashes](../data/industry_portfolios/validation_summary.json)
- [Reproducible preparation script](../prepare_industry_data.py)

Run `python3 prepare_industry_data.py` to regenerate the selected samples from the saved ZIP files. No network access or third-party Python packages are needed. Source return values are percentages: 1.00 means a 1% daily simple return. Preserve the raw archives because source histories may be revised in later releases.

## How the data can answer the research question

For each industry and each window, convert the percentage simple return R_t to a log return with log(1 + R_t/100). Build a log wealth series starting at zero followed by cumulative log returns; pass that series to the existing VarianceRatio implementation. Do not pass the returns directly to an implementation that internally differences its input. These are portfolio total returns, while the original ^GSPC analysis uses a price index, so a direct market-versus-industry comparison needs a consistent return basis.

Use common horizons and robust inference. Adjust for multiple industries and horizons. To assess whether industries differ, add a formal cross-industry comparison that accounts for their common trading dates and correlated returns. Different individual p-values alone do not establish a significant difference. Rejection tests particular random-walk implications, not all forms of market efficiency. The subsequent [US Market Standard Deviation by Sectors](industry_efficiency/us_market_standard_deviation_by_sectors.md) report now provides scores, dispersion, and a conservative uncertainty check.

## Why an all-11-sector ETF panel cannot simply be used for both periods

Yahoo adjusted-close histories of sector ETFs are an alternative for studying tradable sector proxies, but the modern Select Sector SPDR set lacks full common ETF history for 2008. State Street reports inception on October 7, 2015 for [XLRE](https://www.ssga.com/us/en/intermediary/etfs/state-street-real-estate-select-sector-spdr-etf-xlre), and June 18, 2018 for [XLC](https://www.ssga.com/us/en/intermediary/etfs/state-street-communication-services-select-sector-spdr-etf-xlc). These cannot supply actual ETF returns for 2008-2010. A study specifically requiring all 11 GICS sectors would need historical sector-index series and documented classification/back-history treatment. The French portfolios support a clearly labeled industry comparison instead.
