# International GICS industry data

This folder contains frozen files from the [JKP Global Factor Data download](https://www.jkpfactors.com/data), based on Jensen, Kelly, and Pedersen (2023), “Is There a Replication Crisis in Finance?” (*The Journal of Finance*, 78(5), 2465–2518).

The `raw/` archives are country-specific, monthly, value-weighted GICS industry portfolio returns. `jkp_country_classification.xlsx` provides the country-classification snapshot used by the analysis; `jkp_availability.json` records the source's country and data availability. The raw archives were retrieved for this analysis and contain data through December 2025. The GICS sector codes are described in `analyze_global_market_sectors.py` at the project root.

The distributed JKP data are licensed under [Creative Commons Attribution-NonCommercial 4.0](https://creativecommons.org/licenses/by-nc/4.0/). Cite Jensen, Kelly, and Pedersen (2023) and retain the non-commercial restriction when using or redistributing the data. These files are frozen so the report can be reproduced without a network connection; `analyze_global_market_sectors.py` records SHA-256 hashes for its inputs.
