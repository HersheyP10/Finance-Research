"""Validate saved Kenneth French daily TXT archives and select study windows.

Run with Python 3; only the standard library is required. This prepares data,
not statistical tests. Source return values are percentages, not decimals.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import math
import re
import zipfile

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data/industry_portfolios"
WINDOWS = {"2008-2010": ("2008-01-01", "2010-12-31"),
           "2020-2021": ("2020-01-01", "2021-12-31")}


def read_value_weighted(path, number):
    with zipfile.ZipFile(path) as archive:
        names = [n for n in archive.namelist() if n.lower().endswith(".txt")]
        assert len(names) == 1
        text = archive.read(names[0]).decode("utf-8-sig")
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines)
                 if line.strip() == "Average Value Weighted Returns -- Daily")
    columns = lines[start+1].split()
    assert len(columns) == number and len(set(columns)) == number
    records = []
    for line in lines[start+2:]:
        parts = line.split()
        if not parts or not re.fullmatch(r"\d{8}", parts[0]):
            break
        assert len(parts) == number+1
        day = datetime.strptime(parts[0], "%Y%m%d").date().isoformat()
        records.append({"date": day,
                        "returns_pct": dict(zip(columns, map(float, parts[1:])))})
    dates = [r["date"] for r in records]
    assert dates == sorted(set(dates)), "Duplicate or unordered dates"
    return columns, records, lines[0]


def main():
    # The previously downloaded S&P 500 observations independently check the calendar.
    market = json.loads((ROOT / "data/analysis_samples.json").read_text())
    summary = {"prepared_utc": datetime.now(timezone.utc).isoformat(),
               "units": "Daily simple total returns in percent; 1.00 means 1%.",
               "weighting": "Value weighted; equal-weighted alternatives remain in the raw archives.",
               "datasets": []}
    for number in (12, 49):
        path = DATA / f"{number}_Industry_Portfolios_daily_TXT.zip"
        columns, records, header = read_value_weighted(path, number)
        source = f"https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/{path.name}"
        metadata = {"industry_count": number, "source_url": source,
                    "source_header": header, "archive_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                    "first_available_date": records[0]["date"], "last_available_date": records[-1]["date"],
                    "columns": columns, "windows": []}
        selected = {}
        for label, (start, end) in WINDOWS.items():
            panel = [r for r in records if start <= r["date"] <= end]
            expected = sorted({r["date"] for r in market
                               if start <= r["date"] <= end and r["role"] == "in_period_return"})
            assert [r["date"] for r in panel] == expected, "Study-window calendar mismatch"
            missing = sum(v in (-99.99, -999) or not math.isfinite(v)
                          for r in panel for v in r["returns_pct"].values())
            assert missing == 0, "Missing observations must be investigated"
            assert all(v > -100 for r in panel for v in r["returns_pct"].values())
            selected[label] = panel
            metadata["windows"].append({"period": label, "first_return": panel[0]["date"],
                                         "last_return": panel[-1]["date"], "dates_per_industry": len(panel),
                                         "missing_return_cells": missing,
                                         "matches_sp500_trading_dates": True})
        (DATA / f"{number}_industry_value_weighted_study_windows.json").write_text(
            json.dumps({"metadata": metadata, "units": summary["units"], "data": selected}, indent=2)+"\n")
        summary["datasets"].append(metadata)
    (DATA / "validation_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
