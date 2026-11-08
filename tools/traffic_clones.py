#!/usr/bin/env python3
"""Accumulate GitHub clone counts into an all-time total for the README badge.

GitHub's traffic API (GET /repos/{owner}/{repo}/traffic/clones) only reports
the last 14 days, so an all-time count has to be built up by sampling it more
often than that and keeping every day seen.  This script merges one API
response into a history file and writes a shields.io endpoint file:

    traffic_clones.py <api-response.json> <data-dir>

  <data-dir>/clones.json        {"days": {"<ISO day>": {"count": n, "uniques": u}}}
  <data-dir>/clones-badge.json  shields.io endpoint schema (label "downloads")

A day already in the history is overwritten by the API's figure: the current
day's count is partial and only grows, and earlier days do not change.  The
script is idempotent, so running it twice on the same response is harmless.
Run by .github/workflows/traffic.yml, which keeps <data-dir> on the
`traffic-stats` branch.
"""
import json
import os
import sys


def human(n):
    """2512 -> '2.5k', 1234567 -> '1.2M' (the shields.io 'downloads' style)."""
    for div, suffix in ((1_000_000, "M"), (1_000, "k")):
        if n >= div:
            s = f"{n / div:.1f}".rstrip("0").rstrip(".")
            return s + suffix
    return str(n)


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    api_path, data_dir = argv[1], argv[2]
    with open(api_path) as f:
        api = json.load(f)
    if "clones" not in api:  # an error body (bad token, rate limit): keep history
        sys.exit(f"traffic API response has no 'clones': {api}")

    hist_path = os.path.join(data_dir, "clones.json")
    try:
        with open(hist_path) as f:
            hist = json.load(f)
    except FileNotFoundError:
        hist = {"days": {}}

    for day in api["clones"]:
        hist["days"][day["timestamp"]] = {"count": day["count"],
                                          "uniques": day["uniques"]}

    total = sum(d["count"] for d in hist["days"].values())
    hist["total"] = total
    hist["days"] = dict(sorted(hist["days"].items()))

    os.makedirs(data_dir, exist_ok=True)
    with open(hist_path, "w") as f:
        json.dump(hist, f, indent=1)
        f.write("\n")
    with open(os.path.join(data_dir, "clones-badge.json"), "w") as f:
        json.dump({"schemaVersion": 1, "label": "downloads",
                   "message": human(total), "color": "blue"}, f)
        f.write("\n")
    print(f"{len(hist['days'])} days, {total} clones")


if __name__ == "__main__":
    main(sys.argv)
