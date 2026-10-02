"""Check that the 0.143 crossing of the depositor's FSC curve read from a validation record ("load_fsc", as used by
code/sample_unindexed_va.py) equals the index field author_resolution_fsc_0143_value.

Uses the same entries as data/results/index_vs_va_check.json (code/check_index_against_va.py). For each, the record
https://www.ebi.ac.uk/emdb/api/analysis/<id> is downloaded again and the crossing is compared with the index value
under the same agreement rule as that check (difference at most 2% plus 0.01 A). The record itself is not kept; the
URL, retrieval time, size and SHA-256 are.

Output: data/results/own_curve_field_check.json

    python3 code/check_own_curve_field.py data/raw/emdb_index_2026-10-01.csv
"""
import csv, hashlib, json, os, sys, time
from datetime import datetime, timezone

import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_v2 import num  # noqa: E402
from sample_unindexed_va import UA, crossing  # noqa: E402

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
RES = os.path.join(ROOT, "data", "results")


def agree(a, b):
    return abs(a - b) <= 0.02 * b + 0.01


def main(index_csv):
    ids = [e["emdb_id"] for e in json.load(open(os.path.join(RES, "index_vs_va_check.json")))["entries"]]
    idx = {r["emdb_id"]: r for r in csv.DictReader(open(index_csv, newline="", encoding="utf-8"))}
    s, out = requests.Session(), []
    for eid in ids:
        num_id = eid.split("-")[1]
        url = f"https://www.ebi.ac.uk/emdb/api/analysis/{num_id}"
        rec = {"emdb_id": eid, "url": url, "index_value": num(idx[eid]["author_resolution_fsc_0143_value"])}
        for attempt in range(4):
            try:
                resp = s.get(url, headers=UA, timeout=180)
                rec.update({"http_status": resp.status_code, "bytes": len(resp.content),
                            "sha256": hashlib.sha256(resp.content).hexdigest(),
                            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")})
                if resp.status_code == 200:
                    d = resp.json()
                    va = d.get(num_id, d)       # same two record layouts as in sample_unindexed_va.py
                    inner = va.get(num_id) if isinstance(va.get(num_id), dict) and "fsc" in va.get(num_id) else va
                    rec["record_value"] = crossing(inner.get("load_fsc"), "0.143")
                    rec.pop("error", None)
                    break
            except (requests.RequestException, ValueError) as e:
                rec["error"] = repr(e)[:200]
            time.sleep(10 * (attempt + 1))
        iv, rv = rec.get("index_value"), rec.get("record_value")
        rec["outcome"] = ("failed" if rec.get("http_status") != 200 else "both_agree" if iv and rv and agree(rv, iv)
                          else "both_differ" if iv and rv else "index_only" if iv else "record_only" if rv else "neither")
        out.append(rec)
        print(eid, rec["outcome"], iv, rv, flush=True)
        time.sleep(0.3)
    counts = {}
    for r in out:
        counts[r["outcome"]] = counts.get(r["outcome"], 0) + 1
    json.dump({"n": len(out), "counts": counts, "rule": "|record - index| <= 0.02 * index + 0.01 A", "entries": out},
              open(os.path.join(RES, "own_curve_field_check.json"), "w"), indent=1)
    print(counts)


if __name__ == "__main__":
    main(sys.argv[1])
