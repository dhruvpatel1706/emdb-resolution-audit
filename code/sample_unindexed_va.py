"""Random sample of entries whose recalculated FSC values are missing from the EMDB search index, read from their
per-entry validation-analysis (VA) records instead.

Eligible: released, not withdrawn or obsolete; single particle, helical or subtomogram averaging; stated criterion
"FSC 0.143 CUT-OFF"; half-maps listed; no indexed unmasked recalculation. Two strata by release date:
  pre2024  released before 2024-01-01 (500 entries)
  from2024 released 2024-01-01 or later (150 entries)
Drawn with random.Random(20261001) from the sorted entry ids of each stratum.

For each entry the record https://www.ebi.ac.uk/emdb/api/analysis/<id> is downloaded, the 0.143 crossings of
RELION's unmasked, masked and corrected curves ("relion_fsc") and of the depositor's curve ("load_fsc"; checked
against the index field author_resolution_fsc_0143_value on 60 indexed entries: 43 agree, 17 missing in both, 0
differ) are extracted, and one JSON line is appended to data/results/unindexed_sample.jsonl with the URL, retrieval
time, HTTP status, size and SHA-256 of the record. The record itself is not kept. Resumable: entries already in the
output are skipped.

    python3 code/sample_unindexed_va.py data/raw/emdb_index_2026-10-01.csv
"""
import csv, hashlib, json, os, random, sys, time
from datetime import datetime, timezone

import requests

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUT = os.path.join(ROOT, "data", "results", "unindexed_sample.jsonl")
UA = {"User-Agent": "emdb-resolution-audit/2.0 (research metadata audit)"}
METHODS = {"singleParticle", "helical", "subtomogramAveraging"}
SIZES = {"pre2024": 500, "from2024": 150}


def num(x):
    try:
        v = float(str(x).split(",")[0])
        return v if v > 0 else None
    except (TypeError, ValueError):
        return None


def crossing(block, key):
    x = (block or {}).get("intersections", {}).get(key, {}).get("x")
    return None if x in (None, 0) else 1.0 / x


def draw(index_csv):
    rows = list(csv.DictReader(open(index_csv, newline="", encoding="utf-8")))
    strata = {"pre2024": [], "from2024": []}
    for r in rows:
        if not r["header_release_date"] or r["obsolete_date"] or r["withdrawn_date"]:
            continue
        if r["structure_determination_method"] not in METHODS or r["resolution_method"] != "FSC 0.143 CUT-OFF":
            continue
        if not num(r["resolution"]) or not r["half_map_filename"] or num(r["calculated_resolution_fsc_0143_value"]):
            continue
        strata["pre2024" if r["header_release_date"][:10] < "2024-01-01" else "from2024"].append(r)
    chosen = []
    for name, pool in strata.items():
        pool.sort(key=lambda r: r["emdb_id"])
        for r in random.Random(20261001).sample(pool, SIZES[name]):
            chosen.append((name, len(pool), r))
    return chosen


def main(index_csv):
    chosen = draw(index_csv)
    done = set()
    if os.path.exists(OUT):
        done = {r["emdb_id"] for r in map(json.loads, filter(str.strip, open(OUT)))
                if r.get("http_status") == 200 and "has_relion_fsc" in r}      # failed fetches are retried
    s = requests.Session()
    for i, (stratum, pool_n, r) in enumerate(chosen):
        eid = r["emdb_id"]
        if eid in done:
            continue
        num_id = eid.split("-")[1]
        url = f"https://www.ebi.ac.uk/emdb/api/analysis/{num_id}"
        rec = {"emdb_id": eid, "stratum": stratum, "stratum_pool_n": pool_n, "url": url,
               "method": r["structure_determination_method"], "release_date": r["header_release_date"][:10],
               "stated": num(r["resolution"]), "halfmap_pixel": num(r["half_map_pixel_spacing_x"]),
               "map_pixel": num(r["map_pixel_spacing_x"])}
        for attempt in range(4):
            try:
                resp = s.get(url, headers=UA, timeout=180)
                rec.update({"http_status": resp.status_code, "bytes": len(resp.content),
                            "sha256": hashlib.sha256(resp.content).hexdigest(),
                            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")})
                if resp.status_code == 200:
                    rec.pop("error", None)
                    d = resp.json()
                    va = d.get(num_id, d)
                    inner = va.get(num_id) if isinstance(va.get(num_id), dict) and "fsc" in va.get(num_id) else va
                    rec.update({"unmasked": crossing(inner.get("relion_fsc"), "0.143"),
                                "masked": crossing(inner.get("relion_fsc"), "masked"),
                                "corrected": crossing(inner.get("relion_fsc"), "corrected"),
                                "own": crossing(inner.get("load_fsc"), "0.143"),
                                "has_relion_fsc": "relion_fsc" in inner})
                    break
                if resp.status_code < 500:
                    break
            except (requests.RequestException, ValueError) as e:
                rec.update({"error": repr(e)[:200]})
            time.sleep(10 * (attempt + 1))
        with open(OUT, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(f"{i + 1}/{len(chosen)} {eid} {rec.get('http_status')} {rec.get('bytes')} corrected={rec.get('corrected')}", flush=True)
        time.sleep(0.3)


if __name__ == "__main__":
    main(sys.argv[1])
