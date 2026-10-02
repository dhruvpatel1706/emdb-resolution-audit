"""Check that the resolution values in the EMDB search index equal those in EMDB's per-entry validation-analysis
(VA) records, for a reproducible sample of entries.

Sample: 40 entries drawn with a fixed seed from all entries that have an author FSC 0.143 resolution and an unmasked
recalculated value, plus the 20 entries whose masked recalculated resolution is furthest above the author's value
(ratio), so that the tail the paper reports is checked directly.

For each entry the VA record is fetched from https://www.ebi.ac.uk/emdb/api/analysis/<id> and the 0.143 crossings of
the unmasked FSC ("fsc" block) and of the masked and phase-randomization-corrected RELION FSC ("relion_fsc" block) are
converted from spatial frequency (1/A) to resolution (A) and compared with the index values. The index's unmasked
value is RELION's unmasked crossing ("relion_fsc" block, key "0.143"); the separate "fsc" block, a second unmasked
calculation, is recorded for reference only. Raw VA files are kept in
data/raw/va_sample/ with a manifest (URL, time, bytes, SHA-256).

    python3 code/check_index_against_va.py data/raw/emdb_index_2026-10-01.csv
"""
import csv, hashlib, json, os, random, sys, time
from datetime import datetime, timezone

import requests

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUT = os.path.join(ROOT, "data", "raw", "va_sample")
RES = os.path.join(ROOT, "data", "results", "index_vs_va_check.json")
UA = {"User-Agent": "emdb-resolution-audit/2.0 (research metadata audit)"}
METHODS = {"singleParticle", "helical", "subtomogramAveraging"}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def crossing(block, key):
    x = (block or {}).get("intersections", {}).get(key, {}).get("x")
    return None if x in (None, 0) else 1.0 / x


def main(index_csv):
    rows = list(csv.DictReader(open(index_csv, newline="", encoding="utf-8")))
    pool = [r for r in rows if r["structure_determination_method"] in METHODS and r["resolution_method"] == "FSC 0.143 CUT-OFF"
            and num(r["resolution"]) and num(r["calculated_resolution_fsc_0143_value"])]
    pool.sort(key=lambda r: r["emdb_id"])
    rnd = random.Random(20261001).sample(pool, 40)
    tail = sorted([r for r in pool if num(r["calculated_resolution_fsc_masked_0143_value"])],
                  key=lambda r: -num(r["calculated_resolution_fsc_masked_0143_value"]) / num(r["resolution"]))[:20]
    chosen = {r["emdb_id"]: ("random" if r in rnd else "tail", r) for r in rnd + tail}
    os.makedirs(OUT, exist_ok=True)
    manifest, results = [], []
    for eid, (why, r) in sorted(chosen.items()):
        num_id = eid.split("-")[1]
        path = os.path.join(OUT, f"{eid}_va.json")
        url = f"https://www.ebi.ac.uk/emdb/api/analysis/{num_id}"
        if not os.path.exists(path):
            resp = requests.get(url, headers=UA, timeout=120)
            assert resp.status_code == 200, (eid, resp.status_code)
            open(path, "wb").write(resp.content)
            manifest.append({"url": url, "file": os.path.relpath(path, ROOT), "bytes": len(resp.content),
                             "sha256": hashlib.sha256(resp.content).hexdigest(),
                             "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")})
            time.sleep(0.5)
        va = json.load(open(path))[num_id]
        inner = va.get(num_id) if isinstance(va.get(num_id), dict) and "fsc" in va.get(num_id) else va   # two record layouts
        got = {"unmasked": crossing(inner.get("relion_fsc"), "0.143"),      # the index takes RELION's unmasked crossing
               "masked": crossing(inner.get("relion_fsc"), "masked"),
               "corrected": crossing(inner.get("relion_fsc"), "corrected")}
        want = {"unmasked": num(r["calculated_resolution_fsc_0143_value"]),
                "masked": num(r["calculated_resolution_fsc_masked_0143_value"]),
                "corrected": num(r["calculated_resolution_fsc_corrected_0143_value"])}
        agree = {k: (got[k] is None and want[k] is None) or (got[k] is not None and want[k] is not None
                                                             and abs(got[k] - want[k]) <= 0.02 * want[k] + 0.01)
                 for k in got}
        results.append({"emdb_id": eid, "sample": why, "author_A": num(r["resolution"]),
                        "va_fsc_block_unmasked_0143": crossing(inner.get("fsc"), "0.143"),
                        "va_author_resolution": (va.get("resolution") or {}).get("value"),
                        "index": want, "va_record": got, "agree": agree})
    if manifest:
        mpath = os.path.join(OUT, "MANIFEST.json")
        old = json.load(open(mpath)) if os.path.exists(mpath) else []
        json.dump(old + manifest, open(mpath, "w"), indent=1)
    summ = {k: sum(x["agree"][k] for x in results) for k in ("unmasked", "masked", "corrected")}
    json.dump({"n": len(results), "agree_counts": summ, "entries": results}, open(RES, "w"), indent=1)
    print("entries checked:", len(results), "agreement counts:", summ)
    for x in results:
        if not all(x["agree"].values()):
            print("  DISAGREE", x["emdb_id"], x["sample"], "index", x["index"], "va", x["va_record"])


if __name__ == "__main__":
    main(sys.argv[1])
