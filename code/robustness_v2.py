"""Robustness of the tail reported by analyze_v2.py (stated resolution finer than the corrected recalculation by
more than 25%).

1. Method within resolution band: is the higher rate for subtomogram averages explained by their coarser resolution?
2. Clustering by study: many entries come from one paper. Entries are grouped by primary citation (PubMed id, else
   DOI, else the normalized citation title; entries without a citation form their own group). Reported: number of
   studies, studies with at least one tail entry, the largest contributors, and a 95% interval for the entry-level
   tail percentage from a bootstrap that resamples studies rather than entries (fixed seed).

    python3 code/robustness_v2.py data/raw/emdb_index_2026-10-01.csv data/raw/emdb_citations_2026-10-01.csv
"""
import csv, json, os, re, sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_v2 import METHODS, BANDS, num, wilson  # noqa: E402

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def study_key(c):
    if c is None:
        return None
    for k in ("xref_PUBMED", "xref_DOI"):
        v = (c.get(k) or "").split(",")[0].strip()
        if v:
            return f"{k}:{v.lower()}"
    t = re.sub(r"[^a-z0-9]+", " ", (c.get("primary_citation_title") or "").lower()).strip()
    return f"title:{t}" if len(t) > 20 else None


def main(index_csv, cit_csv):
    cit = {r["emdb_id"]: r for r in csv.DictReader(open(cit_csv, newline="", encoding="utf-8"))}
    E = []
    for r in csv.DictReader(open(index_csv, newline="", encoding="utf-8")):
        if not r["header_release_date"] or r["obsolete_date"] or r["withdrawn_date"]:
            continue
        if r["structure_determination_method"] not in METHODS or r["resolution_method"] != "FSC 0.143 CUT-OFF":
            continue
        st, un, co = num(r["resolution"]), num(r["calculated_resolution_fsc_0143_value"]), num(r["calculated_resolution_fsc_corrected_0143_value"])
        if not (st and un and co):
            continue
        px = num(r["half_map_pixel_spacing_x"]) or num(r["map_pixel_spacing_x"])
        ma = num(r["calculated_resolution_fsc_masked_0143_value"])
        if px and any(v is not None and v < 0.999 * 2 * px for v in (un, ma, co)):
            continue                                          # same exclusion as analyze_v2.py
        key = study_key(cit.get(r["emdb_id"])) or f"entry:{r['emdb_id']}"
        E.append({"id": r["emdb_id"], "method": r["structure_determination_method"], "stated": st,
                  "tail": co / st > 1.25, "study": key, "has_citation": not key.startswith("entry:")})
    R = {"n_entries": len(E), "n_tail": sum(x["tail"] for x in E)}

    def band(x):
        for lo, hi in BANDS:
            if lo <= x["stated"] < hi:
                return f"{lo}-{hi if hi < 1e9 else 'up'} A"

    cells = defaultdict(lambda: [0, 0])
    for x in E:
        c = cells[(band(x), x["method"])]
        c[0] += x["tail"]; c[1] += 1
    R["tail_by_band_and_method"] = {f"{b} | {m}": {"tail": k, "n": n, "pct": round(100 * k / n, 2), "ci95": wilson(k, n)}
                                    for (b, m), (k, n) in sorted(cells.items())}

    studies = defaultdict(list)
    for x in E:
        studies[x["study"]].append(x)
    R["n_studies"] = len(studies)
    R["n_entries_with_citation"] = sum(x["has_citation"] for x in E)
    tail_studies = {k: sum(x["tail"] for x in v) for k, v in studies.items() if any(x["tail"] for x in v)}
    R["n_studies_with_tail_entry"] = len(tail_studies)
    top = sorted(tail_studies.items(), key=lambda kv: -kv[1])[:10]
    R["largest_tail_contributors"] = [{"study": k, "tail_entries": n, "entries": len(studies[k])} for k, n in top]
    R["share_of_tail_from_top10_studies"] = round(100 * sum(n for _, n in top) / R["n_tail"], 2)
    keys = sorted(studies)
    k_arr = np.array([sum(x["tail"] for x in studies[k]) for k in keys])
    n_arr = np.array([len(studies[k]) for k in keys])
    rng = np.random.default_rng(20261001)
    pcts = []
    for _ in range(2000):
        idx = rng.integers(0, len(keys), len(keys))
        pcts.append(100 * k_arr[idx].sum() / n_arr[idx].sum())
    R["tail_pct_entry_level"] = round(100 * R["n_tail"] / R["n_entries"], 2)
    R["tail_pct_study_bootstrap_ci95"] = [round(float(np.percentile(pcts, 2.5)), 2), round(float(np.percentile(pcts, 97.5)), 2)]
    R["tail_pct_one_entry_per_study_mean"] = round(100 * float(np.mean(k_arr / n_arr)), 2)
    R["unrounded"] = {"tail_pct_study_bootstrap_ci95": [float(np.percentile(pcts, 2.5)), float(np.percentile(pcts, 97.5))],
                      "tail_pct_one_entry_per_study_mean": 100 * float(np.mean(k_arr / n_arr))}
    json.dump(R, open(os.path.join(ROOT, "data", "results", "robustness_v2.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in R.items() if k != "tail_by_band_and_method"}, indent=1))
    for k, v in R["tail_by_band_and_method"].items():
        print(f"  {k:40s} {v['tail']:4d}/{v['n']:6d} = {v['pct']:6.2f}%  {v['ci95']}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
