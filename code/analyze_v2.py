"""Compare the resolution that depositors state for EMDB entries with (a) the 0.143 crossing of the FSC curve they
deposited and (b) the FSC that EMDB's validation analysis (VA) recalculates from the deposited half-maps, unmasked,
masked with an EMDB-generated mask, and with phase-randomization correction.

Input: data/raw/emdb_index_<date>.csv from code/fetch_emdb_index.py (one row per EMDB entry, CC0).
Output: data/results/analysis_v2.json (every number used in the text), data/results/tail_entries.csv.

Definitions used throughout
  stated     the entry's resolution field, used only when its criterion is "FSC 0.143 CUT-OFF"
  own curve  the 0.143 crossing of the depositor's FSC curve, as extracted by EMDB (author_resolution_fsc_0143_value)
  unmasked   RELION FSC of the deposited half-maps, no mask (calculated_resolution_fsc_0143_value)
  masked     the same with EMDB's automatically generated soft mask (calculated_resolution_fsc_masked_0143_value)
  corrected  masked FSC corrected by phase randomization (calculated_resolution_fsc_corrected_0143_value)
  ratio      recalculated / stated, in angstrom: above 1 means the recalculation gives a coarser (worse) resolution
Tomograms are excluded: a single tomogram has no half-maps, so no FSC can be recalculated for it.
Deterministic: the bootstrap uses a fixed seed.

    python3 code/analyze_v2.py data/raw/emdb_index_2026-10-01.csv
"""
import csv, json, math, os, sys
from collections import Counter

import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUT = os.path.join(ROOT, "data", "results")
METHODS = ("singleParticle", "helical", "subtomogramAveraging")
POLICY = "2022-02-25"        # half-maps mandatory for single-particle, helical and subtomogram-averaging maps
THRESH = (1.10, 1.25, 1.50, 2.00)
BANDS = ((0, 3), (3, 4), (4, 6), (6, 10), (10, 1e9))


def num(x):
    try:
        v = float(str(x).split(",")[0])
        return v if math.isfinite(v) and v > 0 else None
    except (TypeError, ValueError):
        return None


def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(100 * (c - h), 2), round(100 * (c + h), 2)]


def boot_median(x, seed=20261001, b=2000):
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    meds = np.median(x[rng.integers(0, len(x), size=(b, len(x)))], axis=1)
    return [round(float(np.percentile(meds, 2.5)), 4), round(float(np.percentile(meds, 97.5)), 4)]


def summarize(ratio):
    ratio = np.asarray(ratio)
    n = len(ratio)
    out = {"n": n, "median": round(float(np.median(ratio)), 4), "median_ci95": boot_median(ratio),
           "p5_p25_p75_p95": [round(float(v), 4) for v in np.percentile(ratio, [5, 25, 75, 95])]}
    for t in THRESH:
        k_hi, k_lo = int((ratio > t).sum()), int((ratio < 1 / t).sum())
        tag = f"{int(round((t - 1) * 100))}pct"
        out[f"coarser_by_more_than_{tag}"] = {"k": k_hi, "pct": round(100 * k_hi / n, 2), "ci95": wilson(k_hi, n)}
        out[f"finer_by_more_than_{tag}"] = {"k": k_lo, "pct": round(100 * k_lo / n, 2), "ci95": wilson(k_lo, n)}
    return out


def software(r):
    s = (r["final_angle_assignment_software"] or r["image_reconstruction_software"] or r["software"]).upper()
    if not s:
        return "not stated"
    first = s.split(",")[0].strip()
    return "RELION" if "RELION" in first else "cryoSPARC" if "CRYOSPARC" in first else "other"


def main(index_csv):
    rows = list(csv.DictReader(open(index_csv, newline="", encoding="utf-8")))
    R = {"input": os.path.basename(index_csv), "index_rows": len(rows)}
    released = [r for r in rows if r["header_release_date"] and not r["obsolete_date"] and not r["withdrawn_date"]]
    R["released_current"] = len(released)
    R["not_released_or_withdrawn_or_obsolete"] = len(rows) - len(released)
    R["released_by_method"] = dict(Counter(r["structure_determination_method"] for r in released).most_common())
    R["released_by_year"] = dict(sorted(Counter(r["header_release_date"][:4] for r in released).items()))
    R["released_on_or_before_2026-06-30"] = sum(r["header_release_date"][:10] <= "2026-06-30" for r in released)

    pool = [r for r in released if r["structure_determination_method"] in METHODS]
    R["pool_spa_helical_sta"] = len(pool)
    R["pool_with_halfmaps_listed"] = sum(bool(r["half_map_filename"]) for r in pool)
    R["pool_criterion"] = dict(Counter(r["resolution_method"] or "(none)" for r in pool).most_common())

    # entries compared at the 0.143 threshold
    E, nyq_bad = [], []
    for r in pool:
        if r["resolution_method"] != "FSC 0.143 CUT-OFF":
            continue
        st = num(r["resolution"])
        un, ma, co = (num(r[k]) for k in ("calculated_resolution_fsc_0143_value", "calculated_resolution_fsc_masked_0143_value",
                                          "calculated_resolution_fsc_corrected_0143_value"))
        if st is None or un is None:
            continue
        px = num(r["half_map_pixel_spacing_x"]) or num(r["map_pixel_spacing_x"])
        nyq = 2 * px if px else None
        rec = {"id": r["emdb_id"], "method": r["structure_determination_method"], "year": r["header_release_date"][:4],
               "post_policy": r["deposition_date"][:10] >= POLICY, "stated": st, "own": num(r["author_resolution_fsc_0143_value"]),
               "unmasked": un, "masked": ma, "corrected": co, "nyquist": nyq, "software": software(r),
               "point_group": (r["image_reconstruction_point_group"] or "").split(",")[0].strip()}
        bad = nyq is not None and any(v is not None and v < 0.999 * nyq for v in (un, ma, co))
        if bad:
            nyq_bad.append(rec)
            continue
        E.append(rec)
    R["compared_0143"] = {"n": len(E), "excluded_recalculation_finer_than_halfmap_nyquist": len(nyq_bad),
                          "excluded_ids": sorted(x["id"] for x in nyq_bad)[:50],
                          "stated_finer_than_map_nyquist": sum(1 for x in E if x["nyquist"] and x["stated"] < 0.999 * x["nyquist"]),
                          "by_method": dict(Counter(x["method"] for x in E)), "by_year": dict(sorted(Counter(x["year"] for x in E).items())),
                          "deposited_after_policy": sum(x["post_policy"] for x in E)}

    for key in ("unmasked", "masked", "corrected"):
        R[f"stated_vs_{key}"] = summarize([x[key] / x["stated"] for x in E if x[key]])
    own = [x for x in E if x["own"]]
    R["stated_vs_own_curve"] = summarize([x["own"] / x["stated"] for x in own])
    R["stated_vs_own_curve"]["within_1pct"] = round(100 * np.mean([abs(x["own"] / x["stated"] - 1) <= 0.01 for x in own]), 2)
    oc = [x for x in own if x["corrected"]]
    R["own_curve_vs_corrected"] = summarize([x["corrected"] / x["own"] for x in oc])
    um = [x for x in E if x["masked"]]
    R["unmasked_over_masked"] = summarize([x["unmasked"] / x["masked"] for x in um])

    # the tail: stated resolution finer than the corrected recalculation by more than 25%
    T = [x for x in E if x["corrected"] and x["corrected"] / x["stated"] > 1.25]
    base = [x for x in E if x["corrected"]]
    R["tail_definition"] = "corrected / stated > 1.25"
    R["tail_n"], R["tail_base_n"] = len(T), len(base)

    def rate(group_fn):
        g = {}
        for x in base:
            g.setdefault(group_fn(x), [0, 0])[1] += 1
        for x in T:
            g[group_fn(x)][0] += 1
        return {k: {"tail": v[0], "n": v[1], "pct": round(100 * v[0] / v[1], 2), "ci95": wilson(v[0], v[1])}
                for k, v in sorted(g.items(), key=lambda kv: str(kv[0]))}

    def band(x):
        for lo, hi in BANDS:
            if lo <= x["stated"] < hi:
                return f"{lo}-{hi if hi < 1e9 else 'up'} A"

    R["tail_by_band"] = rate(band)
    R["tail_by_method"] = rate(lambda x: x["method"])
    R["tail_by_year"] = rate(lambda x: x["year"])
    R["tail_by_software"] = rate(lambda x: x["software"])
    R["tail_by_symmetry"] = rate(lambda x: "C1" if x["point_group"] == "C1" else "higher or helical" if x["point_group"] else "not stated")
    R["tail_with_own_curve"] = {"n": sum(1 for x in T if x["own"]),
                                "own_curve_agrees_with_stated_within_5pct": sum(1 for x in T if x["own"] and abs(x["own"] / x["stated"] - 1) <= 0.05),
                                "own_curve_agrees_with_corrected_within_10pct": sum(1 for x in T if x["own"] and abs(x["corrected"] / x["own"] - 1) <= 0.10)}
    # the opposite tail: stated resolution coarser than the corrected recalculation by more than 25%
    R["opposite_tail_n"] = sum(1 for x in base if x["corrected"] / x["stated"] < 1 / 1.25)

    os.makedirs(OUT, exist_ok=True)
    json.dump(R, open(os.path.join(OUT, "analysis_v2.json"), "w"), indent=1)
    with open(os.path.join(OUT, "tail_entries.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["emdb_id", "method", "release_year", "stated_A", "own_curve_A", "unmasked_A", "masked_A", "corrected_A",
                    "corrected_over_stated", "halfmap_nyquist_A", "software", "point_group"])
        for x in sorted(T, key=lambda x: -x["corrected"] / x["stated"]):
            w.writerow([x["id"], x["method"], x["year"], x["stated"], x["own"] or "", x["unmasked"], x["masked"] or "", x["corrected"],
                        round(x["corrected"] / x["stated"], 4), x["nyquist"] or "", x["software"], x["point_group"]])
    show = {k: R[k] for k in ("index_rows", "released_current", "released_on_or_before_2026-06-30", "pool_spa_helical_sta", "tail_n", "tail_base_n", "opposite_tail_n")}
    print(json.dumps(show, indent=1))
    print("compared:", json.dumps(R["compared_0143"], default=str)[:600])
    for k in ("stated_vs_unmasked", "stated_vs_masked", "stated_vs_corrected", "stated_vs_own_curve", "own_curve_vs_corrected", "unmasked_over_masked"):
        s = R[k]
        print(f"{k:24s} n={s['n']:6d} median {s['median']} {s['median_ci95']} p5..p95 {s['p5_p25_p75_p95']} "
              f">25% coarser {s['coarser_by_more_than_25pct']['pct']}% >25% finer {s['finer_by_more_than_25pct']['pct']}%")
    for k in ("tail_by_band", "tail_by_method", "tail_by_year", "tail_by_software", "tail_by_symmetry", "tail_with_own_curve"):
        print(k, json.dumps(R[k]))


if __name__ == "__main__":
    main(sys.argv[1])
