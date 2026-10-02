"""Analyse the random sample of entries whose recalculated FSC values are missing from the EMDB search index
(data/results/unindexed_sample.jsonl, written by code/sample_unindexed_va.py), with the definitions of analyze_v2.py:
ratio = recalculated / stated; an entry is compared when its validation record has an unmasked 0.143 crossing; entries
whose recalculation is finer than 0.999 x the half-map Nyquist limit are excluded; the tail is corrected / stated > 1.25.

For each sampled entry the last line with HTTP status 200 and a parsed record is used (the download was resumed after
time-outs, so the file can hold earlier failed lines for the same entry). The sample is compared with the indexed
entries of analysis_v2.json, and a stratified estimate of the tail rate over indexed and unindexed entries together is
given with a bootstrap interval that resamples the sampled entries within each stratum (the indexed entries are a
census and are held fixed).

Output: data/results/unindexed_sample_analysis.json

    python3 code/analyze_unindexed_sample.py
"""
import json, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_v2 import BANDS, wilson  # noqa: E402

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
RES = os.path.join(ROOT, "data", "results")
STRATA = ("pre2024", "from2024")


def load():
    last, n_lines = {}, 0
    for line in open(os.path.join(RES, "unindexed_sample.jsonl")):
        if not line.strip():
            continue
        n_lines += 1
        r = json.loads(line)
        if r.get("http_status") == 200 and "has_relion_fsc" in r:
            last[r["emdb_id"]] = r
        else:
            last.setdefault(r["emdb_id"], r)        # kept only if no successful line exists for this entry
    return list(last.values()), n_lines


def band(stated):
    for lo, hi in BANDS:
        if lo <= stated < hi:
            return f"{lo}-{hi if hi < 1e9 else 'up'} A"


def classify(r):
    """Return (status, record) with status one of: failed, no_relion_fsc, no_unmasked, nyquist_excluded, compared."""
    if r.get("http_status") != 200 or "has_relion_fsc" not in r:
        return "failed"
    if not r["has_relion_fsc"]:
        return "no_relion_fsc"
    if not r.get("unmasked"):
        return "no_unmasked"
    px = r.get("halfmap_pixel") or r.get("map_pixel")
    nyq = 2 * px if px else None
    if nyq is not None and any(v is not None and v < 0.999 * nyq for v in (r["unmasked"], r.get("masked"), r.get("corrected"))):
        return "nyquist_excluded"
    return "compared"


def stratum_summary(recs):
    status = {}
    for r in recs:
        status.setdefault(classify(r), []).append(r)
    E = status.get("compared", [])
    base = [r for r in E if r.get("corrected")]
    tail = [r for r in base if r["corrected"] / r["stated"] > 1.25]
    out = {"sampled": len(recs), "status": {k: len(v) for k, v in sorted(status.items())},
           "compared": len(E), "with_corrected": len(base), "tail": len(tail),
           "tail_pct": round(100 * len(tail) / len(base), 2) if base else None, "tail_ci95": wilson(len(tail), len(base))}
    if base:
        ratio = np.array([r["corrected"] / r["stated"] for r in base])
        out["corrected_over_stated"] = {"median": round(float(np.median(ratio)), 4),
                                        "p25_p75": [round(float(v), 4) for v in np.percentile(ratio, [25, 75])],
                                        "within_10pct": round(100 * float(np.mean(np.abs(ratio - 1) <= 0.10)), 2)}
        out["unrounded"] = {"median": float(np.median(ratio)), "p25_p75": [float(v) for v in np.percentile(ratio, [25, 75])],
                            "within_10pct": 100 * float(np.mean(np.abs(ratio - 1) <= 0.10))}
    un = np.array([r["unmasked"] / r["stated"] for r in E])
    if len(un):
        out["unmasked_over_stated"] = {"median": round(float(np.median(un)), 4),
                                       "coarser_by_more_than_25pct": round(100 * float(np.mean(un > 1.25)), 2)}
        out.setdefault("unrounded", {}).update({"unmasked_median": float(np.median(un)),
                                                "unmasked_coarser_by_more_than_25pct": 100 * float(np.mean(un > 1.25))})
    below4 = [r for r in base if r["stated"] < 4]
    from4 = [r for r in base if r["stated"] >= 4]
    for name, g in (("stated_below_4A", below4), ("stated_4A_or_coarser", from4)):
        k = sum(r["corrected"] / r["stated"] > 1.25 for r in g)
        out[name] = {"tail": k, "n": len(g), "pct": round(100 * k / len(g), 2) if g else None, "ci95": wilson(k, len(g))}
    bands = {}
    for r in base:
        b = bands.setdefault(band(r["stated"]), [0, 0])
        b[1] += 1
        b[0] += r["corrected"] / r["stated"] > 1.25
    out["tail_by_band"] = {k: {"tail": v[0], "n": v[1]} for k, v in sorted(bands.items())}
    own = [r for r in E if r.get("own")]
    out["own_curve"] = {"n": len(own),
                        "agrees_with_stated_within_1pct": sum(abs(r["own"] / r["stated"] - 1) <= 0.01 for r in own)}
    tail_own = [r for r in tail if r.get("own")]
    out["tail_own_curve"] = {"n": len(tail_own),
                             "agrees_with_stated_within_5pct": sum(abs(r["own"] / r["stated"] - 1) <= 0.05 for r in tail_own)}
    out["tail_entries"] = sorted(({"id": r["emdb_id"], "release_date": r["release_date"], "method": r["method"],
                                   "stated": r["stated"], "unmasked": round(r["unmasked"], 3), "corrected": round(r["corrected"], 3),
                                   "own": round(r["own"], 3) if r.get("own") else None} for r in tail), key=lambda d: d["id"])
    return out


def main():
    recs, n_lines = load()
    A = json.load(open(os.path.join(RES, "analysis_v2.json")))
    S = {"input_lines": n_lines, "entries": len(recs), "definition": "tail: corrected / stated > 1.25 (as analysis_v2.json)"}
    by = {s: [r for r in recs if r["stratum"] == s] for s in STRATA}
    for s in STRATA:
        pools = {r["stratum_pool_n"] for r in by[s]}
        assert len(pools) <= 1, (s, pools)
        S[s] = {"pool_n": pools.pop() if pools else None, **stratum_summary(by[s])}

    # stratified estimate over indexed (census) and unindexed (sampled) entries; ratio estimator within each stratum
    ix_tail, ix_base = A["tail_n"], A["tail_base_n"]
    rng = np.random.default_rng(20261001)

    def estimate(samples):
        t, b = float(ix_tail), float(ix_base)
        for s, g in samples.items():
            n = len(g)
            if n == 0:
                return None
            ok = [r for r in g if classify(r) == "compared" and r.get("corrected")]
            t += S[s]["pool_n"] * sum(r["corrected"] / r["stated"] > 1.25 for r in ok) / n
            b += S[s]["pool_n"] * len(ok) / n
        return 100 * t / b, b

    usable = {s: [r for r in by[s] if classify(r) != "failed"] for s in STRATA}
    point = estimate(usable)
    if point:
        boots = []
        for _ in range(2000):
            boots.append(estimate({s: [g[i] for i in rng.integers(0, len(g), len(g))] for s, g in usable.items()})[0])
        S["combined"] = {"indexed_tail": ix_tail, "indexed_base": ix_base, "indexed_pct": round(100 * ix_tail / ix_base, 2),
                         "estimated_base": round(point[1]), "pct": round(point[0], 2),
                         "bootstrap_ci95": [round(float(np.percentile(boots, 2.5)), 2), round(float(np.percentile(boots, 97.5)), 2)],
                         "unrounded": {"pct": point[0], "base": point[1],
                                       "bootstrap_ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]},
                         "note": "failed downloads are left out of each stratum's denominator"}
    out = os.path.join(RES, "unindexed_sample_analysis.json")
    json.dump(S, open(out, "w"), indent=1)
    for s in STRATA:
        d = S[s]
        print(s, {k: d[k] for k in ("pool_n", "sampled", "status", "compared", "with_corrected", "tail", "tail_pct", "tail_ci95")})
        print("   ", {k: d.get(k) for k in ("corrected_over_stated", "unmasked_over_stated", "stated_below_4A", "stated_4A_or_coarser", "own_curve", "tail_own_curve")})
    print("combined", S.get("combined"))


if __name__ == "__main__":
    main()
