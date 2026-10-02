"""
A metadata-scale audit of resolution reporting in the EMDB archive, CPU-only.

From the full census (~46,900 maps) we quantify the resolution distribution over time and by method
(the "resolution revolution"), and the link between reported resolution and whether a fitted atomic
model was deposited. From a stratified per-entry sample we audit how the reported resolution is
substantiated (res_type, the FSC criterion cited) and whether that changed around the February 2022
half-map deposition mandate, and we run a Nyquist sanity check of reported resolution against the
map pixel spacing. Every number traces to the two cached CSVs.
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(__file__)
RAW = os.path.join(HERE, "..", "data", "raw")
RESULTS = os.path.join(HERE, "..", "data", "results")
os.makedirs(RESULTS, exist_ok=True)
MANDATE_YEAR = 2022


def method_group(m):
    if not isinstance(m, str):
        return "other"
    m = m.lower()
    if "single" in m: return "single particle"
    if "subtomogram" in m: return "subtomogram avg"
    if "tomography" in m: return "tomography"
    if "helical" in m: return "helical"
    if "crystal" in m: return "crystallography"
    return "other"


def main():
    c = pd.read_csv(os.path.join(RAW, "emdb_census.csv"), dtype={"deposition": str})
    c["year"] = pd.to_numeric(c["deposition"].str[:4], errors="coerce")
    c["mg"] = c["method"].map(method_group)
    out = {}

    # 1. Resolution revolution: median reported resolution by deposition year
    res = c.dropna(subset=["resolution", "year"])
    by_year = res.groupby("year")["resolution"].agg(["median", "count"])
    yr = {int(y): {"median_res": float(r["median"]), "n": int(r["count"])}
          for y, r in by_year.iterrows() if 2005 <= y <= 2025}
    frac_high = {int(y): float((res[res.year == y]["resolution"] < 4).mean())
                 for y in range(2010, 2026) if (res.year == y).sum() > 30}
    out["resolution_over_time"] = {"median_by_year": yr, "frac_sub4A_by_year": frac_high,
                                   "overall_median": float(res["resolution"].median()),
                                   "n_with_resolution": int(len(res))}

    # 2. By method. We record both the total entry count and the number that actually carry a
    # resolution, because the medians are computed over the with-resolution subset and that
    # denominator differs sharply by method (tomography lists many entries with no resolution).
    out["by_method"] = {g: {"n": int((c.mg == g).sum()),
                            "n_with_res": int(((c.mg == g) & c.resolution.notna()).sum()),
                            "median_res": float(res[res.mg == g]["resolution"].median())
                            if (res.mg == g).any() else None}
                        for g in ["single particle", "subtomogram avg", "tomography", "helical",
                                  "crystallography", "other"]}

    # 3. Fitted-model (PDB) rate vs resolution
    bins = [(0, 3), (3, 4), (4, 6), (6, 10), (10, 1e9)]
    out["pdb_rate_by_resolution"] = {
        f"{lo}-{hi if hi < 1e9 else 'inf'}A": {"n": int(((res.resolution >= lo) & (res.resolution < hi)).sum()),
                                               "pdb_rate": float(res[(res.resolution >= lo) & (res.resolution < hi)]["has_pdb"].mean())}
        for lo, hi in bins}

    # 4. Substantiation (per-entry sample)
    d = pd.read_csv(os.path.join(RAW, "emdb_detail.csv"), dtype={"deposition": str})
    d["year"] = pd.to_numeric(d["deposition"].str[:4], errors="coerce")
    d["mg"] = d["method"].map(method_group)
    fsc143 = d["resolution_method"].fillna("").str.contains("0.143")
    other = d["resolution_method"].fillna("").str.contains("OTHER")
    d["post"] = d["deposition"] >= "20220201"               # the Feb 2022 half-map mandate, matching the sample
    pre_f, post_f = fsc143[~d.post], fsc143[d.post]
    z, p = proportions_ztest(pre_f.sum(), len(pre_f), post_f.sum(), len(post_f))
    # Census-weighted FSC-0.143 estimate: the sample is 50/50 pre/post by construction, but the
    # census is not, so re-weight the sample's pre- and post-mandate rates by the census pre/post
    # split (entries with a deposition date) to get an archive-representative figure.
    cdate = c[c["deposition"].notna()]
    census_pre = int((cdate["deposition"] < "20220201").sum())
    census_post = int((cdate["deposition"] >= "20220201").sum())
    w_pre = census_pre / (census_pre + census_post)
    fsc0143_census_weighted = float(w_pre * pre_f.mean() + (1 - w_pre) * post_f.mean())
    # FSC-0.143 share by deposition-year bucket, to see whether the shift is abrupt or a gradual trend
    # FSC-0.143 adoption fraction per bucket, with a Wilson 95% score interval from the
    # per-bucket counts, so the reported percentages carry sampling uncertainty.
    by_bucket = {}
    for lo, hi in [(2002, 2013), (2014, 2017), (2018, 2021), (2022, 2026)]:
        mask = (d.year >= lo) & (d.year <= hi)
        if mask.sum() > 20:
            x, nb = int(fsc143[mask].sum()), int(mask.sum())
            lo_ci, hi_ci = wilson_ci(x, nb)
            by_bucket[f"{lo}-{hi}"] = {"frac_fsc0143": float(x / nb), "n": nb,
                                       "n_fsc0143": x,
                                       "wilson95_lo": float(lo_ci), "wilson95_hi": float(hi_ci)}
    # Stand-alone artifact for the FSC-0.143 adoption Wilson intervals used in Section 3.3.
    wilson_out = {
        "method": "Wilson score interval, 95% confidence",
        "source": "data/raw/emdb_detail.csv, resolution_method contains '0.143', by deposition-year bucket",
        "buckets": {b: {"n_fsc0143": v["n_fsc0143"], "n": v["n"],
                        "frac": v["frac_fsc0143"], "pct": round(v["frac_fsc0143"] * 100, 1),
                        "wilson95_lo": v["wilson95_lo"], "wilson95_hi": v["wilson95_hi"],
                        "wilson95_lo_pct": round(v["wilson95_lo"] * 100, 1),
                        "wilson95_hi_pct": round(v["wilson95_hi"] * 100, 1)}
                    for b, v in by_bucket.items()},
    }
    with open(os.path.join(RESULTS, "fsc0143_wilson_ci.json"), "w") as fh:
        json.dump(wilson_out, fh, indent=2)
    # New finding: characterize the residual post-mandate entries that still do NOT cite FSC-0.143.
    # The trend is known (Patwardhan 2017); the residual non-adopters are not. Break them down by the
    # criterion they cite instead, by method group, and by year, on the post-2022 half of the sample.
    post = d[d.post].copy()
    post_nonf = post[~fsc143[d.post].values]
    res_lo, res_hi = wilson_ci(len(post_nonf), len(post))
    residual = {
        "n_post": int(len(post)),
        "n_post_fsc0143": int(fsc143[d.post].sum()),
        "frac_post_fsc0143": float(fsc143[d.post].mean()),
        "n_post_non_fsc0143": int(len(post_nonf)),
        "frac_post_non_fsc0143": float(len(post_nonf) / len(post)),
        "frac_post_non_fsc0143_wilson95": [float(res_lo), float(res_hi)],
        "alt_criterion_counts": {str(k): int(v)
                                 for k, v in post_nonf["resolution_method"].fillna("nan").value_counts().items()},
        "non_fsc0143_by_method": {},
        "non_fsc0143_by_year": {},
    }
    for g in ["single particle", "helical", "subtomogram avg", "tomography", "crystallography"]:
        gm = post[post.mg == g]
        if len(gm):
            nf = int((~fsc143[d.post & (d.mg == g)]).sum())
            lo_ci, hi_ci = wilson_ci(nf, len(gm))
            residual["non_fsc0143_by_method"][g] = {"n": int(len(gm)), "n_non_fsc0143": nf,
                                                    "frac_non_fsc0143": float(nf / len(gm)),
                                                    "wilson95_lo": float(lo_ci), "wilson95_hi": float(hi_ci)}
    for y in sorted(post.year.dropna().unique()):
        ym = d.post & (d.year == y)
        nf_y = int((~fsc143[ym]).sum())
        lo_ci, hi_ci = wilson_ci(nf_y, int(ym.sum()))
        residual["non_fsc0143_by_year"][str(int(y))] = {"n": int(ym.sum()),
                                                         "n_non_fsc0143": nf_y,
                                                         "frac_non_fsc0143": float((~fsc143[ym]).mean()),
                                                         "wilson95_lo": float(lo_ci), "wilson95_hi": float(hi_ci)}
    out["substantiation"] = {
        "n": int(len(d)),
        "res_type_counts": d["res_type"].value_counts(dropna=False).to_dict(),
        "method_counts": {str(k): int(v) for k, v in d["resolution_method"].value_counts(dropna=False).head(6).items()},
        "frac_fsc0143": float(fsc143.mean()), "frac_other": float(other.mean()),
        "frac_fsc0143_census_weighted": fsc0143_census_weighted,
        "census_pre_post": {"pre": census_pre, "post": census_post},
        "sample_single_particle_frac": float((d.mg == "single particle").mean()),
        "census_single_particle_frac": float((c.mg == "single particle").mean()),
        "fsc0143_pre_mandate": float(pre_f.mean()), "fsc0143_post_mandate": float(post_f.mean()),
        "n_pre": int(len(pre_f)), "n_post": int(len(post_f)), "mandate_ztest_p": float(p),
        "fsc0143_by_year_bucket": by_bucket,
        "residual_non_fsc0143_post_mandate": residual,
    }

    # 5. Nyquist sanity: reported resolution vs 2 * pixel spacing
    n = d.dropna(subset=["resolution", "pixel_x"]).copy()
    n = n[n.pixel_x > 0]
    n["ratio"] = n["resolution"] / (2.0 * n["pixel_x"])      # < 1 is sub-Nyquist (impossible)
    sub = n[n["ratio"] < 1.0].copy()
    # Of the sub-Nyquist flags, some are not genuine over-claims but a recognizable metadata error:
    # the recorded resolution exactly equals the pixel spacing (the pixel/voxel size landed in the
    # resolution slot), which forces ratio ~ 0.5. Separate those from a substantive sub-Nyquist claim.
    sub["res_eq_pixel"] = (np.abs(sub["resolution"] - sub["pixel_x"]) / sub["pixel_x"]) < 0.001
    substantive = sub[~sub["res_eq_pixel"]]
    out["nyquist"] = {"n": int(len(n)), "median_ratio": float(n["ratio"].median()),
                      "frac_sub_nyquist": float((n["ratio"] < 1.0).mean()),
                      "n_sub_nyquist": int(len(sub)),
                      "n_sub_nyquist_pixel_artifact": int(sub["res_eq_pixel"].sum()),
                      "n_sub_nyquist_substantive": int(len(substantive)),
                      "frac_sub_nyquist_substantive": float(len(substantive) / len(n)),
                      "sub_nyquist_entries": [
                          {"emdb_id": str(r["emdb_id"]), "resolution": float(r["resolution"]),
                           "pixel_x": float(r["pixel_x"]), "ratio": float(r["ratio"]),
                           "method": str(r["method"]), "resolution_method": str(r["resolution_method"]),
                           "res_eq_pixel": bool(r["res_eq_pixel"])}
                          for _, r in sub.sort_values("ratio").iterrows()],
                      "frac_within_1p2x_nyquist": float((n["ratio"] < 1.2).mean()),
                      "median_pixel": float(n["pixel_x"].median())}

    with open(os.path.join(RESULTS, "emdb.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    rev = out["resolution_over_time"]
    print(f"Resolution revolution: median {rev['median_by_year'][2010]['median_res']:.1f}A (2010) -> "
          f"{rev['median_by_year'][2024]['median_res']:.1f}A (2024); sub-4A frac "
          f"{rev['frac_sub4A_by_year'][2012]*100:.0f}% (2012) -> {rev['frac_sub4A_by_year'][2024]*100:.0f}% (2024)")
    print("By method median res:", {g: round(v["median_res"], 1) for g, v in out["by_method"].items() if v["median_res"]})
    print("PDB rate by res:", {k: round(v["pdb_rate"], 2) for k, v in out["pdb_rate_by_resolution"].items()})
    s = out["substantiation"]
    print(f"Substantiation: res_type={s['res_type_counts']}; FSC0.143 {s['frac_fsc0143']*100:.0f}%, OTHER {s['frac_other']*100:.0f}%")
    print(f"  FSC0.143 pre {s['fsc0143_pre_mandate']*100:.0f}% vs post {s['fsc0143_post_mandate']*100:.0f}% mandate (p={s['mandate_ztest_p']:.3f})")
    r = s["residual_non_fsc0143_post_mandate"]
    print(f"  Residual post-mandate NON-FSC0.143: {r['n_post_non_fsc0143']}/{r['n_post']} = {r['frac_post_non_fsc0143']*100:.1f}%")
    print(f"    cite instead: {r['alt_criterion_counts']}")
    print(f"    by method: " + ", ".join(f"{g} {v['frac_non_fsc0143']*100:.0f}% ({v['n_non_fsc0143']}/{v['n']})"
                                          for g, v in r['non_fsc0143_by_method'].items()))
    print(f"    by year: " + ", ".join(f"{y} {v['frac_non_fsc0143']*100:.0f}%"
                                       for y, v in r['non_fsc0143_by_year'].items()))
    ny = out["nyquist"]
    print(f"Nyquist: median ratio res/(2px)={ny['median_ratio']:.2f}, sub-Nyquist {ny['frac_sub_nyquist']*100:.1f}% "
          f"({ny['n_sub_nyquist']}/{ny['n']}); of those {ny['n_sub_nyquist_pixel_artifact']} are res==pixel artifacts, "
          f"{ny['n_sub_nyquist_substantive']} substantive (={ny['frac_sub_nyquist_substantive']*100:.2f}%); "
          f"within 1.2x {ny['frac_within_1p2x_nyquist']*100:.1f}%")
    for e in ny["sub_nyquist_entries"]:
        print(f"    {e['emdb_id']}: res={e['resolution']} px={e['pixel_x']} ratio={e['ratio']:.3f} "
              f"{'[res==px artifact]' if e['res_eq_pixel'] else '[substantive]'} {e['method']} {e['resolution_method']}")
    print("wrote emdb.json")


def proportions_ztest(x1, n1, x2, n2):
    p1, p2 = x1 / n1, x2 / n2
    p = (x1 + x2) / (n1 + n2)
    se = np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    z = (p1 - p2) / se if se > 0 else 0.0
    return z, 2 * (1 - stats.norm.cdf(abs(z)))


def wilson_ci(x, n, conf=0.95):
    """Wilson score interval for a binomial proportion (lo, hi)."""
    if n == 0:
        return (float("nan"), float("nan"))
    z = stats.norm.ppf(1 - (1 - conf) / 2.0)
    phat = x / n
    denom = 1 + z * z / n
    centre = phat + z * z / (2 * n)
    half = z * np.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n))
    return ((centre - half) / denom, (centre + half) / denom)


if __name__ == "__main__":
    main()
