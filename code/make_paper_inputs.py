"""Write every number, table and figure used in the manuscript from the saved results and the index table.

Inputs  data/raw/emdb_index_<date>.csv, data/raw/fetch_log.json, data/results/analysis_v2.json,
        data/results/robustness_v2.json, data/results/index_vs_va_check.json
Outputs paper/numbers.tex (one \\newcommand per number), paper/tab_*.tex, paper/fig*.pdf,
        data/results/paper_extra.json (counts computed here and not stored elsewhere)

    python3 code/make_paper_inputs.py data/raw/emdb_index_2026-10-01.csv
"""
import csv, json, os, sys
from collections import Counter, defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_v2 import METHODS, num, wilson  # noqa: E402
import analyze_unindexed_sample as UQ  # noqa: E402

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
RES = os.path.join(ROOT, "data", "results")
PAPER = os.path.join(ROOT, "paper")
METHOD_NAME = {"singleParticle": "single particle", "helical": "helical", "subtomogramAveraging": "subtomogram averaging"}

plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.linewidth": 0.6, "pdf.fonttype": 42})


REQUIRE_COMPLETE_SAMPLE = True     # the paper is built only from the complete sample (500 + 150 entries)


def fmt_int(n):
    return f"{n:,}".replace(",", "{,}")


def fmt(x, d):
    return f"{x:.{d}f}"


def wil(k, n, z=1.96):
    """Wilson 95% interval in percent, unrounded (analyze_v2.wilson rounds to 2 dp)."""
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * (p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5 / d
    return 100 * (c - h), 100 * (c + h)


def raw_summary(ratio, seed=20261001, b=2000):
    """Unrounded version of analyze_v2.summarize (same bootstrap draws)."""
    x = np.asarray(ratio)
    rng = np.random.default_rng(seed)
    meds = np.median(x[rng.integers(0, len(x), size=(b, len(x)))], axis=1)
    out = {"n": len(x), "median": float(np.median(x)), "median_ci95": [float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))],
           "q": [float(v) for v in np.percentile(x, [5, 25, 75, 95])]}
    for t in (1.10, 1.25, 1.50, 2.00):
        tag = f"{int(round((t - 1) * 100))}pct"
        out["hi" + tag] = int((x > t).sum())
        out["lo" + tag] = int((x < 1 / t).sum())
    return out


def load_entries(index_csv):
    """Same selection as analyze_v2.py, keeping the counts at each step."""
    rows = list(csv.DictReader(open(index_csv, newline="", encoding="utf-8")))
    flow = {"index_rows": len(rows)}
    rel = [r for r in rows if r["header_release_date"] and not r["obsolete_date"] and not r["withdrawn_date"]]
    flow["released_current"] = len(rel)
    flow["released_by_method"] = dict(Counter(r["structure_determination_method"] for r in rel))
    pool = [r for r in rel if r["structure_determination_method"] in METHODS]
    flow["pool"] = len(pool)
    c143 = [r for r in pool if r["resolution_method"] == "FSC 0.143 CUT-OFF"]
    flow["crit_0143"] = len(c143)
    flow["crit_0143_halfmap_listed"] = sum(bool(r["half_map_filename"]) for r in c143)
    withcalc = [r for r in c143 if num(r["resolution"]) and num(r["calculated_resolution_fsc_0143_value"])]
    flow["crit_0143_with_indexed_recalculation"] = len(withcalc)
    flow["with_indexed_recalculation_but_no_halfmap_listed"] = sum(not r["half_map_filename"] for r in withcalc)
    E, nyq = [], 0
    for r in withcalc:
        st = num(r["resolution"])
        un, ma, co = (num(r[k]) for k in ("calculated_resolution_fsc_0143_value", "calculated_resolution_fsc_masked_0143_value",
                                          "calculated_resolution_fsc_corrected_0143_value"))
        hp, mp = num(r["half_map_pixel_spacing_x"]), num(r["map_pixel_spacing_x"])
        px = hp or mp
        if px and any(v is not None and v < 0.999 * 2 * px for v in (un, ma, co)):
            nyq += 1
            continue
        E.append({"id": r["emdb_id"], "method": r["structure_determination_method"], "stated": st, "unmasked": un,
                  "masked": ma, "corrected": co, "own": num(r["author_resolution_fsc_0143_value"]), "hp": hp, "mp": mp,
                  "year": r["header_release_date"][:4]})
    flow["excluded_nyquist"] = nyq
    flow["compared"] = len(E)
    flow["compared_released_2024_on"] = sum(x["year"] >= "2024" for x in E)

    # coverage by release year: pool, half-map listed, indexed recalculation (any criterion)
    cov = defaultdict(lambda: [0, 0, 0])
    for r in pool:
        y = r["header_release_date"][:4]
        y = "2002-2015" if y <= "2015" else y
        cov[y][0] += 1
        cov[y][1] += bool(r["half_map_filename"])
        cov[y][2] += bool(r["half_map_filename"]) and num(r["calculated_resolution_fsc_0143_value"]) is not None
    flow["coverage_by_year"] = {y: {"pool": v[0], "halfmap": v[1], "halfmap_and_indexed": v[2]} for y, v in sorted(cov.items())}
    hm = sum(v[1] for v in cov.values())
    hi = sum(v[2] for v in cov.values())
    hm23 = sum(v[1] for y, v in cov.items() if y <= "2023")
    hi23 = sum(v[2] for y, v in cov.items() if y <= "2023")
    hm24 = sum(v[1] for y, v in cov.items() if y >= "2024")
    hi24 = sum(v[2] for y, v in cov.items() if y >= "2024")
    flow.update({"pool_halfmap": hm, "pool_halfmap_indexed": hi, "halfmap_to2023": hm23, "indexed_to2023": hi23,
                 "halfmap_2024on": hm24, "indexed_2024on": hi24})
    return E, flow


def tail_size_groups(E):
    g = Counter()
    pix = Counter()
    pixn = Counter()
    for x in E:
        if not x["corrected"]:
            continue
        q = x["corrected"] / x["stated"]
        k = "le125" if q <= 1.25 else "125to150" if q <= 1.5 else "150to200" if q <= 2 else "gt200"
        g[k] += 1
        if x["hp"] and x["mp"]:
            pixn[k] += 1
            pix[k] += abs(x["hp"] / x["mp"] - 1) > 0.01
    return {k: {"n": g[k], "pixel_known": pixn[k], "halfmap_pixel_differs_from_map": pix[k]} for k in g}


def figures(E):
    os.makedirs(PAPER, exist_ok=True)
    # Figure 1: share of entries whose 0.143 crossing is more than x times the stated resolution
    fig, ax = plt.subplots(figsize=(3.4, 2.9))
    xs = np.exp(np.linspace(np.log(0.7), np.log(4.0), 400))
    series = [("unmasked", "EMDB, unmasked", "#999999", "-"), ("masked", "EMDB, EMDB mask", "#1f77b4", "-"),
              ("corrected", "EMDB, corrected", "#d62728", "-"), ("own", "depositor's FSC curve", "#000000", "--")]
    for key, lab, col, ls in series:
        v = np.sort([x[key] / x["stated"] for x in E if x[key]])
        share = 1 - np.searchsorted(v, xs, side="right") / len(v)
        share = np.where(share > 0, share, np.nan)
        ax.plot(xs, share, color=col, ls=ls, lw=0.9, label=f"{lab} (n = {len(v):,})")
    ax.axvline(1.25, color="#555555", lw=0.5, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ticks = [0.75, 1, 1.25, 1.5, 2, 3, 4]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.minorticks_off()
    ax.set_xlim(0.7, 4.0)
    ax.set_ylim(1e-4, 1.3)
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.set_xlabel("x = 0.143 crossing / stated resolution")
    ax.set_ylabel("share of entries above x")
    ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.45, -0.24), ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(PAPER, "fig1_ratio_distributions.pdf"), bbox_inches="tight", pad_inches=0.03, metadata={"CreationDate": None})
    plt.close(fig)

    # Figure 2: stated versus corrected recalculation
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    s = np.array([x["stated"] for x in E if x["corrected"]])
    c = np.array([x["corrected"] for x in E if x["corrected"]])
    hb = ax.hexbin(s, c, xscale="log", yscale="log", gridsize=70, bins="log", cmap="Greys", mincnt=1,
                   extent=(np.log10(1.0), np.log10(60), np.log10(1.0), np.log10(200)), linewidths=0.1)
    xx = np.array([1.0, 60])
    ax.plot(xx, xx, color="#1f77b4", lw=0.7, label="equal")
    ax.plot(xx, 1.25 * xx, color="#d62728", lw=0.7, ls="--", label="25% coarser than stated")
    ax.set_xlim(1.0, 60)
    ax.set_ylim(1.0, 200)
    ax.set_xlabel("stated resolution (Å)")
    ax.set_ylabel("EMDB corrected FSC 0.143 crossing (Å)")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    cb = fig.colorbar(hb, ax=ax, pad=0.02)
    cb.set_label("entries per cell", fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PAPER, "fig2_stated_vs_corrected.pdf"), metadata={"CreationDate": None})
    plt.close(fig)


def fig3(rob):
    bands = ["0-3 A", "3-4 A", "4-6 A", "6-10 A", "10-up A"]
    labels = ["< 3", "3 to 4", "4 to 6", "6 to 10", "≥ 10"]
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for j, (m, col, mk) in enumerate([("singleParticle", "#1f77b4", "o"), ("helical", "#2ca02c", "s"),
                                      ("subtomogramAveraging", "#d62728", "^")]):
        xs, ys, lo, hi = [], [], [], []
        for i, b in enumerate(bands):
            cell = rob["tail_by_band_and_method"].get(f"{b} | {m}")
            if not cell or cell["n"] < 20:
                continue
            xs.append(i + (j - 1) * 0.18)
            ys.append(cell["pct"])
            lo.append(cell["pct"] - cell["ci95"][0])
            hi.append(cell["ci95"][1] - cell["pct"])
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt=mk, ms=3, color=col, lw=0.7, capsize=1.5, label=METHOD_NAME[m])
    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels(labels)
    ax.set_xlabel("stated resolution (Å)")
    ax.set_ylabel("entries with\ncorrected / stated > 1.25 (%)")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(PAPER, "fig3_tail_by_band_and_method.pdf"), metadata={"CreationDate": None})
    plt.close(fig)


def main(index_csv):
    A = json.load(open(os.path.join(RES, "analysis_v2.json")))
    B = json.load(open(os.path.join(RES, "robustness_v2.json")))
    C = json.load(open(os.path.join(RES, "index_vs_va_check.json")))
    L = json.load(open(os.path.join(ROOT, "data", "raw", "fetch_log.json")))
    E, flow = load_entries(index_csv)
    assert flow["compared"] == A["compared_0143"]["n"], (flow["compared"], A["compared_0143"]["n"])
    sizes = tail_size_groups(E)
    assert sum(v["n"] for k, v in sizes.items() if k != "le125") == A["tail_n"]
    extra = {"flow": flow, "tail_size_groups": sizes}
    json.dump(extra, open(os.path.join(RES, "paper_extra.json"), "w"), indent=1)

    N = {}

    def put(name, value):
        assert name not in N, name
        N[name] = value

    put("indexRows", fmt_int(flow["index_rows"]))
    put("releasedCurrent", fmt_int(flow["released_current"]))
    put("nTomo", fmt_int(flow["released_by_method"]["tomography"]))
    put("nEC", fmt_int(flow["released_by_method"]["electronCrystallography"]))
    put("poolN", fmt_int(flow["pool"]))
    put("critN", fmt_int(flow["crit_0143"]))
    put("critPct", fmt(100 * flow["crit_0143"] / flow["pool"], 1))
    put("critHalfmap", fmt_int(flow["crit_0143_halfmap_listed"]))
    put("critIndexed", fmt_int(flow["crit_0143_with_indexed_recalculation"]))
    put("nyqExcl", fmt_int(flow["excluded_nyquist"]))
    put("comparedN", fmt_int(flow["compared"]))
    put("comparedRecent", fmt_int(flow["compared_released_2024_on"]))
    put("comparedRecentPct", fmt(100 * flow["compared_released_2024_on"] / flow["compared"], 1))
    put("poolHalfmap", fmt_int(flow["pool_halfmap"]))
    put("poolHalfmapIndexed", fmt_int(flow["pool_halfmap_indexed"]))
    put("hmToTwentyThree", fmt_int(flow["halfmap_to2023"]))
    put("ixToTwentyThree", fmt_int(flow["indexed_to2023"]))
    put("ixToTwentyThreePct", fmt(100 * flow["indexed_to2023"] / flow["halfmap_to2023"], 1))
    put("hmFromTwentyFour", fmt_int(flow["halfmap_2024on"]))
    put("ixFromTwentyFour", fmt_int(flow["indexed_2024on"]))
    put("ixFromTwentyFourPct", fmt(100 * flow["indexed_2024on"] / flow["halfmap_2024on"], 1))
    cy = flow["coverage_by_year"]
    put("ixTwentySixPct", fmt(100 * cy["2026"]["halfmap_and_indexed"] / cy["2026"]["halfmap"], 1))
    put("ixTwentyFivePct", fmt(100 * cy["2025"]["halfmap_and_indexed"] / cy["2025"]["halfmap"], 1))
    put("bySPA", fmt_int(A["compared_0143"]["by_method"]["singleParticle"]))
    put("bySTA", fmt_int(A["compared_0143"]["by_method"]["subtomogramAveraging"]))
    put("byHel", fmt_int(A["compared_0143"]["by_method"]["helical"]))
    put("snapshotUTC", L["retrieved_utc"].replace("T", " ").replace("+00:00", " UTC"))
    put("indexSha", L["sha256"][:16])
    put("indexMB", fmt(L["bytes"] / 1e6, 1))

    RAW = {"stated_vs_unmasked": raw_summary([x["unmasked"] / x["stated"] for x in E]),
           "stated_vs_masked": raw_summary([x["masked"] / x["stated"] for x in E if x["masked"]]),
           "stated_vs_corrected": raw_summary([x["corrected"] / x["stated"] for x in E if x["corrected"]]),
           "stated_vs_own_curve": raw_summary([x["own"] / x["stated"] for x in E if x["own"]]),
           "own_curve_vs_corrected": raw_summary([x["corrected"] / x["own"] for x in E if x["own"] and x["corrected"]]),
           "unmasked_over_masked": raw_summary([x["unmasked"] / x["masked"] for x in E if x["masked"]])}
    for key, r in RAW.items():          # the unrounded values must reproduce the saved, rounded ones
        s = A[key]
        assert r["n"] == s["n"] and round(r["median"], 4) == s["median"], key
        assert [round(v, 4) for v in r["median_ci95"]] == s["median_ci95"], key
        assert [round(v, 4) for v in r["q"]] == s["p5_p25_p75_p95"], key
        for t in ("10pct", "25pct", "50pct", "100pct"):
            assert r["hi" + t] == s[f"coarser_by_more_than_{t}"]["k"] and r["lo" + t] == s[f"finer_by_more_than_{t}"]["k"], key
    extra["unrounded_summaries"] = RAW
    json.dump(extra, open(os.path.join(RES, "paper_extra.json"), "w"), indent=1)

    def ratio_block(prefix, key):
        r = RAW[key]
        n = r["n"]
        put(prefix + "N", fmt_int(n))
        put(prefix + "Med", fmt(r["median"], 3))
        put(prefix + "MedLo", fmt(r["median_ci95"][0], 3))
        put(prefix + "MedHi", fmt(r["median_ci95"][1], 3))
        q = r["q"]
        put(prefix + "Pfive", fmt(q[0], 2)); put(prefix + "Qone", fmt(q[1], 3))
        put(prefix + "Qthree", fmt(q[2], 3)); put(prefix + "Pninetyfive", fmt(q[3], 2))
        for t, w in (("10pct", "Ten"), ("25pct", "TwentyFive"), ("50pct", "Fifty"), ("100pct", "Hundred")):
            put(prefix + "Coarser" + w, fmt(100 * r["hi" + t] / n, 1)); put(prefix + "Coarser" + w + "K", fmt_int(r["hi" + t]))
            put(prefix + "Finer" + w, fmt(100 * r["lo" + t] / n, 1)); put(prefix + "Finer" + w + "K", fmt_int(r["lo" + t]))
        lo, hi = wil(r["hi25pct"], n)
        put(prefix + "CoarserTwentyFiveLo", fmt(lo, 1))
        put(prefix + "CoarserTwentyFiveHi", fmt(hi, 1))

    ratio_block("un", "stated_vs_unmasked")
    ratio_block("ma", "stated_vs_masked")
    ratio_block("co", "stated_vs_corrected")
    ratio_block("own", "stated_vs_own_curve")
    ratio_block("oc", "own_curve_vs_corrected")
    ratio_block("um", "unmasked_over_masked")
    rc = RAW["stated_vs_corrected"]
    put("coWithinTen", fmt(100 * (rc["n"] - rc["hi10pct"] - rc["lo10pct"]) / rc["n"], 1))
    own = [x for x in E if x["own"]]
    within1 = 100 * sum(abs(x["own"] / x["stated"] - 1) <= 0.01 for x in own) / len(own)
    assert round(within1, 2) == A["stated_vs_own_curve"]["within_1pct"]
    put("ownWithinOne", fmt(within1, 1))
    put("tailN", fmt_int(A["tail_n"]))
    put("tailBase", fmt_int(A["tail_base_n"]))
    put("oppN", fmt_int(A["opposite_tail_n"]))
    t = A["tail_with_own_curve"]
    put("tailOwnN", fmt_int(t["n"]))
    put("tailOwnStated", fmt_int(t["own_curve_agrees_with_stated_within_5pct"]))
    put("tailOwnCorr", fmt_int(t["own_curve_agrees_with_corrected_within_10pct"]))
    for k, w in (("0-3 A", "A"), ("3-4 A", "B"), ("4-6 A", "C"), ("6-10 A", "D"), ("10-up A", "E")):
        v = A["tail_by_band"][k]
        lo, hi = wil(v["tail"], v["n"])
        put("band" + w + "Pct", fmt(100 * v["tail"] / v["n"], 1)); put("band" + w + "K", fmt_int(v["tail"])); put("band" + w + "N", fmt_int(v["n"]))
        put("band" + w + "Lo", fmt(lo, 1)); put("band" + w + "Hi", fmt(hi, 1))
    below4 = [A["tail_by_band"][k] for k in ("0-3 A", "3-4 A")]
    above4 = [A["tail_by_band"][k] for k in ("4-6 A", "6-10 A", "10-up A")]
    kb, nb = sum(v["tail"] for v in below4), sum(v["n"] for v in below4)
    ka, na = sum(v["tail"] for v in above4), sum(v["n"] for v in above4)
    put("belowFourPct", fmt(100 * kb / nb, 1)); put("belowFourK", fmt_int(kb)); put("belowFourN", fmt_int(nb))
    put("aboveFourPct", fmt(100 * ka / na, 1)); put("aboveFourK", fmt_int(ka)); put("aboveFourN", fmt_int(na))
    put("aboveFourShare", fmt(100 * ka / A["tail_n"], 0))
    put("aboveFourEntryShare", fmt(100 * na / A["tail_base_n"], 0))
    for m, w in (("singleParticle", "SPA"), ("helical", "Hel"), ("subtomogramAveraging", "STA")):
        v = A["tail_by_method"][m]
        put("meth" + w + "Pct", fmt(100 * v["tail"] / v["n"], 1)); put("meth" + w + "K", fmt_int(v["tail"])); put("meth" + w + "N", fmt_int(v["n"]))
    sta = [B["tail_by_band_and_method"][f"{b} | subtomogramAveraging"] for b in ("0-3 A", "3-4 A", "4-6 A", "6-10 A", "10-up A")]
    put("staAboveFourShare", fmt(100 * sum(v["n"] for v in sta[2:]) / sum(v["n"] for v in sta), 0))
    spa = [B["tail_by_band_and_method"][f"{b} | singleParticle"] for b in ("0-3 A", "3-4 A", "4-6 A", "6-10 A", "10-up A")]
    put("spaAboveFourShare", fmt(100 * sum(v["n"] for v in spa[2:]) / sum(v["n"] for v in spa), 0))
    put("studiesN", fmt_int(B["n_studies"]))
    put("studiesTail", fmt_int(B["n_studies_with_tail_entry"]))
    put("topTenShare", fmt(100 * sum(c["tail_entries"] for c in B["largest_tail_contributors"]) / B["n_tail"], 0))
    put("topStudyTail", fmt_int(B["largest_tail_contributors"][0]["tail_entries"]))
    put("topStudyEntries", fmt_int(B["largest_tail_contributors"][0]["entries"]))
    put("studyCiLo", fmt(B["unrounded"]["tail_pct_study_bootstrap_ci95"][0], 1))
    put("studyCiHi", fmt(B["unrounded"]["tail_pct_study_bootstrap_ci95"][1], 1))
    put("onePerStudy", fmt(B["unrounded"]["tail_pct_one_entry_per_study_mean"], 1))
    put("withCitation", fmt_int(B["n_entries_with_citation"]))
    put("tailPct", fmt(100 * A["tail_n"] / A["tail_base_n"], 1))
    lo, hi = wil(A["tail_n"], A["tail_base_n"])
    put("tailLo", fmt(lo, 1))
    put("tailHi", fmt(hi, 1))
    put("sizeA", fmt_int(sizes["125to150"]["n"])); put("sizeB", fmt_int(sizes["150to200"]["n"])); put("sizeC", fmt_int(sizes["gt200"]["n"]))
    nk = sum(v["pixel_known"] for k, v in sizes.items() if k != "le125")
    nd = sum(v["halfmap_pixel_differs_from_map"] for k, v in sizes.items() if k != "le125")
    put("tailPixDiff", fmt_int(nd)); put("tailPixKnown", fmt_int(nk)); put("tailPixPct", fmt(100 * nd / nk, 1))
    put("restPixDiff", fmt_int(sizes["le125"]["halfmap_pixel_differs_from_map"]))
    put("restPixKnown", fmt_int(sizes["le125"]["pixel_known"]))
    put("restPixPct", fmt(100 * sizes["le125"]["halfmap_pixel_differs_from_map"] / sizes["le125"]["pixel_known"], 1))
    put("checkN", fmt_int(C["n"]))
    put("checkUn", fmt_int(C["agree_counts"]["unmasked"]))
    put("checkMa", fmt_int(C["agree_counts"]["masked"]))
    put("checkCo", fmt_int(C["agree_counts"]["corrected"]))
    put("checkRandom", fmt_int(sum(e["sample"] == "random" for e in C["entries"])))
    put("checkTail", fmt_int(sum(e["sample"] == "tail" for e in C["entries"])))
    off = [e for e in C["entries"] if not all(e["agree"].values())]
    assert len(off) == 1
    put("checkOffId", off[0]["emdb_id"])
    put("checkOffIndex", fmt(off[0]["index"]["masked"], 2))
    put("checkOffRecord", fmt(off[0]["va_record"]["masked"], 2))
    ex = {e["emdb_id"]: e for e in C["entries"]}["EMD-19020"]
    put("exStated", fmt(ex["author_A"], 2)); put("exUn", fmt(ex["index"]["unmasked"], 0)); put("exCo", fmt(ex["index"]["corrected"], 0))

    # random sample of entries whose recalculated values are missing from the index (code/sample_unindexed_va.py,
    # code/analyze_unindexed_sample.py); counts recomputed here from the sample file and checked against the saved JSON
    S = json.load(open(os.path.join(RES, "unindexed_sample_analysis.json")))
    recs, _ = UQ.load()
    for s, p, size in (("pre2024", "sp", 500), ("from2024", "sr", 150)):
        g = [r for r in recs if r["stratum"] == s]
        d = S[s]
        if REQUIRE_COMPLETE_SAMPLE:
            assert len(g) == size == d["sampled"], (s, len(g), d["sampled"])
        st = {}
        for r in g:
            st.setdefault(UQ.classify(r), []).append(r)
        assert {k: len(v) for k, v in st.items()} == d["status"], s
        E_s = st.get("compared", [])
        base_s = [r for r in E_s if r.get("corrected")]
        tail_s = [r for r in base_s if r["corrected"] / r["stated"] > 1.25]
        assert (len(E_s), len(base_s), len(tail_s)) == (d["compared"], d["with_corrected"], d["tail"]), s
        put(p + "Pool", fmt_int(d["pool_n"]))
        put(p + "Sampled", fmt_int(len(g)))
        put(p + "Failed", fmt_int(len(st.get("failed", []))))
        put(p + "WithRecord", fmt_int(len(g) - len(st.get("failed", []))))
        put(p + "NoRelion", fmt_int(len(st.get("no_relion_fsc", []))))
        put(p + "NoUnmasked", fmt_int(len(st.get("no_unmasked", []))))
        put(p + "NyqExcl", fmt_int(len(st.get("nyquist_excluded", []))))
        put(p + "Compared", fmt_int(len(E_s)))
        put(p + "Corr", fmt_int(len(base_s)))
        put(p + "TailK", fmt_int(len(tail_s)))
        put(p + "TailPct", fmt(100 * len(tail_s) / len(base_s), 1))
        lo, hi = wil(len(tail_s), len(base_s))
        put(p + "TailLo", fmt(lo, 1)); put(p + "TailHi", fmt(hi, 1))
        # the abstract and conclusion call the sampled shares similar to the indexed one: each interval contains it
        assert lo <= 100 * A["tail_n"] / A["tail_base_n"] <= hi, (s, lo, hi)
        ratio = np.array([r["corrected"] / r["stated"] for r in base_s])
        u = d["unrounded"]
        assert abs(float(np.median(ratio)) - u["median"]) < 1e-12
        put(p + "Med", fmt(u["median"], 3))
        put(p + "Qone", fmt(u["p25_p75"][0], 3)); put(p + "Qthree", fmt(u["p25_p75"][1], 3))
        put(p + "WithinTen", fmt(u["within_10pct"], 1))
        put(p + "UnMed", fmt(u["unmasked_median"], 3))
        put(p + "UnCoarser", fmt(u["unmasked_coarser_by_more_than_25pct"], 1))
        for key, w in (("stated_below_4A", "BelowFour"), ("stated_4A_or_coarser", "AboveFour")):
            v = d[key]
            put(p + w + "K", fmt_int(v["tail"])); put(p + w + "N", fmt_int(v["n"]))
            put(p + w + "Pct", fmt(100 * v["tail"] / v["n"], 1) if v["n"] else "--")
            if v["n"]:
                lo, hi = wil(v["tail"], v["n"])
                put(p + w + "Lo", fmt(lo, 1)); put(p + w + "Hi", fmt(hi, 1))
        if s == "pre2024":
            # the text says this share is lower than among indexed entries but its interval is wide enough to include it
            v = d["stated_4A_or_coarser"]
            lo, hi = wil(v["tail"], v["n"])
            assert 100 * v["tail"] / v["n"] < 100 * ka / na and lo <= 100 * ka / na <= hi, (v, ka, na)
        put(p + "OwnN", fmt_int(d["own_curve"]["n"]))
        put(p + "OwnWithinOne", fmt(100 * d["own_curve"]["agrees_with_stated_within_1pct"] / d["own_curve"]["n"], 1) if d["own_curve"]["n"] else "--")
        put(p + "TailOwnN", fmt_int(d["tail_own_curve"]["n"]))
        put(p + "TailOwnStated", fmt_int(d["tail_own_curve"]["agrees_with_stated_within_5pct"]))
    cu = S["combined"]["unrounded"]
    put("combPct", fmt(cu["pct"], 1)); put("combLo", fmt(cu["bootstrap_ci95"][0], 1)); put("combHi", fmt(cu["bootstrap_ci95"][1], 1))
    put("combBase", fmt_int(int(round(cu["base"], -2))))
    O = json.load(open(os.path.join(RES, "own_curve_field_check.json")))       # code/check_own_curve_field.py
    oc = O["counts"]
    assert O["n"] == C["n"] and sum(oc.values()) == O["n"] and oc.get("failed", 0) == 0, oc
    put("ofN", fmt_int(O["n"]))
    put("ofBoth", fmt_int(oc.get("both_agree", 0) + oc.get("both_differ", 0)))
    put("ofAgree", fmt_int(oc.get("both_agree", 0)))
    put("ofNeither", fmt_int(oc.get("neither", 0)))
    put("ofOneOnly", fmt_int(oc.get("index_only", 0) + oc.get("record_only", 0)))
    # the text says the two agreed in all entries that had both, and that the others had neither
    assert oc.get("both_differ", 0) == 0 and oc.get("index_only", 0) + oc.get("record_only", 0) == 0, oc

    os.makedirs(PAPER, exist_ok=True)
    with open(os.path.join(PAPER, "numbers.tex"), "w") as fh:
        fh.write("% generated by code/make_paper_inputs.py; do not edit\n")
        for k, v in N.items():
            fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")

    # Table: selection
    rows = [("EMDB index rows (all statuses)", flow["index_rows"]),
            ("Released, not withdrawn or obsolete", flow["released_current"]),
            ("Single particle, helical or subtomogram averaging", flow["pool"]),
            ("Stated criterion FSC 0.143", flow["crit_0143"]),
            ("\\quad with half-maps listed in the index", flow["crit_0143_halfmap_listed"]),
            ("\\quad with an indexed EMDB recalculation", flow["crit_0143_with_indexed_recalculation"]),
            ("Compared, after the sampling-limit exclusion", flow["compared"]),
            ("\\quad with the masked recalculation", A["stated_vs_masked"]["n"]),
            ("\\quad with the corrected recalculation", A["stated_vs_corrected"]["n"]),
            ("\\quad with the depositor's FSC curve", A["stated_vs_own_curve"]["n"])]
    with open(os.path.join(PAPER, "tab_selection.tex"), "w") as fh:
        fh.write("% generated by code/make_paper_inputs.py\n")
        for lab, n in rows:
            fh.write(f"{lab} & {fmt_int(n)} \\\\\n")

    # Table: coverage by release year
    with open(os.path.join(PAPER, "tab_coverage.tex"), "w") as fh:
        fh.write("% generated by code/make_paper_inputs.py\n")
        for y, v in cy.items():
            pct = 100 * v["halfmap_and_indexed"] / v["halfmap"] if v["halfmap"] else 0
            fh.write(f"{y.replace('-', '--')} & {fmt_int(v['pool'])} & {fmt_int(v['halfmap'])} & {fmt_int(v['halfmap_and_indexed'])} & {pct:.1f} \\\\\n")

    # Table: agreement summary
    with open(os.path.join(PAPER, "tab_agreement.tex"), "w") as fh:
        fh.write("% generated by code/make_paper_inputs.py\n")
        for key, lab in (("stated_vs_own_curve", "Depositor's FSC curve"), ("stated_vs_unmasked", "EMDB, unmasked"),
                         ("stated_vs_masked", "EMDB, EMDB mask"), ("stated_vs_corrected", "EMDB, corrected")):
            r = RAW[key]
            n, q = r["n"], r["q"]
            fh.write(f"{lab} & {fmt_int(n)} & {r['median']:.3f} & {q[1]:.3f}--{q[2]:.3f} & "
                     f"{100 * r['hi10pct'] / n:.1f} & {100 * r['hi25pct'] / n:.1f} & "
                     f"{100 * r['hi100pct'] / n:.2f} & {100 * r['lo25pct'] / n:.2f} \\\\\n")

    # Table: tail by band and method
    with open(os.path.join(PAPER, "tab_band_method.tex"), "w") as fh:
        fh.write("% generated by code/make_paper_inputs.py\n")
        for b, lab in (("0-3 A", "$<3$"), ("3-4 A", "3 to $<4$"), ("4-6 A", "4 to $<6$"), ("6-10 A", "6 to $<10$"), ("10-up A", "$\\geq 10$")):
            cells = []
            for m in ("singleParticle", "helical", "subtomogramAveraging"):
                v = B["tail_by_band_and_method"][f"{b} | {m}"]
                cells.append(f"{v['tail']}/{fmt_int(v['n'])} ({100 * v['tail'] / v['n']:.1f})")
            a = A["tail_by_band"][b]
            cells.append(f"{a['tail']}/{fmt_int(a['n'])} ({100 * a['tail'] / a['n']:.1f})")
            fh.write(lab + " & " + " & ".join(cells) + " \\\\\n")

    # Table: indexed entries against the random sample of unindexed entries
    with open(os.path.join(PAPER, "tab_unindexed.tex"), "w") as fh:
        fh.write("% generated by code/make_paper_inputs.py\n")
        rc = RAW["stated_vs_corrected"]
        lo, hi = wil(A["tail_n"], A["tail_base_n"])
        fh.write(f"In the index (all) & {fmt_int(A['tail_base_n'])} & {rc['median']:.3f} & "
                 f"{100 * (rc['n'] - rc['hi10pct'] - rc['lo10pct']) / rc['n']:.1f} & "
                 f"{100 * A['tail_n'] / A['tail_base_n']:.1f} ({lo:.1f}--{hi:.1f}) & {100 * kb / nb:.1f} & {100 * ka / na:.1f} \\\\\n")
        for s, lab in (("pre2024", "Not in the index, released before 2024 (sample)"),
                       ("from2024", "Not in the index, released 2024 or later (sample)")):
            d = S[s]
            u = d["unrounded"]
            lo, hi = wil(d["tail"], d["with_corrected"])
            b4, a4 = d["stated_below_4A"], d["stated_4A_or_coarser"]
            fh.write(f"{lab} & {fmt_int(d['with_corrected'])} & {u['median']:.3f} & {u['within_10pct']:.1f} & "
                     f"{100 * d['tail'] / d['with_corrected']:.1f} ({lo:.1f}--{hi:.1f}) & "
                     f"{100 * b4['tail'] / b4['n']:.1f} & {100 * a4['tail'] / a4['n']:.1f} \\\\\n")

    # Supplementary Table S1: the entries whose corrected recalculation is more than 25% coarser than stated
    cols = ["emdb_id", "method", "release_year", "stated_A", "own_curve_A", "unmasked_A", "masked_A", "corrected_A",
            "corrected_over_stated"]
    rows_t = list(csv.DictReader(open(os.path.join(RES, "tail_entries.csv"), newline="")))
    assert len(rows_t) == A["tail_n"]
    with open(os.path.join(PAPER, "Supplementary_Table_S1.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows_t:
            w.writerow([METHOD_NAME.get(r["method"], r["method"]) if c == "method" else r[c] for c in cols])

    figures(E)
    fig3(B)
    print(f"{len(N)} numbers; flow {json.dumps({k: v for k, v in flow.items() if k != 'coverage_by_year'})}")
    print("tail sizes", json.dumps(sizes))


if __name__ == "__main__":
    main(sys.argv[1])
