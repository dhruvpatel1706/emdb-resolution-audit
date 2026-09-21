#!/usr/bin/env python3
"""Artifact-derived number guard for the EMDB metadata audit (Acta D QE5016 revision candidate).

Re-derives every statistic quoted in the manuscript from data/results/*.json, paper/jac/COHORT_CHECKS.json
and, where no result file carries it, directly from data/raw/emdb_detail.csv (the 1,600-entry sample), with
decimal half-up rounding at the printed precision, and asserts the printed string is in the TeX. The two
supplementary tables are regenerated row by row from the CSV.

  python3 code/check_paper_numbers.py                                  # revision candidate (default)
  EMDB_TEX=paper/acta_d/emdb_ActaD_main.tex python3 code/check_paper_numbers.py   # 1 Sep 2026 submission
Written 21 September 2026. Exit status is non-zero on any failure.
"""
import csv, json, math, os, re, sys
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
TEX = os.environ.get("EMDB_TEX", os.path.join(ROOT, "paper", "acta_d_revision", "emdb_ActaD_main_rev1.tex"))
def J(p): return json.load(open(os.path.join(ROOT, p)))
def R(x, nd): return str(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))
def N(x): return f"{int(x):,}".replace(",", "{,}")
def wilson(k, n, z=1.959963984540054):
    ph = k / n; d = 1 + z * z / n; c = (ph + z * z / (2 * n)) / d
    h = (z / d) * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)); return 100 * (c - h), 100 * (c + h)
def W(k, n): lo, hi = wilson(k, n); return f"[{R(lo,1)}, {R(hi,1)}]" if hi < 99.95 else f"[{R(lo,1)}, 100]"
def norm(t):
    t = re.sub(r"(?<!\\)%.*", "", t)
    t = t.replace(r"\%", "%").replace(r"\,", "").replace(r"\ ", " ").replace("{,}", ",").replace(r"\AA", "A")
    t = re.sub(r"\\(textbf|texttt|emph|text|mathrm)", "", t)
    t = t.replace("$", "").replace("{", "").replace("}", "").replace("~", " ").replace("--", "-")
    t = re.sub(r",\s+", ",", t)
    return re.sub(r"\s+", " ", t)
tex = norm(open(TEX).read())
checks = []
def want(label, s): checks.append((label, norm(s)))

E = J("data/results/emdb.json"); BM = E["by_method"]; RT = E["resolution_over_time"]; PR = E["pdb_rate_by_resolution"]; SB = E["substantiation"]; NY = E["nyquist"]
B = J("data/results/bfactor_law.json"); CF = B["core_fit"]; SC = B["symmetry_control"]; CT = B["counts"]
FW = J("data/results/fsc0143_wilson_ci.json")["buckets"]; NS = J("data/results/nyquist_sensitivity.json"); CH = J("paper/jac/COHORT_CHECKS.json")
rows = list(csv.DictReader(open(os.path.join(ROOT, "data/raw/emdb_detail.csv")))); assert len(rows) == 1600
census = list(csv.DictReader(open(os.path.join(ROOT, "data/raw/emdb_census.csv")))); assert len(census) == 46900
n_res = sum(1 for r in census if r["resolution"] not in ("", "nan")); assert n_res == RT["n_with_resolution"] == CH["sampling_frame_n"] == 45151

# ---- census
want("census n", f"an indexed census of ${N(len(census))}$ entries (${N(n_res)}$ carrying a resolution) and a stratified ${N(len(rows))}$-entry sample")
mb = RT["median_by_year"]
want("year medians", f"in ${2015}$ it drops to ${R(mb['2015']['median_res'],1)}$\\,\\AA, the signature")
want("year medians 2", f"to ${R(mb['2018']['median_res'],1)}$\\,\\AA{{}} by ${2018}$ and ${R(mb['2024']['median_res'],1)}$\\,\\AA{{}} by ${2024}$")
f4 = RT["frac_sub4A_by_year"]
want("sub4 fracs", f"rises from ${R(100*f4['2012'],0)}\\%$ in ${2012}$ to ${R(100*f4['2018'],0)}\\%$ in ${2018}$ to ${R(100*f4['2024'],0)}\\%$ in ${2024}$")
want("overall median", f"The overall median across all ${N(n_res)}$ resolutions is ${R(RT['overall_median'],1)}$\\,\\AA")
assert 8000 < mb["2023"]["n"] < 10000 and mb["2024"]["n"] > 7000     # "about nine thousand"
sp, hl, cr, tm, st = BM["single particle"], BM["helical"], BM["crystallography"], BM["tomography"], BM["subtomogram avg"]
want("method prose", f"the bulk of the archive at ${N(sp['n'])}$ entries, has a median of ${R(sp['median_res'],1)}$\\,\\AA, helical ${R(hl['median_res'],1)}$\\,\\AA, and the small electron (2D) crystallography set ${R(cr['median_res'],1)}$\\,\\AA, while tomography and subtomogram averaging, which image crowded cellular volumes rather than purified particles, sit at ${R(tm['median_res'],1)}$ and ${R(st['median_res'],1)}$\\,\\AA")
want("tomography with res", f"only ${N(tm['n_with_res'])}$ of ${N(tm['n'])}$ tomography entries")
want("T method crystallography", f"electron (2D) crystallography & ${cr['n']}$    & ${cr['n_with_res']}$    & ${R(cr['median_res'],1)}$ \\\\")
want("T method single", f"single particle        & ${N(sp['n'])}$ & ${N(sp['n_with_res'])}$ & ${R(sp['median_res'],1)}$ \\\\")
want("T method helical", f"helical                & ${N(hl['n'])}$  & ${N(hl['n_with_res'])}$  & ${R(hl['median_res'],1)}$ \\\\")
want("T method subtomo", f"subtomogram averaging  & ${N(st['n'])}$  & ${N(st['n_with_res'])}$  & ${R(st['median_res'],1)}$ \\\\")
want("T method tomography", f"tomography             & ${N(tm['n'])}$  & ${tm['n_with_res']}$    & ${R(tm['median_res'],1)}$ \\\\")
want("pdb rates", f"from ${R(100*PR['0-3A']['pdb_rate'],0)}\\%$ below $3$\\,\\AA{{}} to ${R(100*PR['3-4A']['pdb_rate'],0)}\\%$ in the $3$--$4$\\,\\AA{{}} band, ${R(100*PR['4-6A']['pdb_rate'],0)}\\%$ at $4$--$6$\\,\\AA, ${R(100*PR['6-10A']['pdb_rate'],0)}\\%$ at $6$--$10$\\,\\AA, and ${R(100*PR['10-infA']['pdb_rate'],0)}\\%$ beyond $10$\\,\\AA")

# ---- sample: criterion provenance
assert SB["res_type_counts"] == {"BY AUTHOR": 1600}
mc = SB["method_counts"]; n = SB["n"]
want("criterion shares", f"is cited by ${R(100*mc['FSC 0.143 CUT-OFF']/n,0)}\\%$ of the sample overall, the older and looser FSC $0.5$ by ${R(100*mc['FSC 0.5 CUT-OFF']/n,0)}\\%$, and an unspecified ``other'' by only ${R(100*mc['OTHER']/n,0)}\\%$")
want("post-stratified", f"This ${R(100*SB['frac_fsc0143'],1)}\\%$ is computed on the $50/50$ pre/post sample; re-weighting the sample's earlier- and later-stratum rates by the census split at the same boundary (${N(SB['census_pre_post']['pre'])}/{N(SB['census_pre_post']['post'])}$) gives a post-stratified estimate of ${R(100*SB['frac_fsc0143_census_weighted'],1)}\\%$")
assert SB["census_pre_post"]["pre"] + SB["census_pre_post"]["post"] == len(census)
want("method mix", f"the sample's method mix (${R(100*SB['sample_single_particle_frac'],1)}\\%$ single-particle) is close to the census (${R(100*SB['census_single_particle_frac'],1)}\\%$)")
bk = ["2002-2013", "2014-2017", "2018-2021", "2022-2026"]
for b in bk: assert FW[b]["n_fsc0143"] == SB["fsc0143_by_year_bucket"][b]["n_fsc0143"] and FW[b]["n"] == SB["fsc0143_by_year_bucket"][b]["n"]
def bw(b): return f"[{R(FW[b]['wilson95_lo_pct'],1)}, {R(FW[b]['wilson95_hi_pct'],1)}]"
want("fsc buckets", f"rises from ${R(FW[bk[0]]['pct'],1)}\\%$ ($95\\%$ Wilson interval ${bw(bk[0])}$) for entries deposited through ${2013}$ to ${R(FW[bk[1]]['pct'],1)}\\%$ ${bw(bk[1])}$ in ${2014}$--${2017}$, ${R(FW[bk[2]]['pct'],1)}\\%$ ${bw(bk[2])}$ in ${2018}$--${2021}$, and ${R(FW[bk[3]]['pct'],1)}\\%$ ${bw(bk[3])}$ from ${2022}$ on (per-bucket counts ${FW[bk[0]]['n_fsc0143']}/{FW[bk[0]]['n']}$, ${FW[bk[1]]['n_fsc0143']}/{FW[bk[1]]['n']}$, ${FW[bk[2]]['n_fsc0143']}/{FW[bk[2]]['n']}$, ${FW[bk[3]]['n_fsc0143']}/{FW[bk[3]]['n']}$)")
for b in bk:
    lo, hi = wilson(FW[b]["n_fsc0143"], FW[b]["n"]); assert abs(lo - FW[b]["wilson95_lo_pct"]) < 0.051 and abs(hi - FW[b]["wilson95_hi_pct"]) < 0.051
want("abstract buckets", f"(${R(FW[bk[0]]['pct'],0)}\\%$ before ${2014}$, ${R(FW[bk[1]]['pct'],0)}\\%$ in ${2014}$--${2017}$, ${R(FW[bk[2]]['pct'],0)}\\%$ in ${2018}$--${2021}$, ${R(FW[bk[3]]['pct'],0)}\\%$ after)")
want("synopsis buckets", f"show FSC $0.143$ use of ${R(FW[bk[2]]['pct'],0)}$--${R(FW[bk[3]]['pct'],0)}\\%$ after ${2018}$")
want("before after", f"(${R(100*SB['fsc0143_pre_mandate'],0)}\\%$ to ${R(100*SB['fsc0143_post_mandate'],0)}\\%$)")
assert SB["n_pre"] == SB["n_post"] == 800

# ---- later stratum from the CSV
later = [r for r in rows if int(r["deposition"]) >= 20220201]; non = [r for r in later if "0.143" not in r["resolution_method"]]
assert len(later) == CH["later_stratum_n"] == 800 and len(non) == CH["later_stratum_non_fsc0143"] == 34
want("later 34", f"Of the ${len(later)}$ later-stratum entries, ${len(non)}$ (${R(100*len(non)/len(later),1)}\\%$, $95\\%$ Wilson interval ${W(len(non), len(later))}$) do not cite FSC $0.143$")
cc = Counter(r["resolution_method"] for r in non)
want("later criteria", f"${cc['FSC 0.5 CUT-OFF']}$ of the ${len(non)}$ report the looser FSC $0.5$ and ${cc['OTHER']}$ fall in the unspecified ``other'' bin, with single instances of the FSC $3\\sigma$ and FSC $0.33$ cut-offs and one entry recording no criterion at all (${cc['FSC 0.5 CUT-OFF']}+{cc['OTHER']}+{cc['FSC 3 SIGMA CUT-OFF']}+{cc['FSC 0.33 CUT-OFF']}+{cc['']}={len(non)}$)")
assert cc["FSC 3 SIGMA CUT-OFF"] == cc["FSC 0.33 CUT-OFF"] == cc[""] == 1
lm, nm = Counter(r["method"] for r in later), Counter(r["method"] for r in non)
want("later single", f"has a ${R(100*nm['singleParticle']/lm['singleParticle'],1)}\\%$ non-FSC-$0.143$ proportion (${nm['singleParticle']}$ of ${lm['singleParticle']}$, ${W(nm['singleParticle'], lm['singleParticle'])}$)")
want("abstract single", f"single-particle reconstruction accounts for ${nm['singleParticle']}$ of ${lm['singleParticle']}$ (${R(100*nm['singleParticle']/lm['singleParticle'],1)}\\%$, $95\\%$ CI ${W(nm['singleParticle'], lm['singleParticle'])}$)")
want("later helical subtomo", f"(${nm['helical']}$ of ${lm['helical']}$, ${W(nm['helical'], lm['helical'])}$; and ${nm['subtomogramAveraging']}$ of ${lm['subtomogramAveraging']}$, ${W(nm['subtomogramAveraging'], lm['subtomogramAveraging'])}$)")
assert abs(100 * nm["helical"] / lm["helical"] - 8) < 1 and abs(100 * nm["subtomogramAveraging"] / lm["subtomogramAveraging"] - 8) < 1     # "each sit near 8%"
assert nm["tomography"] == lm["tomography"] == 5
want("later tomography", f"all five sampled later-stratum tomography entries report a non-FSC-$0.143$ criterion (${W(5, 5)}$)")
want("tomography three", f"widens it from ${W(5, 5)}$ to ${W(3, 3)}$")
dates = Counter(r["deposition"] for r in non)
want("dates", f"The ${len(non)}$ exceptions occupy only ${len(dates)}$ distinct deposition dates, and {['zero','one','two','three','four','five','six'][sum(1 for v in dates.values() if v > 1)]} dates carry multiple entries (one carries {['zero','one','two','three','four','five'][max(dates.values())]})")
tomo = sorted(r["deposition"] for r in non if r["method"] == "tomography"); assert len(set(tomo)) == 4
yr = {y: [r for r in later if r["deposition"].startswith(y)] for y in ("2022", "2023", "2024", "2025")}
def ys(y): L = yr[y]; k = sum(1 for r in L if "0.143" not in r["resolution_method"]); return f"${R(100*k/len(L),1)}\\%$ in ${y}$ (${W(k, len(L))}$)"
want("later by year", f"The later-stratum entry-level share is {ys('2022')}, {ys('2023')}, and {ys('2024')}, with the ${2025}$ slice too small at ${len(yr['2025'])}$ entries to read")
# policy-date restriction
pol = [r for r in rows if int(r["deposition"]) >= 20220225]; pnon = [r for r in pol if "0.143" not in r["resolution_method"]]
btw = [r for r in later if int(r["deposition"]) < 20220225]; bnon = [r for r in btw if "0.143" not in r["resolution_method"]]
assert len(pol) == CH["on_or_after_policy_n"] == 781 and len(pnon) == CH["on_or_after_policy_non_fsc0143"] == 33 and len(btw) == CH["between_boundaries_n"] == 19
assert [r["emdb_id"] for r in bnon] == CH["between_boundaries_non_fsc0143"] == ["EMD-14430"] and bnon[0]["deposition"] == "20220223" and bnon[0]["resolution_method"] == "OTHER"
psp = [r for r in pol if r["method"] == "singleParticle"]; psn = [r for r in psp if "0.143" not in r["resolution_method"]]
assert len(psp) == CH["date_subset_methods"]["singleParticle"]["n"] and len(psn) == CH["date_subset_methods"]["singleParticle"]["non_fsc0143"]
want("policy date", f"It contains ${len(btw)}$ entries from the intervening dates, one of which, EMD-14430, deposited on ${23}$ February, records the criterion ``OTHER.'' Restricting to dates on or after the effective date leaves ${len(pol)}$ entries and ${len(pnon)}$ non-FSC-$0.143$ records, and the single-particle subset then has ${len(psn)}$ such records among ${len(psp)}$ (${R(100*len(psn)/len(psp),1)}\\%$)")
want("abstract policy", f"leaves ${len(pnon)}$ of ${len(pol)}$")

# ---- Rosenthal-Henderson
assert CT["n_fit_singleParticle"] == 1403 and CT["n_fit_all_methods"] == 1485 and SC["n"] == 602 and B["_meta"]["n_boot"] == 10000
want("rh n", f"across the ${N(CT['n_fit_singleParticle'])}$ single-particle entries in our sample that report both a particle count and a resolution")
want("rh n all", f"would raise the sample to ${N(CT['n_fit_all_methods'])}$")
want("rh slope", f"a positive slope $b = {R(CF['slope_b'],4)}$ ($95\\%$ bootstrap interval $[{R(CF['slope_b_ci95'][0],4)}, {R(CF['slope_b_ci95'][1],4)}]$, ${N(B['_meta']['n_boot'])}$ resamples)")
want("rh B", f"the implied aggregate B-factor is ${R(CF['effective_B_A2'],1)}$\\,\\AA$^2$ ($95\\%$ bootstrap interval $[{R(CF['effective_B_A2_ci95'][0],1)}, {R(CF['effective_B_A2_ci95'][1],1)}]$)")
want("rh sym", f"raises the variance explained from $R^2 = {R(SC['r2_lnN_on_subset'],3)}$ to $R^2 = {R(SC['r2_lnN_times_order_on_subset'],3)}$ on that subset")
want("rh sym again", f"We report this $R^2$ change (${R(SC['r2_lnN_on_subset'],3)}$ to ${R(SC['r2_lnN_times_order_on_subset'],3)}$)")
want("rh r2 r", f"explains only about ${R(100*CF['r2'],0)}\\%$ of the variance in $1/d^2$ across the single-particle sample (Pearson $r = {R(CF['pearson_r'],2)}$)")

# ---- Nyquist
rd = NS["ratio_distribution"]; assert NY["n"] == rd["n"] == 1600 and NY["n_sub_nyquist"] == 3 and NY["n_sub_nyquist_pixel_artifact"] == 2
want("nyquist median", f"Across the sample its median is ${R(NY['median_ratio'],2)}$")
want("abstract nyquist", f"median reported resolution at ${R(NY['median_ratio'],2)}$ times the sampling limit")
want("nyquist three", f"Just three of the ${N(NY['n'])}$ sampled entries (${R(100*NY['frac_sub_nyquist'],1)}\\%$) report a resolution below the Nyquist limit")
se = {e["emdb_id"]: e for e in NY["sub_nyquist_entries"]}
for k in ("EMD-15154", "EMD-45648"): assert se[k]["res_eq_pixel"] and abs(se[k]["resolution"] / se[k]["pixel_x"] - 1) < 0.001
want("nyquist artifacts", f"(EMD-${15154}$ reports ${se['EMD-15154']['resolution']}$\\,\\AA{{}} at a ${se['EMD-15154']['pixel_x']}$\\,\\AA{{}} pixel, EMD-${45648}$ reports ${se['EMD-45648']['resolution']}$\\,\\AA{{}} at a ${se['EMD-45648']['pixel_x']}$\\,\\AA{{}} pixel)")
b44 = NS["borderline_entry"]; assert b44["emdb_id"] == "EMD-44581"
want("nyquist borderline", f"reporting ${b44['resolution']}$\\,\\AA{{}} at a ${b44['pixel_x']}$\\,\\AA{{}} pixel, whose Nyquist limit is ${b44['nyquist_limit']}$\\,\\AA. At a ratio of ${R(b44['ratio'],3)}$ this reported value is only ${R(b44['pct_inside_nyquist'],1)}\\%$ inside its sampling limit")
want("abstract borderline", f"is only ${R(b44['pct_inside_nyquist'],1)}\\%$ inside its limit (${b44['resolution']}$\\,\\AA{{}} at a ${b44['pixel_x']}$\\,\\AA{{}} pixel)")
want("ratio dist", f"the ratio itself has mean ${R(rd['mean'],2)}$, standard deviation ${R(rd['std'],2)}$, and interquartile range ${R(rd['iqr'],2)}$ (its first percentile is already ${R(rd['p01'],2)}$)")
bd = NS["bands"]; assert bd["2pct"]["n_sub_nyquist_within_band"] == 0 and bd["5pct"]["n_sub_nyquist_within_band"] == bd["10pct"]["n_sub_nyquist_within_band"] == 1
want("bands", f"A $\\pm 2\\%$ band ($r \\in [0.98, 1.02]$) covers just ${bd['2pct']['n_within_band']}$ entries and absorbs none of the strictly sub-Nyquist three; a $\\pm 5\\%$ band (${bd['5pct']['n_within_band']}$ entries) and a $\\pm 10\\%$ band (${bd['10pct']['n_within_band']}$ entries) each absorb exactly the one borderline entry")
assert b44["smallest_band_absorbing"] == 0.05 and 2.5 < b44["pct_inside_nyquist"] < 3.5      # "at least about 3%"
want("nyquist rate", f"the substantive rate is ${NY['n_sub_nyquist_substantive']}$ in ${N(NY['n'])}$ (${R(100*NY['frac_sub_nyquist_substantive'],2)}\\%$). In total, ${R(100*NY['frac_within_1p2x_nyquist'],1)}\\%$ come within ${1.2}$ times of the Nyquist limit")

# ---- Table S1 (34 rows, deposition order) and S2
PM = {"singleParticle": "single particle       ", "tomography": "tomography            ", "subtomogramAveraging": "subtomogram averaging ", "helical": "helical               "}
CR = {"FSC 0.5 CUT-OFF": "FSC 0.5", "OTHER": "other", "": "not reported", "FSC 3 SIGMA CUT-OFF": "FSC $3\\sigma$", "FSC 0.33 CUT-OFF": "FSC 0.33"}
for r in sorted(non, key=lambda r: (r["deposition"], r["emdb_id"])):
    d = r["deposition"]; res = r["resolution"]
    want("S1 " + r["emdb_id"], f"{r['emdb_id']} & {PM[r['method']]} & {d[:4]}-{d[4:6]}-{d[6:]} & {res} & {CR[r['resolution_method']]} \\\\")
for e in NY["sub_nyquist_entries"]:
    r = next(x for x in rows if x["emdb_id"] == e["emdb_id"]); d = r["deposition"]
    want("S2 " + e["emdb_id"], f"{e['emdb_id']} & {PM[e['method']].strip()} & {d[:4]}-{d[4:6]}-{d[6:]} & {e['resolution']} & {e['pixel_x']} & {R(2*e['pixel_x'],3).rstrip('0').rstrip('.') if e['emdb_id'] != 'EMD-44581' else R(2*e['pixel_x'],2)} & {R(e['ratio'],3)} &")

fails = [(l, s) for l, s in checks if s not in tex]
out = {"tex": os.path.relpath(TEX, ROOT), "n_checks": len(checks), "n_failed": len(fails), "failed": fails}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(TEX)), "number_check.json"), "w"), indent=1)
print(f"{len(checks)} checks, {len(fails)} failed ({os.path.relpath(TEX, ROOT)})")
for l, s in fails: print("  FAIL", l, "->", s[:200])
sys.exit(1 if fails else 0)
