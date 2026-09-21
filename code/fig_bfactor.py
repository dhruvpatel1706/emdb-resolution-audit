"""Figure: the Rosenthal-Henderson form in pooled EMDB metadata.

Left  : single-particle 1/res^2 vs ln(N_images) with the OLS fit line; the
        predictive slope is positive and significant but R^2 is only ~0.31,
        foregrounding that particle count is a real but minor archive-scale
        predictor.
Right : symmetry control on the point-group subset, R^2 for ln(N) vs
        ln(N * point_group_order); the theoretically predicted variable
        raises R^2 (positive control).

Reads data/raw/emdb_detail.csv and data/results/bfactor_law.json.
Outputs paper/figures/fig_bfactor.png. CPU only, deterministic.
"""
import csv
import json
import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DETAIL = os.path.join(ROOT, "data", "raw", "emdb_detail.csv")
RES = os.path.join(ROOT, "data", "results", "bfactor_law.json")
OUT = os.path.join(ROOT, "paper", "figures", "fig_bfactor.png")


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


d = json.load(open(RES))
cf = d["core_fit"]
sc = d["symmetry_control"]

rows = list(csv.DictReader(open(DETAIL)))
rec = [(num(r["n_images"]), num(r["resolution"]))
       for r in rows if r["method"] == "singleParticle"]
rec = [(n, res) for n, res in rec if n and n > 0 and res and res > 0]
lnN = np.array([math.log(n) for n, res in rec])
inv = np.array([1.0 / res ** 2 for n, res in rec])

fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.2, 3.6))

# Left: scatter + fit.
axL.scatter(lnN, inv, s=6, alpha=0.25, color="#1f4e79", edgecolors="none", zorder=2)
xs = np.linspace(lnN.min(), lnN.max(), 100)
ys = cf["intercept_a"] + cf["slope_b"] * xs
axL.plot(xs, ys, color="#b35900", lw=2.0, zorder=3)
axL.set_xlabel(r"$\ln(N_{\mathrm{images}})$", fontsize=9.5)
axL.set_ylabel(r"$1/d^2$  ($\mathrm{\AA}^{-2}$)", fontsize=9.5)
axL.set_title("single particle (n=%d)" % d["counts"]["n_fit_singleParticle"],
              fontsize=9.5)
txt = (r"$b=%.4f$ [%.4f, %.4f]" % (
           cf["slope_b"], cf["slope_b_ci95"][0], cf["slope_b_ci95"][1]) + "\n" +
       r"$R^2=%.3f$" % cf["r2"] + "\n" +
       r"eff. $B=%.1f\,\mathrm{\AA}^2$" % cf["effective_B_A2"])
axL.text(0.04, 0.96, txt, transform=axL.transAxes, fontsize=8.5,
         va="top", ha="left",
         bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.7", alpha=0.9))
axL.tick_params(labelsize=8.5)
for s in ("top", "right"):
    axL.spines[s].set_visible(False)

# Right: symmetry control bar.
labels = [r"$\ln(N)$", r"$\ln(N\times\,$order$)$"]
vals = [sc["r2_lnN_on_subset"], sc["r2_lnN_times_order_on_subset"]]
colors = ["#7f7f7f", "#1f4e79"]
bars = axR.bar([0, 1], vals, color=colors, width=0.6, zorder=3)
for x, v in zip([0, 1], vals):
    axR.text(x, v + 0.004, "%.3f" % v, ha="center", va="bottom", fontsize=9)
axR.set_xticks([0, 1])
axR.set_xticklabels(labels, fontsize=9)
axR.set_ylabel(r"$R^2$", fontsize=9.5)
axR.set_ylim(0, max(vals) * 1.18)
axR.set_title("symmetry control (point-group subset, n=%d)" % sc["n"],
              fontsize=9.0)
axR.annotate(r"$+%.3f$" % sc["r2_increase"], xy=(0.5, max(vals)),
             ha="center", va="bottom", fontsize=9, color="#1f4e79")
axR.tick_params(labelsize=8.5)
for s in ("top", "right"):
    axR.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig(OUT, dpi=600, bbox_inches="tight")  # 600 d.p.i. per Acta D notes for authors
print("wrote", OUT)
print("b=%.4f CI[%.4f,%.4f] R2=%.3f effB=%.1f" % (
    cf["slope_b"], cf["slope_b_ci95"][0], cf["slope_b_ci95"][1],
    cf["r2"], cf["effective_B_A2"]))
print("symmetry R2 %.3f -> %.3f (+%.3f)" % (
    sc["r2_lnN_on_subset"], sc["r2_lnN_times_order_on_subset"], sc["r2_increase"]))
