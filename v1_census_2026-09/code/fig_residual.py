"""Forest plot of the post-mandate non-FSC-0.143 residual by modality and by year.
All values are read directly from data/results/emdb.json (no new computation).
Outputs paper/figures/fig_residual.png.
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "data", "results", "emdb.json")
OUT = os.path.join(ROOT, "paper", "figures", "fig_residual.png")

d = json.load(open(RES))
r = d["substantiation"]["residual_non_fsc0143_post_mandate"]

# Overall post-mandate residual (well powered reference line).
overall_pct = r["frac_post_non_fsc0143"] * 100.0
overall_lo, overall_hi = [x * 100.0 for x in r["frac_post_non_fsc0143_wilson95"]]

def rows(block, order, label_map):
    out = []
    for key in order:
        v = block[key]
        out.append((
            label_map[key],
            v["n_non_fsc0143"], v["n"],
            v["frac_non_fsc0143"] * 100.0,
            v["wilson95_lo"] * 100.0,
            v["wilson95_hi"] * 100.0,
        ))
    return out

method_order = ["single particle", "helical", "subtomogram avg", "tomography"]
method_labels = {
    "single particle": "single particle",
    "helical": "helical",
    "subtomogram avg": "subtomogram avg",
    "tomography": "tomography",
}
year_order = ["2022", "2023", "2024", "2025"]
year_labels = {y: y for y in year_order}

method_rows = rows(r["non_fsc0143_by_method"], method_order, method_labels)
year_rows = rows(r["non_fsc0143_by_year"], year_order, year_labels)

fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), sharex=True)

def draw(ax, data, title):
    ax.axvspan(overall_lo, overall_hi, color="0.85", zorder=0)
    ax.axvline(overall_pct, color="0.55", lw=1.0, ls="--", zorder=1)
    ys = list(range(len(data)))[::-1]
    labels = []
    for y, (lab, k, n, pct, lo, hi) in zip(ys, data):
        well_powered = (hi - lo) <= 6.0
        color = "#1f4e79" if well_powered else "#b35900"
        ax.errorbar(pct, y, xerr=[[pct - lo], [hi - pct]], fmt="o",
                    color=color, ecolor=color, elinewidth=1.4,
                    capsize=3, ms=5, zorder=3)
        labels.append("%s  (%d/%d)" % (lab, k, n))
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_ylim(-0.6, len(data) - 0.4)
    ax.set_title(title, fontsize=9.5)
    ax.tick_params(axis="x", labelsize=8.5)
    ax.grid(axis="x", color="0.92", zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

draw(axes[0], method_rows, "by modality")
draw(axes[1], year_rows, "by deposition year")
axes[0].set_xlabel("post-mandate non-FSC-0.143 share (%)", fontsize=8.5)
axes[1].set_xlabel("post-mandate non-FSC-0.143 share (%)", fontsize=8.5)
axes[0].set_xlim(-3, 103)

fig.tight_layout()
fig.savefig(OUT, dpi=600, bbox_inches="tight")  # 600 d.p.i. per Acta D notes for authors
print("wrote", OUT)
print("overall %.2f%% [%.1f, %.1f]" % (overall_pct, overall_lo, overall_hi))
for lab, k, n, pct, lo, hi in method_rows + year_rows:
    print("%-18s %2d/%-3d %5.1f%% [%.1f, %.1f]" % (lab, k, n, pct, lo, hi))
