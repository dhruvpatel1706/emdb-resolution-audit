"""Figures for the EMDB resolution-reporting audit. Reads data/results/emdb.json. CPU-only."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
FIG = os.path.join(HERE, "..", "paper", "figures")
os.makedirs(FIG, exist_ok=True)
D = json.load(open(os.path.join(HERE, "..", "data", "results", "emdb.json")))
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})


def main():
    fig, ax = plt.subplots(2, 2, figsize=(9.6, 7.0))

    # A: resolution revolution + archive growth
    my = D["resolution_over_time"]["median_by_year"]
    yrs = np.array(sorted(int(y) for y in my))
    med = np.array([my[str(y)]["median_res"] for y in yrs])
    cnt = np.array([my[str(y)]["count"] if "count" in my[str(y)] else my[str(y)]["n"] for y in yrs])
    a2 = ax[0, 0].twinx()
    a2.bar(yrs, cnt, color="#dadaeb", width=0.8, zorder=0); a2.set_ylabel("maps deposited / year", color="#6a51a3")
    a2.grid(False)
    ax[0, 0].plot(yrs, med, "o-", color="#08519c", zorder=3, lw=1.6)
    ax[0, 0].set_zorder(a2.get_zorder() + 1); ax[0, 0].patch.set_visible(False)
    ax[0, 0].set_yscale("log"); ax[0, 0].set_ylabel("median resolution (Å)", color="#08519c")
    ax[0, 0].set_xlabel("deposition year"); ax[0, 0].set_title("A  the resolution revolution")

    # B: median resolution by method
    bm = {g: v for g, v in D["by_method"].items() if v["median_res"]}
    order = sorted(bm, key=lambda g: bm[g]["median_res"])
    vals = [bm[g]["median_res"] for g in order]
    ax[0, 1].barh(range(len(order)), vals, color="#2c7fb8")
    for i, g in enumerate(order):
        ax[0, 1].text(vals[i] + 0.4, i, f"n={bm[g]['n']:,}", va="center", fontsize=7.5)
    disp = {"crystallography": "electron (2D)\ncrystallography"}
    ax[0, 1].set_yticks(range(len(order))); ax[0, 1].set_yticklabels([disp.get(g, g) for g in order])
    ax[0, 1].set_xlabel("median resolution (Å)"); ax[0, 1].set_title("B  resolution by method")
    ax[0, 1].set_xlim(0, max(vals) * 1.25)

    # C: fitted-model (PDB) rate vs resolution
    pr = D["pdb_rate_by_resolution"]
    labels = list(pr.keys()); rates = [pr[k]["pdb_rate"] for k in labels]
    ax[1, 0].bar(range(len(labels)), rates, color="#41ab5d")
    ax[1, 0].set_xticks(range(len(labels))); ax[1, 0].set_xticklabels([l.replace("A", "Å") for l in labels], fontsize=8)
    ax[1, 0].set_ylabel("fraction with fitted model"); ax[1, 0].set_ylim(0, 1)
    ax[1, 0].set_xlabel("reported resolution"); ax[1, 0].set_title("C  fitted-model rate tracks resolution")

    # D: FSC 0.143 substantiation by era, with the Feb-2022 mandate marked
    bk = D["substantiation"]["fsc0143_by_year_bucket"]
    keys = list(bk.keys()); fr = [bk[k]["frac_fsc0143"] for k in keys]
    cols = ["#fdae6b", "#fd8d3c", "#e6550d", "#a63603"]
    ax[1, 1].bar(range(len(keys)), fr, color=cols)
    for i, k in enumerate(keys):
        ax[1, 1].text(i, fr[i] + 0.02, f"{fr[i]*100:.0f}%", ha="center", fontsize=8)
    ax[1, 1].axvline(2.5, ls="--", c="grey"); ax[1, 1].text(2.55, 0.15, "Feb-2022\nmandate", fontsize=7.5)
    # Mark that the bar immediately left of the mandate line (2018-2021) is already a high plateau,
    # so the line clearly falls on an already-near-universal level rather than on a step up.
    ax[1, 1].annotate("already 91%\npre-mandate", xy=(2, fr[2]), xytext=(0.85, 0.72),
                      fontsize=7, ha="left",
                      arrowprops=dict(arrowstyle="->", color="grey", lw=0.8))
    ax[1, 1].set_xticks(range(len(keys))); ax[1, 1].set_xticklabels(keys, fontsize=8, rotation=15)
    ax[1, 1].set_ylabel("fraction citing FSC 0.143"); ax[1, 1].set_ylim(0, 1.1)
    ax[1, 1].set_title("D  substantiation rose before the mandate")

    fig.tight_layout()
    out = os.path.join(FIG, "fig_emdb.png")
    fig.savefig(out, dpi=600, bbox_inches="tight")  # 600 d.p.i. per Acta D notes for authors
    print("wrote", out)


if __name__ == "__main__":
    main()
