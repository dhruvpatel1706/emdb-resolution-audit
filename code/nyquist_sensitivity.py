"""
Nyquist sensitivity analysis for the EMDB resolution audit (CPU-only, local data).

The main analysis (code/analyze.py) flags entries whose reported resolution sits below
twice the map pixel spacing (the Nyquist limit) and isolates a single borderline candidate,
EMD-44581, at ratio reported_resolution / (2 * pixel_x) = 0.973, that is, only a few percent
inside its sampling limit. The paper calls that "within plausible measurement uncertainty."

This script quantifies that claim. From the same 1,600-entry stratified sample (data/raw/
emdb_detail.csv) it computes the full distribution of the ratio r = resolution / (2 * pixel_x),
its spread (mean, std, quantiles, IQR, MAD), and, for a set of symmetric measurement-uncertainty
bands eps in {2%, 5%, 10%}, how many entries fall (a) within +/- eps of the exact Nyquist limit
(r in [1 - eps, 1 + eps]) and (b) on the sub-Nyquist side of it but still inside the band
(r in [1 - eps, 1)). The point of (b) is to show at which band the lone borderline sub-Nyquist
entry is absorbed: a 2% band does not reach it (it is 2.7% inside), a 5% or 10% band does, so the
"within plausible uncertainty" reading holds only once the band is at least ~3%.

Reads ONLY the local cached CSV. No network, no model. Fixed seed for full determinism even
though the computation is deterministic; output dict keys are written in fixed insertion order
and JSON is emitted with sort_keys for byte-stable artifacts.
"""
import os
import json
import numpy as np
import pandas as pd

SEED = 0
np.random.seed(SEED)

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "raw")
RESULTS = os.path.join(HERE, "..", "data", "results")
os.makedirs(RESULTS, exist_ok=True)

# Symmetric relative bands around the Nyquist limit, as fractions.
BANDS = [0.02, 0.05, 0.10]


def main():
    d = pd.read_csv(os.path.join(RAW, "emdb_detail.csv"), dtype={"deposition": str})
    n = d.dropna(subset=["resolution", "pixel_x"]).copy()
    n = n[n["pixel_x"] > 0]
    # Ratio of reported resolution to the Nyquist limit (2 * pixel spacing).
    # r < 1 is sub-Nyquist (a reported detail finer than the sampling can represent).
    n["ratio"] = n["resolution"] / (2.0 * n["pixel_x"])
    r = n["ratio"].to_numpy()

    # Distribution / spread of the ratio.
    q = np.percentile(r, [1, 5, 10, 25, 50, 75, 90, 95, 99])
    dist = {
        "n": int(len(r)),
        "mean": float(np.mean(r)),
        "std": float(np.std(r, ddof=1)),
        "min": float(np.min(r)),
        "max": float(np.max(r)),
        "median": float(np.median(r)),
        "iqr": float(q[5] - q[3]),
        "mad_about_median": float(np.median(np.abs(r - np.median(r)))),
        "p01": float(q[0]), "p05": float(q[1]), "p10": float(q[2]),
        "p25": float(q[3]), "p50": float(q[4]), "p75": float(q[5]),
        "p90": float(q[6]), "p95": float(q[7]), "p99": float(q[8]),
    }

    # Counts strictly below / at-or-above the Nyquist limit.
    n_sub = int(np.sum(r < 1.0))
    n_at_or_above = int(np.sum(r >= 1.0))

    # Per-band sensitivity. For each band eps:
    #   within_band: entries with r in [1 - eps, 1 + eps] (straddle the limit within eps)
    #   sub_within_band: entries with r in [1 - eps, 1) (below the limit but inside the band)
    #   sub_outside_band: entries with r < 1 - eps (below the limit by more than eps)
    per_band = {}
    for eps in BANDS:
        lo, hi = 1.0 - eps, 1.0 + eps
        within = (r >= lo) & (r <= hi)
        sub_within = (r >= lo) & (r < 1.0)
        sub_outside = r < lo
        key = f"{int(round(eps * 100))}pct"
        per_band[key] = {
            "eps": float(eps),
            "band_lo": float(lo),
            "band_hi": float(hi),
            "n_within_band": int(np.sum(within)),
            "frac_within_band": float(np.mean(within)),
            "n_sub_nyquist_within_band": int(np.sum(sub_within)),
            "n_sub_nyquist_outside_band": int(np.sum(sub_outside)),
        }

    # The borderline candidate the paper discusses, located explicitly so the artifact
    # records which entry each band does or does not absorb. Defined as the smallest-ratio
    # entry whose reported resolution does NOT equal its pixel spacing (i.e. not the two
    # pixel-into-resolution transcription artifacts at ratio ~ 0.5).
    res_eq_pixel = (np.abs(n["resolution"] - n["pixel_x"]) / n["pixel_x"]) < 0.001
    substantive = n[~res_eq_pixel]
    bl = substantive.loc[substantive["ratio"].idxmin()]
    bl_ratio = float(bl["ratio"])
    bl_pct_inside = float((1.0 - bl_ratio) * 100.0)  # how far inside the limit, in percent
    borderline = {
        "emdb_id": str(bl["emdb_id"]),
        "resolution": float(bl["resolution"]),
        "pixel_x": float(bl["pixel_x"]),
        "nyquist_limit": float(2.0 * bl["pixel_x"]),
        "ratio": bl_ratio,
        "pct_inside_nyquist": bl_pct_inside,
        "method": str(bl["method"]),
        "resolution_method": str(bl["resolution_method"]),
        "absorbed_by_band": {
            f"{int(round(eps * 100))}pct": bool(bl_ratio >= 1.0 - eps) for eps in BANDS
        },
        "smallest_band_absorbing": (
            min((eps for eps in BANDS if bl_ratio >= 1.0 - eps), default=None)
        ),
    }

    out = {
        "description": (
            "Sensitivity of the EMDB Nyquist consistency check to measurement-uncertainty "
            "bands around the limit r = reported_resolution / (2 * pixel_spacing) = 1."
        ),
        "source": "data/raw/emdb_detail.csv (1,600-entry stratified per-entry sample)",
        "seed": SEED,
        "ratio_distribution": dist,
        "n_sub_nyquist_strict": n_sub,
        "n_at_or_above_nyquist": n_at_or_above,
        "frac_sub_nyquist_strict": float(n_sub / len(r)),
        "bands": per_band,
        "borderline_entry": borderline,
    }

    with open(os.path.join(RESULTS, "nyquist_sensitivity.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"ratio: n={dist['n']} mean={dist['mean']:.4f} std={dist['std']:.4f} "
          f"median={dist['median']:.4f} IQR={dist['iqr']:.4f} MAD={dist['mad_about_median']:.4f}")
    print(f"min={dist['min']:.4f} p01={dist['p01']:.4f} p05={dist['p05']:.4f} p10={dist['p10']:.4f}")
    print(f"sub-Nyquist (r<1): {n_sub}/{dist['n']} = {out['frac_sub_nyquist_strict']*100:.3f}%")
    for k, v in per_band.items():
        print(f"  band {k}: within [{v['band_lo']:.2f},{v['band_hi']:.2f}] -> {v['n_within_band']} "
              f"({v['frac_within_band']*100:.2f}%); sub-Nyquist inside band "
              f"{v['n_sub_nyquist_within_band']}, sub-Nyquist outside band "
              f"{v['n_sub_nyquist_outside_band']}")
    b = borderline
    print(f"borderline {b['emdb_id']}: ratio={b['ratio']:.4f} ({b['pct_inside_nyquist']:.2f}% inside); "
          f"absorbed by {b['absorbed_by_band']}; smallest band {b['smallest_band_absorbing']}")
    print("wrote nyquist_sensitivity.json")


if __name__ == "__main__":
    main()
