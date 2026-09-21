"""Rosenthal-Henderson form in pooled EMDB metadata (descriptive, not mechanistic).

Fits the Rosenthal & Henderson (2003) ResLog relationship
    1/res^2 = a + b * ln(N_images)
across single-particle entries in the 1,600-entry detail sample, as a
POPULATION regression over heterogeneous deposited metadata, NOT a
per-reconstruction B-factor measurement.

Reported:
  (1) slope b with a nonparametric bootstrap 95% CI, R^2, implied
      "effective B" = 4/b (Angstrom^2) framed explicitly as an aggregate
      regression coefficient over a heterogeneous population.
  (2) symmetry-weighted control: replacing ln(N) by ln(N * point_group_order)
      on the point-group subset, showing R^2 rises (positive control that the
      theoretically predicted variable enters with the right sign).
  (3) R^2 reported explicitly to foreground that particle count explains
      only a minority of resolution variance archive-wide.

Deterministic: fixed RNG seed. CPU only, numpy only.
Reads data/raw/emdb_detail.csv, writes data/results/bfactor_law.json.
"""
import csv
import json
import math
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DETAIL = os.path.join(ROOT, "data", "raw", "emdb_detail.csv")
OUT = os.path.join(ROOT, "data", "results", "bfactor_law.json")

SEED = 20260622
N_BOOT = 10000


def num(x):
    try:
        v = float(x)
        return v
    except (TypeError, ValueError):
        return None


def point_group_order(pg):
    """Order of a Schoenflies point group symbol as used in cryo-EM.

    Returns None if unrecognised. Conventions:
      C1 -> 1, Cn -> n, Dn -> 2n, T -> 12, O -> 24, I -> 60.
    """
    if pg is None:
        return None
    pg = pg.strip().upper()
    if not pg:
        return None
    if pg in ("T", "TET"):
        return 12
    if pg in ("O", "OCT"):
        return 24
    if pg in ("I", "ICO"):
        return 60
    if pg[0] == "C" and pg[1:].isdigit():
        return int(pg[1:])
    if pg[0] == "D" and pg[1:].isdigit():
        return 2 * int(pg[1:])
    return None


def ols(x, y):
    """Ordinary least squares y = a + b x. Returns (a, b, r2)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    X = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b = coef
    yhat = a + b * x
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot
    return float(a), float(b), float(r2)


def bootstrap_slopes(x, y, n_boot, seed):
    """Paired bootstrap of both OLS slopes:
      b_fwd: y = a + b x        (predictive: 1/res^2 ~ ln N)
      b_rev: x = c + d y, d=B/2 (Rosenthal-Henderson ResLog: ln N ~ 1/res^2)
    Returns CIs for b_fwd and for the implied effective B = 2*d.
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    fwd = np.empty(n_boot)
    effB = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        _, b, _ = ols(x[idx], y[idx])
        _, d, _ = ols(y[idx], x[idx])
        fwd[i] = b
        effB[i] = 2.0 * d
    flo, fhi = np.percentile(fwd, [2.5, 97.5])
    blo, bhi = np.percentile(effB, [2.5, 97.5])
    return (float(flo), float(fhi), float(np.mean(fwd)), float(np.std(fwd)),
            float(blo), float(bhi), float(np.mean(effB)), float(np.std(effB)))


def main():
    rows = list(csv.DictReader(open(DETAIL)))

    # Single-particle entries with usable n_images and resolution.
    sp = [r for r in rows if r["method"] == "singleParticle"]
    n_sp_total = len(sp)
    n_sp_missing_nimg = sum(1 for r in sp if num(r["n_images"]) is None)

    rec = []
    for r in sp:
        n_img = num(r["n_images"])
        res = num(r["resolution"])
        if n_img is None or n_img <= 0 or res is None or res <= 0:
            continue
        rec.append({
            "n_images": n_img,
            "res": res,
            "point_group": r["point_group"].strip(),
        })
    n_fit = len(rec)

    # Reconciliation: the 1,485 figure counts ALL methods with usable
    # n_images and res; the defensible single-particle-only subset is 1,403.
    all_methods_fit = sum(
        1 for r in rows
        if num(r["n_images"]) and num(r["n_images"]) > 0
        and num(r["resolution"]) and num(r["resolution"]) > 0
    )

    # (1) Predictive fit: 1/res^2 ~ a + b ln(N) (forward; slope sign + R^2).
    lnN = [math.log(d["n_images"]) for d in rec]
    inv_res2 = [1.0 / (d["res"] ** 2) for d in rec]
    a, b, r2 = ols(lnN, inv_res2)
    # Rosenthal-Henderson ResLog orientation: ln(N) = c + d * (1/res^2),
    # with d = B/2 (Rosenthal & Henderson 2003). Effective B = 2*d, units A^2.
    c_rev, d_rev, r2_rev = ols(inv_res2, lnN)
    eff_B = 2.0 * d_rev
    (flo, fhi, fmean, fsd,
     eff_B_lo, eff_B_hi, eff_B_mean, eff_B_sd) = bootstrap_slopes(
        lnN, inv_res2, N_BOOT, SEED)
    blo, bhi, bmean, bsd = flo, fhi, fmean, fsd
    pearson_r = float(np.corrcoef(lnN, inv_res2)[0, 1])

    # (2) Symmetry-weighted control on the point-group subset.
    sub = []
    for d in rec:
        o = point_group_order(d["point_group"])
        if o is None or o <= 0:
            continue
        sub.append((d, o))
    n_pg = len(sub)

    lnN_pg = [math.log(d["n_images"]) for d, o in sub]
    inv_res2_pg = [1.0 / (d["res"] ** 2) for d, o in sub]
    lnNS_pg = [math.log(d["n_images"] * o) for d, o in sub]

    a_pg, b_pg, r2_pg = ols(lnN_pg, inv_res2_pg)           # ln(N) on subset
    a_sym, b_sym, r2_sym = ols(lnNS_pg, inv_res2_pg)       # ln(N*order) on subset

    out = {
        "_meta": {
            "script": "code/bfactor_law.py",
            "seed": SEED,
            "n_boot": N_BOOT,
            "model": "forward fit 1/res^2 = a + b*ln(n_images); effective B = 2*d, where d is the reverse (ResLog) OLS slope of ln(n_images) on 1/res^2 (Rosenthal-Henderson d = B/2)",
            "framing": (
                "POPULATION regression over heterogeneous deposited metadata; "
                "effective B is an aggregate regression coefficient, NOT a "
                "per-reconstruction B-factor."
            ),
        },
        "counts": {
            "n_singleParticle_total": n_sp_total,
            "n_singleParticle_missing_nimages": n_sp_missing_nimg,
            "n_fit_singleParticle": n_fit,
            "n_fit_all_methods": all_methods_fit,
            "n_pointgroup_subset": n_pg,
            "reconciliation_note": (
                "1485 = all methods with usable n_images and resolution; "
                "1403 = single-particle-only subset used for the fit."
            ),
        },
        "core_fit": {
            "_predictive_orientation": "1/res^2 = a + b*ln(N)",
            "intercept_a": a,
            "slope_b": b,
            "slope_b_bootstrap_mean": bmean,
            "slope_b_bootstrap_sd": bsd,
            "slope_b_ci95": [blo, bhi],
            "r2": r2,
            "pearson_r": pearson_r,
            "slope_positive_and_significant": bool(blo > 0),
            "_reslog_orientation": "ln(N) = c + d*(1/res^2), d = B/2",
            "reslog_slope_d": d_rev,
            "reslog_r2": r2_rev,
            "effective_B_A2": eff_B,
            "effective_B_A2_bootstrap_mean": eff_B_mean,
            "effective_B_A2_bootstrap_sd": eff_B_sd,
            "effective_B_A2_ci95": [eff_B_lo, eff_B_hi],
        },
        "symmetry_control": {
            "n": n_pg,
            "r2_lnN_on_subset": r2_pg,
            "r2_lnN_times_order_on_subset": r2_sym,
            "r2_increase": r2_sym - r2_pg,
            "slope_b_lnN_subset": b_pg,
            "slope_b_sym_subset": b_sym,
            "sym_slope_positive": bool(b_sym > 0),
        },
    }

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
        f.write("\n")

    # Print key numbers.
    print("wrote", OUT)
    print("--- reconciliation ---")
    print("SP total %d, SP missing n_images %d, SP fit %d, all-methods fit %d, pg subset %d"
          % (n_sp_total, n_sp_missing_nimg, n_fit, all_methods_fit, n_pg))
    print("--- predictive fit (1/res^2 ~ a + b ln N), single particle n=%d ---" % n_fit)
    print("a=%.6g  b=%.6g  R2=%.4f  pearson_r=%.4f" % (a, b, r2, pearson_r))
    print("bootstrap b 95%% CI [%.6g, %.6g]  (mean %.6g, sd %.6g)" % (blo, bhi, bmean, bsd))
    print("--- ResLog orientation (ln N ~ 1/res^2, slope=B/2) ---")
    print("d=%.4f  effective B = 2d = %.2f A^2  bootstrap CI [%.2f, %.2f]"
          % (d_rev, eff_B, eff_B_lo, eff_B_hi))
    print("--- symmetry control on point-group subset n=%d ---" % n_pg)
    print("R2 ln(N)            = %.4f" % r2_pg)
    print("R2 ln(N*order)      = %.4f  (increase %.4f)" % (r2_sym, r2_sym - r2_pg))
    print("sym slope b=%.6g (positive=%s)" % (b_sym, b_sym > 0))


if __name__ == "__main__":
    main()
