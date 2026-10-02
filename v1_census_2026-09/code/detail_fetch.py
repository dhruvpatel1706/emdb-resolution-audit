"""Fetch the per-entry "consistency layer" from the EMDB native JSON API for a stratified sample.

For each sampled entry we record how the reported resolution is substantiated: res_type ("BY AUTHOR"
vs other), resolution_method ("FSC 0.143" vs "OTHER" vs ...), the number of particle images, the
applied symmetry, and the map pixel spacing (for a Nyquist sanity check). We stratify the sample
around the February 2022 half-map deposition mandate so the before/after comparison has power.
Cached to CSV as the artifact.
"""
import os, sys, time, random
import requests
import pandas as pd

HERE = os.path.dirname(__file__)
RAW = os.path.join(HERE, "..", "data", "raw")
CENSUS = os.path.join(RAW, "emdb_census.csv")
MANDATE = "20220201"
N_PER_SIDE = int(sys.argv[1]) if len(sys.argv) > 1 else 800
OUT = os.path.join(RAW, "emdb_detail.csv" if N_PER_SIDE >= 200 else "emdb_detail_probe.csv")


def fetch_one(emd, tries=4):
    for k in range(tries):
        try:
            r = requests.get(f"https://www.ebi.ac.uk/emdb/api/entry/{emd}",
                             headers={"User-Agent": "emdb-resolution-audit/1.0 (research metadata audit; contact via repository)",
                                      "Accept": "application/json"}, timeout=40)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except Exception:
            if k == tries - 1:
                return None
            time.sleep(1.5 * (k + 1))


def parse(j):
    try:
        sd = j["structure_determination_list"]["structure_determination"][0]
        method = sd.get("method")
        ip = sd.get("image_processing", [{}])[0]
        fr = ip.get("final_reconstruction", {}) or {}
        res = fr.get("resolution", {}) or {}
        sym = (fr.get("applied_symmetry", {}) or {}).get("point_group")
        m = j.get("map", {}) or {}
        px = (m.get("pixel_spacing", {}) or {}).get("x")
        px = float(px["valueOf_"]) if isinstance(px, dict) else (float(px) if px else None)
        return {"method": method, "res_type": res.get("res_type"),
                "resolution_method": fr.get("resolution_method"),
                "n_images": fr.get("number_images_used"), "point_group": sym, "pixel_x": px}
    except Exception:
        return None


def main():
    df = pd.read_csv(CENSUS, dtype={"deposition": str})
    df = df[df["deposition"].notna() & df["resolution"].notna()].copy()
    pre = df[df["deposition"] < MANDATE]["emdb_id"].tolist()
    post = df[df["deposition"] >= MANDATE]["emdb_id"].tolist()
    random.seed(0)
    sample = random.sample(pre, min(N_PER_SIDE, len(pre))) + random.sample(post, min(N_PER_SIDE, len(post)))
    print(f"sampling {len(sample)} entries ({min(N_PER_SIDE,len(pre))} pre + {min(N_PER_SIDE,len(post))} post mandate)")
    rows = []
    for i, emd in enumerate(sample):
        j = fetch_one(emd)
        if j is None:
            continue
        p = parse(j)
        if p:
            p["emdb_id"] = emd
            rows.append(p)
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(sample)} ({len(rows)} ok)")
    out = pd.DataFrame(rows).merge(df[["emdb_id", "deposition", "resolution", "has_pdb"]], on="emdb_id")
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT}: {len(out)} rows")
    print("res_type counts:\n", out["res_type"].value_counts(dropna=False).to_string())
    print("resolution_method counts:\n", out["resolution_method"].value_counts(dropna=False).head(8).to_string())


if __name__ == "__main__":
    main()
