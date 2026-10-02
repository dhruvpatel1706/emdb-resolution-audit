"""Pull the full EMDB metadata census from the EBI Search REST API and cache it to CSV.

One bulk field set per entry (method, reported resolution, deposition/release dates, whether a fitted
PDB model exists), paginated 100 at a time over all ~46,900 entries. This is the spine of the audit;
it downloads only tiny JSON and is fully reproducible. The cached CSV is the data artifact.
"""
import os, time
import requests
import pandas as pd

RAW = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
os.makedirs(RAW, exist_ok=True)
BASE = "https://www.ebi.ac.uk/ebisearch/ws/rest/emdb"
FIELDS = "method,resolution,depositiondate,mapreleasedate,headerreleasedate,PDB"
HEADERS = {"User-Agent": "emdb-resolution-audit/1.0 (research metadata audit; contact via repository)",
           "Accept": "application/json"}


def get(start, size=100, tries=5):
    for k in range(tries):
        try:
            r = requests.get(BASE, params={"query": "domain_source:emdb", "size": size, "start": start,
                                           "fields": FIELDS, "format": "json"}, headers=HEADERS, timeout=60)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


def first(f, key):
    v = f.get(key) or []
    return v[0] if v else None


def main():
    rows = []
    n = get(0)["hitCount"]
    print(f"hitCount={n}")
    for start in range(0, n, 100):
        j = get(start)
        for e in j["entries"]:
            f = e["fields"]
            res = first(f, "resolution")
            rows.append({
                "emdb_id": e["id"],
                "method": first(f, "method"),
                "resolution": float(res) if res not in (None, "") else None,
                "deposition": first(f, "depositiondate"),
                "map_release": first(f, "mapreleasedate"),
                "header_release": first(f, "headerreleasedate"),
                "has_pdb": bool(f.get("PDB")),
            })
        if start % 5000 == 0:
            print(f"  {start}/{n} ({len(rows)} rows)")
    df = pd.DataFrame(rows)
    out = os.path.join(RAW, "emdb_census.csv")
    df.to_csv(out, index=False)
    print(f"wrote {out}: {len(df)} rows, "
          f"{df['resolution'].notna().sum()} with resolution, {df['method'].notna().sum()} with method")


if __name__ == "__main__":
    main()
