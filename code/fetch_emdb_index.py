"""Fetch one table of EMDB entry metadata, including the EMDB validation-analysis (VA) resolution fields, from the
EMDB search API (Solr index behind https://www.ebi.ac.uk/emdb/search).

One request, all entries, a fixed list of fields, CSV. The query URL, retrieval time (UTC), HTTP status, byte count,
row count and SHA-256 of the response are written to data/raw/fetch_log.json next to the CSV, so the table can be
checked against what was analysed. EMDB archive data are released under CC0 (wwPDB usage policy).

    python3 code/fetch_emdb_index.py            # writes data/raw/emdb_index_<UTC date>.csv and data/raw/fetch_log.json
    python3 code/fetch_emdb_index.py citations  # writes data/raw/emdb_citations_<UTC date>.csv and fetch_log_citations.json
                                                # (entry id and primary citation, used to group entries by study)
"""
import hashlib, json, os, sys, urllib.parse
from datetime import datetime, timezone

import requests

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
RAW = os.path.join(ROOT, "data", "raw")
API = "https://www.ebi.ac.uk/emdb/api/search/"
UA = {"User-Agent": "emdb-resolution-audit/2.0 (research metadata audit)"}

FIELDS = [
    # identity and dates
    "emdb_id", "structure_determination_method", "microscope_type", "deposition_date", "header_release_date",
    "map_release_date", "update_date", "obsolete_date", "withdrawn_date", "deposition_site",
    # the depositor's resolution
    "resolution", "resolution_method", "image_reconstruction_resolution", "image_reconstruction_resolution_type",
    "image_reconstruction_resolution_method", "author_resolution", "author_resolution_fsc_0143_value",
    "author_resolution_fsc_05_value", "author_resolution_fsc_halfbit_value",
    # EMDB validation analysis: FSC recalculated from the deposited half-maps
    "calculated_resolution", "calculated_resolution_fsc_0143_value", "calculated_resolution_fsc_05_value",
    "calculated_resolution_fsc_halfbit_value", "calculated_resolution_fsc_threesig_value",
    "calculated_resolution_fsc_masked_0143_value", "calculated_resolution_fsc_masked_05_value",
    "calculated_resolution_fsc_corrected_0143_value", "calculated_resolution_fsc_corrected_05_value",
    "calculated_resolution_fsc_randomised_0143_value", "calculated_resolution_fsc_randomised_at_resolution_value",
    "calculated_resolution_fsc_overfit_zone_value", "calculated_resolution_fsc_feature_zone_value",
    "calculated_resolution_fsc_masking_area_ratio_value", "relion_mask_coverage_average_value",
    "map_model_resolution_fsc_05_value", "map_model_resolution_fsc_0143_value",
    # what was deposited and how it was made
    "half_map_filename", "map_pixel_spacing_x", "half_map_pixel_spacing_x", "image_reconstruction_point_group",
    "particle_selected_count", "total_number_images", "image_reconstruction_software", "final_angle_assignment_software",
    "final_three_d_classification_software", "software", "acceleration_voltage", "microscope_name", "fitted_pdbs",
    "xref_PDB", "molecular_weight_theoretical_value", "molecular_weight_theoretical_unit", "overall_molecular_weight",
]


CITATION_FIELDS = ["emdb_id", "primary_citation_title", "xref_PUBMED", "xref_DOI", "primary_citation_year"]


def main():
    os.makedirs(RAW, exist_ok=True)
    when = datetime.now(timezone.utc)
    cit = len(sys.argv) > 1 and sys.argv[1] == "citations"
    fields = CITATION_FIELDS if cit else FIELDS
    url = API + urllib.parse.quote("*:*") + "?" + urllib.parse.urlencode(
        {"fl": ",".join(fields), "wt": "csv", "rows": 200000})
    r = requests.get(url, headers=UA, timeout=900)
    body = r.content
    out = os.path.join(RAW, f"emdb_{'citations' if cit else 'index'}_{when:%Y-%m-%d}.csv")
    if os.path.exists(out):
        sys.exit(f"{out} exists; not overwriting")
    open(out, "wb").write(body)
    rows = body.count(b"\n") - 1
    header = body.split(b"\n", 1)[0].decode()
    log = {"url": url, "retrieved_utc": when.isoformat(timespec="seconds"), "http_status": r.status_code,
           "bytes": len(body), "rows_excluding_header": rows, "sha256": hashlib.sha256(body).hexdigest(),
           "fields_requested": fields, "header_returned": header.split(","), "file": os.path.relpath(out, ROOT),
           "licence": "CC0 1.0 (wwPDB usage policy, https://www.wwpdb.org/about/usage-policies)"}
    assert r.status_code == 200, r.status_code
    missing = [f for f in fields if f not in log["header_returned"]]
    log["fields_missing_from_response"] = missing
    json.dump(log, open(os.path.join(RAW, "fetch_log_citations.json" if cit else "fetch_log.json"), "w"), indent=1)
    print(json.dumps({k: log[k] for k in ("http_status", "bytes", "rows_excluding_header", "sha256", "file",
                                           "fields_missing_from_response")}, indent=1))


if __name__ == "__main__":
    main()
