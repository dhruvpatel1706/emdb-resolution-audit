# Stated resolution versus EMDB's FSC recalculation from deposited half-maps

Code and data for a comparison of the resolution stated for each EMDB entry with the resolution EMDB recalculates from
the deposited half-maps: the 0.143 crossing of the Fourier shell correlation without a mask, with EMDB's automatic
mask, and corrected by phase randomization. The manuscript is in preparation and is not included yet. `paper/` holds
the numbers, tables, figures and Supplementary Table S1 that the code writes for it.

Python 3.14 with the packages pinned in `requirements.txt`. EMDB data are CC0 (wwPDB usage policy). The code is under
the MIT license (`LICENSE`). `CHECKSUMS.sha256` lists the SHA-256 of every file shipped for this analysis; check it with
`shasum -a 256 -c CHECKSUMS.sha256`.

Run order (each step writes what the next one reads):

1. `python3 code/fetch_emdb_index.py` and `python3 code/fetch_emdb_index.py citations`: one request each to the EMDB
   search API; writes `data/raw/emdb_index_<date>.csv`, `data/raw/emdb_citations_<date>.csv` and fetch logs with URL,
   time, size and SHA-256. The analysed files are dated 2026-10-01.
2. `python3 code/check_index_against_va.py data/raw/emdb_index_2026-10-01.csv`: downloads the per-entry validation
   records of 60 entries into `data/raw/va_sample/` and compares them with the index
   (`data/results/index_vs_va_check.json`). The 60 records used (about 80 MB) are not stored here;
   `data/raw/va_sample/MANIFEST.json` gives the URL, retrieval time, size and SHA-256 of each copy used. EMDB can update
   a record, so a new download can differ from the copy that was analysed.
3. `python3 code/check_own_curve_field.py data/raw/emdb_index_2026-10-01.csv`: for the same 60 entries, checks that the
   0.143 crossing of the depositor's FSC curve read from the record equals the index field
   `author_resolution_fsc_0143_value` (`data/results/own_curve_field_check.json`).
4. `python3 code/analyze_v2.py data/raw/emdb_index_2026-10-01.csv`: selection and comparison
   (`data/results/analysis_v2.json`, `data/results/tail_entries.csv`).
5. `python3 code/robustness_v2.py data/raw/emdb_index_2026-10-01.csv data/raw/emdb_citations_2026-10-01.csv`:
   resolution band by method and grouping by study (`data/results/robustness_v2.json`).
6. `python3 code/sample_unindexed_va.py data/raw/emdb_index_2026-10-01.csv`: draws fixed-seed random samples of 500
   entries released before 2024 and 150 released later among those with half-maps but no indexed recalculation,
   downloads each validation record and keeps only the 0.143 crossings with the record's URL, time, size and SHA-256
   (`data/results/unindexed_sample.jsonl`). Resumable; failed downloads are retried on a rerun. Four records
   returned HTTP 500 on every attempt on 1 October 2026 and are reported as failed.
7. `python3 code/analyze_unindexed_sample.py`: the same measures for the samples and the combined estimate
   (`data/results/unindexed_sample_analysis.json`).
8. `python3 code/make_paper_inputs.py data/raw/emdb_index_2026-10-01.csv`: every number (`paper/numbers.tex`),
   table, figure and `paper/Supplementary_Table_S1.csv`, recomputed from unrounded values and checked against steps 3-7.
9. `python3 code/check_paper.py`: regenerates step 8, rebuilds the PDF and checks numbers, style and placeholders. It
   needs the manuscript source `paper/main.tex`, which is not in this repository yet.

Supplementary Table S1 columns: EMDB id; method; release year; stated resolution; 0.143 crossing of the depositor's FSC
curve (blank if none deposited); EMDB unmasked, masked and corrected 0.143 crossings (all in angstrom); corrected /
stated.

## Earlier version

The September 2026 census of resolution metadata and criterion audit, a different analysis, is kept unchanged in
`v1_census_2026-09/` and at the tag `v1-actad-2026-09-21`. Its README describes that version as it was then.
