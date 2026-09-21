# Resolution Metadata in the EMDB Cryo-EM Archive: An Indexed Census of 46,900 Entries and a Stratified 1,600-Entry Criterion Audit

Code and data for the paper *"Resolution Metadata in the EMDB Cryo-EM Archive: An Indexed Census of 46,900 Entries and a Stratified 1,600-Entry Criterion Audit"* by Dhruv Patel (Independent Researcher).

A CPU-only, public-API metadata study with two evidence layers: a 46,900-entry EBI-Search-indexed census for historical, method and fitted-model summaries, and a stratified 1,600-entry sample for criterion provenance, particle-count and Nyquist analyses. It does **not** recompute any resolution from density and the sample analyses are not census-wide results. Deterministic; runs on a laptop CPU from cached tables.

The parsed EMDB tables shipped here (a full-census table plus a stratified per-entry sample) are a reusable dataset in their own right: a one-row-per-entry snapshot of EMDB resolution-reporting metadata that others can reanalyze without re-crawling the APIs.

## Repository layout

```
code/                 analysis, figure and number-guard scripts (see RESULTS.md for what each produces)
data/raw/             cached EMDB census and per-entry sample tables (EBI Search and EMDB APIs, June 2026)
data/results/         result JSON files that every number in the paper is bound to
RESULTS.md            per-result provenance
CHECKSUMS.sha256      SHA-256 of every shipped file; verify with: shasum -a 256 -c CHECKSUMS.sha256
requirements.txt      pinned Python environment
```

## Requirements

Python 3 with `requests`, `numpy`, `pandas`, `scipy`, `matplotlib`. CPU-only; sources are public REST/JSON APIs with no authentication and no paid APIs.

## Reproduce

The cached CSVs are included, so the analysis and figures can be reproduced offline. To re-pull the source tables (network required), run the fetch steps first.

```
python code/build_table.py    # -> data/raw/emdb_census.csv (EBI Search census)  [network]
python code/detail_fetch.py   # -> data/raw/emdb_detail.csv (native API sample)   [network]
python code/analyze.py        # -> data/results/emdb.json (all findings)
python code/figures.py        # -> paper/figures/fig_emdb.png
python3 code/check_paper_numbers.py   # guard on paper/acta_d_revision/ (0 failures expected)
cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

Counts reported in the paper reflect the cached snapshot; re-pulling later will reflect archive growth since the snapshot date.

## Data sources

- EBI Search REST API (EMDB census).
- EMDB native entry JSON API (per-entry detail).

Both are public, require no authentication, and the code only parses and tabulates the returned metadata; no reported value is altered.

## License

- Code: MIT (see `LICENSE`).
- Paper text and figures: CC BY 4.0.
- The parsed EMDB metadata tables are derived from public EMDB/EBI data; please also credit EMDB/EBI as the upstream source.

## Citation

If you use this code or the parsed tables, please cite the paper. A Zenodo DOI and public repository URL have not yet been assigned. Replace `[Zenodo DOI: 10.5281/zenodo.XXXXXXX]` and `https://github.com/PLACEHOLDER/emdb-resolution-audit` only after the corresponding deposits exist; do not describe the materials as publicly archived before then.

## 1 September 2026: live Acta D source
The manuscript source that matches the staged `emdb_ActaD.pdf` is `paper/acta_d/emdb_ActaD_main.tex`
(mirrored from `SUBMISSION_READY_2026-07-27/_SUBMIT_PDFS/ACTA_D_SOURCE/`); `paper/main.tex` and
`paper/database_oup/` are older drafts. Per the Acta D Notes for Authors (read 1 Sep 2026): templates
are encouraged, not required; figures must be separate files at 600 d.p.i. or better, so the three
figure scripts now save at 600 d.p.i.; the abstract was shortened toward the requested 250 words and a
synopsis + keywords block added (`paper/acta_d/synopsis_keywords_plaintext.txt` for the form). The
scientific hold in `EMDB_ACTAD_HANDOFF_2026-07-29.md` (sample versus census scope) is unchanged.

## Manuscript

The manuscript is under review (Acta Crystallographica Section D, submitted 1 September 2026). This repository holds the code and the cached data it reads; the paper itself is not included while the review is open.
