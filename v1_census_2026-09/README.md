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

The cached tables are included, so the analysis and figures can be reproduced offline. To re-pull the source tables (network required), run the two fetch steps first; counts will then reflect archive growth since the June 2026 snapshot.

```
pip install -r requirements.txt
python code/build_table.py           # -> data/raw/emdb_census.csv (EBI Search census)   [network, optional]
python code/detail_fetch.py          # -> data/raw/emdb_detail.csv (native API sample)    [network, optional]
python code/analyze.py               # -> data/results/emdb.json (census, substantiation, Nyquist findings)
python code/bfactor_law.py           # -> data/results/bfactor_law.json (Rosenthal-Henderson fit)
python code/nyquist_sensitivity.py   # -> data/results/nyquist_sensitivity.json (tolerance-band sweep)
mkdir -p paper/figures
python code/figures.py               # -> paper/figures/fig_emdb.png
python code/fig_bfactor.py           # -> paper/figures/fig_bfactor.png
python code/fig_residual.py          # -> paper/figures/fig_residual.png
```

`code/check_paper_numbers.py` re-derives every statistic quoted in the manuscript from `data/results/` and compares it with the manuscript source (path given by the `EMDB_TEX` environment variable); the manuscript source is not included in this repository while the paper is under review, so that script is provided for the record.

## Data sources

- EBI Search REST API (EMDB census).
- EMDB native entry JSON API (per-entry detail).

Both are public, require no authentication, and the code only parses and tabulates the returned metadata; no reported value is altered.

## License

- Code: MIT (see `LICENSE`).
- Paper text and figures: CC BY 4.0.
- The parsed EMDB metadata tables are derived from public EMDB/EBI data; please also credit EMDB/EBI as the upstream source.

## Citation

If you use this code or the parsed tables, please cite the paper (under review; the reference will be added here when it is published). The repository URL is https://github.com/dhruvpatel1706/emdb-resolution-audit; a citable archive DOI for the same package has not yet been minted.

## Manuscript

The manuscript is under review at Acta Crystallographica Section D (submitted 30 July 2026). This repository holds the code and the cached data it reads; the paper itself is not included while the review is open.
