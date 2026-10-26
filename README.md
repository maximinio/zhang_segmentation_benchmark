# Zhang segmentation benchmark

Source, segmentation outputs and native RCTD results for the benchmark of RNA
assignment and segmentation geometry in deposited Visium HD profiles of
fibrostenotic Crohn's disease. Authors: Luca Massimino and Federica Ungaro.

## Contents

- `scripts/`: the 20 selected Python/R analytical and final-figure sources.
- `config/` and `protocols/`: exact configurations and scientific protocols.
- `results/segmentations/`: 40 full-section label TIFFs, original-bin label arrays
  and metrics (20 configurations on each of CD1 and CD2).
- `results/rctd/`: receipts for 40 scopes: 38 nonempty native full-mode major-cell
  fits and two explicitly zero-region scopes. The full native RDS containers are
  distributed as release assets, not omitted from the project.
- `results/tables/`: 12 summary and underlying-figure CSVs.
- `figures/`: the three main and 16 supplementary paper figures.
- `metadata/`: file checksums, selection and release-asset inventory.

## Full data

The native RCTD containers and the original counts-only Seurat input are being
published in the project's data release. `metadata/release_assets.json` lists
their logical paths, asset names, sizes and SHA-256 values. The verified release
URL and download URLs will be recorded after the upload completes.

No Git LFS subscription, paid storage or paid download service is required.
Release assets are part of this project even though they do not appear as large
binary blobs in Git history. They are not duplicated in the working checkout
after their uploaded server-side checksums have been verified.

## Scientific provenance and interpretation

Spatial inputs: Zhang et al., *Elucidating a Myofibroblast-dominated Fibrotic Niche
in Crohn's Disease-associated Fibrostenosis Through High-resolution Spatial
Transcriptomics*, GEO [GSE302569](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302569),
PMID [41344440](https://pubmed.ncbi.nlm.nih.gov/41344440/),
DOI [10.1016/j.jcmgh.2025.101701](https://doi.org/10.1016/j.jcmgh.2025.101701).
The native H&E inputs and 8-micrometre bin-level counts come from that deposited
dataset; the histology panels are not copied panels from the published article.

RCTD uses the published Nie et al. (2023) scIBD reference with original integer
`RNAraw` counts and author `major_cluster` labels, as recorded in the protocols
and native bundle provenance. All native values, filtering and explicit excluded
regions are preserved. No coefficient clipping or substitute scoring was used.

RNA recovery is not segmentation accuracy. Dominant coefficients and ranked
relative contributions are not calibrated posterior probabilities or a validated
formal uncertainty score. Under the available deposited bin-level constraints,
no tested configuration was established as a reliable single-cell reconstruction.
The separately generated single-method polygon Seurat was withdrawn and is not
included in this release.

## Reproduction and rights

The sources and configurations document the actual analytical methods. Paths in
the recorded configurations are execution-host paths and must be adapted before
running elsewhere; install the pinned runtimes and obtain the source/reference
inputs identified in the protocols. Publication does not imply that all methods
can be run in a single Python environment. No analysis was recomputed for release.

Third-party datasets, histology and dependency code retain their respective
authors' rights and upstream terms. This repository does not relicense those
materials or claim that the article's license automatically applies to deposited
GEO files. GEO's [data disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)
describes the distinction. No blanket license for all third-party material is
asserted. A project-wide software license has not yet been selected.
