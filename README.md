# DUSP Family Dysregulation in Breast Cancer

Analysis code for a nine-gene DSP expression signature for breast
cancer. The manuscript's central finding is biological: a consistent,
direction-specific expression pattern across the gene set that is discovered
in TCGA-BRCA and reproduced across independent RNA-seq and microarray
cohorts and, in part, other cancer types.


**Manuscript:** Machine Learning Based Identification of Consistently Dysregulated DUSP Family Phosphatases in Breast Cancer
**Status:** Submitted to Computational Biology and Chemistry, [September 2026]

**Authors:** Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran

## Overview

The pipeline first establishes a 9-gene directional expression pattern from
TCGA-BRCA data and tests its statistical significance. A Bayesian-optimised
XGBoost model and a biologically transparent linear directional score (sum of
up-regulated minus down-regulated genes) are both used as quantitative
readouts of this pattern — not as proposed clinical tools — so the finding
can be evaluated for internal consistency, reproduced in independent cohorts,
and checked for specificity to breast cancer versus other tumour types. The
final scripts characterise the biological structure behind the finding: how
the signature genes co-vary with one another and which pathways and Hallmark
gene sets they are associated with.

All stochastic steps use a fixed seed (`42`) for reproducibility.

## Repository structure

```
.
├── 1_training_pipeline.R              # Establishes the directional pattern + readout model
├── 2_permutation_test.R               # Significance of the pattern's discriminative AUC
├── 3_directional_classifier.R         # Linear directional score vs. XGBoost readout
├── 4_pattern_conservation.R           # Pan-cancer conservation of the pattern
├── 5_external_validation.R            # Reproducibility in external GEO RNA-seq cohorts
├── 6_geo_direction_concordance.R      # Per-gene direction check (GEO RNA-seq)
├── 7_crossplatform_validation.R       # Reproducibility in external GEO microarray cohorts
├── 8_microarray_direction_concordance.R  # Per-gene direction check (microarray)
├── 9_intergene_correlation.R          # Gene-gene correlation structure of the pattern
├── 10_enrichr_dotplot.R               # Pathway enrichment behind the pattern
├── 11_hallmark_membership.R           # MSigDB Hallmark membership of the gene set
├── logs/                              # Console logs + session_info per run
├── models/                            # Saved model / signature objects (.rds)
├── data/                              # NOT included — see "Required inputs"
└── output/                            # NOT included — created by the scripts
```

## Pipeline order and dependencies

Run in numeric order; later scripts depend on objects produced by earlier
ones (mainly Script 1).

| # | Script | Depends on | Produces |
|---|--------|-----------|----------|
| 1 | `1_training_pipeline.R` | — (entry point) | directional gene pattern, readout model, directional-score table |
| 2 | `2_permutation_test.R` | 1 | permutation significance of the pattern's discriminative AUC |
| 3 | `3_directional_classifier.R` | 1 | behaviour of the linear directional score vs. the XGBoost readout |
| 4 | `4_pattern_conservation.R` | 1 | evidence for/against the pattern's conservation across cancer types |
| 5 | `5_external_validation.R` | 1 | reproducibility of the pattern in external GEO RNA-seq cohorts |
| 6 | `6_geo_direction_concordance.R` | 1, 5 | per-gene direction concordance (GEO RNA-seq) |
| 7 | `7_crossplatform_validation.R` | 1 | reproducibility of the pattern in external GEO microarray cohorts |
| 8 | `8_microarray_direction_concordance.R` | 1, 7 | per-gene direction concordance (microarray) |
| 9 | `9_intergene_correlation.R` | 1 | inter-gene correlation heatmap |
| 10 | `10_enrichr_dotplot.R` | 1 (gene list); Enrichr TSV exports | pathway enrichment dot plot |
| 11 | `11_hallmark_membership.R` | 1 (gene list) | Hallmark gene-set membership table |

## Required inputs (not included in this repository)

The scripts expect a `data/` folder alongside them, not tracked in git due to
size:

```
data/
├── TCGA-BRCA.star_tpm.tsv.gz                 # TCGA GDC star_tpm, log2(TPM+1)
├── BRCA_gene_list.txt                        # one HGNC gene symbol per line
├── TCGA-[LUAD,KIRC,COAD,GBM,PRAD,SKCM].star_tpm.tsv.gz   # Script 4
├── validation/
│   ├── GSE233242_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz  # Script 5, 6
│   └── GSE58135_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz   # Script 5, 6
│   # GSE70947, GSE15852 microarray series (Scripts 7, 8)
└── enrichr/
    ├── All9_Reactome_Pathways_2024.txt        # Script 10
    ├── All9_BioPlanet_2019.txt
```

Scripts also read/write an `output/tables/` and `output/models/` folder
(intermediate CSVs and `.rds` objects); the `models/` folder in this
repository holds the final locked artifacts (the pattern's readout model and
directional-score table) referenced across scripts:

- `final_xgboost_model_Bayesian.rds`
- `final_preProc_scaler.rds`
- `final_stable_biomarker_signature.rds`
- `Final_Directional_Signature.rds`

## Reproducing a run

1. Install R and the packages each script's header/`library()` calls list
   (not pinned here — add a `renv.lock` or `sessionInfo()` snapshot if you
   need exact package versions; see `logs/*_session_info.txt` for the
   versions originally used).
2. Place the required inputs under `data/` as above.
3. Run scripts in numeric order from the repository root. Each script
   writes its own log to `logs/` and, where applicable, a
   `logs/NN_session_info.txt` capturing the R session used for that run.
4. Script 2 supports crash recovery — set `RECOVERY_MODE = TRUE` inside it to
   resume post-processing from the saved raw permutation results instead of
   re-running the full 1000-permutation loop.

## Notes

- This repository contains only the final scripts used to produce the
  manuscript's results — not the full set of intermediate/exploratory
  analysis scripts run during the project. The `1`–`11` prefixes are simply
  this repo's run order for the curated set; they are not tied to the
  manuscript's figure or table numbering, and are not expected to be
  contiguous with any script numbering used internally during exploratory
  work.
- Filenames are not zero-padded (`1_...` rather than `01_...`), so a plain
  alphabetical listing sorts as 1, 10, 11, 2, 3, ... rather than pipeline
  order; use the table above, not directory sort order, to find the next
  script.
-  internal console title and saved output filenames use an prefix, left over from the script's
  label in the broader internal project. It does not correspond to any
  script in this repository — cosmetic only, outputs are correct.
- `logs/09_InterGene_Corr_Run_20260812_1635.log` has no matching
  `09_session_info.txt` (present for every other script's log).
- Figure title and figure legends are hardcoded, it can vary based on the results  

## Citation

If you use this code, please cite the accompanying manuscript (details to be
added upon publication).
