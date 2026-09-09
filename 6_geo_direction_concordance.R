# =============================================================================
# Script      : 6_geo_direction_concordance.R
# Project     : DSP family dysregulation in breast cancer
#
# PURPOSE
#   Test whether each signature gene individually reproduces its
#   TCGA-BRCA-derived direction of dysregulation in two independent GEO
#   RNA-seq cohorts.
#
#   This is distinct from Script 02, which validated the composite
#   directional SCORE. Here each gene is evaluated separately, so a gene
#   whose contribution is masked within the composite can be identified.
#
# COHORTS
#   GSE58135   84 tumour / 56 normal   (Varley et al. 2014)
#   GSE233242  42 tumour / 42 normal   (Li et al. 2024)
#
# CLASSIFICATION (identical scheme to Script 13, for consistency)
#   Match          concordant direction + significant
#   Weak match     concordant direction + not significant  -> underpowered
#   Divergent      discordant direction + significant      -> real difference
#   Uninformative  discordant direction + not significant
#
#   Direction and significance are scored SEPARATELY. Collapsing them would
#   cause low-effect-size genes to fail automatically in smaller cohorts,
#   confounding statistical power with biological difference.
#
# POWER NOTE — read before interpreting
#   These cohorts are far smaller than TCGA-BRCA (n ~1,222). Genes with
#   modest effect sizes in the discovery cohort (e.g. STYXL1, log2FC +0.20;
#   PTP4A2, marginal) may not reach significance here. Such results are
#   Weak match, not Divergent, and must not be reported as failures to
#   replicate.
#
# OUTPUTS
#   Tables:  14_gene_concordance_full.csv    per gene x dataset
#            14_concordance_summary.csv      per-dataset class counts
#   Figures: Fig_14A_Concordance_Matrix      gene x dataset classification
#            Fig_14B_LogFC_Comparison        TCGA vs GEO effect sizes
#   Report:  14_GEO_CONCORDANCE_REPORT.txt
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026

# =============================================================================

rm(list = ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed            = 42L,
  q_threshold     = 0.05,     # BH-adjusted significance
  logfc_threshold = 0.3,      # matches the screening criterion in Methods 2.2
  min_group_n     = 10L       # minimum samples per group to attempt a test
)

# -----------------------------------------------------------------------------
# DATASET CONFIGURATION
#
# Expression: NCBI-generated normalised count matrices from GEO
#   (*_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz). Row identifiers are ENTREZ
#   GENE IDs, mapped to symbols by this script.
#
# Phenotype: retrieved from GEO via GEOquery. Group is assigned by keyword
#   matching against concatenated title, source, characteristics and
#   description fields.
#
# Keyword precedence: exclude > normal > tumour. Normal is checked before
#   tumour so labels containing both terms (e.g. "uninvolved tissue adjacent
#   to breast tumour") resolve correctly to Normal.
#
# This configuration is IDENTICAL to Script 02, guaranteeing that the score
# validation and this per-gene analysis operate on the same sample sets.
# -----------------------------------------------------------------------------
DATASETS <- list(
  list(
    gse_id           = "GSE233242",
    file_path        = file.path("data", "validation",
                                 "GSE233242_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
    normal_keywords  = c("normal"),
    tumor_keywords   = c("tumor"),
    exclude_keywords = character(0),
    expected_tumour  = 43L,
    expected_normal  = 43L
  ),
  list(
    gse_id           = "GSE58135",
    file_path        = file.path("data", "validation",
                                 "GSE58135_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
    normal_keywords  = c("uninvolved", "mammoplasty", "adjacent",
                         "normal", "control"),
    tumor_keywords   = c("tumor", "breast cancer", "tnbc", "er\\+"),
    exclude_keywords = c("cell line"),
    expected_tumour  = 84L,
    expected_normal  = 56L
  )
)

PATHS <- list(
  sig     = file.path("output", "tables",
                      "final_stable_biomarker_signature.rds"),
  dir_sig = file.path("output", "tables",
                      "Final_Directional_Signature.rds")
)

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, data.table, dplyr, tidyr, tibble, stringr,
               ggplot2, scales, readxl, ggrepel, R.utils)

# Bioconductor packages: GEOquery for phenotype retrieval, org.Hs.eg.db for
# Entrez ID to gene symbol mapping
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("GEOquery", "AnnotationDbi", "org.Hs.eg.db", "Biobase")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
library(GEOquery); library(AnnotationDbi); library(org.Hs.eg.db); library(Biobase)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

dirs <- list(
  out     = here("output", "geo_concordance"),
  figures = here("output", "geo_concordance", "figures"),
  tables  = here("output", "geo_concordance", "tables"),
  logs    = here("output", "geo_concordance", "logs")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(dirs$logs,
                      paste0("14_GEOConcordance_",
                             format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
sink(log_path, append = FALSE, split = TRUE)
t_start <- Sys.time()
set.seed(CONFIG$seed)

cat("=============================================================\n")
cat(" 14_geo_direction_concordance.R\n")
cat(" Start :", as.character(t_start), "\n")
cat("=============================================================\n\n")

# =============================================================================
# AESTHETICS (matching all project scripts)
# =============================================================================
COL <- list(
  match   = "#1B7837",
  weak    = "#A6DBA0",
  diverge = "#B22222",
  uninf   = "grey80",
  tcga    = "#2166AC"
)

theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5,
                                      size = base_size + 2),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey30",
                                      size = base_size - 1,
                                      margin = margin(b = 5)),
      plot.caption     = element_text(size = base_size - 3.5, colour = "grey45",
                                      hjust = 0, lineheight = 1.3,
                                      margin = margin(t = 8)),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.8),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = unit(c(0.8, 1.2, 0.8, 0.8), "cm")
    )
}

save_figure <- function(p, name, w = 10, h = 8) {
  base <- file.path(dirs$figures, name)
  ggplot2::ggsave(paste0(base, ".png"),  p, width = w, height = h, dpi = 300)
  ggplot2::ggsave(paste0(base, ".tiff"), p, width = w, height = h, dpi = 600,
                  device = "tiff", compression = "lzw")
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h)
  print(p); grDevices::dev.off()
  cat("  [SAVED]", name, "\n")
}

read_flexible <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    as.data.frame(readxl::read_excel(path))
  } else {
    as.data.frame(data.table::fread(path, header = TRUE, check.names = FALSE))
  }
}

# =============================================================================
# 1. LOAD SIGNATURE AND ITS TCGA-DERIVED DIRECTIONS
# =============================================================================
cat("[1/5] Loading locked signature...\n")

stable_genes <- readRDS(PATHS$sig)
clean_genes  <- gsub("_", "-", stable_genes)
dir_sig      <- readRDS(PATHS$dir_sig)

up_genes   <- dir_sig$Gene[startsWith(dir_sig$Direction, "\u2191")]
down_genes <- dir_sig$Gene[startsWith(dir_sig$Direction, "\u2193")]
if (length(up_genes) == 0 & length(down_genes) == 0) {
  cat("  [WARN] Arrow parsing failed; falling back to log2FC sign\n")
  up_genes   <- dir_sig$Gene[dir_sig$log2FC > 0]
  down_genes <- dir_sig$Gene[dir_sig$log2FC < 0]
}

PATTERN <- data.frame(
  gene      = c(up_genes, down_genes),
  direction = c(rep("Up", length(up_genes)), rep("Down", length(down_genes))),
  stringsAsFactors = FALSE
)

# Attach the TCGA log2FC for later comparison, if available
if ("log2FC" %in% names(dir_sig)) {
  PATTERN <- PATTERN %>%
    left_join(dir_sig %>% select(gene = Gene, tcga_log2FC = log2FC),
              by = "gene")
} else {
  PATTERN$tcga_log2FC <- NA_real_
}

cat(sprintf("  Signature: %d genes (%d up, %d down)\n",
            nrow(PATTERN), length(up_genes), length(down_genes)))
cat("  Up   :", paste(up_genes,   collapse = ", "), "\n")
cat("  Down :", paste(down_genes, collapse = ", "), "\n\n")

# =============================================================================
# 2. HELPERS — LOADING AND PREPARING A GEO COHORT
# =============================================================================

# Concatenate all descriptive metadata fields into one searchable string
get_meta <- function(gse_obj) {
  meta     <- Biobase::pData(gse_obj)
  tgt_cols <- grep("title|source|characteristic|description",
                   colnames(meta), value = TRUE, ignore.case = TRUE)
  meta$SearchString <- apply(meta[, tgt_cols, drop = FALSE], 1,
                             function(x) tolower(paste(x, collapse = " | ")))
  meta$GSM <- meta$geo_accession
  meta
}

# Precedence: exclude > normal > tumour. Identical to Script 02.
map_group <- function(ss, tkw, nkw, exkw) {
  if (length(exkw) > 0 &&
      any(sapply(exkw, function(k) grepl(k, ss, ignore.case = TRUE))))
    return(NA_character_)
  if (any(sapply(nkw, function(k) grepl(k, ss, ignore.case = TRUE))))
    return("Normal")
  if (any(sapply(tkw, function(k) grepl(k, ss, ignore.case = TRUE))))
    return("Tumour")
  NA_character_
}

load_cohort <- function(ds) {
  
  gse_id <- ds$gse_id
  cat(sprintf("  \u2500\u2500 %s \u2500\u2500\n", gse_id))
  
  if (!file.exists(ds$file_path))
    stop(gse_id, ": expression file not found at ", ds$file_path)
  
  # ── Phenotype from GEO ────────────────────────────────────────────────
  gse_obj  <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)[[1]]
  meta_all <- get_meta(gse_obj) %>%
    mutate(Group = sapply(SearchString, map_group,
                          tkw  = ds$tumor_keywords,
                          nkw  = ds$normal_keywords,
                          exkw = ds$exclude_keywords))
  
  # Audit trail: every sample and how it was classified
  write.csv(meta_all %>% select(GSM, Group, SearchString),
            file.path(dirs$logs, paste0("14_Mapping_Log_", gse_id, ".csv")),
            row.names = FALSE)
  
  n_dropped <- sum(is.na(meta_all$Group))
  if (n_dropped > 0)
    cat(sprintf("    %d samples excluded or unmatched (see mapping log)\n",
                n_dropped))
  
  meta <- meta_all %>% filter(!is.na(Group))
  if (length(unique(meta$Group)) < 2)
    stop(gse_id, ": only one class present after mapping")
  
  # ── Expression matrix: Entrez IDs -> gene symbols ─────────────────────
  expr <- data.table::fread(ds$file_path)
  colnames(expr)[1] <- "EntrezID"
  expr$EntrezID <- as.character(expr$EntrezID)
  
  suppressMessages({
    sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = expr$EntrezID,
                                 column = "SYMBOL", keytype = "ENTREZID",
                                 multiVals = "first")
  })
  expr$gene_symbol <- sym
  
  expr_clean <- expr %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    select(-EntrezID) %>%
    group_by(gene_symbol) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") %>%
    tibble::column_to_rownames("gene_symbol")
  
  cat(sprintf("    Mapped %d Entrez IDs to %d unique gene symbols\n",
              nrow(expr), nrow(expr_clean)))
  
  # ── Align samples ─────────────────────────────────────────────────────
  common <- intersect(meta$GSM, colnames(expr_clean))
  if (length(common) == 0) stop(gse_id, ": no sample overlap")
  meta_al <- meta %>% filter(GSM %in% common)
  expr_al <- as.matrix(expr_clean[, meta_al$GSM, drop = FALSE])
  
  n_t <- sum(meta_al$Group == "Tumour")
  n_n <- sum(meta_al$Group == "Normal")
  cat(sprintf("    Retained: %d tumour, %d normal\n", n_t, n_n))
  
  if (n_t < CONFIG$min_group_n || n_n < CONFIG$min_group_n)
    stop(gse_id, ": insufficient samples (need >= ", CONFIG$min_group_n,
         " per group)")
  
  # ── Log2 guard: NCBI GEO files are linear TPM ─────────────────────────
  mx <- max(expr_al, na.rm = TRUE)
  if (mx > 25) {
    cat(sprintf("    max = %.1f > 25, applying log2(TPM+1)\n", mx))
    expr_al <- log2(expr_al + 1)
  } else {
    cat(sprintf("    max = %.1f, values already log-scaled\n", mx))
  }
  
  list(expr = expr_al, group = meta_al$Group, label = gse_id,
       n_tumour = n_t, n_normal = n_n, n_dropped = n_dropped)
}

# =============================================================================
# 3. LOAD BOTH COHORTS
# =============================================================================
cat("[2/5] Loading GEO cohorts...\n")

cohorts <- lapply(DATASETS, load_cohort)

# ── Verify retained counts against published cohort composition ─────────────
cat("\n  \u2500\u2500 Sample count verification \u2500\u2500\n")
for (i in seq_along(cohorts)) {
  co <- cohorts[[i]]; ds <- DATASETS[[i]]
  ok <- co$n_tumour == ds$expected_tumour && co$n_normal == ds$expected_normal
  cat(sprintf("    %-10s observed %d T / %d N | expected %d T / %d N  %s\n",
              co$label, co$n_tumour, co$n_normal,
              ds$expected_tumour, ds$expected_normal,
              ifelse(ok, "OK",
                     "<-- MISMATCH, check mapping log in output/geo_concordance/logs")))
}
cat("\n")

# =============================================================================
# 4. PER-GENE DIFFERENTIAL EXPRESSION AND CLASSIFICATION
# =============================================================================
cat("[3/5] Computing per-gene differential expression...\n")

test_cohort <- function(co, pattern) {
  
  res <- lapply(seq_len(nrow(pattern)), function(i) {
    
    g   <- pattern$gene[i]
    dir <- pattern$direction[i]
    
    if (!g %in% rownames(co$expr)) {
      return(data.frame(
        dataset = co$label, gene = g, expected = dir,
        n_tumour = co$n_tumour, n_normal = co$n_normal,
        mean_tumour = NA_real_, mean_normal = NA_real_,
        log2FC = NA_real_, pval = NA_real_,
        observed = NA_character_, present = FALSE,
        stringsAsFactors = FALSE))
    }
    
    v_t <- as.numeric(co$expr[g, co$group == "Tumour"])
    v_n <- as.numeric(co$expr[g, co$group == "Normal"])
    
    # Data are already on a log2 scale, so the fold change is a difference
    m_t <- mean(v_t, na.rm = TRUE)
    m_n <- mean(v_n, na.rm = TRUE)
    lfc <- m_t - m_n
    
    p <- tryCatch(
      stats::wilcox.test(v_t, v_n, exact = FALSE)$p.value,
      error = function(e) NA_real_)
    
    data.frame(
      dataset = co$label, gene = g, expected = dir,
      n_tumour = co$n_tumour, n_normal = co$n_normal,
      mean_tumour = m_t, mean_normal = m_n,
      log2FC = lfc, pval = p,
      observed = ifelse(is.na(lfc), NA_character_,
                        ifelse(lfc > 0, "Up", "Down")),
      present = TRUE,
      stringsAsFactors = FALSE)
  })
  
  out <- do.call(rbind, res)
  # BH correction within cohort, across the signature genes only
  out$qval <- p.adjust(out$pval, method = "BH")
  out
}

de_results <- do.call(rbind, lapply(cohorts, test_cohort, pattern = PATTERN))

de_results <- de_results %>%
  left_join(PATTERN %>% select(gene, tcga_log2FC), by = "gene") %>%
  mutate(
    concordant  = !is.na(observed) & observed == expected,
    passes_fc   = !is.na(log2FC) & abs(log2FC) >= CONFIG$logfc_threshold,
    significant = !is.na(qval) & qval < CONFIG$q_threshold & passes_fc,
    class = case_when(
      !present                  ~ "Absent",
      concordant &  significant ~ "Match",
      concordant & !significant ~ "Weak match",
      !concordant &  significant ~ "Divergent",
      TRUE                       ~ "Uninformative"
    )
  )

write.csv(de_results, file.path(dirs$tables, "14_gene_concordance_full.csv"),
          row.names = FALSE)
cat("  [SAVED] 14_gene_concordance_full.csv\n\n")

# ── Console summary per dataset ──────────────────────────────────────────────
for (ds in unique(de_results$dataset)) {
  sub <- de_results %>% filter(dataset == ds)
  cat(sprintf("  ── %s (%d tumour / %d normal) ──\n",
              ds, sub$n_tumour[1], sub$n_normal[1]))
  print(as.data.frame(
    sub %>%
      mutate(log2FC = round(log2FC, 3),
             qval   = ifelse(is.na(qval), NA_character_,
                             ifelse(qval < 0.001, "<0.001",
                                    sprintf("%.4f", qval)))) %>%
      select(gene, expected, observed, log2FC, qval, class)),
    row.names = FALSE)
  cat("\n")
}

# =============================================================================
# 5. SUMMARY AND FIGURES
# =============================================================================
cat("[4/5] Summarising...\n")

n_genes <- nrow(PATTERN)

summary_df <- de_results %>%
  group_by(dataset, n_tumour, n_normal) %>%
  summarise(
    n_match       = sum(class == "Match"),
    n_weak        = sum(class == "Weak match"),
    n_divergent   = sum(class == "Divergent"),
    n_uninf       = sum(class == "Uninformative"),
    n_absent      = sum(class == "Absent"),
    n_concordant  = sum(concordant, na.rm = TRUE),
    n_tested      = sum(present),
    .groups = "drop"
  ) %>%
  mutate(
    full_pattern = n_match == n_genes,
    binom_p = mapply(function(k, n)
      if (n > 0) stats::binom.test(k, n, p = 0.5)$p.value else NA_real_,
      n_concordant, n_tested)
  )

write.csv(summary_df, file.path(dirs$tables, "14_concordance_summary.csv"),
          row.names = FALSE)

cat("\n  ── SUMMARY ──\n")
print(as.data.frame(summary_df %>%
                      mutate(binom_p = signif(binom_p, 3))), row.names = FALSE)
cat("\n  [SAVED] 14_concordance_summary.csv\n\n")

cat("[5/5] Generating figures...\n")

class_levels <- c("Match", "Weak match", "Divergent", "Uninformative", "Absent")
class_cols   <- c("Match" = COL$match, "Weak match" = COL$weak,
                  "Divergent" = COL$diverge, "Uninformative" = COL$uninf,
                  "Absent" = "white")

gene_order <- c(PATTERN$gene[PATTERN$direction == "Up"],
                PATTERN$gene[PATTERN$direction == "Down"])

# ── Fig 14A: classification matrix ──────────────────────────────────────────
p_matrix <- de_results %>%
  mutate(gene    = factor(gene, levels = rev(gene_order)),
         class   = factor(class, levels = class_levels),
         ds_lab  = sprintf("%s\n(%d T / %d N)", dataset, n_tumour, n_normal)) %>%
  ggplot(aes(x = ds_lab, y = gene, fill = class)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(log2FC), "",
                               sprintf("%.2f", log2FC))),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = class_cols, name = NULL, drop = FALSE) +
  labs(
    title    = "Per-gene direction concordance in independent GEO cohorts",
    subtitle = sprintf("Expected direction from TCGA-BRCA  |  q < %.2f, |log2FC| \u2265 %.1f",
                       CONFIG$q_threshold, CONFIG$logfc_threshold),
    x = NULL, y = NULL,
    caption = paste0(
      "Tile labels give the observed log2 fold change in each cohort.\n",
      "Match = concordant direction and significant. ",
      "Weak match = concordant direction, not significant (limited power).\n",
      "Divergent = opposite direction and significant. ",
      "Only Divergent indicates failure to replicate.")
  ) +
  theme_publication() +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(face = "bold"),
        legend.position = "top")

save_figure(p_matrix, "Fig_14A_Concordance_Matrix", w = 8, h = 7)

# ── Fig 14B: TCGA vs GEO effect size comparison ─────────────────────────────
if (any(!is.na(de_results$tcga_log2FC))) {
  
  p_lfc <- de_results %>%
    filter(!is.na(log2FC), !is.na(tcga_log2FC)) %>%
    mutate(class = factor(class, levels = class_levels)) %>%
    ggplot(aes(x = tcga_log2FC, y = log2FC, fill = class)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey40") +
    geom_point(shape = 21, size = 4, colour = "grey20", stroke = 0.4) +
    ggrepel::geom_text_repel(aes(label = gene), size = 3,
                             fontface = "italic", max.overlaps = 20) +
    scale_fill_manual(values = class_cols, name = NULL, drop = FALSE) +
    facet_wrap(~ dataset) +
    labs(
      title    = "Effect size comparison: TCGA-BRCA versus GEO cohorts",
      subtitle = "Points in the lower-left and upper-right quadrants are directionally concordant",
      x = expression("TCGA-BRCA  log"[2]*"FC"),
      y = expression("GEO cohort  log"[2]*"FC"),
      caption = paste0(
        "Dashed line indicates identity. Points below the line show ",
        "attenuated effect sizes in the GEO cohort,\n",
        "which is expected given smaller sample sizes and platform differences.")
    ) +
    theme_publication() +
    theme(legend.position = "top")
  
  save_figure(p_lfc, "Fig_14B_LogFC_Comparison", w = 11, h = 6)
}

# =============================================================================
# REPORT
# =============================================================================
report <- file.path(dirs$out, "14_GEO_CONCORDANCE_REPORT.txt")
con <- file(report, open = "wt")

writeLines(c(
  "=============================================================",
  " PER-GENE DIRECTION CONCORDANCE IN GEO VALIDATION COHORTS",
  paste0(" Generated: ", Sys.time()),
  "=============================================================",
  "",
  " EXPECTED DIRECTIONS (from TCGA-BRCA):",
  paste0("   Up   : ", paste(up_genes,   collapse = ", ")),
  paste0("   Down : ", paste(down_genes, collapse = ", ")),
  "",
  paste0(" q threshold      : ", CONFIG$q_threshold),
  paste0(" |log2FC| minimum : ", CONFIG$logfc_threshold),
  " Test             : Wilcoxon rank-sum, BH-corrected within cohort",
  ""
), con = con)

for (ds in unique(de_results$dataset)) {
  sub <- de_results %>% filter(dataset == ds)
  co  <- Filter(function(z) z$label == ds, cohorts)[[1]]
  dcfg <- Filter(function(z) z$gse_id == ds, DATASETS)[[1]]
  
  writeLines(c("", sprintf("─── %s (%d tumour / %d normal) ───",
                           ds, sub$n_tumour[1], sub$n_normal[1]),
               sprintf("  Samples excluded or unmatched: %d%s", co$n_dropped,
                       if (length(dcfg$exclude_keywords) > 0)
                         paste0(" (exclude keywords: ",
                                paste(dcfg$exclude_keywords, collapse = ", "), ")") else ""),
               ""), con = con)
  capture.output(print(as.data.frame(
    sub %>% mutate(log2FC = round(log2FC, 3),
                   qval = signif(qval, 3)) %>%
      select(gene, expected, observed, log2FC, qval, class)),
    row.names = FALSE), file = con, append = TRUE)
}

writeLines(c("", "─── SUMMARY ───", ""), con = con)
capture.output(print(as.data.frame(summary_df), row.names = FALSE),
               file = con, append = TRUE)

writeLines(c(
  "", "─── INTERPRETATION ───", "",
  "  These cohorts are substantially smaller than TCGA-BRCA. Genes with",
  "  modest discovery effect sizes may fail to reach significance here.",
  "  Such results are Weak match, indicating limited power, and must NOT",
  "  be described as failures to replicate.",
  "",
  "  Only Divergent classifications constitute evidence that a gene does",
  "  not reproduce its assigned direction.",
  "",
  "  Report sample sizes alongside every result.",
  ""
), con = con)
close(con)

cat("  [SAVED] 14_GEO_CONCORDANCE_REPORT.txt\n")
writeLines(capture.output(sessionInfo()),
           file.path(dirs$logs, "14_session_info.txt"))

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 14_geo_direction_concordance.R COMPLETE\n")
cat(" Runtime:", round(difftime(t_end, t_start, units = "mins"), 2), "min\n")
cat("=============================================================\n")
sink()
cat("Done. Report:", report, "\n")