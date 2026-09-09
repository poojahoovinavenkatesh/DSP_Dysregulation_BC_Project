# =============================================================================
# Script      : 8_microarray_direction_concordance.R
# Project     : DSP family dysregulation in breast cancer
#
# PURPOSE
#   Test whether each signature gene individually reproduces its
#   TCGA-BRCA-derived direction of dysregulation on MICROARRAY platforms.
#
#   Companion to Script 14 (RNA-seq GEO cohorts). Together they establish
#   whether directional dysregulation is preserved across both measurement
#   technologies at the level of individual genes, not only the composite
#   score (Scripts 02 and 07).
#
# COHORTS
#   GSE70947   Agilent SurePrint G3 8x60K (GPL13607)
#              148 tumour / 148 paired adjacent normal  (Quigley et al. 2017)
#   GSE15852   Affymetrix U133A (GPL96)
#              43 tumour / 43 paired normal             (Pau Ni et al. 2010)
#
# WHY MICROARRAY MATTERS HERE
#   Direction is scale-free. Unlike the composite score, whose absolute
#   magnitude depends on platform intensity units and is therefore not
#   comparable between RNA-seq and microarray, the SIGN of a fold change
#   is directly comparable across platforms. This analysis is consequently
#   the more rigorous cross-platform test.
#
# CLASSIFICATION (identical to Scripts 13 and 14)
#   Match          concordant direction + significant
#   Weak match     concordant direction + not significant  -> underpowered
#   Divergent      discordant direction + significant      -> real difference
#   Uninformative  discordant direction + not significant
#
# PLATFORM CAVEATS
#   1. Microarray dynamic range is compressed relative to RNA-seq. Absolute
#      log2FC values will be systematically SMALLER than TCGA values; this
#      is a platform property, not attenuated biology. The |log2FC|
#      threshold is therefore configurable and set lower by default.
#   2. Probe coverage is incomplete. Genes absent from a platform are
#      classified "Absent", never imputed as zero, since imputation would
#      fabricate a fold change.
#   3. Multi-probe genes are collapsed by mean.
#
# OUTPUTS
#   Tables:  15_microarray_concordance_full.csv
#            15_microarray_summary.csv
#            15_probe_coverage.csv
#   Figures: Fig_15A_Microarray_Concordance
#            Fig_15B_Platform_LogFC_Comparison
#   Report:  15_MICROARRAY_CONCORDANCE_REPORT.txt
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# =============================================================================

rm(list = ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed            = 42L,
  q_threshold     = 0.05,
  
  # Effect-size threshold(s) applied to |log2FC|.
  #
  # The first is the primary criterion; any further values are reported as
  # sensitivity analyses. Supplying a single value runs one criterion only.
  #
  # Microarray intensity compresses dynamic range relative to RNA-seq, so a
  # lower threshold can be justified. In practice 0.2 and 0.3 produced
  # identical classifications on both cohorts, so the RNA-seq criterion of
  # 0.3 is applied by default for consistency across all validation
  # analyses. To run both, use: logfc_thresholds = c(0.2, 0.3)
  logfc_thresholds = 0.3,
  
  min_group_n      = 10L
)

# Guard against an empty or malformed threshold specification
CONFIG$logfc_thresholds <- unique(sort(as.numeric(
  CONFIG$logfc_thresholds[!is.na(CONFIG$logfc_thresholds)]), decreasing = TRUE))
if (length(CONFIG$logfc_thresholds) == 0)
  stop("CONFIG$logfc_thresholds is empty — supply at least one numeric value")

DATASETS <- list(
  list(
    gse_id           = "GSE70947",
    platform         = "Agilent SurePrint G3 Human GE 8x60K (GPL13607)",
    platform_short   = "Agilent",
    # Sample titles take the form "CM016-normal" / "CM016-tumor";
    # the hyphen prefix prevents substring collisions.
    normal_keywords  = c("-normal"),
    tumor_keywords   = c("-tumor"),
    exclude_keywords = character(0),
    expected_tumour  = 148L,
    expected_normal  = 148L
  ),
  list(
    gse_id           = "GSE15852",
    platform         = "Affymetrix Human Genome U133A (GPL96)",
    platform_short   = "Affymetrix",
    normal_keywords  = c("normal", "non-tumor", "nontumor", "adjacent",
                         "control"),
    tumor_keywords   = c("tumor", "tumour", "cancer", "carcinoma",
                         "breast cancer"),
    exclude_keywords = c("cell line"),
    expected_tumour  = 43L,
    expected_normal  = 43L
  )
)

PATHS <- list(
  sig       = file.path("output", "tables",
                        "final_stable_biomarker_signature.rds"),
  dir_sig   = file.path("output", "tables",
                        "Final_Directional_Signature.rds"),
  geo_cache = here::here("data", "geo_cache")
)

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, data.table, dplyr, tidyr, tibble, stringr,
               ggplot2, scales, ggrepel)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("GEOquery", "Biobase")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
library(GEOquery); library(Biobase)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

dirs <- list(
  out     = here("output", "microarray_concordance"),
  figures = here("output", "microarray_concordance", "figures"),
  tables  = here("output", "microarray_concordance", "tables"),
  logs    = here("output", "microarray_concordance", "logs")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$geo_cache, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(dirs$logs,
                      paste0("15_MicroarrayConcordance_",
                             format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
sink(log_path, append = FALSE, split = TRUE)
t_start <- Sys.time()
set.seed(CONFIG$seed)

cat("=============================================================\n")
cat(" 15_microarray_direction_concordance.R\n")
cat(" Start :", as.character(t_start), "\n")
cat("=============================================================\n\n")

# =============================================================================
# AESTHETICS
# =============================================================================
COL <- list(
  match   = "#1B7837",
  weak    = "#A6DBA0",
  diverge = "#B22222",
  uninf   = "grey80",
  absent  = "white"
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

# =============================================================================
# 1. LOAD SIGNATURE AND TCGA-DERIVED DIRECTIONS
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
# 2. HELPERS
# =============================================================================

# Precedence: exclude > normal > tumour (identical to Scripts 02, 07, 14)
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

# Probe-level matrix -> gene-level matrix.
# Symbol column names differ by platform, hence the wide search list.
probe_to_gene <- function(gse_obj, gse_id) {
  
  expr_raw <- Biobase::exprs(gse_obj)
  cat(sprintf("    Raw matrix: %d probes x %d samples\n",
              nrow(expr_raw), ncol(expr_raw)))
  
  feat <- Biobase::fData(gse_obj)
  if (ncol(feat) == 0)
    stop(gse_id, ": fData() empty — platform annotation did not load. ",
         "Clear the GEO cache and retry with getGPL = TRUE.")
  
  sym_col <- intersect(
    c("Gene Symbol", "GENE_SYMBOL", "gene_symbol", "Symbol", "Gene_Symbol",
      "SYMBOL", "gene symbol", "GeneName", "gene_assignment",
      "mrna_assignment"),
    colnames(feat))[1]
  
  if (is.na(sym_col))
    sym_col <- grep("^symbol$|gene.symbol|gene_sym", colnames(feat),
                    value = TRUE, ignore.case = TRUE)[1]
  
  if (is.na(sym_col))
    stop(gse_id, ": no gene symbol column in fData. Columns present: ",
         paste(colnames(feat), collapse = " | "))
  
  cat(sprintf("    Symbol column: '%s'\n", sym_col))
  
  gene_vec <- as.character(feat[[sym_col]])
  # Multi-gene probe entries e.g. "BRCA1 /// TP53" -> take the first
  gene_vec <- stringr::str_trim(
    stringr::str_split_fixed(gene_vec, " /// |\\|//", 2)[, 1])
  gene_vec[gene_vec %in% c("---", "", "NA", "N/A", "na")] <- NA_character_
  
  cat(sprintf("    Probes with a symbol: %d / %d (%.1f%%)\n",
              sum(!is.na(gene_vec)), length(gene_vec),
              100 * sum(!is.na(gene_vec)) / length(gene_vec)))
  
  expr_df      <- as.data.frame(expr_raw)
  expr_df$Gene <- gene_vec
  
  expr_gene <- expr_df %>%
    filter(!is.na(Gene), nchar(Gene) > 1) %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") %>%
    tibble::column_to_rownames("Gene")
  
  cat(sprintf("    Unique gene symbols: %d\n", nrow(expr_gene)))
  as.matrix(expr_gene)
}

load_microarray_cohort <- function(ds) {
  
  gse_id <- ds$gse_id
  cat(sprintf("  \u2500\u2500 %s (%s) \u2500\u2500\n", gse_id, ds$platform_short))
  
  gse_obj <- tryCatch(
    GEOquery::getGEO(gse_id, destdir = PATHS$geo_cache, GSEMatrix = TRUE,
                     AnnotGPL = TRUE, getGPL = TRUE)[[1]],
    error = function(e) {
      cat("    [ERROR] GEOquery failed:", conditionMessage(e), "\n")
      NULL })
  if (is.null(gse_obj)) return(NULL)
  
  # ── Sample metadata ───────────────────────────────────────────────────
  # Title is the primary search field. Series-level fields such as contact
  # institution contain words like "Cancer Center" for every sample and
  # would misclassify the entire cohort if included.
  meta_raw  <- Biobase::pData(gse_obj)
  title_col <- grep("^title$", colnames(meta_raw), value = TRUE,
                    ignore.case = TRUE)[1]
  char_cols <- grep("characteristics|source_name|description",
                    colnames(meta_raw), value = TRUE, ignore.case = TRUE)
  
  meta_raw$SearchString <- if (!is.na(title_col)) {
    tolower(as.character(meta_raw[[title_col]]))
  } else {
    apply(meta_raw[, char_cols, drop = FALSE], 1,
          function(x) tolower(paste(x, collapse = " | ")))
  }
  meta_raw$GSM <- rownames(meta_raw)
  
  cat("    SearchString preview: ",
      paste(head(meta_raw$SearchString, 2), collapse = " ; "), "\n")
  
  meta_raw$Group <- sapply(meta_raw$SearchString, map_group,
                           tkw = ds$tumor_keywords,
                           nkw = ds$normal_keywords,
                           exkw = ds$exclude_keywords)
  
  # Fallback: widen the search string if title alone finds only one class
  if ((sum(meta_raw$Group == "Tumour", na.rm = TRUE) == 0 ||
       sum(meta_raw$Group == "Normal", na.rm = TRUE) == 0) &&
      length(char_cols) > 0) {
    cat("    [WARN] Title-only mapping incomplete; adding characteristics\n")
    meta_raw$SearchString <- apply(
      meta_raw[, c(title_col, char_cols), drop = FALSE], 1,
      function(x) tolower(paste(x, collapse = " | ")))
    meta_raw$Group <- sapply(meta_raw$SearchString, map_group,
                             tkw = ds$tumor_keywords,
                             nkw = ds$normal_keywords,
                             exkw = ds$exclude_keywords)
  }
  
  write.csv(meta_raw %>% select(GSM, Group, SearchString),
            file.path(dirs$logs, paste0("15_Mapping_Log_", gse_id, ".csv")),
            row.names = FALSE)
  
  n_dropped <- sum(is.na(meta_raw$Group))
  if (n_dropped > 0)
    cat(sprintf("    %d samples excluded or unmatched (see mapping log)\n",
                n_dropped))
  
  meta <- meta_raw %>% filter(!is.na(Group))
  if (length(unique(meta$Group)) < 2)
    stop(gse_id, ": only one class identified")
  
  # ── Expression ────────────────────────────────────────────────────────
  expr_gene <- probe_to_gene(gse_obj, gse_id)
  
  common <- intersect(meta$GSM, colnames(expr_gene))
  if (length(common) == 0) stop(gse_id, ": no sample overlap")
  meta_al <- meta %>% filter(GSM %in% common)
  expr_al <- expr_gene[, meta_al$GSM, drop = FALSE]
  
  n_t <- sum(meta_al$Group == "Tumour")
  n_n <- sum(meta_al$Group == "Normal")
  cat(sprintf("    Retained: %d tumour, %d normal\n", n_t, n_n))
  
  if (n_t < CONFIG$min_group_n || n_n < CONFIG$min_group_n)
    stop(gse_id, ": insufficient samples")
  
  # ── Scale check ───────────────────────────────────────────────────────
  # GEO series matrices for these platforms are normally already log2
  # intensity. Transform only if the range indicates linear values.
  mx <- max(expr_al, na.rm = TRUE)
  if (mx > 100) {
    cat(sprintf("    max = %.1f > 100, applying log2(x+1)\n", mx))
    expr_al <- log2(expr_al + 1)
  } else {
    cat(sprintf("    max = %.1f, already log2 intensity\n", mx))
  }
  
  # ── Probe coverage for signature genes ────────────────────────────────
  found <- intersect(clean_genes, rownames(expr_al))
  miss  <- setdiff(clean_genes, rownames(expr_al))
  cat(sprintf("    Signature coverage: %d / %d genes\n",
              length(found), length(clean_genes)))
  if (length(miss) > 0)
    cat("    [WARN] Not on platform:", paste(miss, collapse = ", "), "\n")
  
  list(expr = expr_al, group = meta_al$Group, label = gse_id,
       platform = ds$platform, platform_short = ds$platform_short,
       n_tumour = n_t, n_normal = n_n, n_dropped = n_dropped,
       found = found, miss = miss)
}

# =============================================================================
# 3. LOAD COHORTS
# =============================================================================
cat("[2/5] Loading microarray cohorts...\n")

cohorts <- Filter(Negate(is.null), lapply(DATASETS, load_microarray_cohort))
if (length(cohorts) == 0) stop("No cohorts loaded successfully")

cat("\n  \u2500\u2500 Sample count verification \u2500\u2500\n")
for (co in cohorts) {
  ds <- Filter(function(z) z$gse_id == co$label, DATASETS)[[1]]
  ok <- co$n_tumour == ds$expected_tumour && co$n_normal == ds$expected_normal
  cat(sprintf("    %-10s observed %d T / %d N | expected %d T / %d N  %s\n",
              co$label, co$n_tumour, co$n_normal,
              ds$expected_tumour, ds$expected_normal,
              ifelse(ok, "OK",
                     "<-- MISMATCH, inspect mapping log")))
}

# Probe coverage record
coverage_df <- do.call(rbind, lapply(cohorts, function(co) data.frame(
  dataset       = co$label,
  platform      = co$platform,
  n_found       = length(co$found),
  n_total       = length(clean_genes),
  pct_covered   = round(100 * length(co$found) / length(clean_genes), 1),
  missing_genes = paste(co$miss, collapse = ";"),
  stringsAsFactors = FALSE)))
write.csv(coverage_df, file.path(dirs$tables, "15_probe_coverage.csv"),
          row.names = FALSE)
cat("\n  [SAVED] 15_probe_coverage.csv\n\n")

# =============================================================================
# 4. PER-GENE DIFFERENTIAL EXPRESSION
# =============================================================================
cat("[3/5] Computing per-gene differential expression...\n")

test_cohort <- function(co, pattern) {
  
  res <- lapply(seq_len(nrow(pattern)), function(i) {
    
    g   <- pattern$gene[i]
    dir <- pattern$direction[i]
    
    # Absent genes are NOT imputed. Imputation would fabricate a fold change.
    if (!g %in% rownames(co$expr)) {
      return(data.frame(
        dataset = co$label, platform = co$platform_short,
        gene = g, expected = dir,
        n_tumour = co$n_tumour, n_normal = co$n_normal,
        mean_tumour = NA_real_, mean_normal = NA_real_,
        log2FC = NA_real_, pval = NA_real_,
        observed = NA_character_, present = FALSE,
        stringsAsFactors = FALSE))
    }
    
    v_t <- as.numeric(co$expr[g, co$group == "Tumour"])
    v_n <- as.numeric(co$expr[g, co$group == "Normal"])
    
    m_t <- mean(v_t, na.rm = TRUE)
    m_n <- mean(v_n, na.rm = TRUE)
    lfc <- m_t - m_n   # values are log2 intensity, so difference = log2FC
    
    p <- tryCatch(stats::wilcox.test(v_t, v_n, exact = FALSE)$p.value,
                  error = function(e) NA_real_)
    
    data.frame(
      dataset = co$label, platform = co$platform_short,
      gene = g, expected = dir,
      n_tumour = co$n_tumour, n_normal = co$n_normal,
      mean_tumour = m_t, mean_normal = m_n,
      log2FC = lfc, pval = p,
      observed = ifelse(is.na(lfc), NA_character_,
                        ifelse(lfc > 0, "Up", "Down")),
      present = TRUE,
      stringsAsFactors = FALSE)
  })
  
  out <- do.call(rbind, res)
  out$qval <- p.adjust(out$pval, method = "BH")
  out
}

de_results <- do.call(rbind, lapply(cohorts, test_cohort, pattern = PATTERN))

classify <- function(df, fc_thresh, tag) {
  df %>%
    mutate(
      threshold   = tag,
      concordant  = !is.na(observed) & observed == expected,
      passes_fc   = !is.na(log2FC) & abs(log2FC) >= fc_thresh,
      significant = !is.na(qval) & qval < CONFIG$q_threshold & passes_fc,
      class = case_when(
        !present                   ~ "Absent",
        concordant &  significant ~ "Match",
        concordant & !significant ~ "Weak match",
        !concordant &  significant ~ "Divergent",
        TRUE                       ~ "Uninformative"
      ))
}

de_all <- bind_rows(lapply(CONFIG$logfc_thresholds, function(fc)
  classify(de_results, fc, sprintf("|log2FC| >= %.1f", fc)))) %>%
  left_join(PATTERN %>% select(gene, tcga_log2FC), by = "gene")

# Primary criterion = first (largest) threshold; any others are sensitivity
primary_tag <- sprintf("|log2FC| >= %.1f", CONFIG$logfc_thresholds[1])

write.csv(de_all,
          file.path(dirs$tables, "15_microarray_concordance_full.csv"),
          row.names = FALSE)
cat("  [SAVED] 15_microarray_concordance_full.csv\n\n")

for (th in unique(de_all$threshold)) {
  cat(sprintf("  \u2500\u2500 %s \u2500\u2500\n", th))
  for (ds in unique(de_all$dataset)) {
    sub <- de_all %>% filter(threshold == th, dataset == ds)
    cat(sprintf("\n  %s (%s, %d T / %d N)\n", ds, sub$platform[1],
                sub$n_tumour[1], sub$n_normal[1]))
    print(as.data.frame(
      sub %>% mutate(
        log2FC = round(log2FC, 3),
        qval   = ifelse(is.na(qval), NA_character_,
                        ifelse(qval < 0.001, "<0.001", sprintf("%.4f", qval)))) %>%
        select(gene, expected, observed, log2FC, qval, class)),
      row.names = FALSE)
  }
  cat("\n")
}

# =============================================================================
# 5. SUMMARY AND FIGURES
# =============================================================================
cat("[4/5] Summarising...\n")

n_genes <- nrow(PATTERN)

summary_df <- de_all %>%
  group_by(threshold, dataset, platform, n_tumour, n_normal) %>%
  summarise(
    n_match      = sum(class == "Match"),
    n_weak       = sum(class == "Weak match"),
    n_divergent  = sum(class == "Divergent"),
    n_uninf      = sum(class == "Uninformative"),
    n_absent     = sum(class == "Absent"),
    n_concordant = sum(concordant, na.rm = TRUE),
    n_tested     = sum(present),
    .groups = "drop"
  ) %>%
  mutate(
    full_pattern = n_match == n_genes,
    binom_p = mapply(function(k, n)
      if (n > 0) stats::binom.test(k, n, p = 0.5)$p.value else NA_real_,
      n_concordant, n_tested)
  )

write.csv(summary_df, file.path(dirs$tables, "15_microarray_summary.csv"),
          row.names = FALSE)

cat("\n  \u2500\u2500 SUMMARY \u2500\u2500\n")
print(as.data.frame(summary_df %>% mutate(binom_p = signif(binom_p, 3))),
      row.names = FALSE)
cat("\n  [SAVED] 15_microarray_summary.csv\n\n")

cat("[5/5] Generating figures...\n")

class_levels <- c("Match", "Weak match", "Divergent", "Uninformative", "Absent")
class_cols   <- c("Match" = COL$match, "Weak match" = COL$weak,
                  "Divergent" = COL$diverge, "Uninformative" = COL$uninf,
                  "Absent" = COL$absent)
gene_order <- c(PATTERN$gene[PATTERN$direction == "Up"],
                PATTERN$gene[PATTERN$direction == "Down"])

# ── Fig 15A: concordance matrix, primary threshold ──────────────────────────
p_matrix <- de_all %>%
  filter(threshold == primary_tag) %>%
  mutate(gene   = factor(gene, levels = rev(gene_order)),
         class  = factor(class, levels = class_levels),
         ds_lab = sprintf("%s\n%s\n(%d T / %d N)",
                          dataset, platform, n_tumour, n_normal)) %>%
  ggplot(aes(x = ds_lab, y = gene, fill = class)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(log2FC), "n/a",
                               sprintf("%.2f", log2FC)),
                colour = class == "Absent"),
            size = 3.2, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = class_cols, name = NULL, drop = FALSE) +
  scale_colour_manual(values = c("TRUE" = "grey40", "FALSE" = "white")) +
  labs(
    title    = "Per-gene direction concordance on microarray platforms",
    subtitle = sprintf("Expected direction from TCGA-BRCA  |  q < %.2f, %s",
                       CONFIG$q_threshold, primary_tag),
    x = NULL, y = NULL,
    caption = paste0(
      "Tile labels give the observed log2 fold change; 'n/a' indicates the ",
      "gene is not represented on that platform.\n",
      "Microarray intensity compresses dynamic range, so effect sizes are ",
      "systematically smaller than RNA-seq values.\n",
      "Absent genes are not imputed. Only Divergent indicates failure to ",
      "replicate.")
  ) +
  theme_publication() +
  theme(axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(face = "bold", size = 9),
        legend.position = "top")

save_figure(p_matrix, "Fig_15A_Microarray_Concordance", w = 8.5, h = 7)

# ── Fig 15B: TCGA vs microarray effect sizes ────────────────────────────────
if (any(!is.na(de_all$tcga_log2FC))) {
  
  p_lfc <- de_all %>%
    filter(threshold == primary_tag, !is.na(log2FC), !is.na(tcga_log2FC)) %>%
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
    facet_wrap(~ paste0(dataset, "\n", platform)) +
    labs(
      title    = "Effect size: TCGA-BRCA RNA-seq versus microarray cohorts",
      subtitle = "Upper-right and lower-left quadrants indicate directional concordance",
      x = expression("TCGA-BRCA  log"[2]*"FC  (RNA-seq)"),
      y = expression("Microarray  log"[2]*"FC  (intensity)"),
      caption = paste0(
        "Dashed line is identity. Points falling below it reflect the ",
        "compressed dynamic range of microarray\n",
        "intensity relative to RNA-seq, not weaker biological effect.")
    ) +
    theme_publication() +
    theme(legend.position = "top")
  
  save_figure(p_lfc, "Fig_15B_Platform_LogFC_Comparison", w = 11, h = 6)
}

# =============================================================================
# REPORT
# =============================================================================
report <- file.path(dirs$out, "15_MICROARRAY_CONCORDANCE_REPORT.txt")
con <- file(report, open = "wt")

writeLines(c(
  "=============================================================",
  " PER-GENE DIRECTION CONCORDANCE — MICROARRAY PLATFORMS",
  paste0(" Generated: ", Sys.time()),
  "=============================================================",
  "",
  " EXPECTED DIRECTIONS (from TCGA-BRCA RNA-seq):",
  paste0("   Up   : ", paste(up_genes,   collapse = ", ")),
  paste0("   Down : ", paste(down_genes, collapse = ", ")),
  "",
  paste0(" q threshold          : ", CONFIG$q_threshold),
  paste0(" |log2FC| threshold(s): ",
         paste(sprintf("%.1f", CONFIG$logfc_thresholds), collapse = ", "),
         if (length(CONFIG$logfc_thresholds) > 1)
           "  (first = primary, remainder = sensitivity)" else ""),
  " Test                 : Wilcoxon rank-sum, BH-corrected within cohort",
  ""
), con = con)

writeLines(c("", "─── PROBE COVERAGE ───", ""), con = con)
capture.output(print(as.data.frame(
  coverage_df %>% select(dataset, platform, n_found, n_total, pct_covered)),
  row.names = FALSE), file = con, append = TRUE)

for (th in unique(de_all$threshold)) {
  writeLines(c("", sprintf("═══ THRESHOLD: %s ═══", th)), con = con)
  for (ds in unique(de_all$dataset)) {
    sub <- de_all %>% filter(threshold == th, dataset == ds)
    co  <- Filter(function(z) z$label == ds, cohorts)[[1]]
    writeLines(c("", sprintf("─── %s | %s (%d T / %d N) ───",
                             ds, sub$platform[1],
                             sub$n_tumour[1], sub$n_normal[1]),
                 sprintf("  Samples excluded or unmatched: %d", co$n_dropped),
                 ""), con = con)
    capture.output(print(as.data.frame(
      sub %>% mutate(log2FC = round(log2FC, 3), qval = signif(qval, 3)) %>%
        select(gene, expected, observed, log2FC, qval, class)),
      row.names = FALSE), file = con, append = TRUE)
  }
}

writeLines(c("", "─── SUMMARY ───", ""), con = con)
capture.output(print(as.data.frame(summary_df), row.names = FALSE),
               file = con, append = TRUE)

writeLines(c(
  "", "─── INTERPRETATION ───", "",
  "  Direction is scale-free and therefore directly comparable between",
  "  RNA-seq and microarray, unlike the composite score, whose absolute",
  "  magnitude depends on platform intensity units.",
  "",
  "  Microarray log2FC values are systematically smaller than RNA-seq",
  "  values because intensity measurement compresses dynamic range. A",
  "  reduced effect size is therefore expected and does not indicate",
  "  weaker biology.",
  "",
  if (length(CONFIG$logfc_thresholds) > 1)
    paste0("  Multiple effect-size thresholds were applied; the first is the\n",
           "  primary criterion and the remainder are sensitivity analyses.")
  else
    paste0("  A single effect-size threshold of ",
           sprintf("%.1f", CONFIG$logfc_thresholds[1]),
           " was applied, identical to\n",
           "  the criterion used in the RNA-seq validation analyses."),
  "",
  "  Genes absent from a platform are reported as Absent and are never",
  "  imputed, since imputation would fabricate a fold change.",
  "",
  "  Only Divergent classifications indicate failure to replicate.",
  "  Uninformative indicates no significant change in either direction,",
  "  which in small cohorts reflects limited power rather than",
  "  contradiction.",
  ""
), con = con)
close(con)

cat("  [SAVED] 15_MICROARRAY_CONCORDANCE_REPORT.txt\n")
writeLines(capture.output(sessionInfo()),
           file.path(dirs$logs, "15_session_info.txt"))

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 15_microarray_direction_concordance.R COMPLETE\n")
cat(" Runtime:", round(difftime(t_end, t_start, units = "mins"), 2), "min\n")
cat("=============================================================\n")
sink()
cat("Done. Report:", report, "\n")