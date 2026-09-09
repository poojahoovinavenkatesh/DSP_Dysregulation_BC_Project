# =============================================================================
# Script      : 4_pattern_conservation.R
# Project     : DSP family dysregulation in breast cancer
#
# PURPOSE
#   Test whether the breast-derived 9-gene DIRECTIONAL PATTERN is reproduced
#   in any other cancer type. The pattern is treated as an indivisible
#   constant: a cancer type "reproduces" it only if ALL nine genes are both
#   significantly changed AND concordant in direction.
#
# WHY THIS REPLACES THE ORGAN-SPECIFICITY SCORE COMPARISON
#   Comparing absolute directional scores across organs confounds
#   cancer-associated dysregulation with baseline tissue-of-origin
#   expression. Direction is scale-free and carries no such confound.
#
# KEY DESIGN DECISION — the four-class scheme
#   Collapsing "significant AND concordant" into a single pass/fail would
#   make every small cohort fail automatically, guaranteeing the specificity
#   result by construction. Direction and significance are therefore scored
#   SEPARATELY:
#
#     Match         concordant direction + significant  -> shared dysregulation
#     Weak match    concordant direction + NS           -> underpowered
#     Divergent     discordant direction + significant  -> real difference
#     Uninformative discordant direction + NS           -> no evidence
#
#   Only Divergent counts as evidence AGAINST pattern sharing.
#
# OUTPUTS
#   Tables:  13_gene_classification_full.csv     per gene x cancer x threshold
#            13_cancer_summary.csv               per-cancer class counts + tests
#            13_gene_conservation.csv            per-gene conservation counts
#            13_excluded_cohorts.csv             what was dropped and why
#   Figures: Fig_13A_Sign_Matrix                 gene x cancer classification
#            Fig_13B_Cancer_Concordance          per-cancer stacked bar
#            Fig_13C_Gene_Conservation           per-gene conservation bar
#   Report:  13_PATTERN_CONSERVATION_REPORT.txt
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# =============================================================================

rm(list = ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed              = 42L,

  # Primary and sensitivity significance thresholds
  q_primary         = 0.05,
  q_sensitivity     = 0.10,

  # Effect size threshold (must match the screening criterion in Methods 2.2)
  logfc_threshold   = 0.3,

  # Minimum total normal samples for a cancer type to be analysed.
  # Cancers below this have no usable reference and are excluded.
  min_normals       = 10L,

  # Minimum tumour samples
  min_tumours       = 30L,

  # Reference cancer — pattern is defined here, so it is 9/9 by construction
  reference_cancer  = "BRCA"
)

PATHS <- list(
  gepia   = file.path("data", "master_gene_expression_gepia_T_PT.csv"),
  samples = file.path("data", "Sample_number.csv")
)

# =============================================================================
# THE CONSTANT — breast-derived directional pattern
# Edit ONLY if the locked signature changes. This is the hypothesis under test.
# =============================================================================
PATTERN <- data.frame(
  gene      = c("CDKN3", "PTP4A2", "PTP4A3", "STYXL1",
                "DUSP1", "DUSP6", "TNS1", "CDC14B", "EPM2A"),
  direction = c("Up", "Up", "Up", "Up",
                "Down", "Down", "Down", "Down", "Down"),
  stringsAsFactors = FALSE
)

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, data.table, dplyr, tidyr, tibble, stringr,
               ggplot2, scales, readxl)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

dirs <- list(
  out     = here("output", "pattern_conservation"),
  figures = here("output", "pattern_conservation", "figures"),
  tables  = here("output", "pattern_conservation", "tables"),
  logs    = here("output", "pattern_conservation", "logs")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(dirs$logs,
                      paste0("13_PatternConservation_",
                             format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
sink(log_path, append = FALSE, split = TRUE)
t_start <- Sys.time()
set.seed(CONFIG$seed)

cat("=============================================================\n")
cat(" 13_pattern_conservation.R\n")
cat(" Start :", as.character(t_start), "\n")
cat("=============================================================\n\n")

# =============================================================================
# AESTHETICS (matching all project scripts)
# =============================================================================
COL <- list(
  match   = "#1B7837",   # concordant + significant
  weak    = "#A6DBA0",   # concordant + NS
  diverge = "#B22222",   # discordant + significant
  uninf   = "grey80"     # discordant + NS
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
# 1. LOAD GEPIA MATRIX
# =============================================================================
cat("[1/6] Loading GEPIA differential expression matrix...\n")

read_flexible <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    as.data.frame(readxl::read_excel(path))
  } else {
    as.data.frame(data.table::fread(path, header = TRUE, check.names = FALSE))
  }
}

gepia <- read_flexible(PATHS$gepia)
cat(sprintf("  Rows: %d | Columns: %s\n", nrow(gepia),
            paste(names(gepia), collapse = ", ")))

# Standardise column names
name_map <- c(
  CancerType = "cancer", cancertype = "cancer", Cancer = "cancer",
  name = "gene", Gene = "gene", gene_name = "gene",
  log2FC = "log2FC", logFC = "log2FC",
  `q-value` = "qval", qvalue = "qval", `q_value` = "qval", adj.P.Val = "qval",
  `p-value` = "pval", pvalue = "pval", `p_value` = "pval"
)
for (old in names(name_map)) {
  if (old %in% names(gepia)) names(gepia)[names(gepia) == old] <- name_map[[old]]
}

req <- c("cancer", "gene", "log2FC", "qval")
missing_cols <- setdiff(req, names(gepia))
if (length(missing_cols) > 0)
  stop("Missing required columns after renaming: ",
       paste(missing_cols, collapse = ", "))

# Coerce numerics defensively — scientific notation stored as text is a
# common failure mode and silently produces NA in every comparison
gepia$log2FC <- suppressWarnings(as.numeric(as.character(gepia$log2FC)))
gepia$qval   <- suppressWarnings(as.numeric(as.character(gepia$qval)))

n_bad <- sum(is.na(gepia$log2FC) | is.na(gepia$qval))
if (n_bad > 0)
  cat(sprintf("  [WARN] %d rows have non-numeric log2FC or q-value\n", n_bad))

cat(sprintf("  Cancer types: %d | Unique genes: %d\n\n",
            length(unique(gepia$cancer)), length(unique(gepia$gene))))

# ── Verify all pattern genes present in every cancer ─────────────────────────
missing_genes <- setdiff(PATTERN$gene, unique(gepia$gene))
if (length(missing_genes) > 0)
  stop("Pattern genes absent from GEPIA matrix: ",
       paste(missing_genes, collapse = ", "))

coverage <- gepia %>%
  filter(gene %in% PATTERN$gene) %>%
  count(cancer, name = "n_genes")

incomplete <- coverage %>% filter(n_genes < nrow(PATTERN))
if (nrow(incomplete) > 0) {
  cat("  [WARN] Cancer types with incomplete gene coverage:\n")
  print(as.data.frame(incomplete), row.names = FALSE)
  cat("\n")
} else {
  cat("  [OK] All 9 pattern genes present in all cancer types\n\n")
}

# =============================================================================
# 2. LOAD SAMPLE SIZES AND APPLY EXCLUSIONS
# =============================================================================
cat("[2/6] Loading sample sizes and applying cohort exclusions...\n")

samples_raw <- read_flexible(PATHS$samples)

# The sample file may arrive unparsed (all fields in one column). Handle both.
if (ncol(samples_raw) == 1) {
  cat("  [INFO] Single-column file detected; re-parsing as CSV text\n")
  txt <- c(names(samples_raw)[1], samples_raw[[1]])
  samples_raw <- as.data.frame(
    data.table::fread(text = paste(txt, collapse = "\n"),
                      header = TRUE, check.names = FALSE))
}

names(samples_raw) <- names(samples_raw) %>%
  str_trim() %>% str_replace_all("\\s+", "_")

sn <- names(samples_raw)
col_abbr   <- sn[grepl("^Abbr",   sn, ignore.case = TRUE)][1]
col_tumour <- sn[grepl("^Tumor|^Tumour", sn, ignore.case = TRUE)][1]
col_normal <- sn[grepl("^Normal", sn, ignore.case = TRUE)][1]
col_gtex   <- sn[grepl("GTEx_Num|GTEx_num", sn, ignore.case = TRUE)][1]

num_or_zero <- function(x) {
  x <- as.character(x); x[x %in% c("-", "", "NA", "N/A")] <- "0"
  suppressWarnings(as.numeric(x))
}

samples <- data.frame(
  cancer       = as.character(samples_raw[[col_abbr]]),
  n_tumour     = num_or_zero(samples_raw[[col_tumour]]),
  n_tcga_norm  = num_or_zero(samples_raw[[col_normal]]),
  n_gtex_norm  = if (!is.na(col_gtex)) num_or_zero(samples_raw[[col_gtex]]) else 0,
  stringsAsFactors = FALSE
) %>%
  mutate(n_normal_total = n_tcga_norm + n_gtex_norm)

cat(sprintf("  Sample metadata for %d cancer types\n", nrow(samples)))

# ── Exclusions ───────────────────────────────────────────────────────────────
samples <- samples %>%
  mutate(
    exclude_reason = case_when(
      n_normal_total < CONFIG$min_normals ~
        sprintf("insufficient normals (%d < %d)", n_normal_total, CONFIG$min_normals),
      n_tumour < CONFIG$min_tumours ~
        sprintf("insufficient tumours (%d < %d)", n_tumour, CONFIG$min_tumours),
      TRUE ~ NA_character_
    ),
    included = is.na(exclude_reason)
  )

excluded <- samples %>% filter(!included)
if (nrow(excluded) > 0) {
  cat("\n  ── EXCLUDED COHORTS ──\n")
  print(as.data.frame(excluded %>%
          select(cancer, n_tumour, n_tcga_norm, n_gtex_norm,
                 n_normal_total, exclude_reason)), row.names = FALSE)
  write.csv(excluded, file.path(dirs$tables, "13_excluded_cohorts.csv"),
            row.names = FALSE)
}

included_cancers <- samples$cancer[samples$included]
cat(sprintf("\n  Cancers retained: %d of %d\n\n",
            length(included_cancers), nrow(samples)))

if (!CONFIG$reference_cancer %in% included_cancers)
  stop("Reference cancer ", CONFIG$reference_cancer, " was excluded — check inputs")

# =============================================================================
# 3. CLASSIFY EACH GENE IN EACH CANCER
# =============================================================================
cat("[3/6] Classifying gene-by-cancer directional concordance...\n")

classify <- function(df, q_thresh, label) {
  df %>%
    mutate(
      observed_dir = ifelse(log2FC > 0, "Up", "Down"),
      concordant   = observed_dir == direction,
      passes_fc    = abs(log2FC) >= CONFIG$logfc_threshold,
      significant  = !is.na(qval) & qval < q_thresh & passes_fc,
      class = case_when(
         concordant &  significant ~ "Match",
         concordant & !significant ~ "Weak match",
        !concordant &  significant ~ "Divergent",
        TRUE                       ~ "Uninformative"
      ),
      threshold = label
    )
}

base_df <- gepia %>%
  filter(gene %in% PATTERN$gene, cancer %in% included_cancers) %>%
  left_join(PATTERN, by = "gene") %>%
  left_join(samples %>% select(cancer, n_tumour, n_normal_total), by = "cancer")

classified <- bind_rows(
  classify(base_df, CONFIG$q_primary,
           sprintf("q < %.2f", CONFIG$q_primary)),
  classify(base_df, CONFIG$q_sensitivity,
           sprintf("q < %.2f", CONFIG$q_sensitivity))
)

write.csv(
  classified %>% select(threshold, cancer, gene, direction, observed_dir,
                        log2FC, qval, concordant, significant, class,
                        n_tumour, n_normal_total),
  file.path(dirs$tables, "13_gene_classification_full.csv"), row.names = FALSE)
cat("  [SAVED] 13_gene_classification_full.csv\n\n")

# =============================================================================
# 4. PER-CANCER SUMMARY AND PATTERN CONSERVATION TEST
# =============================================================================
cat("[4/6] Summarising per cancer type...\n")

n_genes <- nrow(PATTERN)

cancer_summary <- classified %>%
  group_by(threshold, cancer, n_tumour, n_normal_total) %>%
  summarise(
    n_match       = sum(class == "Match"),
    n_weak        = sum(class == "Weak match"),
    n_divergent   = sum(class == "Divergent"),
    n_uninf       = sum(class == "Uninformative"),
    n_concordant  = sum(concordant),
    .groups = "drop"
  ) %>%
  mutate(
    # PRIMARY OUTCOME: strict reproduction of the full pattern
    full_pattern  = n_match == n_genes,
    # Direction-only concordance vs chance (0.5 per gene)
    binom_p       = mapply(function(k) {
                      stats::binom.test(k, n_genes, p = 0.5,
                                        alternative = "two.sided")$p.value
                    }, n_concordant)
  ) %>%
  group_by(threshold) %>%
  mutate(binom_q = p.adjust(binom_p, method = "BH")) %>%
  ungroup() %>%
  arrange(threshold, desc(n_match), desc(n_concordant))

write.csv(cancer_summary,
          file.path(dirs$tables, "13_cancer_summary.csv"), row.names = FALSE)

for (th in unique(cancer_summary$threshold)) {
  cat(sprintf("\n  ── %s ──\n", th))
  sub <- cancer_summary %>% filter(threshold == th)
  full <- sub$cancer[sub$full_pattern]
  cat(sprintf("  Cancers reproducing the FULL 9-gene pattern: %s\n",
              if (length(full) == 0) "none" else paste(full, collapse = ", ")))
  cat("\n  Top 10 by Match count:\n")
  print(as.data.frame(sub %>%
          select(cancer, n_tumour, n_normal_total, n_match, n_weak,
                 n_divergent, n_concordant, binom_q) %>%
          mutate(binom_q = signif(binom_q, 3)) %>%
          head(10)), row.names = FALSE)
}
cat("\n  [SAVED] 13_cancer_summary.csv\n\n")

# =============================================================================
# 5. PER-GENE CONSERVATION
# =============================================================================
cat("[5/6] Per-gene conservation across cancer types...\n")

gene_conservation <- classified %>%
  filter(cancer != CONFIG$reference_cancer) %>%   # exclude reference: 9/9 by construction
  group_by(threshold, gene, direction) %>%
  summarise(
    n_cancers     = n(),
    n_match       = sum(class == "Match"),
    n_weak        = sum(class == "Weak match"),
    n_divergent   = sum(class == "Divergent"),
    pct_match     = round(100 * n_match / n(), 1),
    .groups = "drop"
  ) %>%
  arrange(threshold, desc(n_match))

write.csv(gene_conservation,
          file.path(dirs$tables, "13_gene_conservation.csv"), row.names = FALSE)

cat(sprintf("\n  (reference cancer %s excluded — 9/9 by construction)\n",
            CONFIG$reference_cancer))
for (th in unique(gene_conservation$threshold)) {
  cat(sprintf("\n  ── %s ──\n", th))
  print(as.data.frame(gene_conservation %>% filter(threshold == th) %>%
          select(gene, direction, n_match, n_weak, n_divergent, n_cancers)),
        row.names = FALSE)
}
cat("\n  [SAVED] 13_gene_conservation.csv\n\n")

# =============================================================================
# 6. FIGURES
# =============================================================================
cat("[6/6] Generating figures...\n")

class_levels <- c("Match", "Weak match", "Divergent", "Uninformative")
class_cols   <- c("Match" = COL$match, "Weak match" = COL$weak,
                  "Divergent" = COL$diverge, "Uninformative" = COL$uninf)

for (th in unique(classified$threshold)) {

  th_tag <- gsub("[^0-9]", "", th)
  sub_cl <- classified %>% filter(threshold == th)
  sub_su <- cancer_summary %>% filter(threshold == th)

  cancer_order <- sub_su$cancer[order(-sub_su$n_match, -sub_su$n_concordant)]
  gene_order   <- c(PATTERN$gene[PATTERN$direction == "Up"],
                    PATTERN$gene[PATTERN$direction == "Down"])

  # ── Fig 13A: sign / classification matrix ──────────────────────────────
  p_matrix <- sub_cl %>%
    mutate(cancer = factor(cancer, levels = cancer_order),
           gene   = factor(gene, levels = rev(gene_order)),
           class  = factor(class, levels = class_levels)) %>%
    ggplot(aes(x = cancer, y = gene, fill = class)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = class_cols, name = NULL, drop = FALSE) +
    labs(
      title    = "Conservation of the breast-derived DSP directional pattern",
      subtitle = sprintf("Significance threshold: %s  |  |log2FC| \u2265 %.1f",
                         th, CONFIG$logfc_threshold),
      x = NULL, y = NULL,
      caption = paste0(
        "Match = concordant direction and significant change. ",
        "Weak match = concordant direction, not significant (underpowered).\n",
        "Divergent = opposite direction and significant. ",
        "Only Divergent is evidence against pattern sharing.\n",
        "Cancer types ordered by number of Matches. ",
        "Reference cancer is concordant by construction.")
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 8.5),
          axis.text.y = element_text(face = "italic"),
          legend.position = "top")

  save_figure(p_matrix, sprintf("Fig_13A_Sign_Matrix_q%s", th_tag),
              w = 13, h = 6.5)

  # ── Fig 13B: per-cancer stacked classification ─────────────────────────
  p_bar <- sub_cl %>%
    mutate(cancer = factor(cancer, levels = cancer_order),
           class  = factor(class, levels = class_levels)) %>%
    count(cancer, class) %>%
    ggplot(aes(x = cancer, y = n, fill = class)) +
    geom_col(colour = "black", linewidth = 0.25) +
    geom_hline(yintercept = n_genes, linetype = "dashed", colour = "grey30") +
    scale_fill_manual(values = class_cols, name = NULL, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, n_genes, 1), expand = c(0, 0)) +
    labs(
      title    = "Per-cancer classification of the nine-gene pattern",
      subtitle = sprintf("Significance threshold: %s", th),
      x = NULL, y = "Number of genes",
      caption = paste0(
        "Full pattern conservation requires all nine genes classified as Match.\n",
        "Cancer types with many Weak matches are underpowered rather than divergent.")
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 8.5),
          legend.position = "top")

  save_figure(p_bar, sprintf("Fig_13B_Cancer_Concordance_q%s", th_tag),
              w = 13, h = 6.5)

  # ── Fig 13C: per-gene conservation ─────────────────────────────────────
  gc_sub <- gene_conservation %>% filter(threshold == th)
  p_gene <- gc_sub %>%
    mutate(gene = factor(gene, levels = gc_sub$gene[order(gc_sub$n_match)])) %>%
    select(gene, direction, Match = n_match, `Weak match` = n_weak,
           Divergent = n_divergent) %>%
    pivot_longer(c(Match, `Weak match`, Divergent),
                 names_to = "class", values_to = "n") %>%
    mutate(class = factor(class, levels = class_levels)) %>%
    ggplot(aes(x = n, y = gene, fill = class)) +
    geom_col(colour = "black", linewidth = 0.25) +
    scale_fill_manual(values = class_cols, name = NULL, drop = FALSE) +
    labs(
      title    = "Per-gene conservation across non-reference cancer types",
      subtitle = sprintf("Significance threshold: %s  |  reference cancer excluded",
                         th),
      x = "Number of cancer types", y = NULL,
      caption = paste0(
        "Genes with high Match counts represent pan-cancer DSP dysregulation.\n",
        "Genes matching in few cancer types are candidate breast-restricted events.")
    ) +
    theme_publication() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "top")

  save_figure(p_gene, sprintf("Fig_13C_Gene_Conservation_q%s", th_tag),
              w = 10, h = 7)
}

# =============================================================================
# REPORT
# =============================================================================
report <- file.path(dirs$out, "13_PATTERN_CONSERVATION_REPORT.txt")
con <- file(report, open = "wt")

writeLines(c(
  "=============================================================",
  " PAN-CANCER DIRECTIONAL PATTERN CONSERVATION",
  paste0(" Generated: ", Sys.time()),
  "=============================================================",
  "",
  " PATTERN UNDER TEST (breast-derived, treated as indivisible):",
  paste0("   Up in tumour   : ",
         paste(PATTERN$gene[PATTERN$direction == "Up"], collapse = ", ")),
  paste0("   Down in tumour : ",
         paste(PATTERN$gene[PATTERN$direction == "Down"], collapse = ", ")),
  "",
  paste0(" |log2FC| threshold : ", CONFIG$logfc_threshold),
  paste0(" Primary q          : ", CONFIG$q_primary),
  paste0(" Sensitivity q      : ", CONFIG$q_sensitivity),
  paste0(" Min normals        : ", CONFIG$min_normals),
  paste0(" Min tumours        : ", CONFIG$min_tumours),
  paste0(" Cancers analysed   : ", length(included_cancers)),
  paste0(" Cancers excluded   : ", nrow(excluded)),
  ""
), con = con)

if (nrow(excluded) > 0) {
  writeLines(c("─── EXCLUDED COHORTS ───", ""), con = con)
  capture.output(print(as.data.frame(
    excluded %>% select(cancer, n_tumour, n_normal_total, exclude_reason)),
    row.names = FALSE), file = con, append = TRUE)
}

for (th in unique(cancer_summary$threshold)) {
  sub <- cancer_summary %>% filter(threshold == th)
  full <- sub$cancer[sub$full_pattern]
  writeLines(c("", sprintf("─── THRESHOLD %s ───", th), "",
    sprintf("  PRIMARY OUTCOME — cancers reproducing all 9 genes as Match: %s",
            if (length(full) == 0) "none" else paste(full, collapse = ", ")),
    ""), con = con)
  capture.output(print(as.data.frame(
    sub %>% select(cancer, n_tumour, n_normal_total, n_match, n_weak,
                   n_divergent, n_uninf, n_concordant, binom_q) %>%
      mutate(binom_q = signif(binom_q, 3))), row.names = FALSE),
    file = con, append = TRUE)
}

writeLines(c("", "─── PER-GENE CONSERVATION ───",
             sprintf("  (reference cancer %s excluded)", CONFIG$reference_cancer),
             ""), con = con)
capture.output(print(as.data.frame(gene_conservation), row.names = FALSE),
               file = con, append = TRUE)

writeLines(c(
  "", "─── HOW TO INTERPRET ───", "",
  "  Only BRCA reproduces the full pattern",
  "    -> the nine-gene directional architecture is breast-restricted.",
  "",
  "  Other cancers reproduce it too",
  "    -> the pattern is shared; reframe as conserved across those types.",
  "",
  "  Many Weak matches, few Divergent",
  "    -> those cohorts are underpowered, NOT biologically different.",
  "       Do not claim specificity on the basis of Weak matches.",
  "",
  "  A cancer with high Divergent count",
  "    -> genuine biological difference. This is real evidence.",
  ""
), con = con)
close(con)

cat("\n  [SAVED] 13_PATTERN_CONSERVATION_REPORT.txt\n")

writeLines(capture.output(sessionInfo()),
           file.path(dirs$logs, "13_session_info.txt"))

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 13_pattern_conservation.R COMPLETE\n")
cat(" Runtime:", round(difftime(t_end, t_start, units = "mins"), 2), "min\n")
cat("=============================================================\n")
sink()
cat("Done. Report:", report, "\n")
