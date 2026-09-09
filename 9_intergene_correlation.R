# =============================================================================
# Script      : 9_intergene_correlation.R
# Project     : DSP family dysregulation in breast cancer
# Description : Generates a publication-ready HALF-heatmap (lower triangle) 
#               for inter-gene correlation of the DSP signature in TCGA-BRCA.
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# Seed        : 42
# Depends on  : 01_training_pipeline.R (signature genes), TCGA-BRCA TPM matrix
# =============================================================================

rm(list = ls()); gc()

# =============================================================================
# CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed = 42L,
  
  # ---------------------------------------------------------------------------
  # Magnitude bands for interpreting correlation strength.
  # Defined a priori, BEFORE inspecting results, so that thresholds are not
  # fitted to the observed data.
  # ---------------------------------------------------------------------------
  rho_strong   = 0.5,   # |rho| >= 0.50        -> strong
  rho_moderate = 0.3,   # 0.30 <= |rho| < 0.50 -> moderate
  # |rho| <  0.30        -> weak
  
  # ---------------------------------------------------------------------------
  # Significance markers on heatmap tiles.
  # Default FALSE. With n ~ 1,100 tumour samples, nominal significance is
  # reached at |rho| ~ 0.06, so asterisks would mark nearly every pair and
  # would misleadingly imply widespread co-regulation. Magnitude is the
  # interpretable quantity here; p and q values are reported in the
  # supplementary table instead.
  # ---------------------------------------------------------------------------
  show_significance = FALSE,
  q_threshold       = 0.05
)

PATHS <- list(
  sig     = file.path("data", "directional_score", "final_stable_biomarker_signature.rds"),
  dir_sig = file.path("data", "directional_score", "Final_Directional_Signature.rds"),
  brca    = file.path("data",   "TCGA-BRCA.star_tpm.tsv.gz")
)

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here, data.table, dplyr, tibble, stringr, tidyr,
  ggplot2, ggcorrplot, scales,
  AnnotationDbi, org.Hs.eg.db
)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# =============================================================================
# DIRECTORIES
# =============================================================================
dirs <- list(
  out     = here("output_intergene_correlation"),
  figures = here("output_intergene_correlation", "figures"),
  tables  = here("output_intergene_correlation", "tables"),
  logs    = here("output_intergene_correlation", "logs")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Verify every directory exists and is writable; fail early with a clear
# message rather than silently losing outputs later.
for (nm in names(dirs)) {
  d <- dirs[[nm]]
  if (!dir.exists(d))
    stop("Could not create directory '", nm, "' at: ", d,
         "\nCheck working directory and write permissions.")
}
cat("Output directories:\n")
for (nm in names(dirs))
  cat(sprintf("  %-8s %s\n", nm, normalizePath(dirs[[nm]], winslash = "/")))
cat("\n")

log_path <- file.path(dirs$logs,
                      paste0("11_InterGene_Corr_Run_",
                             format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
sink(log_path, append = FALSE, split = TRUE)
cat("============================================================\n")
cat(" 11_intergene_correlation.R\n")
cat(" Start :", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

set.seed(CONFIG$seed)

# =============================================================================
# AESTHETICS 
# =============================================================================
COL <- list(
  tumor    = "#B22222",  normal   = "#2166AC",
  up       = "#A32D2D",  down     = "#185FA5",
  high     = "#D53E4F",  low      = "#3288BD",
  mid      = "#FFFFBF"
)

theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5,
                                      size = base_size + 2),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey30",
                                      size = base_size - 1,
                                      margin = margin(b = 5)),
      plot.caption     = element_text(size = base_size - 3.5,
                                      colour = "grey45", hjust = 0,
                                      lineheight = 1.3,
                                      margin = margin(t = 8)),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), # Keeps empty space clean
      panel.border     = element_blank(), # Removed border for a floating triangle look
      axis.line        = element_blank(), 
      axis.ticks       = element_blank(), # Removes tick marks
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = unit(c(0.8, 1.2, 0.8, 0.8), "cm")
    )
}

# Write a table and confirm it exists, reporting the absolute path.
# Silent write failures are the main cause of "no output" reports.
save_table <- function(x, filename, dir = dirs$tables) {
  fp <- file.path(dir, filename)
  ok <- tryCatch({
    utils::write.csv(x, fp, row.names = FALSE)
    TRUE
  }, error = function(e) {
    cat("  [ERROR] Could not write", filename, ":", conditionMessage(e), "\n")
    FALSE
  })
  if (ok && file.exists(fp)) {
    cat(sprintf("  [SAVED] %-45s (%d rows) -> %s\n",
                filename, nrow(as.data.frame(x)),
                normalizePath(fp, winslash = "/")))
  } else if (ok) {
    cat("  [WARN] write.csv reported success but file not found:", fp, "\n")
  }
  invisible(fp)
}

# Matrix variant: preserves row names, which write.csv drops when
# row.names = FALSE.
save_matrix <- function(m, filename, dir = dirs$tables) {
  fp <- file.path(dir, filename)
  df <- data.frame(Gene = rownames(m), as.data.frame(m),
                   check.names = FALSE, stringsAsFactors = FALSE)
  ok <- tryCatch({
    utils::write.csv(df, fp, row.names = FALSE)
    TRUE
  }, error = function(e) {
    cat("  [ERROR] Could not write", filename, ":", conditionMessage(e), "\n")
    FALSE
  })
  if (ok && file.exists(fp))
    cat(sprintf("  [SAVED] %-45s (%dx%d) -> %s\n",
                filename, nrow(m), ncol(m),
                normalizePath(fp, winslash = "/")))
  invisible(fp)
}

save_figure <- function(p, name, w = 9, h = 8) {
  base <- file.path(dirs$figures, name)
  ggplot2::ggsave(paste0(base, ".png"),  plot = p, width = w, height = h, dpi = 300)
  ggplot2::ggsave(paste0(base, ".tiff"), plot = p, width = w, height = h, dpi = 600,
                  device = "tiff", compression = "lzw")
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h)
  print(p); grDevices::dev.off()
  cat("  [SAVED]", name, "\n")
}

# =============================================================================
# 1. INTER-GENE CORRELATION HEATMAP
# =============================================================================
cat("[1/1] Inter-gene correlation heatmap\n")

# Load signature
stable_genes <- readRDS(PATHS$sig)
clean_genes  <- gsub("_", "-", stable_genes)
cat("  Signature genes:", paste(clean_genes, collapse = ", "), "\n")

# Load directional info to colour gene labels
dir_sig <- readRDS(PATHS$dir_sig)
up_genes   <- dir_sig$Gene[startsWith(dir_sig$Direction, "\u2191")]
down_genes <- dir_sig$Gene[startsWith(dir_sig$Direction, "\u2193")]
if (length(up_genes) == 0 & length(down_genes) == 0) {
  up_genes   <- dir_sig$Gene[dir_sig$log2FC > 0]
  down_genes <- dir_sig$Gene[dir_sig$log2FC < 0]
}

# Load TCGA-BRCA, extract signature genes, filter to TUMOUR samples
df_brca <- data.table::fread(PATHS$brca, header = TRUE, sep = "\t")
if (!"gene_id" %in% names(df_brca)) names(df_brca)[1] <- "gene_id"
df_brca$gene_id_clean <- str_split_fixed(df_brca$gene_id, "\\.", 2)[, 1]
suppressMessages({
  df_brca$gene_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = df_brca$gene_id_clean,
    column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
})

expr_mat <- df_brca %>%
  filter(!is.na(gene_symbol), gene_symbol %in% clean_genes) %>%
  group_by(gene_symbol) %>%
  summarise(across(starts_with("TCGA"), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix() %>% t() %>% as.data.frame() %>%
  mutate(sample_type = substr(rownames(.), 14, 15)) %>%
  filter(sample_type == "01") %>%   # TUMOUR ONLY
  select(-sample_type) %>%
  na.omit()

cat(sprintf("  Tumour samples: %d × %d genes\n",
            nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# CORRELATION WITH SIGNIFICANCE TESTING
# =============================================================================
# cor() alone returns no p-values. Pairwise cor.test() is used so that
# p-values and BH-adjusted q-values can be reported in the supplementary
# table, even though interpretation is by magnitude (see CONFIG note).
# =============================================================================
X <- as.matrix(expr_mat[, clean_genes])
g <- colnames(X)
k <- length(g)

cor_mat  <- matrix(NA_real_, k, k, dimnames = list(g, g))
pval_mat <- matrix(NA_real_, k, k, dimnames = list(g, g))
diag(cor_mat)  <- 1
diag(pval_mat) <- 0

for (i in seq_len(k - 1)) {
  for (j in (i + 1):k) {
    ct <- suppressWarnings(
      stats::cor.test(X[, i], X[, j], method = "spearman", exact = FALSE))
    cor_mat[i, j]  <- cor_mat[j, i]  <- as.numeric(ct$estimate)
    pval_mat[i, j] <- pval_mat[j, i] <- ct$p.value
  }
}

n_samp <- nrow(X)

# ── Long-format pairwise table: rho, p, BH q, magnitude band ───────────────
pair_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)

pair_df <- data.frame(
  Gene1     = rownames(cor_mat)[pair_idx[, "row"]],
  Gene2     = colnames(cor_mat)[pair_idx[, "col"]],
  rho       = cor_mat[pair_idx],
  p_value   = pval_mat[pair_idx],
  stringsAsFactors = FALSE
)
pair_df$q_value <- p.adjust(pair_df$p_value, method = "BH")

# =============================================================================
# CRITICAL |rho| AT NOMINAL AND BH-ADJUSTED SIGNIFICANCE
# =============================================================================
# Two distinct quantities:
#
#   rho_crit_p  the |rho| at which an UNADJUSTED two-sided test reaches
#               p < alpha. Fixed by sample size alone, computed analytically
#               via the Fisher z transformation.
#
#   rho_crit_q  the |rho| at which the BH-ADJUSTED q falls below alpha.
#               NOT a fixed formula: Benjamini-Hochberg rejects when
#               p(i) <= (i/m) * alpha, so the effective p-cutoff depends on
#               the observed distribution of p-values. It is therefore
#               derived empirically from the data, as the largest p-value
#               among pairs that survive correction.
#
# rho_crit_q >= rho_crit_p always, since correction is conservative.
# =============================================================================
alpha  <- CONFIG$q_threshold
m_test <- nrow(pair_df)

# Fisher z back-transform: rho = tanh(z / sqrt(n - 3))
p_to_rho <- function(p, n) {
  tanh(stats::qnorm(1 - p / 2) / sqrt(n - 3))
}

# ── Nominal threshold ─────────────────────────────────────────────────────
rho_crit_p <- p_to_rho(alpha, n_samp)

# ── BH-adjusted threshold, derived from the observed p-values ─────────────
sig_idx <- which(pair_df$q_value < alpha)

if (length(sig_idx) > 0) {
  # Largest p-value still surviving BH correction = the effective cutoff
  p_bh_cutoff <- max(pair_df$p_value[sig_idx])
  rho_crit_q  <- p_to_rho(p_bh_cutoff, n_samp)
  bh_note     <- sprintf("empirical, from %d of %d pairs surviving correction",
                         length(sig_idx), m_test)
} else {
  # No pair survives: the most stringent BH cutoff applies (rank 1)
  p_bh_cutoff <- alpha / m_test
  rho_crit_q  <- p_to_rho(p_bh_cutoff, n_samp)
  bh_note     <- "no pair survived correction; rank-1 BH cutoff shown"
}

# Retained for backward compatibility with downstream captions
rho_crit <- rho_crit_p

cat("\n  ── Critical |rho| thresholds ──\n")
cat(sprintf("    n = %d samples, %d pairwise tests, alpha = %.2f\n",
            n_samp, m_test, alpha))
cat(sprintf("    Nominal   p < %.2f  reached at |rho| >= %.4f\n",
            alpha, rho_crit_p))
cat(sprintf("    BH-adj    q < %.2f  reached at |rho| >= %.4f  (p <= %.3g)\n",
            alpha, rho_crit_q, p_bh_cutoff))
cat(sprintf("    %s\n", bh_note))
cat("    NOTE: both thresholds are far below the pre-specified magnitude\n")
cat(sprintf("          bands (moderate >= %.1f, strong >= %.1f). Statistical\n",
            CONFIG$rho_moderate, CONFIG$rho_strong))
cat("          significance is therefore not the limiting criterion here.\n")

# Spearman note: the Fisher z approximation assumes Pearson correlation.
# For Spearman the standard error is inflated by roughly 1.03 (Bonett &
# Wright 2000), so the true thresholds are marginally higher. The
# difference is negligible at this sample size but is reported for accuracy.
rho_crit_p_spearman <- tanh(stats::qnorm(1 - alpha / 2) *
                              sqrt(1.06) / sqrt(n_samp - 3))
cat(sprintf("    Spearman-corrected nominal threshold: |rho| >= %.4f\n\n",
            rho_crit_p_spearman))

pair_df <- pair_df %>%
  mutate(
    Gene1_direction = ifelse(Gene1 %in% up_genes, "Up", "Down"),
    Gene2_direction = ifelse(Gene2 %in% up_genes, "Up", "Down"),
    magnitude = dplyr::case_when(
      abs(rho) >= CONFIG$rho_strong   ~ "Strong",
      abs(rho) >= CONFIG$rho_moderate ~ "Moderate",
      TRUE                            ~ "Weak"
    ),
    significant = q_value < CONFIG$q_threshold,
    n_samples   = n_samp
  ) %>%
  arrange(desc(abs(rho)))

save_table(
  pair_df %>%
    mutate(rho     = round(rho, 4),
           p_value = signif(p_value, 3),
           q_value = signif(q_value, 3)),
  "11A_inter_gene_correlation_pairwise.csv")

# ── Magnitude band summary ─────────────────────────────────────────────────
band_tab <- pair_df %>%
  count(magnitude, name = "n_pairs") %>%
  mutate(pct = round(100 * n_pairs / nrow(pair_df), 1))

cat("\n  ── Magnitude bands across all", nrow(pair_df), "pairs ──\n")
print(as.data.frame(band_tab), row.names = FALSE)

cat(sprintf("\n  Pairs with q < %.2f: %d of %d (%.1f%%)\n",
            CONFIG$q_threshold, sum(pair_df$significant), nrow(pair_df),
            100 * sum(pair_df$significant) / nrow(pair_df)))
cat("  NOTE: high significance rate reflects sample size, not co-regulation.\n")
cat("        Interpret by magnitude, not by significance.\n\n")

cat("  ── Strongest pairs (|rho| >= ", CONFIG$rho_moderate, ") ──\n", sep = "")
strong_pairs <- pair_df %>% filter(abs(rho) >= CONFIG$rho_moderate)
if (nrow(strong_pairs) == 0) {
  cat("    none\n\n")
} else {
  print(as.data.frame(strong_pairs %>%
                        mutate(rho = round(rho, 3), q_value = signif(q_value, 3)) %>%
                        select(Gene1, Gene2, rho, q_value, magnitude)), row.names = FALSE)
  cat("\n")
}

# Save square matrices
save_matrix(round(cor_mat, 4),  "11A_inter_gene_correlation.csv")
save_matrix(signif(pval_mat, 3), "11A_inter_gene_pvalues.csv")

# BH-adjusted q-value matrix, reconstructed from the pairwise table so that
# correction is applied once across all unique pairs, not per column.
qval_mat <- matrix(NA_real_, k, k, dimnames = list(g, g))
diag(qval_mat) <- 0
for (r in seq_len(nrow(pair_df))) {
  i <- pair_df$Gene1[r]; j <- pair_df$Gene2[r]
  qval_mat[i, j] <- qval_mat[j, i] <- pair_df$q_value[r]
}
save_matrix(signif(qval_mat, 3), "11A_inter_gene_qvalues.csv")

# Magnitude band summary
save_table(band_tab, "11A_correlation_magnitude_summary.csv")

# Analysis parameters, for reproducibility
save_table(
  data.frame(
    parameter = c("n_tumour_samples", "n_genes", "n_unique_pairs",
                  "correlation_method", "rho_strong_threshold",
                  "rho_moderate_threshold", "alpha",
                  "critical_rho_nominal_p",
                  "critical_rho_nominal_p_spearman_corrected",
                  "critical_rho_BH_adjusted_q",
                  "effective_BH_p_cutoff",
                  "n_pairs_significant_BH",
                  "multiple_testing", "seed"),
    value = c(as.character(n_samp), as.character(k),
              as.character(nrow(pair_df)), "Spearman",
              as.character(CONFIG$rho_strong),
              as.character(CONFIG$rho_moderate),
              as.character(CONFIG$q_threshold),
              sprintf("%.4f", rho_crit_p),
              sprintf("%.4f", rho_crit_p_spearman),
              sprintf("%.4f", rho_crit_q),
              sprintf("%.3g", p_bh_cutoff),
              as.character(sum(pair_df$significant)),
              "Benjamini-Hochberg", as.character(CONFIG$seed)),
    stringsAsFactors = FALSE),
  "11A_analysis_parameters.csv")

# Build heatmap with hierarchical clustering
hc <- hclust(as.dist(1 - abs(cor_mat)), method = "average")
ord <- hc$order
cor_mat_ord <- cor_mat[ord, ord]

# =============================================================================
# MASK UPPER TRIANGLE FOR HALF-HEATMAP
# Change `diag = FALSE` to `diag = TRUE` if you want to also hide the 1.00 diagonal.
# =============================================================================
cor_mat_ord[upper.tri(cor_mat_ord, diag = FALSE)] <- NA

# Convert to long-format for ggplot
cor_long <- as.data.frame(cor_mat_ord) %>%
  tibble::rownames_to_column("RowGene") %>%
  pivot_longer(-RowGene, names_to = "ColGene", values_to = "rho") %>%
  filter(!is.na(rho)) %>%   # Drop the NA values so they don't plot as empty tiles
  mutate(
    # Set X-axis left-to-right
    ColGene = factor(ColGene, levels = colnames(cor_mat_ord)),
    # Reverse Y-axis top-to-bottom so the triangle aligns perfectly to the bottom-left
    RowGene = factor(RowGene, levels = rev(rownames(cor_mat_ord)))
  )

# ── Attach BH q-values so tiles can optionally carry significance markers ──
# Join is symmetric: the pairwise table holds each pair once, so match on
# both orderings.
qlook <- rbind(
  pair_df %>% select(RowGene = Gene1, ColGene = Gene2, q_value),
  pair_df %>% select(RowGene = Gene2, ColGene = Gene1, q_value)
) %>% distinct()

cor_long <- cor_long %>%
  mutate(RowGene_chr = as.character(RowGene),
         ColGene_chr = as.character(ColGene)) %>%
  left_join(qlook, by = c("RowGene_chr" = "RowGene",
                          "ColGene_chr" = "ColGene")) %>%
  mutate(
    # Diagonal has no test
    q_value = ifelse(RowGene_chr == ColGene_chr, NA_real_, q_value),
    stars = dplyr::case_when(
      is.na(q_value)   ~ "",
      q_value < 0.001  ~ "***",
      q_value < 0.01   ~ "**",
      q_value < 0.05   ~ "*",
      TRUE             ~ ""
    ),
    tile_label = if (isTRUE(CONFIG$show_significance)) {
      paste0(sprintf("%.2f", rho), ifelse(stars == "", "", paste0("\n", stars)))
    } else {
      sprintf("%.2f", rho)
    }
  )

# Extract levels so we can assign axis text colors perfectly
genes_x <- levels(cor_long$ColGene)
genes_y <- levels(cor_long$RowGene)

col_x <- ifelse(genes_x %in% up_genes, COL$up, COL$down)
col_y <- ifelse(genes_y %in% up_genes, COL$up, COL$down)

# Plotting the Half Heatmap
p_corr <- ggplot(cor_long, aes(x = ColGene, y = RowGene, fill = rho)) +
  geom_tile(colour = "white", linewidth = 1) +  # Thicker white line for cleaner grid
  geom_text(aes(label = tile_label,
                colour = abs(rho) >= CONFIG$rho_strong),
            size = 3.5, fontface = "bold", lineheight = 0.9) +
  scale_fill_gradient2(
    low = COL$low, mid = "white", high = COL$high,
    midpoint = 0, limits = c(-1, 1),
    name = expression(Spearman ~ rho),
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white"),
                      guide = "none") +
  # Move X axis to the bottom to properly frame the triangle
  scale_x_discrete(position = "bottom") + 
  coord_equal() +
  labs(
    title    = "Inter-gene correlation across the DSP signature",
    subtitle = paste0("Spearman correlation in TCGA-BRCA tumour samples  |  n = ", nrow(expr_mat)),
    x = NULL, y = NULL,
    caption  = paste0(
      "Block structure reveals subfamily-specific co-regulation. ",
      "Up-regulated genes (red labels), down-regulated genes (blue labels).\n",
      "Correlations are interpreted by magnitude: |rho| >= ",
      sprintf("%.1f", CONFIG$rho_strong), " strong, ",
      sprintf("%.1f", CONFIG$rho_moderate), "-",
      sprintf("%.1f", CONFIG$rho_strong), " moderate, < ",
      sprintf("%.1f", CONFIG$rho_moderate), " weak.\n",
      "At n = ", n_samp, ", BH-adjusted significance (q < ",
      sprintf("%.2f", CONFIG$q_threshold), ") is reached at |rho| >= ",
      sprintf("%.2f", rho_crit_q), "; p and q values are given in\n",
      "Supplementary Table SX. Significance is therefore not the limiting ",
      "criterion at this sample size.\n",
      "Low cross-cluster correlation supports independent pathway membership.")
  ) +
  theme_publication() +
  theme(
    # Angle text correctly for bottom axis
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               colour = col_x, face = "bold", size = 12),
    axis.text.y = element_text(colour = col_y,
                               face = "bold", size = 12),
    legend.position = "right",
    legend.justification = c(0, 1) # Align legend to the top-right
  )

save_figure(p_corr, "Fig_InterGene_Correlation_Half", w = 9, h = 8)

# =============================================================================
# WRAP UP
# =============================================================================
# =============================================================================
# OUTPUT MANIFEST
# Lists what is actually present on disk, so a missing file is visible
# immediately rather than being discovered later.
# =============================================================================
cat("\n=============================================================\n")
cat(" OUTPUT MANIFEST\n")
cat("=============================================================\n")

expected <- list(
  tables = c("11A_inter_gene_correlation.csv",
             "11A_inter_gene_pvalues.csv",
             "11A_inter_gene_qvalues.csv",
             "11A_inter_gene_correlation_pairwise.csv",
             "11A_correlation_magnitude_summary.csv",
             "11A_analysis_parameters.csv"),
  figures = c("Fig_InterGene_Correlation_Half.png",
              "Fig_InterGene_Correlation_Half.tiff",
              "Fig_InterGene_Correlation_Half.pdf")
)

n_missing <- 0L
for (grp in names(expected)) {
  cat(sprintf("\n %s  (%s)\n", toupper(grp),
              normalizePath(dirs[[grp]], winslash = "/")))
  for (f in expected[[grp]]) {
    fp  <- file.path(dirs[[grp]], f)
    ok  <- file.exists(fp)
    sz  <- if (ok) file.info(fp)$size else NA
    if (!ok) n_missing <- n_missing + 1L
    cat(sprintf("   %s %-45s %s\n",
                ifelse(ok, "[OK]     ", "[MISSING]"), f,
                ifelse(ok, sprintf("%.1f KB", sz / 1024), "")))
  }
}

if (n_missing > 0) {
  cat(sprintf("\n [WARN] %d expected file(s) missing.\n", n_missing))
  cat(" Check the console above for [ERROR] lines and confirm write\n")
  cat(" permissions on the output directory.\n")
} else {
  cat("\n All expected outputs written successfully.\n")
}

cat("\n=============================================================\n")
cat(" 11_intergene_correlation.R  COMPLETE\n")
cat(" End:", as.character(Sys.time()), "\n")
cat("=============================================================\n")

sink()
cat("Done.\n")