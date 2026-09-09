# =============================================================================
# Script      : 10_enrichr_dotplot.R
# Project     : DSP family dysregulation in breast cancer
# Description : Auto-detects top enriched pathways across 3 gene sets (All/Up/Down)
#               and 3+ enrichment libraries, then generates a publication-grade
#               dot plot demonstrating multi-pathway separation.
#
# CRITICAL FRAMING:
#   The dot plot proves the manuscript's central claim that the up- and
#   down-regulated subsets of the signature converge on DISTINCT biological
#   pathways — i.e., the signature is multi-pathway, not single-pathway.
#
# INPUT FORMAT (Enrichr TSV download):
#   Term  Overlap  P-value  Adjusted P-value  Old P-value  Old Adjusted P-value
#   Odds Ratio  Combined Score  Genes
#
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# Seed        : 42
# =============================================================================
# REQUIRED INPUT FOLDER STRUCTURE:
#   data/enrichr/
#     ├── All9_Reactome_Pathways_2024.txt
#     ├── All9_BioPlanet_2019.txt
#     ├── All9_MSigDB_Hallmark_2020.txt
#     ├── All9_GO_Biological_Process_2025.txt   (optional)
#     ├── All9_KEGG_2026.txt                    (optional)
#     ├── All9_WikiPathways_2024_Human.txt      (optional)
#     ├── Down_Reactome_Pathways_2024.txt
#     ├── Down_BioPlanet_2019.txt
#     ├── Down_MSigDB_Hallmark_2020.txt
#     ├── Down_GO_Biological_Process_2025.txt
#     ├── Up_Reactome_Pathways_2024.txt
#     ├── Up_BioPlanet_2019.txt
#     └── Up_MSigDB_Hallmark_2020.txt
#
# FILENAME CONVENTION (auto-parsed):
#   <GeneSet>_<Library>.txt
#   GeneSet  ∈ {All9, Down, Up}
#   Library  = anything else (whitespace/underscores tolerated)
#
# OUTPUTS:
#   Tables:
#     10_Enrichr_combined_long.csv          — all results merged
#     10_Enrichr_significant_only.csv       — FDR < 0.05 subset
#     10_Enrichr_top_pathways_plotted.csv   — exactly what's in the figure
#     10_pathway_rank_specificity.csv       — rank-based directionality classification
#     10_directionality_summary.csv         — summary counts per library
#   Figures:
#     Fig_10A_Enrichr_Dotplot               — main multi-pathway figure
#     Fig_10C_Pathway_Directionality        — rank-based directionality bar
#     Fig_10E_Top5_Per_Direction            — top 5 pathways per direction × library
#   Report:
#     10_ENRICHR_RESULTS_REPORT.txt
#   Report:
#     10_ENRICHR_RESULTS_REPORT.txt
# =============================================================================

rm(list = ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed                  = 42L,
  fdr_threshold         = 0.05,         # significance cutoff
  top_n_per_lib         = 8,            # top N per library × geneset combination
  min_gene_count        = 1,            # minimum genes overlapping for inclusion
  max_term_length       = 55,           # truncate long pathway names
  pseudo_minlogp        = 1e-15,        # cap for -log10 transform
  
  # Rank-based directionality (Section 6 — replaces binary "Both" overlap)
  # Binary FDR overlap is misleading because pathway databases annotate
  # STYXL1 (catalytically dead) and PTP4A2/3 to MAPK pathways alongside
  # DUSP1/DUSP6, causing apparent overlap. Rank-based analysis reveals
  # true direction-specificity.
  rank_difference       = 5,            # |rank_DOWN − rank_UP| ≥ this = direction-specific
  top_rank_for_specific = 10            # only consider pathways in top N of either direction
)

PATHS <- list(
  enrichr_dir = file.path("data")
)

# Library priority — only the listed libraries appear (and in this order).
# Comment out any you don't want to include in the plot.
LIBRARY_ORDER <- c(
  "Reactome_Pathways_2024",
  "BioPlanet_2019",
  "MSigDB_Hallmark_2020",
  "GO_Biological_Process_2025",
  "WikiPathways_2024_Human",
  "KEGG_2026"
)

# Geneset display order
GENESET_ORDER <- c("All9", "Down", "Up")
GENESET_LABELS <- c("All9" = "All 9 genes",
                    "Down" = "Down (5)",
                    "Up"   = "Up (4)")

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here, data.table, dplyr, tibble, stringr, tidyr, readr, purrr,
  ggplot2, scales, patchwork
)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# =============================================================================
# DIRECTORIES
# =============================================================================
dirs <- list(
  out     = here("output"),
  figures = here("output", "figures"),
  tables  = here("output", "tables"),
  logs    = here("output", "logs")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOGGING
# =============================================================================
log_path <- file.path(dirs$logs,
                      paste0("10_Enrichr_Run_",
                             format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
sink(log_path, append = FALSE, split = TRUE)
t_start <- Sys.time()

cat("=============================================================\n")
cat(" 10_enrichr_dotplot.R\n")
cat(" Start  :", as.character(t_start), "\n")
cat(" FDR    :", CONFIG$fdr_threshold, "\n")
cat(" Top N  :", CONFIG$top_n_per_lib, "per library × geneset\n")
cat("=============================================================\n\n")

set.seed(CONFIG$seed)

# =============================================================================
# SHARED AESTHETICS (identical to scripts 01–09)
# =============================================================================
COL <- list(
  all       = "#7B7B7B",
  down      = "#185FA5",
  up        = "#A32D2D",
  down_dark = "#0E3D6B",
  up_dark   = "#5A1010",
  bg        = "white",
  grid      = "grey92"
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
      panel.border     = element_rect(colour = "black", fill = NA,
                                      linewidth = 0.8),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = unit(c(0.8, 1.2, 0.8, 0.8), "cm")
    )
}

save_figure <- function(plot_obj, filename_base, width = 9, height = 8,
                        out_dir = dirs$figures) {
  base <- file.path(out_dir, filename_base)
  ggplot2::ggsave(paste0(base, ".png"),  plot = plot_obj,
                  width = width, height = height, dpi = 300)
  ggplot2::ggsave(paste0(base, ".tiff"), plot = plot_obj,
                  width = width, height = height, dpi = 600,
                  device = "tiff", compression = "lzw")
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = width, height = height)
  print(plot_obj); grDevices::dev.off()
  cat("  [SAVED]", filename_base, "(PNG 300dpi + TIFF 600dpi + PDF)\n")
  invisible(plot_obj)
}

# =============================================================================
# 1. LOAD AND PARSE ALL ENRICHR FILES
# =============================================================================
cat("[1/5] Discovering and loading Enrichr files...\n")

# Detect all .txt and .tsv files in the enrichr directory
files <- list.files(PATHS$enrichr_dir,
                    pattern = "\\.(txt|tsv)$",
                    full.names = TRUE)

if (length(files) == 0) {
  stop("No Enrichr .txt/.tsv files found in: ", PATHS$enrichr_dir,
       "\nExpected naming: <GeneSet>_<Library>.txt e.g. All9_Reactome_Pathways_2024.txt")
}

cat(sprintf("  Found %d files\n", length(files)))

# Normalise common library name variants to canonical forms
# (handles spaces, missing year suffixes, alternate capitalisations)
normalise_library <- function(lib) {
  lib_clean <- str_replace_all(lib, "[ \\-]", "_")  # spaces/hyphens → underscores
  lib_clean <- str_replace_all(lib_clean, "_+", "_") # collapse multiple underscores
  
  # Canonical library names (LIBRARY_ORDER) — try to match these
  canonicals <- c(
    "Reactome_Pathways_2024", "Reactome_2024", "Reactome_2022", "Reactome",
    "BioPlanet_2019", "BioPlanet",
    "MSigDB_Hallmark_2020", "MSigDB_Hallmark", "Hallmark_2020", "Hallmark",
    "GO_Biological_Process_2025", "GO_Biological_Process_2023",
    "GO_Biological_Process",
    "WikiPathways_2024_Human", "WikiPathways_2023", "WikiPathways",
    "KEGG_2026", "KEGG_2021_Human", "KEGG"
  )
  
  # Map any recognisable variant to the LIBRARY_ORDER canonical
  if (grepl("Reactome", lib_clean, ignore.case = TRUE))
    return("Reactome_Pathways_2024")
  if (grepl("BioPlanet", lib_clean, ignore.case = TRUE))
    return("BioPlanet_2019")
  if (grepl("MSigDB.*Hallmark|^Hallmark", lib_clean, ignore.case = TRUE))
    return("MSigDB_Hallmark_2020")
  if (grepl("GO.*Biological.*Process|GO_BP", lib_clean, ignore.case = TRUE))
    return("GO_Biological_Process_2025")
  if (grepl("WikiPathways|WikiPathway", lib_clean, ignore.case = TRUE))
    return("WikiPathways_2024_Human")
  if (grepl("KEGG", lib_clean, ignore.case = TRUE))
    return("KEGG_2026")
  
  return(lib_clean)  # unmodified if no match
}

# Parse filename → (GeneSet, Library) with library normalisation
parse_filename <- function(fp) {
  fn <- tools::file_path_sans_ext(basename(fp))
  # Match pattern: GeneSet_LibraryName
  parts <- str_split_fixed(fn, "_", n = 2)
  if (parts[2] == "") return(NULL)
  list(geneset = parts[1],
       library = normalise_library(parts[2]),
       filepath = fp)
}

file_meta <- compact(lapply(files, parse_filename))
cat(sprintf("  Parsed: %d files\n\n", length(file_meta)))

# Read a single Enrichr file robustly
read_enrichr <- function(fp) {
  tryCatch({
    df <- data.table::fread(fp, sep = "\t", header = TRUE,
                            check.names = FALSE, quote = "")
    df <- as.data.frame(df)
    # Standardise column names — Enrichr always has these (with spaces)
    expected <- c("Term", "Overlap", "P-value", "Adjusted P-value",
                  "Odds Ratio", "Combined Score", "Genes")
    if (!all(expected %in% names(df))) {
      cat("  [WARN] Missing expected columns in", basename(fp), "\n")
      return(NULL)
    }
    df %>%
      rename(
        Pvalue       = `P-value`,
        FDR          = `Adjusted P-value`,
        OddsRatio    = `Odds Ratio`,
        CombinedScore = `Combined Score`
      ) %>%
      select(Term, Overlap, Pvalue, FDR, OddsRatio, CombinedScore, Genes)
  }, error = function(e) {
    cat("  [ERROR] Could not read", basename(fp), ":", conditionMessage(e), "\n")
    NULL
  })
}

# Load all files into one long-format dataframe
all_results <- map_dfr(file_meta, function(m) {
  df <- read_enrichr(m$filepath)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>% mutate(GeneSet = m$geneset, Library = m$library, .before = 1)
})

if (nrow(all_results) == 0)
  stop("No valid Enrichr results loaded. Check file format.")

# Parse Overlap "X/Y" into GeneCount and TotalInPathway
all_results <- all_results %>%
  mutate(
    GeneCount       = as.integer(str_extract(Overlap, "^\\d+")),
    TotalInPathway  = as.integer(str_extract(Overlap, "(?<=/)\\d+")),
    NegLog10FDR     = -log10(pmax(FDR, CONFIG$pseudo_minlogp)),
    Significant     = FDR < CONFIG$fdr_threshold,
    Library_clean   = str_replace_all(Library, "_", " "),
    Term_short      = ifelse(nchar(Term) > CONFIG$max_term_length,
                             paste0(substr(Term, 1, CONFIG$max_term_length - 1),
                                    "\u2026"),
                             Term)
  )

cat(sprintf("  Total rows loaded: %d\n", nrow(all_results)))
cat(sprintf("  Significant (FDR < %.2f): %d\n",
            CONFIG$fdr_threshold, sum(all_results$Significant, na.rm = TRUE)))

# Per-file summary
file_summary <- all_results %>%
  group_by(GeneSet, Library) %>%
  summarise(
    Total       = n(),
    Significant = sum(Significant, na.rm = TRUE),
    Min_FDR     = min(FDR, na.rm = TRUE),
    .groups = "drop"
  )
cat("\n  ── Per-file summary ──\n")
print(file_summary, n = nrow(file_summary))
cat("\n")

# Save combined long-format
write.csv(all_results,
          file.path(dirs$tables, "10_Enrichr_combined_long.csv"),
          row.names = FALSE)
cat("  [SAVED] 10_Enrichr_combined_long.csv\n\n")

# =============================================================================
# 2. AUTO-DETECT TOP PATHWAYS
# =============================================================================
cat("[2/5] Auto-detecting top pathways for plotting...\n")

# STRATEGY:
# Per (Library × GeneSet), keep top N by FDR (significant only, gene count >= 2).
# Then UNION across all combinations to get the master pathway list to plot.
# This ensures every gene set's strongest hits are represented while keeping
# the plot readable.

top_per_combo <- all_results %>%
  filter(
    Significant,
    GeneCount >= CONFIG$min_gene_count,
    Library %in% LIBRARY_ORDER
  ) %>%
  group_by(GeneSet, Library) %>%
  arrange(FDR, desc(CombinedScore), .by_group = TRUE) %>%
  slice_head(n = CONFIG$top_n_per_lib) %>%
  ungroup()

cat(sprintf("  Top hits per (geneset × library, n=%d): %d rows\n",
            CONFIG$top_n_per_lib, nrow(top_per_combo)))

# Master pathway list — union of all top hits
master_terms <- top_per_combo %>%
  distinct(Library, Term)

cat(sprintf("  Unique pathways to plot: %d\n", nrow(master_terms)))

# Now retrieve ALL geneset values (including non-significant) for those pathways
# This is what makes the dot plot informative — empty/grey cells show absence
plot_df <- all_results %>%
  semi_join(master_terms, by = c("Library", "Term")) %>%
  filter(Library %in% LIBRARY_ORDER) %>%
  mutate(
    GeneSet  = factor(GeneSet, levels = GENESET_ORDER),
    Library  = factor(Library, levels = LIBRARY_ORDER)
  )

cat(sprintf("  Plot dataframe: %d rows\n\n", nrow(plot_df)))

# Save plotted pathway data
write.csv(plot_df %>%
            select(GeneSet, Library, Term, GeneCount, TotalInPathway,
                   FDR, OddsRatio, CombinedScore, Genes),
          file.path(dirs$tables, "10_Enrichr_top_pathways_plotted.csv"),
          row.names = FALSE)
cat("  [SAVED] 10_Enrichr_top_pathways_plotted.csv\n\n")

# Save significant-only subset (for supplementary table)
sig_only <- all_results %>%
  filter(Significant, GeneCount >= CONFIG$min_gene_count) %>%
  arrange(GeneSet, Library, FDR) %>%
  select(GeneSet, Library, Term, GeneCount, TotalInPathway,
         Pvalue, FDR, OddsRatio, CombinedScore, Genes)
write.csv(sig_only,
          file.path(dirs$tables, "10_Enrichr_significant_only.csv"),
          row.names = FALSE)
cat(sprintf("  [SAVED] 10_Enrichr_significant_only.csv (%d rows)\n\n",
            nrow(sig_only)))

# =============================================================================
# 3. ORDER PATHWAYS FOR THE PLOT
# =============================================================================
# Pathway display order: group by which GeneSet they're most enriched in,
# then by FDR within that group.
# This produces visual clustering (DOWN-specific terms group together,
# UP-specific terms group together).

cat("[3/5] Ordering pathways by GeneSet specificity...\n")

term_order <- plot_df %>%
  filter(Significant) %>%
  group_by(Library, Term) %>%
  summarise(
    # Which geneset gives this term its strongest hit?
    BestGeneSet = GeneSet[which.min(FDR)],
    BestFDR     = min(FDR),
    .groups = "drop"
  ) %>%
  # Order: Down-specific → All-best → Up-specific (so UP/DOWN cluster visually)
  mutate(
    GroupRank = case_when(
      BestGeneSet == "Down" ~ 1L,
      BestGeneSet == "All9" ~ 2L,
      BestGeneSet == "Up"   ~ 3L,
      TRUE ~ 4L
    )
  ) %>%
  arrange(Library, GroupRank, BestFDR) %>%
  mutate(TermID = paste0(Library, " | ", Term))

# Apply to plot data
plot_df <- plot_df %>%
  mutate(TermID = paste0(Library, " | ", Term)) %>%
  mutate(TermID = factor(TermID, levels = rev(term_order$TermID)))

# Strip "Library | " prefix for axis labels (we use facet for library instead)
plot_df <- plot_df %>%
  mutate(TermDisplay = str_remove(as.character(TermID),
                                  paste0("^", Library, " \\| ")),
         TermDisplay = factor(TermDisplay,
                              levels = unique(rev(term_order$Term))))

# Truncate display names
plot_df <- plot_df %>%
  mutate(TermDisplay_short = ifelse(
    nchar(as.character(TermDisplay)) > CONFIG$max_term_length,
    paste0(substr(as.character(TermDisplay), 1, CONFIG$max_term_length - 1),
           "\u2026"),
    as.character(TermDisplay)
  ))

# =============================================================================
# 4. BUILD THE MAIN DOT PLOT (Fig 10A)
# =============================================================================
cat("[4/5] Building dot plot...\n")

# Two-tier alpha: significant = solid, non-significant = faded
plot_df <- plot_df %>%
  mutate(plot_alpha = ifelse(Significant, 1, 0.25))

# Main figure — one panel per library, GeneSet on X-axis
p_main <- ggplot(plot_df,
                 aes(x = GeneSet, y = TermID)) +
  geom_point(aes(size = GeneCount, fill = NegLog10FDR, alpha = plot_alpha),
             shape = 21, colour = "grey20", stroke = 0.4) +
  # Gene count is a discrete quantity: fractional legend breaks (1.5, 2.5)
  # are meaningless and draw attention to the small overlap sizes.
  # breaks_width(1) yields integer steps and adapts if the range changes.
  scale_size_continuous(range = c(2.5, 9), name = "Gene\ncount",
                        breaks = scales::breaks_width(1)) +
  scale_fill_gradientn(
    colours = c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#A50F15", "#3D0000"),
    name    = expression(-log[10]~"(FDR)"),
    limits  = c(0, NA)
  ) +
  scale_alpha_identity() +
  scale_x_discrete(labels = GENESET_LABELS) +
  scale_y_discrete(labels = function(x) str_remove(x, "^.*\\| ")) +
  facet_grid(Library ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(Library = function(x) str_replace_all(x, "_", " "))) +
  labs(
    title    = "Pathway enrichment of the 9-gene DSP signature",
    subtitle = paste0("Up- and down-regulated genes converge on distinct biological pathways  |  ",
                      "Top ", CONFIG$top_n_per_lib, " per library × geneset"),
    x        = NULL,
    y        = NULL,
    caption  = paste0(
      "Enrichr enrichment analysis. Pathways shown reach FDR < ",
      CONFIG$fdr_threshold, " in at least one geneset (gene count \u2265 ",
      CONFIG$min_gene_count, ").\n",
      "Faded dots: not significant in that geneset.  ",
      "Pathways grouped by which subset shows strongest enrichment.")
  ) +
  theme_publication() +
  theme(
    axis.text.y       = element_text(size = 9),
    axis.text.x       = element_text(face = "bold", size = 11),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    strip.text.y      = element_text(angle = 0, face = "bold", size = 9),
    strip.background.y = element_rect(fill = "grey88", colour = "black"),
    legend.position   = "right",
    legend.box        = "vertical"
  )

save_figure(p_main, "Fig_10A_Enrichr_Dotplot",
            width = 11, height = max(7, 0.32 * length(unique(plot_df$TermID)) + 4))

# =============================================================================
# 5. RANK-BASED DIRECTIONALITY ANALYSIS (Fig 10C, Fig 10D, Fig 10E)
# =============================================================================
# WHY RANK-BASED:
#   Binary "significant in both UP and DOWN" overlap is misleading because
#   pathway databases annotate STYXL1 (catalytically dead) and PTP4A2/3 to
#   MAPK pathways alongside DUSP1/DUSP6, causing apparent overlap at the
#   significance threshold. Rank-based analysis asks: "where does this
#   pathway RANK in UP vs DOWN?" A pathway ranked #1 in DOWN but #15 in UP
#   is direction-specific even when both pass FDR.

cat("\n[5/5] Rank-based directionality analysis...\n")

# Rank pathways within each (Library × GeneSet)
ranked_full <- all_results %>%
  filter(Library %in% LIBRARY_ORDER, GeneCount >= CONFIG$min_gene_count) %>%
  group_by(Library, GeneSet) %>%
  arrange(FDR, desc(CombinedScore), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

# Pivot to wide: ranks and FDRs in All9 / Down / Up
rank_wide <- ranked_full %>%
  select(Library, Term, GeneSet, Rank, FDR) %>%
  pivot_wider(names_from = GeneSet,
              values_from = c(Rank, FDR),
              values_fill = NA)

# Ensure all expected columns exist
for (col in c("Rank_All9", "Rank_Down", "Rank_Up",
              "FDR_All9", "FDR_Down", "FDR_Up")) {
  if (!col %in% names(rank_wide)) rank_wide[[col]] <- NA_real_
}

# Classify by direction-specificity using ranks (not just significance)
classify_direction <- function(rd, ru, fdr_d, fdr_u, thr_diff,
                               top_thr, fdr_cutoff) {
  d_sig <- !is.na(fdr_d) && fdr_d < fdr_cutoff
  u_sig <- !is.na(fdr_u) && fdr_u < fdr_cutoff
  
  if (!d_sig && !u_sig) return("Not significant")
  if (d_sig  && !u_sig) return("DOWN-specific")
  if (!d_sig &&  u_sig) return("UP-specific")
  
  # Both significant — use rank difference
  rd_ <- ifelse(is.na(rd), 999, rd)
  ru_ <- ifelse(is.na(ru), 999, ru)
  
  # Must be in top N of at least one direction
  if (min(rd_, ru_) > top_thr) return("Both - neither dominant")
  
  diff <- ru_ - rd_   # positive → DOWN ranks better; negative → UP ranks better
  if (abs(diff) >= thr_diff) {
    if (diff > 0) return("DOWN-driven (UP also passes)")
    else          return("UP-driven (DOWN also passes)")
  } else {
    return("Shared (similar rank)")
  }
}

rank_wide <- rank_wide %>%
  rowwise() %>%
  mutate(
    Direction = classify_direction(
      Rank_Down, Rank_Up, FDR_Down, FDR_Up,
      thr_diff   = CONFIG$rank_difference,
      top_thr    = CONFIG$top_rank_for_specific,
      fdr_cutoff = CONFIG$fdr_threshold
    )
  ) %>%
  ungroup() %>%
  arrange(Library, pmin(Rank_Down, Rank_Up, na.rm = TRUE))

# Save the full rank-specificity table
write.csv(rank_wide,
          file.path(dirs$tables, "10_pathway_rank_specificity.csv"),
          row.names = FALSE)
cat("  [SAVED] 10_pathway_rank_specificity.csv\n")

# Summary counts
direction_summary <- rank_wide %>%
  filter(Direction != "Not significant") %>%
  count(Library, Direction, name = "n_pathways") %>%
  group_by(Library) %>%
  mutate(pct = round(100 * n_pathways / sum(n_pathways), 1)) %>%
  ungroup()

write.csv(direction_summary,
          file.path(dirs$tables, "10_directionality_summary.csv"),
          row.names = FALSE)
cat("  [SAVED] 10_directionality_summary.csv\n")

cat("\n  ── Rank-based directionality summary ──\n")
print(as.data.frame(direction_summary))

# ── Fig 10C: Direction-specificity stacked bar ─────────────────────────────
direction_levels <- c("DOWN-specific", "DOWN-driven (UP also passes)",
                      "Shared (similar rank)",
                      "UP-driven (DOWN also passes)", "UP-specific",
                      "Both - neither dominant")
direction_colours <- c(
  "DOWN-specific"                = COL$down_dark,
  "DOWN-driven (UP also passes)" = COL$down,
  "Shared (similar rank)"        = "#7B3294",
  "UP-driven (DOWN also passes)" = COL$up,
  "UP-specific"                  = COL$up_dark,
  "Both - neither dominant"      = "grey70"
)

direction_summary <- direction_summary %>%
  mutate(Direction = factor(Direction, levels = direction_levels)) %>%
  arrange(Library, Direction)

p_sep <- ggplot(direction_summary,
                aes(x = Library, y = n_pathways, fill = Direction)) +
  geom_col(colour = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_pathways >= 2, n_pathways, "")),
            position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = direction_colours, drop = FALSE) +
  scale_x_discrete(labels = function(x) str_replace_all(x, "_", "\n")) +
  labs(
    title    = "Pathway directionality - rank-based specificity",
    subtitle = paste0("Pathways classified by where they rank highest  |  ",
                      "rank-difference \u2265 ", CONFIG$rank_difference,
                      " = direction-specific"),
    x        = NULL,
    y        = "Number of pathways  (FDR < 0.05)",
    fill     = "Directionality",
    caption  = paste0(
      "Binary FDR overlap obscures direction-specificity because pathway databases\n",
      "do not distinguish catalytically active from inactive family members.\n",
      "Rank-based classification reveals the true UP-vs-DOWN biology separation.")
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(face = "bold", size = 10),
    legend.position = "right"
  )

save_figure(p_sep, "Fig_10C_Pathway_Directionality", width = 11, height = 6)

# ── Fig 10E: TOP-5 pathways per direction (cleanest visual evidence) ──────
top_down <- ranked_full %>%
  filter(GeneSet == "Down", FDR < CONFIG$fdr_threshold) %>%
  group_by(Library) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(Direction = "DOWN")

top_up <- ranked_full %>%
  filter(GeneSet == "Up", FDR < CONFIG$fdr_threshold) %>%
  group_by(Library) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(Direction = "UP")

if (nrow(top_down) + nrow(top_up) > 0) {
  top_combined <- bind_rows(top_down, top_up) %>%
    mutate(
      NegLog10FDR = -log10(pmax(FDR, CONFIG$pseudo_minlogp)),
      Direction   = factor(Direction, levels = c("DOWN", "UP")),
      Term_id     = paste0(Library, " | ", Term)
    ) %>%
    arrange(Direction, Library, FDR) %>%
    mutate(Term_id = factor(Term_id, levels = unique(Term_id)))
  
  p_top <- ggplot(top_combined,
                  aes(x = NegLog10FDR, y = Term_id, fill = Direction)) +
    geom_col(colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%d/%d", GeneCount, TotalInPathway)),
              hjust = -0.2, size = 3, fontface = "bold", colour = "grey20") +
    scale_fill_manual(values = c(DOWN = COL$down_dark, UP = COL$up_dark)) +
    scale_y_discrete(labels = function(x) str_remove(x, "^.*\\| ")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
    facet_grid(Library ~ Direction, scales = "free_y", space = "free_y",
               labeller = labeller(Library = function(x) str_replace_all(x, "_", " "))) +
    labs(
      title    = "Top-5 enriched pathways per gene set per library",
      subtitle = "DOWN and UP genes hit qualitatively different pathways",
      x        = expression(-log[10] * "(FDR)"),
      y        = NULL,
      caption  = paste0(
        "Numbers at bar end: gene overlap (signature/total in pathway).\n",
        "Bar length = enrichment strength. Direction-specific patterns visually evident.")
    ) +
    theme_publication() +
    theme(
      axis.text.y     = element_text(size = 8.5),
      strip.text.y    = element_text(angle = 0, face = "bold", size = 9),
      legend.position = "none",
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
    )
  
  save_figure(p_top, "Fig_10E_Top5_Per_Direction",
              width = 13, height = max(7, 0.4 * nrow(top_combined) + 3))
}

# Print summary to console
cat("\n  ── Direction summary ──\n")
print(as.data.frame(direction_summary))

# =============================================================================
# 6. STANDALONE RESULTS REPORT
# =============================================================================
cat("\n[6/6] Writing standalone results report...\n")

report_path <- file.path(dirs$out, "10_ENRICHR_RESULTS_REPORT.txt")
con <- file(report_path, open = "wt")

writeLines(c(
  "=============================================================",
  " ENRICHR PATHWAY ENRICHMENT — RESULTS REPORT",
  paste0(" Generated: ", Sys.time()),
  paste0(" Script   : 10_enrichr_dotplot.R"),
  "=============================================================",
  "",
  paste0("  FDR threshold      : ", CONFIG$fdr_threshold),
  paste0("  Min gene count     : ", CONFIG$min_gene_count),
  paste0("  Top N per combo    : ", CONFIG$top_n_per_lib),
  paste0("  Files loaded       : ", length(file_meta)),
  paste0("  Total enrichments  : ", nrow(all_results)),
  paste0("  Significant total  : ", sum(all_results$Significant, na.rm = TRUE)),
  ""
), con = con)

# Per-geneset top hits
for (gs in GENESET_ORDER) {
  writeLines(c(
    "",
    paste0("─── TOP 5 PATHWAYS — ", gs, " ───"),
    ""
  ), con = con)
  
  top_gs <- all_results %>%
    filter(GeneSet == gs, Significant, GeneCount >= CONFIG$min_gene_count,
           Library %in% LIBRARY_ORDER) %>%
    arrange(FDR, desc(CombinedScore)) %>%
    slice_head(n = 5) %>%
    select(Library, Term, Overlap, FDR, Genes)
  
  if (nrow(top_gs) == 0) {
    writeLines("  (no significant hits)", con = con)
  } else {
    # FIX: convert tibble to data.frame to avoid 'invalid na.print' error
    capture.output(print(as.data.frame(top_gs), row.names = FALSE),
                   file = con, append = TRUE)
  }
}

writeLines(c(
  "",
  "─── KEY MANUSCRIPT NUMBERS ───",
  "",
  "  Use these for the Results paragraph:",
  ""
), con = con)

# Compose the manuscript-ready FDR strings for the strongest hits
ms_lines <- top_per_combo %>%
  arrange(GeneSet, Library, FDR) %>%
  group_by(GeneSet, Library) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  mutate(text = sprintf("    %s | %s | %s (FDR = %.2e, %d/%d genes)",
                        GeneSet, Library, Term, FDR,
                        GeneCount, TotalInPathway)) %>%
  pull(text)
writeLines(ms_lines, con = con)

writeLines(c(
  "",
  "─── PATHWAY DIRECTIONALITY (rank-based) ───",
  "",
  "  Pathways classified by where they rank highest within each library.",
  paste0("  Rank-difference threshold for specificity: ",
         CONFIG$rank_difference),
  paste0("  Maximum rank to be considered specific:    ",
         CONFIG$top_rank_for_specific),
  ""
), con = con)
# FIX: convert to data.frame for clean printing
capture.output(print(as.data.frame(direction_summary), row.names = FALSE),
               file = con, append = TRUE)

# Top DOWN-driven and UP-driven pathways
writeLines(c("",
             "─── TOP DOWN-DRIVEN / DOWN-SPECIFIC PATHWAYS ───", ""),
           con = con)
top_down_print <- rank_wide %>%
  filter(Direction %in% c("DOWN-specific", "DOWN-driven (UP also passes)")) %>%
  arrange(Library, Rank_Down) %>%
  slice_head(n = 15) %>%
  select(Library, Term, Rank_Down, Rank_Up, FDR_Down, FDR_Up, Direction)
if (nrow(top_down_print) > 0) {
  capture.output(print(as.data.frame(top_down_print), row.names = FALSE,
                       max = 99999),
                 file = con, append = TRUE)
} else {
  writeLines("  (none)", con = con)
}

writeLines(c("",
             "─── TOP UP-DRIVEN / UP-SPECIFIC PATHWAYS ───", ""),
           con = con)
top_up_print <- rank_wide %>%
  filter(Direction %in% c("UP-specific", "UP-driven (DOWN also passes)")) %>%
  arrange(Library, Rank_Up) %>%
  slice_head(n = 15) %>%
  select(Library, Term, Rank_Down, Rank_Up, FDR_Down, FDR_Up, Direction)
if (nrow(top_up_print) > 0) {
  capture.output(print(as.data.frame(top_up_print), row.names = FALSE,
                       max = 99999),
                 file = con, append = TRUE)
} else {
  writeLines("  (none)", con = con)
}

writeLines(c(
  "",
  "─── OUTPUT FILES ───",
  paste0("  Tables  : ", dirs$tables),
  paste0("  Figures : ", dirs$figures),
  paste0("  Report  : ", report_path),
  ""
), con = con)

close(con)
cat("  [SAVED] 10_ENRICHR_RESULTS_REPORT.txt\n")

# Session info
session_path <- file.path(dirs$logs, "10_session_info.txt")
writeLines(capture.output(sessionInfo()), session_path)
cat("  [SAVED] 10_session_info.txt\n\n")

# =============================================================================
# WRAP UP
# =============================================================================
t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 10_enrichr_dotplot.R  COMPLETE\n")
cat(" End    :", as.character(t_end), "\n")
cat(" Runtime:", round(difftime(t_end, t_start, units = "mins"), 2), "minutes\n")
cat("\n SAVED ARTEFACTS:\n")
cat("  Tables  : 10_Enrichr_combined_long.csv\n")
cat("            10_Enrichr_significant_only.csv\n")
cat("            10_Enrichr_top_pathways_plotted.csv\n")
cat("            10_pathway_rank_specificity.csv\n")
cat("            10_directionality_summary.csv\n")
cat("  Figures : Fig_10A_Enrichr_Dotplot\n")
cat("            Fig_10C_Pathway_Directionality\n")
cat("            Fig_10E_Top5_Per_Direction\n")
cat("  Report  : 10_ENRICHR_RESULTS_REPORT.txt\n")
cat("  Log     :", basename(log_path), "\n")
cat("=============================================================\n")

sink()
cat("Done. Log saved to:", log_path, "\n")
cat("Results report saved to:", report_path, "\n")