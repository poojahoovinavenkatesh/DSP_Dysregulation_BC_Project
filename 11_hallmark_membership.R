# =============================================================================
# Script      : 11_hallmark_membership.R
# Project     : DSP family dysregulation in breast cancer
#
# PURPOSE
#   Document which MSigDB Hallmark gene sets contain each signature gene.
#
# WHY THIS IS NOT AN ENRICHMENT ANALYSIS
#   Over-representation analysis (Enrichr, Script 10) asks a SET-level
#   question: is the overlap between the gene set and a pathway larger than
#   expected by chance? With a nine-gene input split into subsets of four and
#   five, most overlaps are one or two genes, and the resulting FDR is close
#   to a statement about relative set sizes rather than about biology.
#
#   This script asks a GENE-level question instead: which Hallmark processes
#   is each gene annotated to? That is a membership lookup, not a statistical
#   inference. No p-values are computed and none are appropriate — the answer
#   is a fact about the MSigDB annotation, not an estimate.
#
# INTERPRETIVE LIMIT — state this in the manuscript
#   Hallmark sets were derived by identifying genes co-expressed across many
#   tissues and conditions (Liberzon et al., 2015). Membership therefore
#   indicates that a gene's expression covaries with a hallmark programme in
#   the source data. It does NOT establish that the gene drives that process,
#   nor that the association holds in breast cancer specifically. Membership
#   counts describe annotation breadth, not causal involvement.
#
# OUTPUTS
#   Tables:  16_hallmark_membership_matrix.csv    gene x hallmark, 0/1
#            16_hallmark_membership_long.csv      one row per gene-hallmark
#            16_gene_hallmark_counts.csv          hallmarks per gene
#            16_hallmark_gene_counts.csv          signature genes per hallmark
#   Figures: Fig_16A_Hallmark_Membership_Tile
#            Fig_16B_Hallmarks_Per_Gene
#   Report:  16_HALLMARK_MEMBERSHIP_REPORT.txt
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# =============================================================================

rm(list = ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed = 42L,
  
  # Show only hallmarks containing at least one signature gene.
  # FALSE would render all 50, mostly empty.
  drop_empty_hallmarks = TRUE,
  
  # Order hallmarks on the tile plot by how many signature genes they
  # contain. FALSE keeps alphabetical order.
  order_by_count = TRUE
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
pacman::p_load(here, dplyr, tidyr, tibble, stringr, ggplot2, scales)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
library(msigdbr)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

dirs <- list(
  out     = here("output_hallmark_membership"),
  figures = here("output_hallmark_membership", "figures"),
  tables  = here("output_hallmark_membership", "tables"),
  logs    = here("output_hallmark_membership", "logs")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
for (nm in names(dirs))
  if (!dir.exists(dirs[[nm]]))
    stop("Could not create directory '", nm, "' at: ", dirs[[nm]])

cat("Output directories:\n")
for (nm in names(dirs))
  cat(sprintf("  %-8s %s\n", nm, normalizePath(dirs[[nm]], winslash = "/")))
cat("\n")

log_path <- file.path(dirs$logs,
                      paste0("16_HallmarkMembership_",
                             format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
sink(log_path, append = FALSE, split = TRUE)
t_start <- Sys.time()
set.seed(CONFIG$seed)

cat("=============================================================\n")
cat(" 16_hallmark_membership.R\n")
cat(" Start :", as.character(t_start), "\n")
cat("=============================================================\n\n")

# =============================================================================
# AESTHETICS (matching all project scripts)
# =============================================================================
COL <- list(
  up      = "#A32D2D",
  down    = "#185FA5",
  present = "#4D4D4D",
  absent  = "grey94"
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

save_table <- function(x, filename, dir = dirs$tables) {
  fp <- file.path(dir, filename)
  ok <- tryCatch({ utils::write.csv(x, fp, row.names = FALSE); TRUE },
                 error = function(e) {
                   cat("  [ERROR] Could not write", filename, ":",
                       conditionMessage(e), "\n"); FALSE })
  if (ok && file.exists(fp))
    cat(sprintf("  [SAVED] %-42s (%d rows) -> %s\n", filename,
                nrow(as.data.frame(x)), normalizePath(fp, winslash = "/")))
  invisible(fp)
}

# =============================================================================
# 1. LOAD SIGNATURE AND DIRECTIONS
# =============================================================================
cat("[1/4] Loading locked signature...\n")

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

gene_info <- data.frame(
  gene      = c(up_genes, down_genes),
  direction = c(rep("Up", length(up_genes)), rep("Down", length(down_genes))),
  stringsAsFactors = FALSE
)

cat(sprintf("  Signature: %d genes (%d up, %d down)\n",
            nrow(gene_info), length(up_genes), length(down_genes)))
cat("  Up   :", paste(up_genes,   collapse = ", "), "\n")
cat("  Down :", paste(down_genes, collapse = ", "), "\n\n")

# =============================================================================
# 2. RETRIEVE HALLMARK GENE SETS
# =============================================================================
cat("[2/4] Retrieving MSigDB Hallmark collection...\n")

hallmark <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "H"),
  error = function(e) {
    # msigdbr >= 10 renamed the argument
    msigdbr(species = "Homo sapiens", collection = "H")
  })

hallmark <- as.data.frame(hallmark)

# Column names differ between msigdbr versions
sym_col <- intersect(c("gene_symbol", "human_gene_symbol"),
                     colnames(hallmark))[1]
set_col <- intersect(c("gs_name", "gs_id"), colnames(hallmark))[1]
if (is.na(sym_col) || is.na(set_col))
  stop("Unexpected msigdbr columns: ", paste(colnames(hallmark), collapse = ", "))

hallmark <- hallmark %>%
  select(hallmark = all_of(set_col), symbol = all_of(sym_col)) %>%
  distinct()

set_sizes <- hallmark %>% count(hallmark, name = "set_size")

cat(sprintf("  Hallmark sets: %d\n", nrow(set_sizes)))
cat(sprintf("  Set size: median %d, range %d-%d\n",
            median(set_sizes$set_size), min(set_sizes$set_size),
            max(set_sizes$set_size)))
cat(sprintf("  msigdbr version: %s\n\n",
            as.character(utils::packageVersion("msigdbr"))))

# =============================================================================
# 3. MEMBERSHIP LOOKUP
# =============================================================================
cat("[3/4] Building membership matrix...\n")

membership_long <- hallmark %>%
  filter(symbol %in% clean_genes) %>%
  rename(gene = symbol) %>%
  left_join(gene_info, by = "gene") %>%
  left_join(set_sizes, by = "hallmark") %>%
  mutate(hallmark_label = str_to_title(str_replace_all(
    str_remove(hallmark, "^HALLMARK_"), "_", " ")))

if (nrow(membership_long) == 0)
  stop("No signature gene appears in any Hallmark set — check gene symbols")

save_table(membership_long %>%
             select(gene, direction, hallmark, hallmark_label, set_size) %>%
             arrange(gene, hallmark),
           "16_hallmark_membership_long.csv")

# ── Per-gene counts ────────────────────────────────────────────────────────
gene_counts <- gene_info %>%
  left_join(membership_long %>% count(gene, name = "n_hallmarks"),
            by = "gene") %>%
  mutate(n_hallmarks = tidyr::replace_na(n_hallmarks, 0L)) %>%
  arrange(desc(n_hallmarks))

save_table(gene_counts, "16_gene_hallmark_counts.csv")

cat("\n  ── Hallmark sets containing each signature gene ──\n")
print(as.data.frame(gene_counts), row.names = FALSE)

n_multi <- sum(gene_counts$n_hallmarks >= 2)
cat(sprintf("\n  Genes in >= 2 hallmark sets: %d of %d\n",
            n_multi, nrow(gene_counts)))
cat(sprintf("  Genes in 0 hallmark sets   : %d\n",
            sum(gene_counts$n_hallmarks == 0)))

# ── Per-hallmark counts ────────────────────────────────────────────────────
hallmark_counts <- membership_long %>%
  count(hallmark, hallmark_label, set_size, name = "n_signature_genes") %>%
  mutate(pct_of_set = round(100 * n_signature_genes / set_size, 2)) %>%
  arrange(desc(n_signature_genes))

save_table(hallmark_counts, "16_hallmark_gene_counts.csv")

cat("\n  ── Hallmark sets containing the most signature genes ──\n")
print(as.data.frame(head(hallmark_counts %>%
                           select(hallmark_label, n_signature_genes, set_size, pct_of_set), 12)),
      row.names = FALSE)

# ── Wide 0/1 matrix ────────────────────────────────────────────────────────
membership_wide <- membership_long %>%
  select(gene, hallmark_label) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = gene, values_from = present, values_fill = 0L) %>%
  arrange(hallmark_label)

save_table(membership_wide, "16_hallmark_membership_matrix.csv")
cat("\n")

# =============================================================================
# 4. FIGURES
# =============================================================================
cat("[4/4] Generating figures...\n")

gene_order <- c(gene_info$gene[gene_info$direction == "Up"],
                gene_info$gene[gene_info$direction == "Down"])
gene_cols  <- ifelse(gene_order %in% up_genes, COL$up, COL$down)

hm_order <- if (isTRUE(CONFIG$order_by_count)) {
  hallmark_counts$hallmark_label[order(hallmark_counts$n_signature_genes,
                                       decreasing = TRUE)]
} else {
  sort(unique(membership_long$hallmark_label))
}

# Complete grid so absent combinations render as empty tiles
tile_df <- expand.grid(
  gene           = gene_order,
  hallmark_label = hm_order,
  stringsAsFactors = FALSE) %>%
  left_join(membership_long %>%
              select(gene, hallmark_label) %>%
              mutate(present = TRUE),
            by = c("gene", "hallmark_label")) %>%
  mutate(
    present   = tidyr::replace_na(present, FALSE),
    direction = ifelse(gene %in% up_genes, "Up", "Down"),
    fill_val  = dplyr::case_when(
      !present            ~ "Absent",
      direction == "Up"   ~ "Up-regulated",
      TRUE                ~ "Down-regulated"),
    gene           = factor(gene, levels = gene_order),
    hallmark_label = factor(hallmark_label, levels = rev(hm_order)))

fill_cols <- c("Up-regulated"   = COL$up,
               "Down-regulated" = COL$down,
               "Absent"         = COL$absent)

# ── Fig 16A: membership tile ───────────────────────────────────────────────
p_tile <- ggplot(tile_df, aes(x = gene, y = hallmark_label, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  scale_fill_manual(values = fill_cols, name = NULL,
                    breaks = c("Up-regulated", "Down-regulated")) +
  scale_x_discrete(position = "top") +
  labs(
    title    = "Hallmark process annotation of the DSP signature genes",
    subtitle = "Individual genes are annotated to multiple cancer hallmark processes",
    x = NULL, y = NULL,
    caption = paste0(
      "Filled tiles indicate that the gene is a member of that MSigDB ",
      "Hallmark gene set. This is an annotation\n",
      "lookup, not an enrichment test; no significance values apply. ",
      "Hallmark sets were defined by\n",
      "co-expression across tissues and conditions, so membership indicates ",
      "coordinated expression with a\n",
      "hallmark programme rather than causal involvement in that process.")
  ) +
  theme_publication() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 0, face = "bold.italic",
                                   colour = gene_cols, size = 11),
    axis.text.y     = element_text(size = 9),
    panel.border    = element_blank(),
    legend.position = "bottom")

save_figure(p_tile, "Fig_16A_Hallmark_Membership_Tile",
            w = 9, h = max(6, 0.28 * length(hm_order) + 3))

# ── Fig 16B: hallmarks per gene ────────────────────────────────────────────
p_count <- gene_counts %>%
  mutate(gene = factor(gene, levels = gene_counts$gene[
    order(gene_counts$n_hallmarks)])) %>%
  ggplot(aes(x = n_hallmarks, y = gene, fill = direction)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.72) +
  geom_text(aes(label = n_hallmarks), hjust = -0.4, size = 3.6,
            fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = c(Up = COL$up, Down = COL$down),
                    name = "Direction\nin tumour",
                    labels = c(Up = "Up-regulated", Down = "Down-regulated")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)),
                     breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title    = "Number of cancer hallmark processes per signature gene",
    subtitle = sprintf("%d of %d genes are annotated to two or more hallmark sets",
                       n_multi, nrow(gene_counts)),
    x = "Number of MSigDB Hallmark gene sets", y = NULL,
    caption = paste0(
      "Counts reflect annotation breadth in MSigDB, not the number of ",
      "processes a gene has been shown to drive.\n",
      "Genes belonging to many sets tend to be broadly co-expressed rather ",
      "than functionally pleiotropic.")
  ) +
  theme_publication() +
  theme(axis.text.y = element_text(face = "bold.italic"))

save_figure(p_count, "Fig_16B_Hallmarks_Per_Gene", w = 9, h = 6)

# =============================================================================
# REPORT
# =============================================================================
report <- file.path(dirs$out, "16_HALLMARK_MEMBERSHIP_REPORT.txt")
con <- file(report, open = "wt")

writeLines(c(
  "=============================================================",
  " MSigDB HALLMARK MEMBERSHIP OF THE DSP SIGNATURE",
  paste0(" Generated: ", Sys.time()),
  "=============================================================",
  "",
  " This is an ANNOTATION LOOKUP, not an enrichment analysis.",
  " No p-values or FDR values are computed, and none apply.",
  "",
  paste0(" Hallmark sets in collection : ", nrow(set_sizes)),
  paste0(" msigdbr version             : ",
         as.character(utils::packageVersion("msigdbr"))),
  paste0(" Signature genes             : ", nrow(gene_info)),
  ""
), con = con)

writeLines(c("", "─── HALLMARK SETS PER GENE ───", ""), con = con)
capture.output(print(as.data.frame(gene_counts), row.names = FALSE),
               file = con, append = TRUE)

writeLines(c("",
             sprintf("  Genes annotated to >= 2 hallmark sets: %d of %d",
                     n_multi, nrow(gene_counts)), ""), con = con)

writeLines(c("", "─── SIGNATURE GENES PER HALLMARK SET ───", ""), con = con)
capture.output(print(as.data.frame(
  hallmark_counts %>% select(hallmark_label, n_signature_genes,
                             set_size, pct_of_set)),
  row.names = FALSE), file = con, append = TRUE)

writeLines(c("", "─── FULL GENE-HALLMARK PAIRS ───", ""), con = con)
capture.output(print(as.data.frame(
  membership_long %>% select(gene, direction, hallmark_label, set_size) %>%
    arrange(gene, hallmark_label)), row.names = FALSE),
  file = con, append = TRUE)

writeLines(c(
  "", "─── INTERPRETIVE LIMITS ───", "",
  "  1. Membership is not enrichment. This analysis reports which sets",
  "     contain each gene. It does not test whether the signature is",
  "     over-represented in any set, and no significance value applies.",
  "",
  "  2. Membership is not causal involvement. Hallmark sets were derived",
  "     by identifying genes co-expressed across many tissues and",
  "     conditions. A gene's presence indicates that its expression",
  "     covaries with the hallmark programme in those source data, not",
  "     that it drives the process.",
  "",
  "  3. Membership is not breast-cancer-specific. Hallmark annotation is",
  "     pan-tissue. Involvement in breast cancer specifically requires",
  "     separate evidence.",
  "",
  "  4. High counts may reflect broad co-expression. Genes appearing in",
  "     many sets are often those with wide expression variation across",
  "     conditions rather than genuinely pleiotropic function.",
  ""
), con = con)
close(con)

cat("  [SAVED] 16_HALLMARK_MEMBERSHIP_REPORT.txt\n")
writeLines(capture.output(sessionInfo()),
           file.path(dirs$logs, "16_session_info.txt"))

# ── Manifest ───────────────────────────────────────────────────────────────
cat("\n=============================================================\n")
cat(" OUTPUT MANIFEST\n")
cat("=============================================================\n")
expected <- list(
  tables  = c("16_hallmark_membership_matrix.csv",
              "16_hallmark_membership_long.csv",
              "16_gene_hallmark_counts.csv",
              "16_hallmark_gene_counts.csv"),
  figures = c("Fig_16A_Hallmark_Membership_Tile.png",
              "Fig_16B_Hallmarks_Per_Gene.png"))
n_missing <- 0L
for (grp in names(expected)) {
  cat(sprintf("\n %s  (%s)\n", toupper(grp),
              normalizePath(dirs[[grp]], winslash = "/")))
  for (f in expected[[grp]]) {
    fp <- file.path(dirs[[grp]], f); ok <- file.exists(fp)
    if (!ok) n_missing <- n_missing + 1L
    cat(sprintf("   %s %s\n", ifelse(ok, "[OK]     ", "[MISSING]"), f))
  }
}
cat(ifelse(n_missing > 0,
           sprintf("\n [WARN] %d file(s) missing.\n", n_missing),
           "\n All expected outputs written.\n"))

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 16_hallmark_membership.R COMPLETE\n")
cat(" Runtime:", round(difftime(t_end, t_start, units = "mins"), 2), "min\n")
cat("=============================================================\n")
sink()
cat("Done. Report:", report, "\n")