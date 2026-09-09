# =============================================================================
# Script      : 3_directional_classifier.R
# Project     : DSP family dysregulation in breast cancer
# Description : Clinical translation analysis of the locked directional
#               diagnostic engine. Answers: can a biologically transparent
#               linear score (Σ Up − Σ Down) match or complement a Bayesian
#               XGBoost model? Is it organ-specific? What are its clinical
#               operating characteristics at each decision threshold?
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# Seed        : 42
# Depends on  : 1_training_pipeline.R  
# =============================================================================
# REQUIRED INPUTS:
#   output/models/final_xgboost_model_Bayesian.rds
#   output/models/final_preProc_scaler.rds
#   output/tables/final_stable_biomarker_signature.rds
#   output/tables/Final_Directional_Signature.rds
#     └─ columns: Gene (hyphen) | Direction (↑/↓ Tumor) | log2FC | SHAP_importance
#   output/tables/Biological_Validation_DE_Results.csv
#     └─ columns: Gene | Mean_Normal | Mean_Tumor | Log2FC | P_Value | FDR | DE_Direction
#   data/TCGA-BRCA.star_tpm.tsv.gz
#   data/TCGA-[LUAD,KIRC,COAD,GBM,PRAD,SKCM].star_tpm.tsv.gz  (Section D)
#
# DIRECTIONAL SCORE FORMULA:
#   Score = Σ(expression of ↑ Tumor genes) − Σ(expression of ↓ Tumor genes)
#   Applied to raw log2(TPM+1) values (NOT scaled) to preserve biological units.
#   Analogous to Oncotype DX — transparent and platform-robust.
#
# KEY OUTPUT FOR 05b:
#   output/directional/B_sample_scores.csv
#     └─ columns: SampleID | Group | XGB_Prob | Dir_Score_Raw | Hybrid_Prob
# =============================================================================

rm(list=ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed        = 42L,
  cv_folds    = 5L,    # folds for hybrid classifier OOF estimation
  bootstrap_n = 500L,
  conf_level  = 0.95
)

COHORTS_D <- list(
  list(label="Lung",         file=file.path("data","TCGA-LUAD.star_tpm.tsv.gz")),
  list(label="Kidney",       file=file.path("data","TCGA-KIRC.star_tpm.tsv.gz")),
  list(label="Colon",        file=file.path("data","TCGA-COAD.star_tpm.tsv.gz")),
  list(label="Glioblastoma", file=file.path("data","TCGA-GBM.star_tpm.tsv.gz")),
  list(label="Prostate",     file=file.path("data","TCGA-PRAD.star_tpm.tsv.gz")),
  list(label="Melanoma",     file=file.path("data","TCGA-SKCM.star_tpm.tsv.gz"))
)

PATHS <- list(
  model   = file.path("output","models","final_xgboost_model_Bayesian.rds"),
  scaler  = file.path("output","models","final_preProc_scaler.rds"),
  sig     = file.path("output","tables","final_stable_biomarker_signature.rds"),
  dir_sig = file.path("output","tables","Final_Directional_Signature.rds"),
  de_res  = file.path("output","tables","Biological_Validation_DE_Results.csv"),
  brca    = file.path("data","TCGA-BRCA.star_tpm.tsv.gz")
)

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman",quietly=TRUE)) install.packages("pacman")
pacman::p_load(
  here, data.table, dplyr, tibble, stringr, tidyr,
  xgboost, caret, pROC, ggplot2, scales,
  AnnotationDbi, org.Hs.eg.db
)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# =============================================================================
# DIRECTORIES
# =============================================================================
dirs <- list(
  out     = here("output","directional"),
  figures = here("output","directional","figures"),
  logs    = here("output","directional","logs")
)
lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# LOGGING
# =============================================================================
log_path <- file.path(dirs$logs,
                      paste0("05a_Directional_Run_",format(Sys.time(),"%Y%m%d_%H%M"),".log"))
sink(log_path, append=FALSE, split=TRUE)
t_start <- Sys.time()
cat("=============================================================\n")
cat(" 05a_directional_classifier.R\n")
cat(" Start:", as.character(t_start), "\n")
cat("=============================================================\n\n")
set.seed(CONFIG$seed)

# =============================================================================
# SHARED AESTHETICS (identical across all scripts)
# =============================================================================
COL <- list(
  tumor    = "#B22222", normal   = "#2166AC",
  low_expr = "#3288BD", high_expr= "#D53E4F", mid_expr = "#FFFFBF"
)

CANCER_COLOURS <- c(
  Breast="#B22222", Lung="#E69F00", Kidney="#56B4E9",
  Colon="#009E73", Glioblastoma="#CC79A7",
  Prostate="#D55E00", Melanoma="#0072B2"
)

theme_publication <- function(base_size=12) {
  theme_bw(base_size=base_size) +
    theme(
      plot.title       = element_text(face="bold",hjust=0.5,size=base_size+2),
      plot.subtitle    = element_text(hjust=0.5,colour="grey30",
                                      size=base_size-1,margin=margin(b=5)),
      plot.caption     = element_text(size=base_size-3.5,colour="grey45",
                                      hjust=0,lineheight=1.3,margin=margin(t=8)),
      axis.title       = element_text(face="bold",size=base_size),
      axis.text        = element_text(size=base_size-1,colour="black"),
      legend.title     = element_text(face="bold",size=base_size-1),
      legend.text      = element_text(size=base_size-2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour="black",fill=NA,linewidth=0.8),
      strip.background = element_rect(fill="grey90",colour="black"),
      strip.text       = element_text(face="bold",size=base_size-1),
      plot.margin      = unit(c(0.8,1.2,0.8,0.8),"cm")
    )
}

save_figure <- function(plot_obj, filename_base, width=9, height=8,
                        out_dir=dirs$figures) {
  base <- file.path(out_dir, filename_base)
  ggplot2::ggsave(paste0(base,".png"),  plot=plot_obj,
                  width=width, height=height, dpi=300)
  ggplot2::ggsave(paste0(base,".tiff"), plot=plot_obj,
                  width=width, height=height, dpi=600,
                  device="tiff", compression="lzw")
  grDevices::cairo_pdf(paste0(base,".pdf"), width=width, height=height)
  print(plot_obj); grDevices::dev.off()
  cat("  [SAVED]", filename_base, "(PNG 300dpi + TIFF 600dpi + PDF)\n")
  invisible(plot_obj)
}

# =============================================================================
# DIRECTIONAL SCORE FUNCTION
# Score = Σ(Up-gene expression) − Σ(Down-gene expression)
# Applied to RAW log2(TPM+1) values; hyphen column names required.
# =============================================================================
compute_dir_score <- function(X_raw_df, up_g, down_g) {
  up_p  <- intersect(up_g,   colnames(X_raw_df))
  dn_p  <- intersect(down_g, colnames(X_raw_df))
  score <- rep(0, nrow(X_raw_df))
  if (length(up_p) > 0) score <- score + rowSums(X_raw_df[, up_p, drop=FALSE])
  if (length(dn_p) > 0) score <- score - rowSums(X_raw_df[, dn_p, drop=FALSE])
  score
}

# =============================================================================
# 1. LOAD LOCKED DIAGNOSTIC ENGINE
# =============================================================================
cat("[INIT] Loading locked diagnostic engine...\n")

req <- c(PATHS$model, PATHS$scaler, PATHS$sig,
         PATHS$dir_sig, PATHS$de_res, PATHS$brca)
miss <- req[!file.exists(req)]
if (length(miss)>0)
  stop("Missing artefacts:\n", paste(miss, collapse="\n"),
       "\nRun 01_training_pipeline.R first.")

final_model    <- readRDS(PATHS$model)
preProc_scaler <- readRDS(PATHS$scaler)
stable_genes   <- readRDS(PATHS$sig)         # underscore names (XGBoost convention)
clean_genes    <- gsub("_","-",stable_genes) # hyphen names for display
n_genes        <- length(stable_genes)
dir_sig        <- readRDS(PATHS$dir_sig)     # Gene | Direction | log2FC | SHAP_importance
de_results     <- read.csv(PATHS$de_res, stringsAsFactors=FALSE)

# ── Parse directional gene vectors ───────────────────────────────────────────
# Dir_Symbol written by Script 01: "\u2191 Tumor" | "\u2193 Tumor" | "-"
# startsWith on the arrow character avoids cross-matching.
up_genes   <- dir_sig$Gene[startsWith(dir_sig$Direction, "\u2191")]
down_genes <- dir_sig$Gene[startsWith(dir_sig$Direction, "\u2193")]
if (length(up_genes) == 0 & length(down_genes) == 0) {
  cat("  [WARN] Arrow parsing returned 0 genes — falling back to log2FC sign\n")
  up_genes   <- dir_sig$Gene[dir_sig$log2FC > 0]
  down_genes <- dir_sig$Gene[dir_sig$log2FC < 0]
}
if (length(up_genes) + length(down_genes) == 0)
  stop("Could not parse directional genes from Final_Directional_Signature.rds")

cat("  Genes       :", paste(clean_genes, collapse=", "), "\n")
cat("  Up in Tumour:", paste(up_genes,    collapse=", "), "\n")
cat("  Down        :", paste(down_genes,  collapse=", "), "\n\n")

# ── Load TCGA-BRCA ────────────────────────────────────────────────────────────
cat("  Loading TCGA-BRCA...\n")
df_brca <- data.table::fread(PATHS$brca, header=TRUE, sep="\t")
if (!"gene_id" %in% names(df_brca)) names(df_brca)[1] <- "gene_id"
df_brca$gene_id_clean <- str_split_fixed(df_brca$gene_id, "\\.", 2)[, 1]
suppressMessages({
  df_brca$gene_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys=df_brca$gene_id_clean,
    column="SYMBOL", keytype="ENSEMBL", multiVals="first")
})

master_data <- df_brca %>%
  filter(!is.na(gene_symbol), gene_symbol %in% clean_genes) %>%
  group_by(gene_symbol) %>%
  summarise(across(starts_with("TCGA"), \(x) mean(x,na.rm=TRUE)),
            .groups="drop") %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix() %>% t() %>% as.data.frame() %>%
  mutate(sample_type_code=substr(rownames(.),14,15)) %>%
  filter(sample_type_code %in% c("01","11")) %>%
  mutate(Group=factor(ifelse(sample_type_code=="01","Tumor","Normal"),
                      levels=c("Normal","Tumor"))) %>%
  select(-sample_type_code) %>% na.omit()

colnames(master_data) <- gsub("-","_",colnames(master_data))

y_full     <- master_data$Group
y_full_num <- as.numeric(y_full=="Tumor")
X_full     <- master_data %>% select(-Group)
X_sc       <- predict(preProc_scaler, X_full[, stable_genes])

cat("  TCGA-BRCA:", nrow(master_data), "samples",
    "| Tumour:", sum(y_full=="Tumor"),
    "| Normal:", sum(y_full=="Normal"), "\n\n")

# Compute both scores on all TCGA-BRCA samples
xgb_probs <- predict(final_model,
                     xgboost::xgb.DMatrix(data=as.matrix(X_sc)))

X_full_hyp           <- X_full[, stable_genes]
colnames(X_full_hyp) <- clean_genes
dir_scores           <- compute_dir_score(X_full_hyp, up_genes, down_genes)

cat("  Directional Score range :",
    round(min(dir_scores),3), "to", round(max(dir_scores),3), "\n")
cat("  Tumour mean:", round(mean(dir_scores[y_full=="Tumor"]),3),
    "| Normal mean:", round(mean(dir_scores[y_full=="Normal"]),3), "\n\n")

# =============================================================================
# SECTION A — DIRECTIONAL PROFILE DEFINITION
# =============================================================================
cat("─────────────────────────────────────────────────────────────\n")
cat("SECTION A: Directional profile\n")
cat("─────────────────────────────────────────────────────────────\n\n")

dir_profile <- de_results %>%
  filter(Gene %in% clean_genes) %>%
  select(Gene, Mean_Normal, Mean_Tumor, Log2FC, FDR, DE_Direction) %>%
  left_join(dir_sig %>% select(Gene, Direction, SHAP_importance), by="Gene") %>%
  arrange(desc(abs(Log2FC)))

cat("  Directional profile (TCGA-BRCA training — pipeline-derived only):\n\n")
for (i in seq_len(nrow(dir_profile))) {
  cat(sprintf("  %-14s  %s  Normal=%.3f  Tumour=%.3f  Log2FC=%+.3f  FDR=%.2e  SHAP=%s\n",
              dir_profile$Gene[i],     dir_profile$Direction[i],
              dir_profile$Mean_Normal[i], dir_profile$Mean_Tumor[i],
              dir_profile$Log2FC[i],   dir_profile$FDR[i],
              dir_profile$SHAP_importance[i]))
}

write.csv(dir_profile, file.path(dirs$out,"A_directional_profile.csv"),
          row.names=FALSE)
cat("\n  [SAVED] A_directional_profile.csv\n\n")

# ── Figure A1: Log2FC bar coloured by biological direction ────────────────────
dir_col_map <- c("Up in Tumor"=COL$tumor, "Down in Tumor"=COL$normal,
                 "No Sig. Change"="grey60")

p_dirbar <- ggplot(dir_profile,
                   aes(x=reorder(Gene,-Log2FC), y=Log2FC, fill=DE_Direction)) +
  geom_col(colour="black", linewidth=0.35, width=0.72) +
  geom_hline(yintercept=0, colour="black", linewidth=0.5) +
  geom_text(aes(label=sprintf("%+.2f",Log2FC),
                vjust=ifelse(Log2FC>=0,-0.45,1.3)),
            size=3.4, colour="grey20") +
  scale_fill_manual(values=dir_col_map,
                    name="Biological Direction\n(DEG Analysis)") +
  labs(
    title    = "Locked Directional Signature: Log\u2082 Fold Change",
    subtitle = paste0("Tumour vs Normal  |  TCGA-BRCA training data  |  ",
                      n_genes," stable genes"),
    x = NULL,
    y = "Log\u2082 Fold Change  (Tumour / Normal)",
    caption  = paste0(
      "Direction from Wilcoxon DEG analysis (FDR < 0.05).\n",
      "Pipeline-derived from TCGA-BRCA training data only — ",
      "no GEPIA3 or external database used.")
  ) +
  theme_publication() +
  theme(axis.text.x=element_text(face="italic", size=11, angle=30, hjust=1),
        legend.position="top")
save_figure(p_dirbar, "Fig_05a_A1_Profile", width=10, height=6)

# =============================================================================
# SECTION B — CLASSIFIER COMPARISON
# =============================================================================
cat("─────────────────────────────────────────────────────────────\n")
cat("SECTION B: Classifier comparison (XGBoost vs Directional vs Hybrid)\n")
cat("─────────────────────────────────────────────────────────────\n\n")

roc_xgb <- pROC::roc(y_full_num, xgb_probs,  direction="<", quiet=TRUE)
roc_dir <- pROC::roc(y_full_num, dir_scores,  direction="<", quiet=TRUE)
ci_xgb  <- pROC::ci.auc(roc_xgb, conf.level=CONFIG$conf_level, method="delong")
ci_dir  <- pROC::ci.auc(roc_dir, conf.level=CONFIG$conf_level, method="delong")
dl_dx   <- pROC::roc.test(roc_dir, roc_xgb, method="delong")

cat(sprintf("  XGBoost AUC     : %.4f [%.4f\u2013%.4f]\n",
            as.numeric(ci_xgb[2]), as.numeric(ci_xgb[1]), as.numeric(ci_xgb[3])))
cat(sprintf("  Directional AUC : %.4f [%.4f\u2013%.4f]  DeLong p vs XGB: %s\n",
            as.numeric(ci_dir[2]), as.numeric(ci_dir[1]), as.numeric(ci_dir[3]),
            ifelse(dl_dx$p.value<0.001,"<0.001",sprintf("%.4f",dl_dx$p.value))))

# Hybrid: 5-fold CV logistic combination (out-of-fold predictions)
cat("  Fitting 5-fold CV hybrid classifier...\n")
set.seed(CONFIG$seed)
hyb_folds <- caret::createFolds(y_full, k=CONFIG$cv_folds,
                                list=TRUE, returnTrain=FALSE)
hyb_oof <- rep(NA_real_, length(y_full_num))

for (fi in seq_along(hyb_folds)) {
  te <- hyb_folds[[fi]]; tr <- setdiff(seq_along(y_full_num), te)
  glm_h       <- glm(y ~ xgb + dir, family=binomial(),
                     data=data.frame(y=y_full_num[tr],
                                     xgb=xgb_probs[tr], dir=dir_scores[tr]))
  hyb_oof[te] <- predict(glm_h,
                         newdata=data.frame(xgb=xgb_probs[te], dir=dir_scores[te]),
                         type="response")
}

roc_hyb <- pROC::roc(y_full_num, hyb_oof, direction="<", quiet=TRUE)
ci_hyb  <- pROC::ci.auc(roc_hyb, conf.level=CONFIG$conf_level, method="delong")
dl_hx   <- pROC::roc.test(roc_hyb, roc_xgb, method="delong")

cat(sprintf("  Hybrid AUC      : %.4f [%.4f\u2013%.4f]  DeLong p vs XGB: %s\n\n",
            as.numeric(ci_hyb[2]), as.numeric(ci_hyb[1]), as.numeric(ci_hyb[3]),
            ifelse(dl_hx$p.value<0.001,"<0.001",sprintf("%.4f",dl_hx$p.value))))

get_youden <- function(roc_obj)
  as.data.frame(pROC::coords(roc_obj,"best",best.method="youden",
                             ret=c("specificity","sensitivity","threshold"),
                             transpose=FALSE))[1,]
y_xgb <- get_youden(roc_xgb)
y_dir <- get_youden(roc_dir)
y_hyb <- get_youden(roc_hyb)

clf_sum <- data.frame(
  Classifier      = c("XGBoost (Bayesian)","Directional Score","Hybrid (XGB+Dir)"),
  AUC             = round(c(as.numeric(ci_xgb[2]),as.numeric(ci_dir[2]),
                            as.numeric(ci_hyb[2])),4),
  AUC_CI_Lo       = round(c(as.numeric(ci_xgb[1]),as.numeric(ci_dir[1]),
                            as.numeric(ci_hyb[1])),4),
  AUC_CI_Hi       = round(c(as.numeric(ci_xgb[3]),as.numeric(ci_dir[3]),
                            as.numeric(ci_hyb[3])),4),
  Sensitivity     = round(c(y_xgb$sensitivity,y_dir$sensitivity,y_hyb$sensitivity),4),
  Specificity     = round(c(y_xgb$specificity,y_dir$specificity,y_hyb$specificity),4),
  DeLong_p_vs_XGB = c(NA, round(dl_dx$p.value,4), round(dl_hx$p.value,4)),
  stringsAsFactors=FALSE)

write.csv(clf_sum, file.path(dirs$out,"B_classifier_comparison.csv"), row.names=FALSE)
cat("  [SAVED] B_classifier_comparison.csv\n")
print(clf_sum)

# B_sample_scores.csv — KEY OUTPUT consumed by 05b_subtype_specificity.R
scores_export <- data.frame(
  SampleID      = rownames(X_full),
  Group         = as.character(y_full),
  XGB_Prob      = round(xgb_probs, 5),
  Dir_Score_Raw = round(dir_scores, 4),
  Hybrid_Prob   = round(hyb_oof,   5)
)
write.csv(scores_export, file.path(dirs$out,"B_sample_scores.csv"), row.names=FALSE)
cat("  [SAVED] B_sample_scores.csv  ← input for 05b_subtype_specificity.R\n\n")

# ── Figure B1: Triple ROC ─────────────────────────────────────────────────────
get_roc_df <- function(roc_obj, lbl)
  pROC::coords(roc_obj,"all",ret=c("specificity","sensitivity"),
               transpose=FALSE) %>%
  as.data.frame() %>% mutate(fpr=1-specificity, Label=lbl)

lbl_xgb <- sprintf("XGBoost  (AUC=%.4f [%.4f\u2013%.4f])",
                   as.numeric(ci_xgb[2]),as.numeric(ci_xgb[1]),as.numeric(ci_xgb[3]))
lbl_dir <- sprintf("Directional  (AUC=%.4f [%.4f\u2013%.4f])",
                   as.numeric(ci_dir[2]),as.numeric(ci_dir[1]),as.numeric(ci_dir[3]))
lbl_hyb <- sprintf("Hybrid  (AUC=%.4f [%.4f\u2013%.4f])",
                   as.numeric(ci_hyb[2]),as.numeric(ci_hyb[1]),as.numeric(ci_hyb[3]))

roc_all <- rbind(get_roc_df(roc_xgb,lbl_xgb),
                 get_roc_df(roc_dir,lbl_dir),
                 get_roc_df(roc_hyb,lbl_hyb))
roc_all$Label <- factor(roc_all$Label, levels=c(lbl_xgb,lbl_dir,lbl_hyb))
clf_cols <- setNames(c(COL$tumor,"#1B7837","#762A83"),
                     c(lbl_xgb,lbl_dir,lbl_hyb))

p_triple <- ggplot(roc_all, aes(x=fpr,y=sensitivity,colour=Label)) +
  geom_line(linewidth=1.15) +
  geom_abline(slope=1,intercept=0,linetype="dashed",
              colour="grey50",linewidth=0.6) +
  geom_point(aes(x=1-y_xgb$specificity,y=y_xgb$sensitivity),
             colour=COL$tumor, size=3.5,shape=18,inherit.aes=FALSE) +
  geom_point(aes(x=1-y_dir$specificity,y=y_dir$sensitivity),
             colour="#1B7837",size=3.5,shape=18,inherit.aes=FALSE) +
  geom_point(aes(x=1-y_hyb$specificity,y=y_hyb$sensitivity),
             colour="#762A83",size=3.5,shape=18,inherit.aes=FALSE) +
  scale_colour_manual(values=clf_cols,
                      name="Classifier  (AUC [95% CI])") +
  scale_x_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
  scale_y_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
  labs(
    title    = "Classifier Comparison: XGBoost vs Directional vs Hybrid",
    subtitle = paste0("TCGA-BRCA Tumour vs Normal  |  ",n_genes,"-gene locked signature"),
    x = "1 \u2212 Specificity  (False Positive Rate)",
    y = "Sensitivity  (True Positive Rate)",
    caption  = paste0(
      "\u25c6 = Youden optimal threshold.  95% CI by DeLong.\n",
      "Directional Score = \u03a3(",paste(up_genes,collapse="+"),
      ") \u2212 \u03a3(",paste(down_genes,collapse="+"),") in log\u2082(TPM+1).\n",
      "Hybrid = 5-fold CV logistic combination of both scores (out-of-fold).")
  ) +
  theme_publication() +
  theme(aspect.ratio=1,
        legend.position=c(0.97,0.05),legend.justification=c(1,0),
        legend.background=element_rect(fill="white",colour="black",linewidth=0.5),
        legend.text=element_text(size=8.5))
save_figure(p_triple,"Fig_05a_B1_Triple_ROC",width=10,height=9.5)

# ── Figure B2: Score violin (Tumour vs Normal) ────────────────────────────────
score_vln <- data.frame(Group=y_full, XGB=xgb_probs, Dir=dir_scores) %>%
  tidyr::pivot_longer(cols=c(XGB,Dir),names_to="Score",values_to="Value") %>%
  mutate(Score=recode(Score,
                      XGB="XGBoost Probability",
                      Dir="Directional Score (raw)"))

wt_xgb <- wilcox.test(xgb_probs[y_full=="Tumor"],xgb_probs[y_full=="Normal"])
wt_dir <- wilcox.test(dir_scores[y_full=="Tumor"],dir_scores[y_full=="Normal"])

p_vln_b <- ggplot(score_vln, aes(x=Group,y=Value,fill=Group)) +
  geom_violin(trim=FALSE,alpha=0.75,colour="black",linewidth=0.4) +
  geom_jitter(width=0.15,size=0.45,alpha=0.25,colour="black") +
  # FIX: Changed colour="white" to "black" and width=0.35 to 0.45 for contrast
  stat_summary(fun=median,geom="crossbar",
               width=0.45,colour="black",linewidth=0.8) +
  facet_wrap(~Score,scales="free_y",ncol=2) +
  scale_fill_manual(values=c(Normal=COL$normal,Tumor=COL$tumor),
                    labels=c("Normal","Tumour")) +
  labs(
    title    = "Score Distributions: Tumour vs Normal (TCGA-BRCA)",
    subtitle = sprintf(
      "XGBoost Wilcoxon p %s  |  Directional Score Wilcoxon p %s",
      ifelse(wt_xgb$p.value<0.001,"< 0.001",sprintf("%.4f",wt_xgb$p.value)),
      ifelse(wt_dir$p.value<0.001,"< 0.001",sprintf("%.4f",wt_dir$p.value))),
    x=NULL, y="Score",
    caption="Black crossbar = median."
  ) +
  theme_publication() +
  theme(legend.position="none",
        axis.text.x=element_text(face="bold",size=11))
save_figure(p_vln_b,"Fig_05a_B2_ScoreViolin",width=10,height=6.5)

# =============================================================================
# SECTION D — ORGAN SPECIFICITY OF DIRECTIONAL SCORE
# =============================================================================
cat("\n─────────────────────────────────────────────────────────────\n")
cat("SECTION D: Organ specificity\n")
cat("─────────────────────────────────────────────────────────────\n\n")

load_cohort_dir <- function(co_label, co_file, stable_genes_us,
                            clean_genes_hyp, scaler, model, up_g, down_g) {
  if (!file.exists(co_file)) {
    cat("  [SKIP]",co_label,"— file not found\n"); return(NULL)
  }
  cat("  Loading",co_label,"...\n")
  df2 <- data.table::fread(co_file,header=TRUE,sep="\t")
  if (!"gene_id" %in% names(df2)) names(df2)[1] <- "gene_id"
  df2$gene_id_clean <- str_split_fixed(df2$gene_id,"\\.",2)[,1]
  suppressMessages({
    df2$gene_symbol <- AnnotationDbi::mapIds(
      org.Hs.eg.db, keys=df2$gene_id_clean,
      column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  })
  expr <- df2 %>%
    filter(!is.na(gene_symbol), gene_symbol %in% clean_genes_hyp) %>%
    group_by(gene_symbol) %>%
    summarise(across(starts_with("TCGA"),\(x) mean(x,na.rm=TRUE)),.groups="drop") %>%
    tibble::column_to_rownames("gene_symbol")
  keep <- colnames(expr)[substr(colnames(expr),14,15)=="01"]
  if (length(keep)==0) keep <- colnames(expr)
  expr <- expr[,keep,drop=FALSE]
  rownames(expr) <- gsub("-","_",rownames(expr))
  miss_g <- setdiff(stable_genes_us,rownames(expr))
  if (length(miss_g)>0) {
    z <- matrix(0,nrow=length(miss_g),ncol=ncol(expr),
                dimnames=list(miss_g,colnames(expr)))
    expr <- rbind(expr,z)
  }
  X_raw <- t(expr[stable_genes_us,,drop=FALSE])
  if (max(X_raw,na.rm=TRUE)>25) X_raw <- log2(X_raw+1)
  X_raw[is.na(X_raw)] <- 0
  X_sc2   <- predict(scaler,as.data.frame(X_raw))
  xgb_p2  <- predict(model,xgboost::xgb.DMatrix(data=as.matrix(X_sc2)))
  X_hyp2  <- as.data.frame(X_raw); colnames(X_hyp2) <- gsub("_","-",colnames(X_hyp2))
  dir_sc2 <- compute_dir_score(X_hyp2,up_g,down_g)
  data.frame(SampleID=rownames(X_raw),Cohort=co_label,Is_Breast=0L,
             XGB_Prob=round(xgb_p2,5),Dir_Score=round(dir_sc2,4),
             stringsAsFactors=FALSE)
}

breast_d_df <- data.frame(
  SampleID=rownames(X_full), Cohort="Breast", Is_Breast=1L,
  XGB_Prob=round(xgb_probs,5), Dir_Score=round(dir_scores,4),
  stringsAsFactors=FALSE
)
orgspec_d <- breast_d_df
for (co in COHORTS_D) {
  res_d <- tryCatch(
    load_cohort_dir(co$label,co$file,stable_genes,clean_genes,
                    preProc_scaler,final_model,up_genes,down_genes),
    error=function(e){cat("  [ERROR]",co$label,":",conditionMessage(e),"\n");NULL})
  if (!is.null(res_d)) orgspec_d <- rbind(orgspec_d,res_d)
}
orgspec_d$Cohort <- factor(orgspec_d$Cohort,
                           levels=c("Breast",sapply(COHORTS_D,`[[`,"label")))
non_breast_labels <- levels(orgspec_d$Cohort)
non_breast_labels <- non_breast_labels[non_breast_labels != "Breast"]

auc_d_rows <- lapply(non_breast_labels, function(nm) {
  sub  <- orgspec_d %>% filter(Cohort %in% c("Breast",nm))
  r_d  <- pROC::roc(sub$Is_Breast,sub$Dir_Score,levels=c(0,1),direction="<",quiet=TRUE)
  r_x  <- pROC::roc(sub$Is_Breast,sub$XGB_Prob, levels=c(0,1),direction="<",quiet=TRUE)
  ci_d <- pROC::ci.auc(r_d,conf.level=CONFIG$conf_level,method="delong")
  ci_x <- pROC::ci.auc(r_x,conf.level=CONFIG$conf_level,method="delong")
  cat(sprintf("  Breast vs %-14s  Dir=%.4f  XGB=%.4f\n",
              nm,as.numeric(ci_d[2]),as.numeric(ci_x[2])))
  data.frame(Cancer_Type=nm,
             Dir_AUC=round(as.numeric(ci_d[2]),4),
             Dir_CI_Lo=round(as.numeric(ci_d[1]),4),
             Dir_CI_Hi=round(as.numeric(ci_d[3]),4),
             XGB_AUC=round(as.numeric(ci_x[2]),4),
             XGB_CI_Lo=round(as.numeric(ci_x[1]),4),
             XGB_CI_Hi=round(as.numeric(ci_x[3]),4),
             stringsAsFactors=FALSE)
})
auc_d_df <- do.call(rbind, auc_d_rows)
write.csv(auc_d_df, file.path(dirs$out,"D_orgspec_AUC.csv"), row.names=FALSE)
cat("  [SAVED] D_orgspec_AUC.csv\n\n")

auc_d_long <- auc_d_df %>%
  mutate(Cancer_Type=factor(Cancer_Type,
                            levels=auc_d_df$Cancer_Type[order(auc_d_df$Dir_AUC)])) %>%
  tidyr::pivot_longer(cols=c(Dir_AUC,XGB_AUC),
                      names_to="Score",values_to="AUC") %>%
  mutate(CI_Lo=ifelse(Score=="Dir_AUC",Dir_CI_Lo,XGB_CI_Lo),
         CI_Hi=ifelse(Score=="Dir_AUC",Dir_CI_Hi,XGB_CI_Hi),
         Score=recode(Score,Dir_AUC="Directional Score",XGB_AUC="XGBoost"))

# FIX: Dynamically calculate the minimum x limit based on your actual data
min_x <- min(auc_d_long$CI_Lo, na.rm = TRUE) - 0.05

p_forest_d <- ggplot(auc_d_long,
                     aes(x=AUC,y=Cancer_Type,colour=Score,
                         xmin=CI_Lo,xmax=CI_Hi)) +
  geom_errorbarh(height=0.28,linewidth=0.85,position=position_dodge(0.5)) +
  geom_point(size=3.5,shape=18,position=position_dodge(0.5)) +
  geom_vline(xintercept=0.5,linetype="dashed",colour="grey40",linewidth=0.6) +
  scale_colour_manual(values=c("Directional Score"="#1B7837","XGBoost"=COL$tumor),
                      name="Classifier") +
  scale_x_continuous(breaks=seq(0.1, 1.0, by=0.1)) +
  # FIX: Apply dynamic limits via coord_cartesian to prevent clipping of low values
  coord_cartesian(xlim = c(min_x, 1.05)) +
  labs(
    title    = "Organ Specificity: Directional Score vs XGBoost AUC",
    subtitle = "Breast vs each non-breast TCGA cohort  |  95% CI by DeLong",
    x="AUC  (Breast vs Cancer Type)", y=NULL,
    caption  = paste0(
      "Dashed = AUC 0.5 (no discrimination).\n",
      "A breast-specific signature should score HIGH for Breast and LOW for all others.")
  ) +
  theme_publication() +
  theme(legend.position="top",axis.text.y=element_text(face="bold",size=12))
save_figure(p_forest_d,"Fig_05a_D1_OrgSpec_Forest",width=10,height=6.5)

# =============================================================================
# SECTION E — CLINICAL RELEVANCE THRESHOLD ANALYSIS
# =============================================================================
cat("\n─────────────────────────────────────────────────────────────\n")
cat("SECTION E: Clinical relevance (threshold analysis)\n")
cat("─────────────────────────────────────────────────────────────\n\n")

compute_clinical_metrics <- function(y_true, score_vec,
                                     comparison_label, n_steps=20) {
  thresholds <- unique(quantile(score_vec, probs=seq(0.02,0.98,length.out=n_steps)))
  do.call(rbind, lapply(thresholds, function(thr) {
    TP <- sum(score_vec>=thr & y_true==1); TN <- sum(score_vec<thr & y_true==0)
    FP <- sum(score_vec>=thr & y_true==0); FN <- sum(score_vec<thr & y_true==1)
    N  <- length(y_true)
    Se <- TP/(TP+FN); Sp <- TN/(TN+FP)
    PPV <- TP/(TP+FP); NPV <- TN/(TN+FN)
    LRp <- ifelse(Sp<1, Se/(1-Sp), NA); LRn <- ifelse(Se>0,(1-Se)/Sp,NA)
    mcc_d <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    data.frame(
      Comparison=comparison_label, Threshold=round(thr,4),
      N_Total=N, N_Positive=TP+FN, N_Negative=TN+FP,
      TP=TP,TN=TN,FP=FP,FN=FN,
      Sensitivity=round(Se,4), Specificity=round(Sp,4),
      PPV=round(PPV,4), NPV=round(NPV,4),
      LR_Positive=round(LRp,3), LR_Negative=round(LRn,4),
      Accuracy=round((TP+TN)/N,4), Youden_J=round(Se+Sp-1,4),
      MCC=round(ifelse(mcc_d==0,0,((TP*TN)-(FP*FN))/mcc_d),4),
      F1=round(2*TP/(2*TP+FP+FN),4),
      stringsAsFactors=FALSE)
  }))
}

clin_brca <- compute_clinical_metrics(y_full_num, dir_scores,
                                      "Breast_Tumour_vs_Normal")
y_pooled  <- orgspec_d$Is_Breast; s_pooled <- orgspec_d$Dir_Score
clin_pool <- compute_clinical_metrics(y_pooled, s_pooled,
                                      "Breast_vs_All_NonBreast_Pooled")
clin_per  <- do.call(rbind, lapply(non_breast_labels, function(nm) {
  sub <- orgspec_d %>% filter(Cohort %in% c("Breast",nm))
  compute_clinical_metrics(sub$Is_Breast,sub$Dir_Score,paste0("Breast_vs_",nm))
}))

clin_all <- rbind(clin_brca,clin_pool,clin_per)
write.csv(clin_all, file.path(dirs$out,"E_clinical_threshold_table.csv"),
          row.names=FALSE)

best_thresh <- clin_all %>%
  group_by(Comparison) %>%
  slice_max(order_by=Youden_J, n=1, with_ties=FALSE) %>%
  ungroup()
write.csv(best_thresh, file.path(dirs$out,"E_best_threshold_summary.csv"),
          row.names=FALSE)
cat("  [SAVED] E_clinical_threshold_table.csv\n")
cat("  [SAVED] E_best_threshold_summary.csv\n\n")
cat("  ── Youden-optimal thresholds ──\n")
print(best_thresh[,c("Comparison","Threshold","Sensitivity","Specificity",
                     "PPV","NPV","LR_Positive","LR_Negative","Youden_J","MCC")])

# ── Figure E1: Metric curves (Tumour vs Normal) ───────────────────────────────
clin_long_e <- clin_brca %>%
  select(Threshold,Sensitivity,Specificity,PPV,NPV,Youden_J) %>%
  tidyr::pivot_longer(cols=-Threshold,names_to="Metric",values_to="Value") %>%
  mutate(Metric=factor(Metric,
                       levels=c("Sensitivity","Specificity","PPV","NPV","Youden_J"),
                       labels=c("Sensitivity","Specificity","PPV","NPV","Youden J")))

p_clin_e1 <- ggplot(clin_long_e, aes(x=Threshold,y=Value,
                                     colour=Metric,group=Metric)) +
  geom_line(linewidth=1.1) + geom_point(size=2.5) +
  scale_colour_manual(
    values=c(Sensitivity="#B22222",Specificity="#2166AC",
             PPV="#1B7837",NPV="#762A83","Youden J"="darkorange2"),
    name="Metric") +
  scale_y_continuous(limits=c(0,1.05),breaks=seq(0,1,0.2)) +
  labs(
    title    = "Clinical Relevance: Directional Score Threshold Analysis",
    subtitle = "Breast Tumour vs Normal (TCGA-BRCA)",
    x        = "Directional Score Threshold",  # <-- FIX: Removed the giant string of gene names
    y        = "Metric Value",
    caption  = "Score = \u03a3(Up genes) \u2212 \u03a3(Down genes) in log\u2082(TPM+1) units."
  ) +
  theme_publication() + theme(legend.position="right")
save_figure(p_clin_e1,"Fig_05a_E1_Clinical_Thresholds",width=10,height=6.5)

# ── Figure E2: LR+ / LR− bar chart across all comparisons ────────────────────
lr_df <- best_thresh %>%
  filter(!is.na(LR_Positive),!is.infinite(LR_Positive),
         !is.na(LR_Negative),!is.infinite(LR_Negative)) %>%
  mutate(Comparison=stringr::str_replace_all(Comparison,"_"," ")) %>%
  arrange(desc(LR_Positive))

# Calculate dynamic offset for text to look good regardless of range
max_lr <- max(lr_df$LR_Positive, na.rm=TRUE)

p_lr_e2 <- ggplot(lr_df, aes(x=reorder(Comparison,LR_Positive))) +
  geom_col(aes(y=LR_Positive),  fill=COL$tumor,
           alpha=0.85,colour="black",linewidth=0.3,width=0.6) +
  geom_col(aes(y=-LR_Negative), fill=COL$normal,
           alpha=0.85,colour="black",linewidth=0.3,width=0.6) +
  geom_hline(yintercept=0,colour="black",linewidth=0.5) +
  geom_text(aes(y=LR_Positive + (max_lr * 0.05),
                label=sprintf("LR+ = %.1f",LR_Positive)),
            size=3.5,colour=COL$tumor,hjust=0, fontface="bold") +
  geom_text(aes(y=-LR_Negative - (max_lr * 0.05),
                label=sprintf("LR\u2212 = %.3f",LR_Negative)),
            size=3.5,colour=COL$normal,hjust=1, fontface="bold") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) + # <-- FIX: Adds padding for text
  labs(
    title    = "Likelihood Ratios at Youden-Optimal Threshold",
    subtitle = "Directional Score  |  All comparisons",
    x=NULL, y="LR+  (red, upward)  /  LR\u2212  (blue, downward)",
    caption  = paste0(
      "LR+ > 10 = strong evidence for disease.  ",
      "LR\u2212 < 0.1 = strong evidence against.\n",
      "Threshold selected by Youden J index per comparison.")
  ) +
  theme_publication() +
  theme(axis.text.y=element_text(size=10, face="bold"))
save_figure(p_lr_e2,"Fig_05a_E2_LikelihoodRatios",width=11,height=6.5)

# =============================================================================
# SAVE RESULTS OBJECT
# =============================================================================
saveRDS(
  list(dir_profile    = dir_profile,
       clf_comparison = clf_sum,
       scores_brca    = scores_export,
       orgspec_scores = orgspec_d,
       orgspec_auc    = auc_d_df,
       clinical_table = clin_all,
       best_thresholds= best_thresh),
  file.path(dirs$out,"05a_Directional_Results.rds")
)
cat("\n  [SAVED] 05a_Directional_Results.rds\n")

# =============================================================================
# WRAP UP
# =============================================================================
sink(file.path(dirs$logs,"05a_session_info.txt"))
print(sessionInfo()); sink()
sink(log_path, append=TRUE, split=TRUE)

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 05a_directional_classifier.R  COMPLETE\n")
cat(" Runtime:", round(difftime(t_end,t_start,units="mins"),2),"mins\n")
cat("\n SAVED ARTEFACTS:\n")
cat("  Tables  : A_directional_profile.csv\n")
cat("            B_classifier_comparison.csv\n")
cat("            B_sample_scores.csv          \u2190 KEY INPUT for 05b\n")
cat("            D_orgspec_AUC.csv\n")
cat("            E_clinical_threshold_table.csv\n")
cat("            E_best_threshold_summary.csv\n")
cat("            05a_Directional_Results.rds\n")
cat("  Figures : Fig_05a_A1_Profile.*\n")
cat("            Fig_05a_B1_Triple_ROC.*\n")
cat("            Fig_05a_B2_ScoreViolin.*\n")
cat("            Fig_05a_D1_OrgSpec_Forest.*\n")
cat("            Fig_05a_E1_Clinical_Thresholds.*\n")
cat("            Fig_05a_E2_LikelihoodRatios.*\n")
cat("            (PNG 300dpi + TIFF 600dpi + PDF)\n")
cat("  Logs    : 05a_Directional_Run_[timestamp].log\n")
cat("            05a_session_info.txt\n")
cat("=============================================================\n")
sink()