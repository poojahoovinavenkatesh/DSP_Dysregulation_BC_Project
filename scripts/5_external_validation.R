# =============================================================================
# Script      : 5_external_validation.R
# Project     : DSP family dysregulation in breast cancer
# Description : External validation of the locked diagnostic engine on independent
#               GEO datasets. Computes BOTH:
#               (A) ML Probability Score  — Bayesian XGBoost model output (0–1)
#               (B) Linear Directional Score — sum(Up genes) − sum(Down genes)
#                   This is platform-robust, transparent, and clinically
#                   interpretable (analogous to Oncotype DX scoring logic).
#               Both scores are validated together, allowing direct comparison.
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# Seed        : 42
# Depends on  : 01_training_pipeline.R (must be run first)
# =============================================================================
# REQUIRED INPUTS:
#   output/models/final_xgboost_model_Bayesian.rds
#   output/models/final_preProc_scaler.rds
#   output/tables/final_stable_biomarker_signature.rds
#   output/tables/Final_Directional_Signature.rds     ← Gene | Direction (↑/↓) | log2FC | SHAP_importance
#   data/validation/GSE233242_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
#   data/validation/GSE58135_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
#
# ADDING MORE DATASETS: edit ONLY the DATASETS list in the CONFIG block.
# =============================================================================

rm(list=ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(seed=42L, bootstrap_n=500L, conf_level=0.95)

DATASETS <- list(
  list(
    gse_id           = "GSE233242",
    file_path        = file.path("data","validation",
                                 "GSE233242_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
    normal_keywords  = c("normal"),
    tumor_keywords   = c("tumor"),
    exclude_keywords = character(0)
  ),
  list(
    gse_id           = "GSE58135",
    file_path        = file.path("data","validation",
                                 "GSE58135_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz"),
    normal_keywords  = c("uninvolved","mammoplasty","adjacent","normal","control"),
    tumor_keywords   = c("tumor","breast cancer","tnbc","er\\+"),
    exclude_keywords = c("cell line")
  )
)

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman",quietly=TRUE)) install.packages("pacman")
pacman::p_load(
  here, GEOquery, xgboost, dplyr, tibble, pROC,
  ggplot2, data.table, stringr, tidyr,
  AnnotationDbi, org.Hs.eg.db, pheatmap, caret, scales
)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# =============================================================================
# DIRECTORIES
# =============================================================================
dirs <- list(
  models  = here("output","models"),
  tables  = here("output","tables"),
  val     = here("output","validation"),
  figures = here("output","validation","figures"),
  logs    = here("output","validation","logs")
)
lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# LOGGING
# =============================================================================
log_path <- file.path(dirs$logs,
                      paste0("02_Validation_Run_",format(Sys.time(),"%Y%m%d_%H%M"),".log"))
sink(log_path, append=FALSE, split=TRUE)
t_start <- Sys.time()
cat("=============================================================\n")
cat(" 02_external_validation.R\n")
cat(" Start   :", as.character(t_start), "\n")
cat(" Datasets:", length(DATASETS), "\n")
cat("=============================================================\n\n")
set.seed(CONFIG$seed)

# =============================================================================
# SHARED AESTHETICS (identical across all 5 scripts)
# =============================================================================
COL <- list(
  tumor    = "#B22222", normal   = "#2166AC",
  low_expr = "#3288BD", high_expr= "#D53E4F", mid_expr = "#FFFFBF",
  score_up = "#B22222", score_dn = "#2166AC"
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
  ggplot2::ggsave(paste0(base,".png"), plot=plot_obj,
                  width=width, height=height, dpi=300)
  ggplot2::ggsave(paste0(base,".tiff"), plot=plot_obj,
                  width=width, height=height, dpi=600,
                  device="tiff", compression="lzw")
  grDevices::cairo_pdf(paste0(base,".pdf"), width=width, height=height)
  print(plot_obj); grDevices::dev.off()
  cat("  [SAVED]", filename_base, "(PNG 300dpi + TIFF 600dpi + PDF)\n")
  invisible(plot_obj)
}

DATASET_COLOURS <- c("#E69F00","#56B4E9","#009E73",
                     "#CC79A7","#D55E00","#0072B2","#F0E442")

# =============================================================================
# 1. LOAD LOCKED DIAGNOSTIC ENGINE
# =============================================================================
cat("[1/5] Loading locked diagnostic engine...\n")

req_files <- c(
  file.path(dirs$models,"final_xgboost_model_Bayesian.rds"),
  file.path(dirs$models,"final_preProc_scaler.rds"),
  file.path(dirs$tables,"final_stable_biomarker_signature.rds"),
  file.path(dirs$tables,"Final_Directional_Signature.rds")
)
miss <- req_files[!file.exists(req_files)]
if (length(miss)>0) stop("Missing artefacts:\n",paste(miss,collapse="\n"),
                         "\nRun 01_training_pipeline.R first.")

final_model    <- readRDS(req_files[1])
preProc_scaler <- readRDS(req_files[2])
stable_genes   <- readRDS(req_files[3])         # underscore format
clean_genes    <- gsub("_","-",stable_genes)    # hyphen format for display
dir_sig        <- readRDS(req_files[4])         # Final_Directional_Signature.rds

# Extract directional vectors from the locked signature
# Parse directional vectors from Dir_Symbol column.
# Dir_Symbol values set by Script 01:
#   "\u2191 Tumor"  (↑ Tumor) — up-regulated in tumour
#   "\u2193 Tumor"  (↓ Tumor) — down-regulated in tumour
#   "-"             — no significant change
# Match on the arrow character ONLY to avoid cross-matching.
up_genes   <- dir_sig$Gene[startsWith(dir_sig$Direction, "\u2191")]
down_genes <- dir_sig$Gene[startsWith(dir_sig$Direction, "\u2193")]
# Fallback: use Log2FC sign (handles any future encoding variants)
if (length(up_genes) == 0 & length(down_genes) == 0) {
  cat("  [WARN] Arrow parsing returned 0 genes — falling back to log2FC sign\n")
  up_genes   <- dir_sig$Gene[dir_sig$log2FC > 0]
  down_genes <- dir_sig$Gene[dir_sig$log2FC < 0]
}
if (length(up_genes) + length(down_genes) == 0)
  stop("Could not parse any directional genes from Final_Directional_Signature.rds")

cat("  Genes:", paste(clean_genes, collapse=", "), "\n")
cat("  Up in Tumour   :", paste(up_genes, collapse=", "), "\n")
cat("  Down in Tumour :", paste(down_genes, collapse=", "), "\n\n")

# =============================================================================
# DIRECTIONAL SCORE FUNCTION
# Computes: Score = sum(Up genes expression) - sum(Down genes expression)
# This is the REMARK/TRIPOD-compliant clinically transparent score.
# Applied to RAW (unscaled) log2(TPM+1) values to maintain biological units.
# =============================================================================
compute_dir_score <- function(X_raw_df, up_g, down_g) {
  # X_raw_df: samples × genes data frame, HYPHEN gene names
  up_g_present   <- intersect(up_g,   colnames(X_raw_df))
  down_g_present <- intersect(down_g, colnames(X_raw_df))
  score <- 0
  if (length(up_g_present)   > 0)
    score <- score + rowSums(X_raw_df[, up_g_present, drop=FALSE])
  if (length(down_g_present) > 0)
    score <- score - rowSums(X_raw_df[, down_g_present, drop=FALSE])
  score
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
get_meta <- function(gse_obj) {
  meta     <- pData(gse_obj)
  tgt_cols <- grep("title|source|characteristic|description",
                   colnames(meta),value=TRUE,ignore.case=TRUE)
  meta$SearchString <- apply(meta[,tgt_cols,drop=FALSE],1,
                             function(x) tolower(paste(x,collapse=" | ")))
  meta$GSM <- meta$geo_accession
  meta
}

map_group <- function(ss, tkw, nkw, exkw) {
  if (length(exkw)>0 &&
      any(sapply(exkw,function(k) grepl(k,ss,ignore.case=TRUE))))
    return(NA_character_)
  if (any(sapply(nkw, function(k) grepl(k,ss,ignore.case=TRUE))))
    return("Normal")
  if (any(sapply(tkw, function(k) grepl(k,ss,ignore.case=TRUE))))
    return("Tumor")
  NA_character_
}

# =============================================================================
# 2. VALIDATION ENGINE
# =============================================================================
run_validation <- function(ds, model, genes_us, genes_disp,
                           scaler, up_g, down_g, seed) {
  gse_id <- ds$gse_id
  cat("──────────────────────────────────────────────\n")
  cat("Dataset:", gse_id, "\n")
  
  if (!file.exists(ds$file_path)) {
    cat("  [SKIP] File not found:", ds$file_path, "\n\n"); return(NULL)
  }
  
  # Metadata
  gse_obj  <- GEOquery::getGEO(gse_id, GSEMatrix=TRUE, AnnotGPL=FALSE)[[1]]
  meta_all <- get_meta(gse_obj) %>%
    mutate(Group=sapply(SearchString,map_group,
                        tkw=ds$tumor_keywords,nkw=ds$normal_keywords,
                        exkw=ds$exclude_keywords))
  write.csv(meta_all %>% select(GSM,Group,SearchString),
            file.path(dirs$logs,paste0("Mapping_Log_",gse_id,".csv")),
            row.names=FALSE)
  meta <- meta_all %>% filter(!is.na(Group))
  cat("  Tumor:", sum(meta$Group=="Tumor"),
      "| Normal:", sum(meta$Group=="Normal"), "\n")
  if (length(unique(meta$Group))<2) stop("Only one class in ",gse_id)
  
  # Expression matrix
  expr <- data.table::fread(ds$file_path)
  colnames(expr)[1] <- "EntrezID"
  expr$EntrezID <- as.character(expr$EntrezID)
  suppressMessages({
    sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=expr$EntrezID,
                                 column="SYMBOL",keytype="ENTREZID",
                                 multiVals="first")
  })
  expr$gene_symbol <- sym
  expr_clean <- expr %>%
    filter(!is.na(gene_symbol) & gene_symbol!="") %>%
    select(-EntrezID) %>%
    group_by(gene_symbol) %>%
    summarise(across(everything(),\(x) mean(x,na.rm=TRUE)),.groups="drop") %>%
    tibble::column_to_rownames("gene_symbol")
  
  common_sam <- intersect(meta$GSM, colnames(expr_clean))
  if (length(common_sam)==0) stop("No sample overlap for ",gse_id)
  meta_al  <- meta %>% filter(GSM %in% common_sam)
  expr_al  <- expr_clean[, meta_al$GSM, drop=FALSE]
  cat("  Matched samples:", length(common_sam), "\n")
  
  # Gene alignment
  rownames(expr_al) <- gsub("-","_",rownames(expr_al))
  missing_g <- setdiff(genes_us, rownames(expr_al))
  present_g <- intersect(genes_us, rownames(expr_al))
  
  write.csv(
    data.frame(Dataset=gse_id,N_Present=length(present_g),
               N_Imputed=length(missing_g),
               Pct=round(length(present_g)/length(genes_us)*100,1),
               Missing_Genes=paste(gsub("_","-",missing_g),collapse=";"),
               stringsAsFactors=FALSE),
    file.path(dirs$logs,paste0("Gene_Coverage_",gse_id,".csv")),
    row.names=FALSE)
  
  if (length(missing_g)>0) {
    cat("  [WARN] Imputing",length(missing_g),"absent genes with 0\n")
    z <- matrix(0,nrow=length(missing_g),ncol=ncol(expr_al),
                dimnames=list(missing_g,colnames(expr_al)))
    expr_al <- rbind(expr_al,z)
  }
  
  X_raw <- t(expr_al[genes_us,,drop=FALSE])   # samples × genes (underscore names)
  
  # Log2 guard: TCGA star_tpm already log2(TPM+1); GEO NCBI files = linear TPM
  if (max(X_raw,na.rm=TRUE) > 25) {
    cat("  Max =",round(max(X_raw,na.rm=TRUE),1),"→ applying log2(TPM+1)\n")
    X_raw <- log2(X_raw + 1)
  }
  
  # ── (A) ML Probability Score — XGBoost ──────────────────────────────────
  X_sc     <- predict(scaler, as.data.frame(X_raw))
  xgb_prob <- predict(model, xgboost::xgb.DMatrix(data=as.matrix(X_sc)))
  
  # ── (B) Linear Directional Score — biology-first ────────────────────────
  # Use HYPHEN names for directional score (matches dir_sig$Gene format)
  X_raw_hyp           <- as.data.frame(X_raw)
  colnames(X_raw_hyp) <- gsub("_","-",colnames(X_raw_hyp))
  dir_score           <- compute_dir_score(X_raw_hyp, up_g, down_g)
  
  y_true <- as.numeric(meta_al$Group=="Tumor")
  cat("  Gene coverage:", length(present_g), "/", length(genes_us), "\n")
  
  # ROC for both scores
  roc_xgb <- pROC::roc(y_true, xgb_prob,   levels=c(0,1),direction="<",quiet=TRUE)
  roc_dir <- pROC::roc(y_true, dir_score,   levels=c(0,1),direction="<",quiet=TRUE)
  ci_xgb  <- pROC::ci.auc(roc_xgb, conf.level=CONFIG$conf_level, method="delong")
  ci_dir  <- pROC::ci.auc(roc_dir, conf.level=CONFIG$conf_level, method="delong")
  delong  <- pROC::roc.test(roc_xgb, roc_dir, method="delong")
  
  cat(sprintf("  XGBoost AUC        : %.4f [%.4f–%.4f]\n",
              as.numeric(ci_xgb[2]),as.numeric(ci_xgb[1]),as.numeric(ci_xgb[3])))
  cat(sprintf("  Directional AUC    : %.4f [%.4f–%.4f]\n",
              as.numeric(ci_dir[2]),as.numeric(ci_dir[1]),as.numeric(ci_dir[3])))
  cat(sprintf("  DeLong (XGB vs Dir): p = %s\n",
              ifelse(delong$p.value<0.001,"<0.001",sprintf("%.4f",delong$p.value))))
  
  # Bootstrap CI band on XGBoost sensitivity (for individual plot)
  set.seed(seed)
  ci_se <- pROC::ci.se(roc_xgb, specificities=seq(0,1,by=0.04),
                       conf.level=CONFIG$conf_level,method="bootstrap",
                       boot.n=CONFIG$bootstrap_n,quiet=TRUE)
  ci_band <- data.frame(fpr=1-seq(0,1,by=0.04),
                        lower=as.numeric(ci_se[,1]),
                        upper=as.numeric(ci_se[,3]))
  
  youden_xgb <- as.data.frame(
    pROC::coords(roc_xgb,"best",best.method="youden",
                 ret=c("specificity","sensitivity","threshold"),
                 transpose=FALSE))[1,]
  youden_dir <- as.data.frame(
    pROC::coords(roc_dir,"best",best.method="youden",
                 ret=c("specificity","sensitivity","threshold"),
                 transpose=FALSE))[1,]
  
  # Confusion matrices
  get_cm_metrics <- function(probs, y, thr) {
    pred_cl <- factor(ifelse(probs>=thr,"Tumor","Normal"),levels=c("Normal","Tumor"))
    true_cl <- factor(ifelse(y==1,"Tumor","Normal"),levels=c("Normal","Tumor"))
    cm <- caret::confusionMatrix(pred_cl,true_cl,positive="Tumor")
    TP<-cm$table[2,2];TN<-cm$table[1,1];FP<-cm$table[2,1];FN<-cm$table[1,2]
    mcc_d <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    c(Sensitivity=unname(cm$byClass["Sensitivity"]),
      Specificity=unname(cm$byClass["Specificity"]),
      Precision  =unname(cm$byClass["Pos Pred Value"]),
      F1         =unname(cm$byClass["F1"]),
      MCC        =ifelse(mcc_d==0,0,((TP*TN)-(FP*FN))/mcc_d),
      Accuracy   =unname(cm$overall["Accuracy"]))
  }
  m_xgb <- get_cm_metrics(xgb_prob, y_true, youden_xgb$threshold)
  m_dir <- get_cm_metrics(dir_score, y_true, youden_dir$threshold)
  
  sum_row <- data.frame(
    Dataset=gse_id, N_Total=nrow(meta_al),
    N_Tumor=sum(meta_al$Group=="Tumor"),
    N_Normal=sum(meta_al$Group=="Normal"),
    Genes_Present=length(present_g), Genes_Imputed=length(missing_g),
    XGB_AUC=round(as.numeric(ci_xgb[2]),4),
    XGB_CI_Lo=round(as.numeric(ci_xgb[1]),4),
    XGB_CI_Hi=round(as.numeric(ci_xgb[3]),4),
    XGB_Sensitivity=round(m_xgb["Sensitivity"],4),
    XGB_Specificity=round(m_xgb["Specificity"],4),
    XGB_F1=round(m_xgb["F1"],4), XGB_MCC=round(m_xgb["MCC"],4),
    Dir_AUC=round(as.numeric(ci_dir[2]),4),
    Dir_CI_Lo=round(as.numeric(ci_dir[1]),4),
    Dir_CI_Hi=round(as.numeric(ci_dir[3]),4),
    Dir_Sensitivity=round(m_dir["Sensitivity"],4),
    Dir_Specificity=round(m_dir["Specificity"],4),
    Dir_F1=round(m_dir["F1"],4), Dir_MCC=round(m_dir["MCC"],4),
    DeLong_P_XGBvsDir=round(delong$p.value,4),
    stringsAsFactors=FALSE)
  
  write.csv(sum_row,
            file.path(dirs$val,paste0("Validation_Summary_",gse_id,".csv")),
            row.names=FALSE)
  cat("  [SAVED] Validation_Summary_",gse_id,".csv\n",sep="")
  
  # ── ROC curves: XGBoost + Directional overlaid ───────────────────────────
  roc_xgb_df <- pROC::coords(roc_xgb,"all",ret=c("specificity","sensitivity"),
                             transpose=FALSE) %>%
    as.data.frame() %>% mutate(fpr=1-specificity,
                               Score="XGBoost Probability")
  roc_dir_df <- pROC::coords(roc_dir,"all",ret=c("specificity","sensitivity"),
                             transpose=FALSE) %>%
    as.data.frame() %>% mutate(fpr=1-specificity,
                               Score="Linear Directional Score")
  
  lbl_xgb <- sprintf("XGBoost (AUC=%.4f [%.4f–%.4f])",
                     as.numeric(ci_xgb[2]),as.numeric(ci_xgb[1]),as.numeric(ci_xgb[3]))
  lbl_dir <- sprintf("Directional (AUC=%.4f [%.4f–%.4f])",
                     as.numeric(ci_dir[2]),as.numeric(ci_dir[1]),as.numeric(ci_dir[3]))
  
  roc_df_both <- rbind(
    roc_xgb_df %>% mutate(Label=lbl_xgb),
    roc_dir_df %>% mutate(Label=lbl_dir)
  )
  roc_df_both$Label <- factor(roc_df_both$Label,
                              levels=c(lbl_xgb,lbl_dir))
  roc_colours <- setNames(c(COL$tumor,"#1B7837"),c(lbl_xgb,lbl_dir))
  
  delong_lbl <- ifelse(delong$p.value<0.001,"< 0.001",
                       sprintf("= %.4f",delong$p.value))
  
  p_roc <- ggplot() +
    geom_ribbon(data=ci_band,
                aes(x=fpr,ymin=lower,ymax=upper),
                fill=COL$tumor,alpha=0.12) +
    geom_line(data=roc_df_both,
              aes(x=fpr,y=sensitivity,colour=Label),linewidth=1.15) +
    geom_abline(slope=1,intercept=0,linetype="dashed",
                colour="grey50",linewidth=0.6) +
    geom_point(aes(x=1-youden_xgb$specificity,y=youden_xgb$sensitivity),
               colour=COL$tumor,size=3.5,shape=18) +
    geom_point(aes(x=1-youden_dir$specificity,y=youden_dir$sensitivity),
               colour="#1B7837",size=3.5,shape=18) +
    scale_colour_manual(values=roc_colours,
                        name="Score  (AUC [95% CI])") +
    scale_x_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
    scale_y_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
    labs(
      title    = paste0(gse_id," \u2014 External Validation ROC"),
      subtitle = sprintf("XGBoost AUC = %.4f  |  Directional AUC = %.4f  |  DeLong p %s",
                         as.numeric(ci_xgb[2]),as.numeric(ci_dir[2]),delong_lbl),
      x = "1 \u2212 Specificity  (False Positive Rate)",
      y = "Sensitivity  (True Positive Rate)",
      caption = paste0(
        "\u25c6 = Youden optimal threshold.  ",
        "Shaded = 95% bootstrap CI on XGBoost sensitivity (",CONFIG$bootstrap_n," reps).\n",
        "Genes present: ",length(present_g),"/",length(genes_us),".  ",
        "Stable genes: ",paste(genes_disp,collapse=", "))
    ) +
    theme_publication() +
    theme(aspect.ratio=1,
          legend.position=c(0.97,0.05),legend.justification=c(1,0),
          legend.background=element_rect(fill="white",colour="black",
                                         linewidth=0.5),
          legend.text=element_text(size=8.5))
  save_figure(p_roc, paste0("Fig_Val_ROC_",gse_id), width=9, height=8.5)
  
  # ── Directional Score Distribution (Violin + Box) ────────────────────────
  score_df <- data.frame(
    Group     = meta_al$Group,
    XGB_Score = xgb_prob,
    Dir_Score = dir_score
  )
  score_long <- score_df %>%
    tidyr::pivot_longer(cols=c(XGB_Score,Dir_Score),
                        names_to="Score_Type",values_to="Score") %>%
    mutate(Score_Type=recode(Score_Type,
                             XGB_Score="XGBoost Probability",
                             Dir_Score="Linear Directional Score"))
  
  wt_xgb <- wilcox.test(XGB_Score~Group, data=score_df)
  wt_dir <- wilcox.test(Dir_Score~Group, data=score_df)
  
  p_violin <- ggplot(score_long,
                     aes(x=Group,y=Score,fill=Group)) +
    geom_violin(trim=FALSE,alpha=0.75,colour="black",linewidth=0.4) +
    geom_jitter(width=0.15,size=0.6,alpha=0.3,colour="black") +
    stat_summary(fun=median,geom="crossbar",
                 width=0.35,colour="white",linewidth=0.7) +
    facet_wrap(~Score_Type, scales="free_y", ncol=2) +
    scale_fill_manual(values=c(Normal=COL$normal,Tumor=COL$tumor),
                      labels=c("Normal","Tumour")) +
    labs(
      title    = paste0(gse_id," \u2014 Score Distribution: Tumour vs Normal"),
      subtitle = sprintf(
        "XGBoost Wilcoxon p %s  |  Directional Score Wilcoxon p %s",
        ifelse(wt_xgb$p.value<0.001,"< 0.001",sprintf("%.4f",wt_xgb$p.value)),
        ifelse(wt_dir$p.value<0.001,"< 0.001",sprintf("%.4f",wt_dir$p.value))),
      x = NULL, y = "Score",
      caption = paste0(
        "White crossbar = median.  ",
        "Directional Score = \u03a3(Up genes) \u2212 \u03a3(Down genes) in log\u2082(TPM+1) units.\n",
        "Both scores derived from the locked signature trained on TCGA-BRCA only.")
    ) +
    theme_publication() +
    theme(legend.position="none",
          axis.text.x=element_text(face="bold",size=11))
  save_figure(p_violin, paste0("Fig_Val_ScoreDist_",gse_id), width=10, height=6.5)
  
  # ── Heatmap ──────────────────────────────────────────────────────────────
  heat_dat  <- t(X_sc)
  rownames(heat_dat) <- genes_disp
  # Add directional annotation to rows
  # Heatmap row annotation: derive direction from locked signature
  # Convert Dir_Symbol arrows to readable labels for pheatmap
  hm_dir_labels <- dplyr::case_when(
    startsWith(dir_sig$Direction[match(genes_disp, dir_sig$Gene)], "\u2191") ~ "Up in Tumour",
    startsWith(dir_sig$Direction[match(genes_disp, dir_sig$Gene)], "\u2193") ~ "Down in Tumour",
    TRUE ~ "No Sig. Change"
  )
  row_anno <- data.frame(Direction=hm_dir_labels, row.names=genes_disp)
  anno_row  <- list(Direction=c("Up in Tumour"=COL$tumor,
                                "Down in Tumour"=COL$normal,
                                "No Sig. Change"="grey70"))
  anno_col  <- data.frame(Group=meta_al$Group, row.names=meta_al$GSM)
  anno_cols <- list(Group=c(Normal=COL$normal, Tumor=COL$tumor))
  heat_pal  <- colorRampPalette(c("navy","white","firebrick3"))(100)
  for (ext in c("png","tiff")) {
    pheatmap::pheatmap(
      mat=heat_dat, color=heat_pal, scale="row",
      annotation_col=anno_col, annotation_colors=c(anno_cols,anno_row),
      annotation_row=row_anno,
      cluster_cols=TRUE, cluster_rows=TRUE,
      show_colnames=FALSE, show_rownames=TRUE, fontsize_row=10,
      filename=file.path(dirs$figures,
                         paste0("Fig_Val_Heatmap_",gse_id,".",ext)),
      width=8, height=max(5,length(genes_disp)*0.4+2)
    )
  }
  cat("  [SAVED] Fig_Val_Heatmap_",gse_id," (PNG + TIFF)\n",sep="")
  
  list(AUC_XGB=as.numeric(ci_xgb[2]),CI_Lo_XGB=as.numeric(ci_xgb[1]),
       CI_Hi_XGB=as.numeric(ci_xgb[3]),
       AUC_Dir=as.numeric(ci_dir[2]),CI_Lo_Dir=as.numeric(ci_dir[1]),
       CI_Hi_Dir=as.numeric(ci_dir[3]),
       ROC_XGB_df=roc_xgb_df, ROC_Dir_df=roc_dir_df,
       Summary=sum_row)
}

# =============================================================================
# 3. RUN ALL DATASETS
# =============================================================================
cat("[2/5] Validating all datasets...\n\n")
all_results <- list()
for (ds in DATASETS) {
  res <- tryCatch(
    run_validation(ds,final_model,stable_genes,clean_genes,
                   preProc_scaler,up_genes,down_genes,CONFIG$seed),
    error=function(e){cat("  [ERROR]",ds$gse_id,":",conditionMessage(e),"\n\n");NULL})
  if (!is.null(res)) all_results[[ds$gse_id]] <- res
}
cat("\n  Completed:",length(all_results),"/",length(DATASETS),"datasets\n\n")

# =============================================================================
# 4. COMBINED ROC — BOTH SCORES ACROSS ALL DATASETS
# =============================================================================
cat("[3/5] Combined ROC plots...\n")
if (length(all_results) > 0) {
  
  ds_names <- names(all_results)
  
  # Combined XGBoost ROC
  comb_xgb <- do.call(rbind, lapply(ds_names, function(nm) {
    r <- all_results[[nm]]
    r$ROC_XGB_df %>% mutate(Dataset=nm, Score="XGBoost")
  }))
  
  # Combined Directional ROC
  comb_dir <- do.call(rbind, lapply(ds_names, function(nm) {
    r <- all_results[[nm]]
    r$ROC_Dir_df %>% mutate(Dataset=nm, Score="Directional")
  }))
  
  build_combined_roc_plot <- function(roc_data, score_label,
                                      ci_lo_col, ci_hi_col, auc_col,
                                      line_col, title_suffix) {
    leg_lbl <- sapply(ds_names, function(nm) {
      r <- all_results[[nm]]
      sprintf("%s\nAUC=%.3f [%.3f\u2013%.3f]",
              nm, r[[auc_col]], r[[ci_lo_col]], r[[ci_hi_col]])
    })
    roc_data$Label <- factor(roc_data$Dataset,levels=ds_names,labels=leg_lbl)
    col_map <- setNames(line_col[seq_along(ds_names)], leg_lbl)
    
    ggplot(roc_data, aes(x=fpr,y=sensitivity,colour=Label)) +
      geom_line(linewidth=1.15) +
      geom_abline(slope=1,intercept=0,linetype="dashed",
                  colour="grey50",linewidth=0.6) +
      scale_colour_manual(values=col_map, name=paste0(score_label,"\n(AUC [95% CI])")) +
      scale_x_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
      scale_y_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
      labs(
        title    = paste0("External Validation \u2014 ",score_label),
        subtitle = paste0("Locked diagnostic engine  |  ",
                          length(stable_genes)," stable genes  |  ",
                          length(all_results)," independent cohort(s)"),
        x = "1 \u2212 Specificity  (False Positive Rate)",
        y = "Sensitivity  (True Positive Rate)",
        caption  = paste0("95% CI by DeLong method.  ",
                          "Individual CI bands in per-dataset figures.\n",
                          "Stable genes: ",
                          stringr::str_wrap(paste(clean_genes,collapse=", "),width=90))
      ) +
      theme_publication() +
      theme(aspect.ratio=1,
            legend.position=c(0.97,0.05),legend.justification=c(1,0),
            legend.background=element_rect(fill="white",colour="black",
                                           linewidth=0.5),
            legend.key.size=unit(0.8,"cm"),legend.text=element_text(size=9))
  }
  
  p_comb_xgb <- build_combined_roc_plot(
    comb_xgb, "XGBoost Probability Score",
    "CI_Lo_XGB","CI_Hi_XGB","AUC_XGB",
    DATASET_COLOURS, "XGBoost")
  save_figure(p_comb_xgb,"Fig_Val_Combined_ROC_XGBoost",width=10,height=9.5)
  
  p_comb_dir <- build_combined_roc_plot(
    comb_dir, "Linear Directional Score",
    "CI_Lo_Dir","CI_Hi_Dir","AUC_Dir",
    c("#2D6A4F","#52B788","#95D5B2","#1B4332","#40916C"),
    "Directional")
  save_figure(p_comb_dir,"Fig_Val_Combined_ROC_Directional",width=10,height=9.5)
  
  # ── Combined summary table ────────────────────────────────────────────────
  cat("[4/5] Saving combined summary...\n")
  comb_sum <- do.call(rbind,lapply(all_results,`[[`,"Summary"))
  write.csv(comb_sum,
            file.path(dirs$val,"Validation_Summary_ALL_Datasets.csv"),
            row.names=FALSE)
  saveRDS(all_results,
          file.path(dirs$val,"02_all_validation_results.rds"))
  cat("  [SAVED] Validation_Summary_ALL_Datasets.csv\n")
  cat("  [SAVED] 02_all_validation_results.rds\n\n")
  
  # ── Score comparison bar chart: XGBoost vs Directional AUC per dataset ───
  auc_compare_df <- do.call(rbind, lapply(ds_names, function(nm) {
    r <- all_results[[nm]]
    data.frame(
      Dataset = rep(nm,2),
      Score   = c("XGBoost Probability","Linear Directional Score"),
      AUC     = c(r$AUC_XGB, r$AUC_Dir),
      CI_Lo   = c(r$CI_Lo_XGB, r$CI_Lo_Dir),
      CI_Hi   = c(r$CI_Hi_XGB, r$CI_Hi_Dir)
    )
  }))
  auc_compare_df$Score <- factor(auc_compare_df$Score,
                                 levels=c("XGBoost Probability",
                                          "Linear Directional Score"))
  
  p_auc_bar <- ggplot(auc_compare_df,
                      aes(x=Dataset,y=AUC,fill=Score,
                          ymin=CI_Lo,ymax=CI_Hi)) +
    geom_col(position=position_dodge(0.75),
             colour="black",linewidth=0.35,width=0.7) +
    geom_errorbar(position=position_dodge(0.75),width=0.2,linewidth=0.7) +
    geom_text(aes(label=sprintf("%.3f",AUC),y=AUC+0.005),
              position=position_dodge(0.75),vjust=-0.5,size=3.3) +
    scale_fill_manual(values=c("XGBoost Probability"=COL$tumor,
                               "Linear Directional Score"="#1B7837"),
                      name="Score Type") +
    scale_y_continuous(limits=c(0,1.1),breaks=seq(0,1,0.2)) +
    geom_hline(yintercept=0.9,linetype="dashed",colour="grey40",linewidth=0.6) +
    annotate("text",x=0.6,y=0.905,label="AUC = 0.90",
             size=3.2,colour="grey40",hjust=0) +
    labs(
      title    = "External Validation: XGBoost vs Directional Score AUC",
      subtitle = "Both scores derived from the locked TCGA-BRCA diagnostic engine",
      x = "Validation Dataset", y = "AUC (95% CI by DeLong)",
      caption  = paste0(
        "Error bars = 95% DeLong CI.  Dashed = AUC threshold of 0.90.\n",
        "Linear Directional Score = \u03a3(Up genes) \u2212 \u03a3(Down genes) ",
        "in log\u2082(TPM+1) units.")
    ) +
    theme_publication()
  save_figure(p_auc_bar,"Fig_Val_AUC_Comparison",width=9,height=6.5)
  
  cat("\n  ── Combined Validation Performance ──\n")
  print(comb_sum[,c("Dataset","N_Total","XGB_AUC","XGB_CI_Lo","XGB_CI_Hi",
                    "Dir_AUC","Dir_CI_Lo","Dir_CI_Hi","DeLong_P_XGBvsDir")])
}

# =============================================================================
# WRAP UP
# =============================================================================
cat("[5/5] Session info...\n")
sink(file.path(dirs$logs,"02_session_info.txt"))
print(sessionInfo()); sink()
sink(log_path, append=TRUE, split=TRUE)

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 02_external_validation.R  COMPLETE\n")
cat(" Runtime:", round(difftime(t_end,t_start,units="mins"),2), "mins\n")
cat("\n SAVED ARTEFACTS:\n")
cat("  Tables  : Validation_Summary_[dataset].csv (one per dataset)\n")
cat("            Validation_Summary_ALL_Datasets.csv\n")
cat("            02_all_validation_results.rds\n")
cat("  Figures : Fig_Val_ROC_[dataset].*         (XGBoost + Directional overlay)\n")
cat("            Fig_Val_ScoreDist_[dataset].*   (violin + score distributions)\n")
cat("            Fig_Val_Heatmap_[dataset].*     (with directional row annotation)\n")
cat("            Fig_Val_Combined_ROC_XGBoost.*\n")
cat("            Fig_Val_Combined_ROC_Directional.*\n")
cat("            Fig_Val_AUC_Comparison.*        (XGBoost vs Directional bar)\n")
cat("            (PNG 300dpi + TIFF 600dpi + PDF for each)\n")
cat("  Logs    : 02_Validation_Run_[timestamp].log\n")
cat("            02_session_info.txt\n")
cat("            Mapping_Log_[dataset].csv\n")
cat("            Gene_Coverage_[dataset].csv\n")
cat("=============================================================\n")
sink()