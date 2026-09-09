# =============================================================================
# Script      : 7_crossplatform_validation.R
# Project     : DSP family dysregulation in breast cancer
# Description : Cross-platform validation of the locked diagnostic engine on
#               independent GEO MICROARRAY datasets. Demonstrates that the
#               Directional Score direction is preserved when moving from the
#               RNA-seq platform used for discovery (TCGA-BRCA, log2(TPM+1))
#               to Affymetrix and Agilent microarray platforms (log2 intensity).
#
#               Datasets validated:
#               (A) GSE70947 — Agilent SurePrint G3 Human GE 8x60K microarray
#                   148 breast adenocarcinomas + 148 paired adjacent normal
#                   Source: UCSF Ashworth Lab, Quigley & Kristensen 2016
#               (B) GSE15852 — Affymetrix Human Genome U133A Array (GPL96)
#                   43 breast tumours + 43 paired normal breast tissues
#                   Source: Institute for Medical Research, Malaysia
#
#               WHY MICROARRAY DATASETS:
#               These datasets test whether the Directional Score — derived
#               from RNA-seq data — retains its discriminative power across a
#               completely different assay technology. The Directional Score is
#               computed in raw log2 microarray intensity units, with NO
#               threshold or scaler transferred from TCGA. The XGBoost model is
#               also applied but results are labelled as cross-platform
#               extrapolation and interpreted with appropriate caution.
#
#               NOTE ON XGBoost IN CROSS-PLATFORM CONTEXT:
#               The preProc_scaler was fitted on TCGA RNA-seq log2(TPM+1).
#               Microarray log2 intensity occupies a different range (~3-15 vs
#               ~0-20 for RNA-seq). The scaler is applied for completeness but
#               XGBoost cross-platform results are interpreted as secondary to
#               the Directional Score results, which require no scaler.
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# Seed        : 42
# Depends on  : 01_training_pipeline.R (must be run first)
# =============================================================================

rm(list=ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed        = 42L,
  bootstrap_n = 500L,
  conf_level  = 0.95
)

DATASETS <- list(
  list(
    gse_id           = "GSE70947",
    platform         = "Agilent SurePrint G3 Human GE 8x60K (GPL13607)",
    platform_short   = "Agilent microarray",
    # Sample titles are "CM016-normal", "CM016-tumor" — match suffix exactly
    normal_keywords  = c("-normal"),
    tumor_keywords   = c("-tumor"),
    exclude_keywords = character(0),
    platform_note    = "Agilent log\u2082 intensity.  148 tumours + 148 paired adjacent normals."
  ),
  list(
    gse_id           = "GSE15852",
    platform         = "Affymetrix Human Genome U133A Array (GPL96)",
    platform_short   = "Affymetrix U133A microarray",
    normal_keywords  = c("normal","non-tumor","nontumor","adjacent","control"),
    tumor_keywords   = c("tumor","tumour","cancer","carcinoma","breast cancer"),
    exclude_keywords = c("cell line"),
    platform_note    = "Affymetrix log\u2082 intensity.  43 breast tumours + 43 paired normals."
  )
)

PATHS <- list(
  models    = here::here("output","models"),
  tables    = here::here("output","tables"),
  geo_cache = here::here("data","geo_cache")
)

# =============================================================================
# 0. PACKAGES
# FIX 1: Biobase added explicitly. fData(), pData(), exprs() are Biobase
#         functions. GEOquery imports Biobase but does not always attach it,
#         causing "could not find function 'fData'" at runtime.
# =============================================================================
if (!requireNamespace("pacman",quietly=TRUE)) install.packages("pacman")
pacman::p_load(
  here, GEOquery, Biobase,
  xgboost, dplyr, tibble, pROC,
  ggplot2, data.table, stringr, tidyr,
  pheatmap, caret, scales
)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

dir.create(PATHS$geo_cache, recursive=TRUE, showWarnings=FALSE)
options(GEOquery.inmemory.gpl=FALSE)

# =============================================================================
# DIRECTORIES & LOGGING
# =============================================================================
dirs <- list(
  models  = PATHS$models,
  tables  = PATHS$tables,
  val     = here::here("output","crossplatform_validation"),
  figures = here::here("output","crossplatform_validation","figures"),
  logs    = here::here("output","crossplatform_validation","logs")
)
lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

log_path <- file.path(dirs$logs,
                      paste0("07_CrossPlatform_Run_",format(Sys.time(),"%Y%m%d_%H%M"),".log"))
sink(log_path, append=FALSE, split=TRUE)
t_start <- Sys.time()
cat("=============================================================\n")
cat(" 07_crossplatform_validation.R\n")
cat(" Start   :", as.character(t_start), "\n")
cat(" Datasets:", length(DATASETS), "\n")
cat("=============================================================\n\n")
set.seed(CONFIG$seed)

# =============================================================================
# SHARED AESTHETICS
# =============================================================================
COL <- list(
  tumor    = "#B22222", normal   = "#2166AC",
  low_expr = "#3288BD", high_expr= "#D53E4F", mid_expr = "#FFFFBF"
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

DATASET_COLOURS <- c("#E69F00","#56B4E9","#009E73",
                     "#CC79A7","#D55E00","#0072B2","#F0E442")

# =============================================================================
# DIRECTIONAL SCORE FUNCTION
# FIX 2: Column names in X_raw_df must be HYPHEN format (HGNC symbols).
#         Underscore conversion happens ONLY at the XGBoost/scaler step.
# =============================================================================
compute_dir_score <- function(X_raw_df, up_g, down_g) {
  up_p  <- intersect(up_g,   colnames(X_raw_df))
  dn_p  <- intersect(down_g, colnames(X_raw_df))
  score <- rep(0, nrow(X_raw_df))
  if (length(up_p) > 0) score <- score + rowSums(X_raw_df[,up_p,drop=FALSE])
  if (length(dn_p) > 0) score <- score - rowSums(X_raw_df[,dn_p,drop=FALSE])
  score
}

map_group <- function(ss, tkw, nkw, exkw) {
  # Exclusion check runs first (cell lines etc.)
  if (length(exkw) > 0 &&
      any(sapply(exkw, function(k) grepl(k,ss,ignore.case=TRUE))))
    return(NA_character_)
  # TUMOUR checked BEFORE NORMAL — a sample matching both is always Tumour.
  # This prevents "adjacent normal tissue from adenocarcinoma patient" from
  # being classified as Normal when it contains tumour keywords too.
  if (any(sapply(tkw, function(k) grepl(k,ss,ignore.case=TRUE))))
    return("Tumor")
  if (any(sapply(nkw, function(k) grepl(k,ss,ignore.case=TRUE))))
    return("Normal")
  NA_character_
}

# =============================================================================
# HELPER: extract_expression_matrix
#
# FIX 2: Receives target_genes_hyp in HYPHEN format.
#         Returns expr_gene with HYPHEN rownames.
#
# FIX 3: Expanded gene symbol column detection:
#         - Affymetrix GPL96/GPL570: "Gene Symbol" (space, capitalised)
#         - Agilent GPL13607:        "GENE_SYMBOL"  (underscore, uppercase)
#         - Older platforms:         various mixed cases
#         Also strips " /// " multi-gene entries and "---" placeholders.
# =============================================================================
extract_expression_matrix <- function(gse_obj, gse_id, target_genes_hyp) {
  
  cat("  Extracting expression matrix from ExpressionSet...\n")
  expr_raw <- Biobase::exprs(gse_obj)
  cat(sprintf("  Raw matrix: %d probes x %d samples\n",
              nrow(expr_raw), ncol(expr_raw)))
  
  feat <- Biobase::fData(gse_obj)
  
  if (ncol(feat) == 0) {
    cat("  [WARN] fData() is empty — AnnotGPL=TRUE may not have loaded.\n")
    cat("  Gene coverage will be 0. Check GEO cache and retry.\n")
    gene_vec <- rownames(expr_raw)
  } else {
    cat("  fData columns:", paste(colnames(feat), collapse=", "), "\n")
    
    # FIX 3: comprehensive column name search
    sym_col <- intersect(
      c("Gene Symbol",       # Affymetrix GPL96, GPL570
        "GENE_SYMBOL",       # Agilent GPL13607
        "gene_symbol",
        "Symbol",
        "Gene_Symbol",
        "SYMBOL",
        "gene symbol",
        "GeneName",
        "gene_assignment",   # some Affymetrix platforms
        "mrna_assignment"),
      colnames(feat))[1]
    
    if (is.na(sym_col))
      sym_col <- grep("^symbol$|gene.symbol|gene_sym",
                      colnames(feat), value=TRUE, ignore.case=TRUE)[1]
    
    if (is.na(sym_col)) {
      cat("  [WARN] No gene symbol column found in fData.\n")
      cat("  ACTION: Add correct column name to sym_col intersect() list.\n")
      cat("  All fData columns:", paste(colnames(feat), collapse=" | "), "\n")
      gene_vec <- rownames(expr_raw)
    } else {
      cat(sprintf("  Gene symbol column: '%s'\n", sym_col))
      gene_vec <- as.character(feat[[sym_col]])
      
      # FIX 3: Handle multi-gene entries e.g. "BRCA1 /// TP53"
      gene_vec <- stringr::str_trim(
        stringr::str_split_fixed(gene_vec, " /// |\\|//", 2)[,1])
      
      # FIX 3: Remove uninformative placeholders
      gene_vec[gene_vec %in% c("---","","NA","N/A","na")] <- NA_character_
      
      n_mapped <- sum(!is.na(gene_vec))
      cat(sprintf("  Probes with symbol: %d / %d (%.1f%%)\n",
                  n_mapped, length(gene_vec),
                  100*n_mapped/length(gene_vec)))
    }
  }
  
  # Collapse duplicate probes → mean per gene symbol
  expr_df      <- as.data.frame(expr_raw)
  expr_df$Gene <- gene_vec
  expr_gene <- expr_df %>%
    filter(!is.na(Gene), nchar(Gene) > 1) %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm=TRUE)),
              .groups="drop") %>%
    tibble::column_to_rownames("Gene")
  
  cat(sprintf("  Unique gene symbols: %d\n", nrow(expr_gene)))
  
  found <- intersect(target_genes_hyp, rownames(expr_gene))
  miss  <- setdiff(target_genes_hyp, rownames(expr_gene))
  cat(sprintf("  Signature coverage: %d / %d genes\n",
              length(found), length(target_genes_hyp)))
  if (length(miss) > 0)
    cat("  [WARN] Missing:", paste(miss, collapse=", "), "\n")
  
  list(expr_gene=expr_gene, found=found, miss=miss)
}

# =============================================================================
# 1. LOAD LOCKED DIAGNOSTIC ENGINE
# =============================================================================
cat("[1/5] Loading locked diagnostic engine...\n")

req_files <- c(
  file.path(dirs$models, "final_xgboost_model_Bayesian.rds"),
  file.path(dirs$models, "final_preProc_scaler.rds"),
  file.path(dirs$tables, "final_stable_biomarker_signature.rds"),
  file.path(dirs$tables, "Final_Directional_Signature.rds")
)
miss_req <- req_files[!file.exists(req_files)]
if (length(miss_req) > 0)
  stop("Missing artefacts:\n", paste(miss_req,collapse="\n"),
       "\nRun 01_training_pipeline.R first.")

final_model    <- readRDS(req_files[1])
preProc_scaler <- readRDS(req_files[2])
stable_genes   <- readRDS(req_files[3])       # underscore (XGBoost convention)
clean_genes    <- gsub("_","-", stable_genes) # hyphen (HGNC display + biology)
n_genes        <- length(stable_genes)
dir_sig        <- readRDS(req_files[4])

up_genes   <- dir_sig$Gene[startsWith(dir_sig$Direction,"\u2191")]
down_genes <- dir_sig$Gene[startsWith(dir_sig$Direction,"\u2193")]
if (length(up_genes)==0 & length(down_genes)==0) {
  cat("  [WARN] Arrow parsing 0 genes — fallback to log2FC sign\n")
  up_genes   <- dir_sig$Gene[dir_sig$log2FC > 0]
  down_genes <- dir_sig$Gene[dir_sig$log2FC < 0]
}
if (length(up_genes)+length(down_genes)==0)
  stop("Could not parse directional genes.")

cat("  Signature genes :", paste(clean_genes, collapse=", "), "\n")
cat("  Up in Tumour    :", paste(up_genes,    collapse=", "), "\n")
cat("  Down in Tumour  :", paste(down_genes,  collapse=", "), "\n\n")

# =============================================================================
# 2. CROSS-PLATFORM VALIDATION ENGINE
# =============================================================================
run_crossplatform_validation <- function(ds, model, genes_us, genes_disp,
                                         scaler, up_g, down_g, seed) {
  gse_id <- ds$gse_id
  cat("──────────────────────────────────────────────\n")
  cat("Dataset:", gse_id, "\n")
  cat("Platform:", ds$platform, "\n")
  
  # Download via GEOquery (cached after first run)
  cat("  Downloading from GEO (cached after first run)...\n")
  gse_obj <- tryCatch(
    GEOquery::getGEO(gse_id,
                     destdir   = PATHS$geo_cache,
                     GSEMatrix = TRUE,
                     AnnotGPL  = TRUE,   # provides gene symbols in fData
                     getGPL    = TRUE)[[1]],
    error=function(e){
      cat("  [ERROR] GEOquery failed:", conditionMessage(e), "\n")
      NULL})
  if (is.null(gse_obj)) return(NULL)
  
  # -- Sample metadata -------------------------------------------------------
  meta_raw <- Biobase::pData(gse_obj)
  
  # Build SearchString from sample-level annotation columns ONLY.
  # CRITICAL: Do NOT include institution/contact/series-level fields that
  # contain words like "Cancer Center" for every sample regardless of group.
  # Priority order: title > source_name > characteristics > description
  # The title column alone is usually sufficient and most reliable.
  title_col <- grep("^title$", colnames(meta_raw),
                    value=TRUE, ignore.case=TRUE)[1]
  char_cols  <- grep("characteristics|source_name|description",
                     colnames(meta_raw), value=TRUE, ignore.case=TRUE)
  
  # Use title as PRIMARY search string — most specific, least contaminated
  if (!is.na(title_col)) {
    meta_raw$SearchString <- tolower(as.character(meta_raw[[title_col]]))
  } else {
    meta_raw$SearchString <- apply(
      meta_raw[, char_cols, drop=FALSE], 1,
      function(x) tolower(paste(x, collapse=" | ")))
  }
  meta_raw$GSM <- rownames(meta_raw)
  
  # Diagnostic: show first 3 SearchStrings to confirm correct column
  cat("  SearchString preview (first 3):\n")
  cat("   ", paste(head(meta_raw$SearchString, 3), collapse="\n    "), "\n")
  
  meta_raw$Group <- sapply(meta_raw$SearchString, map_group,
                           tkw  = ds$tumor_keywords,
                           nkw  = ds$normal_keywords,
                           exkw = ds$exclude_keywords)
  
  # Fallback: if title-only fails to find both groups, try adding characteristics
  n_tumor_raw  <- sum(meta_raw$Group=="Tumor",  na.rm=TRUE)
  n_normal_raw <- sum(meta_raw$Group=="Normal", na.rm=TRUE)
  if ((n_tumor_raw == 0 || n_normal_raw == 0) && length(char_cols) > 0) {
    cat("  [WARN] Title-only mapping failed (T=", n_tumor_raw,
        "N=", n_normal_raw, ") — adding characteristics columns\n")
    meta_raw$SearchString <- apply(
      meta_raw[, c(title_col, char_cols), drop=FALSE], 1,
      function(x) tolower(paste(x, collapse=" | ")))
    meta_raw$Group <- sapply(meta_raw$SearchString, map_group,
                             tkw  = ds$tumor_keywords,
                             nkw  = ds$normal_keywords,
                             exkw = ds$exclude_keywords)
  }
  
  write.csv(meta_raw %>% select(GSM, Group, SearchString),
            file.path(dirs$logs, paste0("Mapping_Log_",gse_id,".csv")),
            row.names=FALSE)
  
  meta <- meta_raw %>% filter(!is.na(Group))
  cat(sprintf("  Tumour: %d  |  Normal: %d\n",
              sum(meta$Group=="Tumor"), sum(meta$Group=="Normal")))
  if (length(unique(meta$Group)) < 2) {
    cat("  [SKIP] Could not identify both groups.\n\n")
    return(NULL)
  }
  
  # -- Expression matrix (FIX 2: pass hyphen genes) --------------------------
  expr_info <- extract_expression_matrix(gse_obj, gse_id, genes_disp)
  expr_gene <- expr_info$expr_gene
  found_g   <- expr_info$found
  miss_g    <- expr_info$miss
  
  if (length(found_g) == 0) {
    cat("  [ERROR] Zero signature genes found. Check fData gene symbol column.\n")
    cat("  First 5 gene names in matrix:",
        paste(head(rownames(expr_gene),5),collapse=", "),"\n\n")
    return(NULL)
  }
  
  # -- Sample alignment ------------------------------------------------------
  common_sam <- intersect(meta$GSM, colnames(expr_gene))
  if (length(common_sam)==0) {
    cat("  [ERROR] No sample overlap.\n\n"); return(NULL)
  }
  meta_al <- meta %>% filter(GSM %in% common_sam)
  expr_al  <- expr_gene[, meta_al$GSM, drop=FALSE]
  cat(sprintf("  Matched samples: %d\n", nrow(meta_al)))
  
  # -- Impute missing genes --------------------------------------------------
  write.csv(
    data.frame(Dataset=gse_id, Platform=ds$platform_short,
               N_Present=length(found_g), N_Imputed=length(miss_g),
               Coverage_Pct=round(length(found_g)/length(genes_disp)*100,1),
               Missing_Genes=paste(miss_g,collapse=";"),
               stringsAsFactors=FALSE),
    file.path(dirs$logs,paste0("Gene_Coverage_",gse_id,".csv")),
    row.names=FALSE)
  
  if (length(miss_g) > 0) {
    cat(sprintf("  [WARN] Imputing %d gene(s) with 0\n",length(miss_g)))
    z <- matrix(0, nrow=length(miss_g), ncol=ncol(expr_al),
                dimnames=list(miss_g,colnames(expr_al)))
    expr_al <- rbind(expr_al, z)
  }
  
  # FIX 5: Align to genes_disp AFTER imputation (miss_g now present in rows)
  X_raw_hyp <- as.data.frame(t(expr_al[genes_disp,,drop=FALSE]))
  # samples x genes, HYPHEN column names throughout
  
  # FIX 4: Log2 guard on FULL GEO expression matrix range
  expr_max_full <- max(Biobase::exprs(gse_obj), na.rm=TRUE)
  if (expr_max_full > 100) {
    cat(sprintf("  [WARN] Max value in full matrix = %.1f > 100 — applying log2(x+1)\n",
                expr_max_full))
    X_raw_hyp <- as.data.frame(log2(as.matrix(X_raw_hyp)+1))
  } else if (expr_max_full > 25) {
    cat(sprintf("  [WARN] Max = %.1f > 25 — applying log2(x+1)\n", expr_max_full))
    X_raw_hyp <- as.data.frame(log2(as.matrix(X_raw_hyp)+1))
  } else {
    cat(sprintf("  Expression range: %.2f to %.2f (log\u2082 confirmed)\n",
                min(X_raw_hyp,na.rm=TRUE), max(X_raw_hyp,na.rm=TRUE)))
  }
  
  y_true <- as.numeric(meta_al$Group=="Tumor")
  
  # -- (A) Directional Score — PRIMARY (FIX 2: hyphen names) ----------------
  dir_score <- compute_dir_score(X_raw_hyp, up_g, down_g)
  cat(sprintf("  Dir Score: Tumour mean=%.3f  Normal mean=%.3f\n",
              mean(dir_score[y_true==1]), mean(dir_score[y_true==0])))
  
  # -- (B) XGBoost — SECONDARY (FIX 2: underscore conversion HERE only) -----
  X_raw_us           <- X_raw_hyp
  colnames(X_raw_us) <- gsub("-","_",colnames(X_raw_hyp))
  X_sc <- tryCatch(
    predict(scaler, as.data.frame(X_raw_us)),
    error=function(e){ cat("  [WARN] Scaler failed.\n"); as.data.frame(X_raw_us) })
  xgb_prob <- tryCatch(
    predict(model, xgboost::xgb.DMatrix(data=as.matrix(X_sc))),
    error=function(e){ cat("  [WARN] XGBoost failed.\n"); rep(NA_real_,nrow(X_sc)) })
  
  # -- ROC -------------------------------------------------------------------
  roc_dir <- pROC::roc(y_true, dir_score,
                       levels=c(0,1), direction="<", quiet=TRUE)
  ci_dir  <- pROC::ci.auc(roc_dir, conf.level=CONFIG$conf_level,
                          method="delong")
  cat(sprintf("  Directional AUC: %.4f [%.4f\u2013%.4f]\n",
              as.numeric(ci_dir[2]),as.numeric(ci_dir[1]),as.numeric(ci_dir[3])))
  
  # FIX 8: Initialise all XGBoost objects as NULL before conditional block
  roc_xgb <- NULL; ci_xgb <- NULL; delong <- NULL; youden_xgb <- NULL; lbl_xgb <- NULL
  if (!all(is.na(xgb_prob))) {
    roc_xgb    <- pROC::roc(y_true, xgb_prob, levels=c(0,1),
                            direction="<", quiet=TRUE)
    ci_xgb     <- pROC::ci.auc(roc_xgb, conf.level=CONFIG$conf_level,
                               method="delong")
    delong     <- pROC::roc.test(roc_dir, roc_xgb, method="delong")
    youden_xgb <- as.data.frame(
      pROC::coords(roc_xgb,"best",best.method="youden",
                   ret=c("specificity","sensitivity","threshold"),
                   transpose=FALSE))[1,]
    lbl_xgb    <- sprintf("XGBoost* (AUC=%.4f [%.4f\u2013%.4f])",
                          as.numeric(ci_xgb[2]),
                          as.numeric(ci_xgb[1]),
                          as.numeric(ci_xgb[3]))
    cat(sprintf("  XGBoost AUC (extrap.): %.4f [%.4f\u2013%.4f]  DeLong p=%s\n",
                as.numeric(ci_xgb[2]),as.numeric(ci_xgb[1]),
                as.numeric(ci_xgb[3]),
                ifelse(delong$p.value<0.001,"<0.001",
                       sprintf("%.4f",delong$p.value))))
  }
  
  # Bootstrap CI band (Directional Score)
  set.seed(seed)
  ci_se_dir <- pROC::ci.se(roc_dir, specificities=seq(0,1,by=0.04),
                           conf.level=CONFIG$conf_level, method="bootstrap",
                           boot.n=CONFIG$bootstrap_n, quiet=TRUE)
  ci_band_dir <- data.frame(fpr=1-seq(0,1,by=0.04),
                            lower=as.numeric(ci_se_dir[,1]),
                            upper=as.numeric(ci_se_dir[,3]))
  
  youden_dir <- as.data.frame(
    pROC::coords(roc_dir,"best",best.method="youden",
                 ret=c("specificity","sensitivity","threshold"),
                 transpose=FALSE))[1,]
  
  # -- Confusion matrix metrics ----------------------------------------------
  get_cm_metrics <- function(probs,y,thr){
    pred_cl <- factor(ifelse(probs>=thr,"Tumor","Normal"),
                      levels=c("Normal","Tumor"))
    true_cl <- factor(ifelse(y==1,"Tumor","Normal"),
                      levels=c("Normal","Tumor"))
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
  m_dir <- get_cm_metrics(dir_score, y_true, youden_dir$threshold)
  
  # -- Summary ---------------------------------------------------------------
  sum_row <- data.frame(
    Dataset=gse_id, Platform=ds$platform_short,
    N_Total=nrow(meta_al), N_Tumor=sum(meta_al$Group=="Tumor"),
    N_Normal=sum(meta_al$Group=="Normal"),
    Genes_Present=length(found_g), Genes_Imputed=length(miss_g),
    Dir_AUC    =round(as.numeric(ci_dir[2]),4),
    Dir_CI_Lo  =round(as.numeric(ci_dir[1]),4),
    Dir_CI_Hi  =round(as.numeric(ci_dir[3]),4),
    Dir_Sensitivity=round(m_dir["Sensitivity"],4),
    Dir_Specificity=round(m_dir["Specificity"],4),
    Dir_F1     =round(m_dir["F1"],4),
    Dir_MCC    =round(m_dir["MCC"],4),
    XGB_AUC    =if(!is.null(ci_xgb)) round(as.numeric(ci_xgb[2]),4) else NA,
    XGB_CI_Lo  =if(!is.null(ci_xgb)) round(as.numeric(ci_xgb[1]),4) else NA,
    XGB_CI_Hi  =if(!is.null(ci_xgb)) round(as.numeric(ci_xgb[3]),4) else NA,
    DeLong_P   =if(!is.null(delong))  round(delong$p.value,4)        else NA,
    XGB_Note   ="Platform extrapolation — TCGA-fitted scaler on microarray",
    stringsAsFactors=FALSE)
  write.csv(sum_row,
            file.path(dirs$val,paste0("CrossPlatform_Summary_",gse_id,".csv")),
            row.names=FALSE)
  cat("  [SAVED] CrossPlatform_Summary_",gse_id,".csv\n",sep="")
  
  # ==========================================================================
  # FIGURES
  # ==========================================================================
  
  # -- Figure 1: ROC ---------------------------------------------------------
  roc_dir_df <- pROC::coords(roc_dir,"all",ret=c("specificity","sensitivity"),
                             transpose=FALSE) %>%
    as.data.frame() %>% mutate(fpr=1-specificity)
  lbl_dir <- sprintf("Directional Score (AUC=%.4f [%.4f\u2013%.4f])",
                     as.numeric(ci_dir[2]),as.numeric(ci_dir[1]),
                     as.numeric(ci_dir[3]))
  roc_plot_df <- roc_dir_df %>% mutate(Label=lbl_dir)
  roc_colours <- setNames("#1B7837", lbl_dir)
  
  # FIX 8: lbl_xgb is NULL unless XGBoost succeeded — safe conditional
  if (!is.null(roc_xgb) && !is.null(lbl_xgb)) {
    roc_xgb_df  <- pROC::coords(roc_xgb,"all",ret=c("specificity","sensitivity"),
                                transpose=FALSE) %>%
      as.data.frame() %>% mutate(fpr=1-specificity)
    roc_plot_df <- rbind(roc_plot_df, roc_xgb_df %>% mutate(Label=lbl_xgb))
    roc_colours <- c(roc_colours, setNames(COL$tumor, lbl_xgb))
  }
  roc_plot_df$Label <- factor(roc_plot_df$Label, levels=names(roc_colours))
  
  p_roc <- ggplot() +
    geom_ribbon(data=ci_band_dir, aes(x=fpr,ymin=lower,ymax=upper),
                fill="#1B7837", alpha=0.12) +
    geom_line(data=roc_plot_df, aes(x=fpr,y=sensitivity,colour=Label),
              linewidth=1.15) +
    geom_abline(slope=1,intercept=0,linetype="dashed",
                colour="grey50",linewidth=0.6) +
    geom_point(aes(x=1-youden_dir$specificity, y=youden_dir$sensitivity),
               colour="#1B7837", size=3.5, shape=18) +
    # FIX 8: youden_xgb checked for NULL before reference
    {if (!is.null(youden_xgb))
      geom_point(aes(x=1-youden_xgb$specificity, y=youden_xgb$sensitivity),
                 colour=COL$tumor, size=3.5, shape=18)} +
    scale_colour_manual(values=roc_colours, name="Score  (AUC [95% CI])") +
    scale_x_continuous(limits=c(0,1), expand=c(0.01,0.01)) +
    scale_y_continuous(limits=c(0,1), expand=c(0.01,0.01)) +
    labs(
      title    = paste0(gse_id," \u2014 Cross-Platform Validation ROC"),
      subtitle = sprintf(
        "Directional AUC = %.4f  |  Platform: %s  |  n=%d tumours, %d normals",
        as.numeric(ci_dir[2]), ds$platform_short,
        sum(meta_al$Group=="Tumor"), sum(meta_al$Group=="Normal")),
      x = "1 \u2212 Specificity  (False Positive Rate)",
      y = "Sensitivity  (True Positive Rate)",
      caption = paste0(
        "\u25c6 = Youden optimal threshold.  Shaded = 95% bootstrap CI (",
        CONFIG$bootstrap_n," reps).\n",
        "Genes present: ",length(found_g),"/",length(genes_disp),".\n",
        "Directional Score = \u03a3(Up) \u2212 \u03a3(Down) in log\u2082 microarray intensity.\n",
        if(!is.null(roc_xgb))
          "* XGBoost: TCGA-fitted scaler on microarray \u2014 platform extrapolation."
        else "")
    ) +
    theme_publication() +
    theme(aspect.ratio=1,
          legend.position=c(0.97,0.05), legend.justification=c(1,0),
          legend.background=element_rect(fill="white",colour="black",
                                         linewidth=0.5),
          legend.text=element_text(size=8.5))
  save_figure(p_roc, paste0("Fig_CP_ROC_",gse_id), width=9, height=8.5)
  
  # -- Figure 2: Score violin ------------------------------------------------
  score_df <- data.frame(Group=meta_al$Group, Dir_Score=dir_score)
  wt_dir   <- wilcox.test(Dir_Score~Group, data=score_df)
  
  p_violin <- ggplot(score_df, aes(x=Group, y=Dir_Score, fill=Group)) +
    geom_violin(trim=FALSE, alpha=0.75, colour="black", linewidth=0.4) +
    geom_jitter(width=0.15, size=0.6, alpha=0.3, colour="black") +
    stat_summary(fun=median, geom="crossbar",
                 width=0.35, colour="white", linewidth=0.7) +
    scale_fill_manual(values=c(Normal=COL$normal, Tumor=COL$tumor),
                      labels=c("Normal","Tumour")) +
    labs(
      title    = paste0(gse_id," \u2014 Directional Score: Tumour vs Normal"),
      subtitle = sprintf(
        "Wilcoxon p %s  |  Platform: %s  |  n=%d tumours, %d normals",
        ifelse(wt_dir$p.value<0.001,"< 0.001",sprintf("= %.4f",wt_dir$p.value)),
        ds$platform_short,
        sum(meta_al$Group=="Tumor"), sum(meta_al$Group=="Normal")),
      x=NULL, y="Directional Score  [\u03a3(Up) \u2212 \u03a3(Down)]",
      caption=paste0("White crossbar = median.  ",
                     "Score in raw log\u2082 microarray intensity.  ",
                     "No TCGA scaler.\n", ds$platform_note)
    ) +
    theme_publication() +
    theme(legend.position="none",
          axis.text.x=element_text(face="bold",size=12))
  save_figure(p_violin, paste0("Fig_CP_ScoreDist_",gse_id), width=8, height=6.5)
  
  # -- Figure 3: Heatmap (FIX 6: fallback to X_raw_us if X_sc has NA) -------
  heat_input <- if (anyNA(X_sc)) X_raw_us else X_sc
  heat_dat   <- t(heat_input)
  rownames(heat_dat) <- genes_disp
  
  hm_dir_labels <- dplyr::case_when(
    startsWith(dir_sig$Direction[match(genes_disp,dir_sig$Gene)],"\u2191") ~
      "Up in Tumour",
    startsWith(dir_sig$Direction[match(genes_disp,dir_sig$Gene)],"\u2193") ~
      "Down in Tumour",
    TRUE ~ "No Sig. Change")
  
  row_anno  <- data.frame(Direction=hm_dir_labels, row.names=genes_disp)
  anno_row  <- list(Direction=c("Up in Tumour"=COL$tumor,
                                "Down in Tumour"=COL$normal,
                                "No Sig. Change"="grey70"))
  anno_col  <- data.frame(Group=meta_al$Group, row.names=meta_al$GSM)
  anno_cols <- list(Group=c(Normal=COL$normal, Tumor=COL$tumor))
  heat_pal  <- colorRampPalette(c("navy","white","firebrick3"))(100)
  
  for (ext in c("png","tiff")) {
    pheatmap::pheatmap(
      mat=heat_dat, color=heat_pal, scale="row",
      annotation_col=anno_col, annotation_row=row_anno,
      annotation_colors=c(anno_cols,anno_row),
      cluster_cols=TRUE, cluster_rows=TRUE,
      show_colnames=FALSE, show_rownames=TRUE, fontsize_row=10,
      filename=file.path(dirs$figures,
                         paste0("Fig_CP_Heatmap_",gse_id,".",ext)),
      width=8, height=max(5,n_genes*0.4+2),
      main=paste0(gse_id," \u2014 Signature Heatmap\n",ds$platform_short))
  }
  cat("  [SAVED] Fig_CP_Heatmap_",gse_id," (PNG + TIFF)\n",sep="")
  
  # -- Return ----------------------------------------------------------------
  list(
    AUC_Dir    = as.numeric(ci_dir[2]),
    CI_Lo_Dir  = as.numeric(ci_dir[1]),
    CI_Hi_Dir  = as.numeric(ci_dir[3]),
    AUC_XGB    = if(!is.null(ci_xgb)) as.numeric(ci_xgb[2]) else NA,
    CI_Lo_XGB  = if(!is.null(ci_xgb)) as.numeric(ci_xgb[1]) else NA,
    CI_Hi_XGB  = if(!is.null(ci_xgb)) as.numeric(ci_xgb[3]) else NA,
    ROC_Dir_df = roc_dir_df,
    ROC_XGB_df = if(!is.null(roc_xgb))
      pROC::coords(roc_xgb,"all",ret=c("specificity","sensitivity"),
                   transpose=FALSE) %>%
      as.data.frame() %>% mutate(fpr=1-specificity)
    else NULL,
    Summary    = sum_row,
    Platform   = ds$platform_short
  )
}

# =============================================================================
# 3. RUN ALL DATASETS
# =============================================================================
cat("[2/5] Validating all microarray datasets...\n\n")
all_results <- list()
for (ds in DATASETS) {
  res <- tryCatch(
    run_crossplatform_validation(
      ds, final_model, stable_genes, clean_genes,
      preProc_scaler, up_genes, down_genes, CONFIG$seed),
    error=function(e){
      cat("  [ERROR]",ds$gse_id,":",conditionMessage(e),"\n\n"); NULL})
  if (!is.null(res)) all_results[[ds$gse_id]] <- res
}
cat("\n  Completed:",length(all_results),"/",length(DATASETS),"datasets\n\n")

# =============================================================================
# 4. COMBINED FIGURES AND SUMMARY
# =============================================================================
cat("[3/5] Combined figures...\n")

if (length(all_results) > 0) {
  
  ds_names <- names(all_results)
  
  # -- Combined ROC ----------------------------------------------------------
  comb_dir <- do.call(rbind, lapply(ds_names, function(nm) {
    r <- all_results[[nm]]
    r$ROC_Dir_df %>% mutate(Dataset=nm, Platform=r$Platform)
  }))
  
  # FIX 7: Build label lookup by name, map via Dataset column
  leg_lbl_dir <- setNames(
    sapply(ds_names, function(nm) {
      r <- all_results[[nm]]
      sprintf("%s [%s]\nAUC=%.3f [%.3f\u2013%.3f]",
              nm, r$Platform, r$AUC_Dir, r$CI_Lo_Dir, r$CI_Hi_Dir)
    }), ds_names)
  
  comb_dir$Label <- factor(leg_lbl_dir[comb_dir$Dataset],
                           levels=unname(leg_lbl_dir))
  col_map_dir    <- setNames(DATASET_COLOURS[seq_along(ds_names)],
                             unname(leg_lbl_dir))
  
  p_comb_dir <- ggplot(comb_dir, aes(x=fpr, y=sensitivity, colour=Label)) +
    geom_line(linewidth=1.15) +
    geom_abline(slope=1,intercept=0,linetype="dashed",
                colour="grey50",linewidth=0.6) +
    scale_colour_manual(values=col_map_dir,
                        name="Dataset  (AUC [95% CI])") +
    scale_x_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
    scale_y_continuous(limits=c(0,1),expand=c(0.01,0.01)) +
    labs(
      title    = "Cross-Platform Validation \u2014 Directional Score ROC",
      subtitle = paste0("Locked engine  |  ",n_genes," stable genes  |  ",
                        length(all_results)," microarray cohort(s)"),
      x = "1 \u2212 Specificity  (False Positive Rate)",
      y = "Sensitivity  (True Positive Rate)",
      caption = paste0(
        "Directional Score = \u03a3(Up) \u2212 \u03a3(Down) in log\u2082 microarray intensity.\n",
        "No TCGA scaler applied.  Platform-invariant score.  95% CI by DeLong.\n",
        "Stable genes: ",
        stringr::str_wrap(paste(clean_genes,collapse=", "),width=90))
    ) +
    theme_publication() +
    theme(aspect.ratio=1,
          legend.position=c(0.97,0.05), legend.justification=c(1,0),
          legend.background=element_rect(fill="white",colour="black",
                                         linewidth=0.5),
          legend.key.size=unit(0.8,"cm"),
          legend.text=element_text(size=9))
  save_figure(p_comb_dir,"Fig_CP_Combined_ROC_Directional",width=10,height=9.5)
  
  # -- AUC bar chart (publication quality) ----------------------------------
  auc_bar_df <- do.call(rbind, lapply(ds_names, function(nm) {
    r <- all_results[[nm]]
    rbind(
      data.frame(Dataset=nm, Platform=r$Platform,
                 Score="Directional Score",
                 AUC=r$AUC_Dir, CI_Lo=r$CI_Lo_Dir, CI_Hi=r$CI_Hi_Dir,
                 stringsAsFactors=FALSE),
      if (!is.na(r$AUC_XGB))
        data.frame(Dataset=nm, Platform=r$Platform,
                   Score="XGBoost (platform extrap.)",
                   AUC=r$AUC_XGB, CI_Lo=r$CI_Lo_XGB, CI_Hi=r$CI_Hi_XGB,
                   stringsAsFactors=FALSE)
      else NULL)
  }))
  auc_bar_df$Score <- factor(auc_bar_df$Score,
                             levels=c("Directional Score",
                                      "XGBoost (platform extrap.)"))
  # Concise x-axis label: dataset ID only, platform in subtitle
  auc_bar_df$Dataset_Label <- auc_bar_df$Dataset
  
  # Compute y-ceiling dynamically so labels always have breathing room
  y_ceil <- min(1.05, max(auc_bar_df$CI_Hi, na.rm=TRUE) + 0.12)
  
  # Colour alpha: Directional Score = full; XGBoost = slightly muted
  bar_colours <- c("Directional Score"       = "#1B7837",
                   "XGBoost (platform extrap.)" = "#B22222")
  bar_alphas  <- c("Directional Score"       = 1.0,
                   "XGBoost (platform extrap.)" = 0.75)
  
  p_auc_bar <- ggplot(auc_bar_df,
                      aes(x=Dataset_Label, y=AUC,
                          fill=Score, alpha=Score,
                          ymin=CI_Lo, ymax=CI_Hi)) +
    # Reference band for AUC >= 0.90 (shaded, not just a line)
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.90, ymax=y_ceil,
             fill="grey92", alpha=0.6) +
    geom_col(position=position_dodge(0.72),
             colour="grey20", linewidth=0.4, width=0.68) +
    # Error bars in matching dark shade, slightly narrower than bar
    geom_errorbar(aes(colour=Score),
                  position=position_dodge(0.72),
                  width=0.18, linewidth=0.9) +
    # AUC value labels: placed ABOVE the CI upper bound, not above bar top
    geom_text(aes(y=CI_Hi + 0.022,
                  label=sprintf("%.3f", AUC),
                  colour=Score),
              position=position_dodge(0.72),
              size=3.6, fontface="bold", vjust=0) +
    # Reference line at AUC = 0.90
    geom_hline(yintercept=0.90, linetype="dashed",
               colour="grey45", linewidth=0.7) +
    # Label for the reference line — placed at right edge, BELOW the line
    annotate("text",
             x=length(ds_names) + 0.52,
             y=0.885,
             label="AUC = 0.90",
             size=3.1, colour="grey40", hjust=1, fontface="italic") +
    scale_fill_manual(values=bar_colours,
                      name=NULL,
                      guide=guide_legend(override.aes=list(alpha=1))) +
    scale_colour_manual(values=c("Directional Score"       = "#145228",
                                 "XGBoost (platform extrap.)" = "#7B1010"),
                        guide="none") +
    scale_alpha_manual(values=bar_alphas, guide="none") +
    scale_y_continuous(
      limits  = c(0, y_ceil),
      breaks  = seq(0, 1, 0.2),
      expand  = c(0, 0),
      labels  = scales::number_format(accuracy=0.1)) +
    scale_x_discrete(expand=expansion(add=0.6)) +
    labs(
      title    = "Cross-Platform Validation \u2014 AUC by Dataset",
      subtitle = paste0(
        "Directional Score (primary) vs XGBoost (secondary)  |  ",
        paste(sapply(ds_names, function(nm) all_results[[nm]]$Platform),
              collapse=" and ")),
      x = NULL,
      y = "AUC  (95% CI by DeLong)",
      caption = paste0(
        "Shaded region = AUC \u2265 0.90.  Error bars = 95% DeLong CI.\n",
        "XGBoost uses TCGA-fitted scaler applied to microarray \u2014 ",
        "platform extrapolation; interpret with caution.\n",
        "Directional Score = \u03a3(Up genes) \u2212 \u03a3(Down genes) ",
        "in log\u2082 intensity \u2014 no scaler required.")
    ) +
    theme_publication(base_size=13) +
    theme(
      legend.position      = "top",
      legend.justification = "left",
      legend.key.size      = unit(0.55,"cm"),
      legend.text          = element_text(size=11),
      legend.margin        = margin(b=4),
      axis.text.x          = element_text(face="bold", size=12,
                                          colour="grey10"),
      axis.ticks.x         = element_blank(),
      panel.grid.major.y   = element_line(colour="grey90", linewidth=0.4),
      panel.grid.major.x   = element_blank(),
      plot.subtitle        = element_text(size=10)
    )
  save_figure(p_auc_bar,"Fig_CP_AUC_Comparison",width=9,height=6.5)
  
  # -- Summary table ---------------------------------------------------------
  cat("[4/5] Saving combined summary...\n")
  comb_sum <- do.call(rbind,lapply(all_results,`[[`,"Summary"))
  write.csv(comb_sum,
            file.path(dirs$val,"CrossPlatform_Summary_ALL_Datasets.csv"),
            row.names=FALSE)
  saveRDS(all_results,
          file.path(dirs$val,"07_all_crossplatform_results.rds"))
  cat("  [SAVED] CrossPlatform_Summary_ALL_Datasets.csv\n")
  cat("  [SAVED] 07_all_crossplatform_results.rds\n\n")
  
  cat("\n  ── Cross-Platform Performance Summary ──\n")
  print(comb_sum[,c("Dataset","Platform","N_Total","N_Tumor","N_Normal",
                    "Genes_Present","Dir_AUC","Dir_CI_Lo","Dir_CI_Hi",
                    "Dir_Sensitivity","Dir_Specificity","Dir_F1")])
}

# =============================================================================
# WRAP UP
# =============================================================================
cat("[5/5] Session info...\n")
sink(file.path(dirs$logs,"07_session_info.txt")); print(sessionInfo()); sink()
sink(log_path,append=TRUE,split=TRUE)

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 07_crossplatform_validation.R  COMPLETE\n")
cat(" Runtime:",round(difftime(t_end,t_start,units="mins"),2),"mins\n")
cat("\n SCIENTIFIC NOTE:\n")
cat("   PRIMARY metric: Directional Score — no scaler, platform-invariant.\n")
cat("   SECONDARY metric: XGBoost — TCGA scaler on microarray = extrapolation.\n")
cat("\n SAVED ARTEFACTS:\n")
cat("  Tables  : CrossPlatform_Summary_[dataset].csv\n")
cat("            CrossPlatform_Summary_ALL_Datasets.csv\n")
cat("            07_all_crossplatform_results.rds\n")
cat("  Figures : Fig_CP_ROC_[dataset].*\n")
cat("            Fig_CP_ScoreDist_[dataset].*\n")
cat("            Fig_CP_Heatmap_[dataset].*\n")
cat("            Fig_CP_Combined_ROC_Directional.*\n")
cat("            Fig_CP_AUC_Comparison.*\n")
cat("            (PNG 300dpi + TIFF 600dpi + PDF)\n")
cat("  Logs    : 07_CrossPlatform_Run_[timestamp].log\n")
cat("            Mapping_Log_[dataset].csv\n")
cat("            Gene_Coverage_[dataset].csv\n")
cat("=============================================================\n")
sink()