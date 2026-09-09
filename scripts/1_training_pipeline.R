# =============================================================================
# Script      : 1_training_pipeline.R
# Project     : DSP family dysregulation in breast cancer
# Description : Repeated Nested CV (100 folds), 5-method ensemble feature
#               selection, 4-tier consensus, 4-model soft-vote ensemble,
#               Bayesian-optimised XGBoost final model, full SHAP suite.
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# Seed        : 42
# =============================================================================
# REQUIRED INPUTS (place in data/):
#   data/TCGA-BRCA.star_tpm.tsv.gz   — TCGA GDC star_tpm (already log2(TPM+1))
#   data/BRCA_gene_list.txt           — one HGNC gene symbol per line
#
# KEY DESIGN NOTES:
#   - TCGA GDC star_tpm files are log2(TPM+1); NO additional log2 transform.
#   - Gene names stored with UNDERSCORES internally (XGBoost requirement).
#     Hyphens restored for all display labels.
#   - preProc_scaler fitted on TCGA-BRCA training data ONLY; applied
#     identically to all external datasets in downstream scripts.
#   - BEST PRACTICE: Biological Direction (DEG) dictates the signature direction;
#     XGBoost SHAP determines feature importance magnitude.
# =============================================================================

rm(list = ls()); gc()

# =============================================================================
# USER CONFIGURATION — edit this block only
# =============================================================================
CONFIG <- list(
  seed             = 42L,
  outer_folds      = 5L,
  repeats          = 20L,      # 5 × 20 = 100 total evaluation folds
  inner_folds      = 3L,
  tune_length      = 5L,
  stability_cut    = 0.80,     # gene must appear in ≥80% of folds
  consensus_votes  = 4L,       # primary threshold (≥4/5 methods)
  bayes_init       = 5L,       # Bayesian optimisation initial points
  bayes_iter       = 10L,      # Bayesian optimisation iterations
  n_cores_leave    = 2L        # cores to leave free for OS
)

PATHS <- list(
  data_file  = file.path("data", "TCGA-BRCA.star_tpm.tsv.gz"),
  gene_list  = file.path("data", "BRCA_gene_list.txt")
)

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here, data.table, dplyr, stringr, tibble,
  caret, glmnet, ranger, xgboost, e1071, kernlab,
  Boruta, pROC, ggplot2, ggbeeswarm, uwot, Rtsne,
  SHAPforxgboost, doParallel, foreach, parallel,
  ParBayesianOptimization, UpSetR,
  AnnotationDbi, org.Hs.eg.db,
  scales, tools
)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# =============================================================================
# DIRECTORIES
# =============================================================================
dirs <- list(
  data    = here("data"),
  output  = here("output"),
  figures = here("output", "figures"),
  tables  = here("output", "tables"),
  models  = here("output", "models"),
  logs    = here("output", "logs")
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOGGING
# =============================================================================
log_path <- file.path(dirs$logs,
                      paste0("01_Pipeline_Run_", format(Sys.time(), "%Y%m%d_%H%M"), ".log"))
sink(log_path, append = FALSE, split = TRUE)
t_start <- Sys.time()

cat("=============================================================\n")
cat(" 01_training_pipeline.R\n")
cat(" Start  :", as.character(t_start), "\n")
cat(" Seed   :", CONFIG$seed, "\n")
cat(" Folds  :", CONFIG$outer_folds, "× repeats", CONFIG$repeats,
    "=", CONFIG$outer_folds * CONFIG$repeats, "total folds\n")
cat("=============================================================\n\n")

set.seed(CONFIG$seed)

# =============================================================================
# SHARED AESTHETICS — identical across all 5 scripts
# =============================================================================

# Colour palette
COL <- list(
  tumor    = "#B22222",
  normal   = "#2166AC",
  stable   = "#B22222",
  unstable = "grey70",
  low_expr = "#3288BD",
  high_expr= "#D53E4F",
  mid_expr = "#FFFFBF"
)

# Publication ggplot2 theme
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey30", size = base_size - 1, margin = margin(b = 5)),
      plot.caption     = element_text(size = base_size - 3.5, colour = "grey45", hjust = 0, lineheight = 1.3, margin = margin(t = 8)),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      plot.margin      = unit(c(0.8, 1.2, 0.8, 0.8), "cm")
    )
}

# Multi-format figure save: PNG (300dpi) + TIFF (600dpi LZW) + PDF (vector)
save_figure <- function(plot_obj, filename_base, width = 9, height = 8, out_dir = dirs$figures) {
  base <- file.path(out_dir, filename_base)
  ggplot2::ggsave(paste0(base, ".png"),  plot = plot_obj, width = width, height = height, dpi = 300)
  ggplot2::ggsave(paste0(base, ".tiff"), plot = plot_obj, width = width, height = height, dpi = 600, device = "tiff", compression = "lzw")
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = width, height = height)
  print(plot_obj); grDevices::dev.off()
  cat("  [SAVED]", filename_base, "(PNG 300dpi + TIFF 600dpi + PDF vector)\n")
  invisible(plot_obj)
}

# =============================================================================
# 1. DATA LOADING & PREPROCESSING
# =============================================================================
cat("[1/9] Loading and formatting data...\n")

target_genes <- readLines(PATHS$gene_list)
cat("  Gene list loaded:", length(target_genes), "candidate genes\n")

df <- data.table::fread(PATHS$data_file, header = TRUE, sep = "\t")
if (!"gene_id" %in% names(df)) names(df)[1] <- "gene_id"
df$gene_id_clean <- str_split_fixed(df$gene_id, "\\.", 2)[, 1]

suppressMessages({
  df$gene_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = df$gene_id_clean, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
  )
})

master_data <- df %>%
  filter(!is.na(gene_symbol), gene_symbol %in% target_genes) %>%
  group_by(gene_symbol) %>%
  summarise(across(starts_with("TCGA"), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames("gene_symbol") %>%
  as.matrix() %>% t() %>% as.data.frame() %>%
  mutate(sample_type_code = substr(rownames(.), 14, 15)) %>%
  filter(sample_type_code %in% c("01", "11")) %>%
  mutate(Group = factor(ifelse(sample_type_code == "01", "Tumor", "Normal"), levels = c("Normal", "Tumor"))) %>%
  select(-sample_type_code) %>%
  na.omit()

# CRITICAL: underscore names match XGBoost feature naming convention
colnames(master_data) <- gsub("-", "_", colnames(master_data))

n_tumor  <- sum(master_data$Group == "Tumor")
n_normal <- sum(master_data$Group == "Normal")

cat("  Samples loaded:", nrow(master_data), "| Tumor:", n_tumor, "| Normal:", n_normal, "| Genes:", ncol(master_data) - 1, "\n")
cat("  Class ratio (T/N):", round(n_tumor / n_normal, 3), "\n\n")

# Save a snapshot of the raw feature matrix for reproducibility audit
saveRDS(master_data, file.path(dirs$tables, "01_master_data_snapshot.rds"))
cat("  [SAVED] output/tables/01_master_data_snapshot.rds\n\n")

# =============================================================================
# 2. EXPLORATORY VISUALISATION (UMAP + t-SNE)
# =============================================================================
cat("[2/9] Generating UMAP and t-SNE projections...\n")

X_full <- master_data %>% select(-Group)
y_full <- master_data$Group

set.seed(CONFIG$seed)
umap_res <- uwot::umap(X_full, n_neighbors = 15, min_dist = 0.1, n_threads = 1)
p_umap <- ggplot(data.frame(U1 = umap_res[, 1], U2 = umap_res[, 2], Group = y_full), aes(U1, U2, colour = Group)) +
  geom_point(alpha = 0.75, size = 1.4) +
  scale_colour_manual(values = c(Normal = COL$normal, Tumor = COL$tumor), labels = c("Normal", "Tumour")) +
  labs(title = "UMAP Projection: Tumour vs Normal", subtitle = paste0("TCGA-BRCA  |  n = ", nrow(master_data), "  (", n_tumor, " Tumour, ", n_normal, " Normal)"), x = "UMAP 1", y = "UMAP 2", caption = paste0("log\u2082(TPM+1) expression of ", ncol(X_full), " candidate genes.  UMAP: n_neighbors=15, min_dist=0.1.")) +
  theme_publication() + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
save_figure(p_umap, "Fig01A_UMAP", width = 7, height = 6)

set.seed(CONFIG$seed)
tsne_res <- Rtsne::Rtsne(as.matrix(X_full), dims = 2, perplexity = min(30, floor(nrow(X_full) / 3)), check_duplicates = FALSE)
p_tsne <- ggplot(data.frame(T1 = tsne_res$Y[, 1], T2 = tsne_res$Y[, 2], Group = y_full), aes(T1, T2, colour = Group)) +
  geom_point(alpha = 0.75, size = 1.4) +
  scale_colour_manual(values = c(Normal = COL$normal, Tumor = COL$tumor), labels = c("Normal", "Tumour")) +
  labs(title = "t-SNE Projection: Tumour vs Normal", subtitle = paste0("TCGA-BRCA  |  n = ", nrow(master_data)), x = "t-SNE 1", y = "t-SNE 2", caption = paste0("t-SNE: perplexity = ", min(30, floor(nrow(X_full) / 3)), ".")) +
  theme_publication() + guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
save_figure(p_tsne, "Fig01B_tSNE", width = 7, height = 6)

# =============================================================================
# 3. PARALLEL REPEATED NESTED CV
# =============================================================================
cat("\n[3/9] Running Repeated Nested CV (", CONFIG$outer_folds * CONFIG$repeats, "folds)...\n")

n_cores <- max(1L, parallel::detectCores() - CONFIG$n_cores_leave)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)
cat("  Using", n_cores, "cores of", parallel::detectCores(), "detected\n")

multi_folds <- caret::createMultiFolds(y = master_data$Group, k = CONFIG$outer_folds, times = CONFIG$repeats)

cv_results <- foreach(
  f = seq_along(multi_folds),
  .packages = c("caret","glmnet","xgboost","ranger","dplyr","pROC","Boruta","kernlab"),
  .errorhandling = "pass"
) %dopar% {
  
  set.seed(CONFIG$seed + f)
  train_idx <- multi_folds[[f]]; test_idx <- setdiff(seq_len(nrow(master_data)), train_idx)
  
  X_tr <- master_data[train_idx, ] %>% select(-Group); y_tr <- master_data$Group[train_idx]
  X_te <- master_data[test_idx,  ] %>% select(-Group); y_te <- master_data$Group[test_idx]
  if (length(unique(y_te)) < 2) return(list(Status = "Skipped_SingleClass"))
  
  nzv <- caret::nearZeroVar(X_tr)
  if (length(nzv) > 0) { X_tr <- X_tr[, -nzv]; X_te <- X_te[, -nzv] }
  pp   <- caret::preProcess(X_tr, method = c("center", "scale"))
  X_tr <- predict(pp, X_tr); X_te <- predict(pp, X_te)
  
  obs_w   <- ifelse(y_tr == "Tumor", sum(y_tr == "Normal") / sum(y_tr == "Tumor"), 1)
  y_tr_n  <- as.numeric(y_tr == "Tumor"); fold_seed <- CONFIG$seed + f
  
  cv_el  <- tryCatch(glmnet::cv.glmnet(as.matrix(X_tr), y_tr_n, family = "binomial", alpha = 0.5, weights = obs_w, nfolds = 3), error = function(e) NULL)
  if (!is.null(cv_el)) {
    el_coef <- as.matrix(coef(cv_el, s = "lambda.min"))
    pre_fs  <- rownames(el_coef)[el_coef[, 1] != 0 & rownames(el_coef) != "(Intercept)"]
    pre_fs  <- pre_fs[pre_fs %in% colnames(X_tr)]
  } else { pre_fs <- character(0) }
  if (length(pre_fs) < 5) pre_fs <- names(sort(apply(X_tr, 2, sd), decreasing = TRUE)[1:min(30, ncol(X_tr))])
  X_tr_fs <- X_tr[, pre_fs, drop = FALSE]
  
  set.seed(fold_seed)
  lasso_g <- tryCatch({ cv_l <- glmnet::cv.glmnet(as.matrix(X_tr_fs), y_tr_n, family = "binomial", alpha = 1, weights = obs_w, nfolds = 3); rownames(coef(cv_l, "lambda.min"))[coef(cv_l,"lambda.min")[,1]!=0][-1] }, error = function(e) character(0))
  set.seed(fold_seed + 1)
  rf_g <- tryCatch({ rf_m <- ranger::ranger(x = X_tr_fs, y = y_tr, importance = "impurity", case.weights = obs_w, seed = fold_seed + 1); names(sort(rf_m$variable.importance, decreasing = TRUE)[1:min(10, ncol(X_tr_fs))]) }, error = function(e) character(0))
  set.seed(fold_seed + 2)
  xgb_g <- tryCatch({ dtr <- xgboost::xgb.DMatrix(data = as.matrix(X_tr_fs), label = y_tr_n, weight = obs_w); xm <- xgboost::xgb.train(params = list(objective = "binary:logistic", eta = 0.1, max_depth = 3, nthread = 1, seed = fold_seed + 2), data = dtr, nrounds = 50, verbose = 0); xgboost::xgb.importance(model = xm)$Feature[1:min(10, ncol(X_tr_fs))] }, error = function(e) character(0))
  set.seed(fold_seed + 3)
  bor_g <- tryCatch({ br <- Boruta::Boruta(x = X_tr_fs, y = y_tr, maxRuns = 30, doTrace = 0); Boruta::getSelectedAttributes(br, withTentative = TRUE) }, error = function(e) character(0))
  set.seed(fold_seed + 4)
  svm_g <- tryCatch({ rc <- caret::rfeControl(functions = caretFuncs, method = "cv", number = 3, allowParallel = FALSE, verbose = FALSE); sz <- unique(c(min(5, ncol(X_tr_fs)), min(10, ncol(X_tr_fs)))); predictors(caret::rfe(X_tr_fs, y_tr, sizes = sz, rfeControl = rc, method = "svmRadial")) }, error = function(e) character(0))
  
  all_fs   <- c(lasso_g, rf_g, xgb_g, bor_g, svm_g)
  gcounts  <- table(all_fs)
  con_g    <- names(gcounts[gcounts >= CONFIG$consensus_votes])
  if (length(con_g) < 2) con_g <- names(gcounts[gcounts >= 3])
  if (length(con_g) < 2) con_g <- names(gcounts[gcounts >= 2])
  if (length(con_g) < 2) con_g <- names(sort(gcounts, decreasing=TRUE)[1:min(5,length(gcounts))])
  if (length(con_g) < 1) return(list(Status = "Skipped_NoGenes"))
  
  X_tr_f <- X_tr[, con_g, drop = FALSE]; X_te_f <- X_te[, con_g, drop = FALSE]
  ctrl_in <- caret::trainControl(method = "cv", number = CONFIG$inner_folds, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = FALSE)
  
  p_xgb <- tryCatch({ di <- xgboost::xgb.DMatrix(data=as.matrix(X_tr_f), label=y_tr_n, weight=obs_w); mx <- xgboost::xgb.train(params=list(objective="binary:logistic", eval_metric="auc", max_depth=3, eta=0.1, nthread=1), data=di, nrounds=50, verbose=0); predict(mx, xgboost::xgb.DMatrix(data=as.matrix(X_te_f))) }, error = function(e) rep(0.5, nrow(X_te_f)))
  p_rf <- tryCatch({ m <- caret::train(x=X_tr_f, y=y_tr, method="ranger", trControl=ctrl_in, tuneLength=CONFIG$tune_length, metric="ROC", weights=obs_w); predict(m, X_te_f, type="prob")[,"Tumor"] }, error = function(e) rep(0.5, nrow(X_te_f)))
  p_en <- tryCatch({ m <- caret::train(x=X_tr_f, y=y_tr, method="glmnet", trControl=ctrl_in, tuneLength=CONFIG$tune_length, metric="ROC", weights=obs_w); predict(m, X_te_f, type="prob")[,"Tumor"] }, error = function(e) rep(0.5, nrow(X_te_f)))
  p_svm <- tryCatch({ m <- caret::train(x=X_tr_f, y=y_tr, method="svmRadial", trControl=ctrl_in, tuneLength=CONFIG$tune_length, metric="ROC", weights=obs_w); predict(m, X_te_f, type="prob")[,"Tumor"] }, error = function(e) rep(0.5, nrow(X_te_f)))
  
  prob_ens <- (p_xgb + p_rf + p_en + p_svm) / 4
  
  roc_f  <- pROC::roc(as.numeric(y_te=="Tumor"), prob_ens, direction="<", quiet=TRUE)
  ci_f   <- pROC::ci.auc(roc_f, conf.level=0.95)
  brier  <- mean((prob_ens - as.numeric(y_te=="Tumor"))^2)
  cm     <- caret::confusionMatrix(factor(ifelse(prob_ens>0.5,"Tumor","Normal"), levels=c("Normal","Tumor")), y_te, positive="Tumor")
  TP <- cm$table[2,2]; TN <- cm$table[1,1]; FP <- cm$table[2,1]; FN <- cm$table[1,2]
  mcc_d <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc   <- ifelse(mcc_d==0, 0, ((TP*TN)-(FP*FN))/mcc_d)
  
  list(Status="Success", Metrics=c(AUC=as.numeric(ci_f[2]), AUC_Lower=as.numeric(ci_f[1]), AUC_Upper=as.numeric(ci_f[3]), Brier=brier, Accuracy=cm$overall["Accuracy"], Sensitivity=cm$byClass["Sensitivity"], Specificity=cm$byClass["Specificity"], Precision=cm$byClass["Pos Pred Value"], F1=cm$byClass["F1"], MCC=mcc), Genes=con_g, Preds=data.frame(Observed=y_te, Predicted=prob_ens), FSLists=list(LASSO=lasso_g, RF=rf_g, XGB=xgb_g, Boruta=bor_g, SVM=svm_g))
}

parallel::stopCluster(cl)
cat("  Parallel CV complete.\n")

err_idx  <- which(sapply(cv_results, inherits, "error"))
ok_idx   <- which(sapply(cv_results, function(x) is.list(x) && identical(x$Status, "Success")))

cat("  Results: Success =", length(ok_idx), "| Errors =", length(err_idx), "\n")
valid_results <- cv_results[ok_idx]
if (length(valid_results) == 0) stop("All folds failed.")

# =============================================================================
# 4. PERFORMANCE METRICS & NESTED ROC
# =============================================================================
cat("\n[4/9] Computing performance metrics and gene stability...\n")

metrics_mat <- do.call(rbind, lapply(valid_results, `[[`, "Metrics"))
n_v <- nrow(metrics_mat)
col_mean <- colMeans(metrics_mat, na.rm = TRUE); col_sd <- apply(metrics_mat, 2, sd, na.rm = TRUE); col_se <- col_sd / sqrt(n_v)
t_crit   <- qt(0.975, df = n_v - 1)

perf_df <- data.frame(Metric=colnames(metrics_mat), Mean=round(col_mean,4), SD=round(col_sd,4), CI95_Lower=round(col_mean-t_crit*col_se,4), CI95_Upper=round(col_mean+t_crit*col_se,4))
write.csv(perf_df, file.path(dirs$tables, "Model_Performance_Summary.csv"), row.names=FALSE)

auc_mean <- col_mean["AUC"]; auc_sd <- col_sd["AUC"]
cat(sprintf("  Nested CV AUC: %.4f ± %.4f SD\n", auc_mean, auc_sd))

gene_freq    <- table(unlist(lapply(valid_results, `[[`, "Genes")))
stability_df <- data.frame(Gene=names(gene_freq), Frequency=as.numeric(gene_freq), SelectionRate=as.numeric(gene_freq)/length(valid_results)) %>% arrange(desc(SelectionRate))
write.csv(stability_df, file.path(dirs$tables, "Gene_Bootstrap_Stability.csv"), row.names=FALSE)

stable_genes       <- stability_df$Gene[stability_df$SelectionRate >= CONFIG$stability_cut]
clean_stable_genes <- gsub("_", "-", stable_genes)
n_stable           <- length(stable_genes)

cat("  Stable genes (≥", CONFIG$stability_cut, "):", n_stable, "\n")
writeLines(stable_genes, file.path(dirs$tables, "final_stable_biomarker_signature.txt"))
saveRDS(stable_genes, file.path(dirs$tables, "final_stable_biomarker_signature.rds"))

p_stab <- ggplot(stability_df[1:min(25, nrow(stability_df)),] %>% mutate(Gene=gsub("_","-",Gene), Stable=SelectionRate>=CONFIG$stability_cut), aes(x=reorder(Gene,SelectionRate), y=SelectionRate, fill=Stable)) + geom_col(colour="black", linewidth=0.35, width=0.75) + geom_hline(yintercept=CONFIG$stability_cut, linetype="dashed", colour=COL$tumor, linewidth=0.8) + coord_flip() + scale_fill_manual(values=c("FALSE"=COL$unstable, "TRUE"=COL$stable)) + labs(title="Feature Stability", x=NULL, y="Selection Rate") + theme_publication() + theme(axis.text.y=element_text(face="italic"))
save_figure(p_stab, "Fig01D_Feature_Stability", width=8, height=max(6, min(25, nrow(stability_df))*0.35+2))

oof_preds <- do.call(rbind, lapply(valid_results, `[[`, "Preds"))
brier_oof <- mean((oof_preds$Predicted - as.numeric(oof_preds$Observed=="Tumor"))^2)

calib_df <- oof_preds %>% mutate(Bin=ntile(Predicted,10)) %>% group_by(Bin) %>% summarise(Mean_Pred=mean(Predicted), Obs_Prob=mean(as.numeric(Observed=="Tumor")), .groups="drop")
p_calib <- ggplot(calib_df, aes(x=Mean_Pred, y=Obs_Prob)) + geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") + geom_line(colour=COL$tumor, linewidth=0.9) + geom_point(size=3, colour=COL$tumor) + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + labs(title="Calibration Curve", subtitle=paste0("Brier = ", round(brier_oof,4)), x="Mean Predicted Prob", y="Observed Prob") + theme_publication() + theme(aspect.ratio=1)
save_figure(p_calib, "Fig01C_Calibration", width=6, height=6)

roc_cv <- pROC::roc(as.numeric(oof_preds$Observed=="Tumor"), oof_preds$Predicted, direction="<", quiet=TRUE)
roc_coords <- pROC::coords(roc_cv, "all", ret=c("specificity","sensitivity"), transpose=FALSE) %>% as.data.frame() %>% mutate(fpr=1-specificity)
set.seed(CONFIG$seed)
ci_se <- pROC::ci.se(roc_cv, specificities=seq(0,1,by=0.04), conf.level=0.95, method="bootstrap", boot.n=500, quiet=TRUE)
ci_df <- data.frame(fpr=1-seq(0,1,by=0.04), se_lower=as.numeric(ci_se[,1]), se_upper=as.numeric(ci_se[,3]))

p_roc_cv <- ggplot() + geom_ribbon(data=ci_df, aes(x=fpr, ymin=se_lower, ymax=se_upper), fill=COL$tumor, alpha=0.15) + geom_line(data=roc_coords, aes(x=fpr, y=sensitivity), colour=COL$tumor, linewidth=1.1) + geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)) + labs(title="Nested CV Ensemble ROC", subtitle=sprintf("AUC = %.4f ± %.4f SD (%d folds)", auc_mean, auc_sd, n_v), x="1 - Specificity", y="Sensitivity") + theme_publication() + theme(aspect.ratio=1)
save_figure(p_roc_cv, "Fig01E_Nested_CV_ROC", width=8, height=7.5)

# =============================================================================
# 5. UPSET PLOT & FEATURE SELECTION SUMMARY
# =============================================================================
cat("\n[5/9] Generating UpSet plot...\n")

X_up <- master_data %>% select(-Group); y_up <- master_data$Group
nzv_up <- caret::nearZeroVar(X_up)
if (length(nzv_up) > 0) X_up <- X_up[, -nzv_up]

set.seed(CONFIG$seed)
up_data <- caret::upSample(x=X_up, y=y_up, yname="Group")
Xb <- up_data %>% select(-Group); ybn <- as.numeric(up_data$Group == "Tumor")

set.seed(CONFIG$seed+10)
u_lasso <- tryCatch({ cv_ul <- glmnet::cv.glmnet(as.matrix(Xb), ybn, family="binomial", alpha=1, nfolds=5); rownames(coef(cv_ul,"lambda.min"))[coef(cv_ul,"lambda.min")[,1]!=0 & rownames(coef(cv_ul,"lambda.min"))!="(Intercept)"] }, error=function(e) character(0))
set.seed(CONFIG$seed+11)
u_rf <- tryCatch({ rfm <- randomForest::randomForest(x=Xb, y=up_data$Group, ntree=300, importance=TRUE); rownames(randomForest::importance(rfm))[order(randomForest::importance(rfm)[,"MeanDecreaseGini"], decreasing=TRUE)][1:min(10,ncol(Xb))] }, error=function(e) character(0))
set.seed(CONFIG$seed+12)
u_xgb <- tryCatch({ du <- xgboost::xgb.DMatrix(data=as.matrix(Xb), label=ybn); xu <- xgboost::xgb.train(params=list(objective="binary:logistic",eta=0.1,max_depth=3,nthread=1), data=du, nrounds=50, verbose=0); xgboost::xgb.importance(model=xu)$Feature[1:min(10,ncol(Xb))] }, error=function(e) character(0))
set.seed(CONFIG$seed+13)
u_bor <- tryCatch({ Boruta::getSelectedAttributes(Boruta::Boruta(x=Xb, y=up_data$Group, maxRuns=30, doTrace=0), withTentative=TRUE) }, error=function(e) character(0))
set.seed(CONFIG$seed+14)
u_svm <- tryCatch({ rc <- caret::rfeControl(functions=caretFuncs, method="cv", number=3, verbose=FALSE); predictors(caret::rfe(Xb, up_data$Group, sizes=c(5,10,15), rfeControl=rc, method="svmRadial")) }, error=function(e) character(0))

upset_list <- Filter(function(x) length(x)>0, list(LASSO=u_lasso, RF=u_rf, XGB=u_xgb, Boruta=u_bor, SVM=u_svm))
saveRDS(upset_list, file.path(dirs$tables, "01_upset_gene_lists.rds"))

fs_txt_path <- file.path(dirs$tables, "Feature_Selection_Summary.txt")
sink(fs_txt_path)
cat("=============================================================\n 5-METHOD FEATURE SELECTION SUMMARY\n=============================================================\n\n")
for (m in names(upset_list)) {
  genes_m <- upset_list[[m]]
  cat(sprintf("Method : %s\nCount  : %d genes selected\nGenes  : %s\n-------------------------------------------------------------\n", m, length(genes_m), paste(gsub("_", "-", genes_m), collapse=", ")))
}
sink()
cat("  [SAVED] Feature_Selection_Summary.txt\n")

nm <- length(upset_list)
if (nm >= 2) {
  for (ext in c("png","tiff","pdf")) {
    if (ext == "pdf") { grDevices::cairo_pdf(file.path(dirs$figures, paste0("Fig02_UpSet.", ext)), width=12, height=max(5, nm*1.2+2))
    } else { grDevices::png(file.path(dirs$figures, paste0("Fig02_UpSet.", ext)), width=12, height=max(5, nm*1.2+2), units="in", res=ifelse(ext=="tiff",600,300)) }
    p_up <- UpSetR::upset(UpSetR::fromList(upset_list), order.by="freq", main.bar.color=COL$tumor, sets.bar.color=COL$normal, text.scale=c(1.8, 1.4, 1.3, 1.3, 1.5, 1.3), point.size=2.5, line.size=0.8)
    print(p_up)
    grDevices::dev.off()
  }
  cat("  [SAVED] Fig02_UpSet (PNG + TIFF + PDF)\n")
}

# =============================================================================
# 6. FINAL MODEL — BAYESIAN OPTIMISATION
# =============================================================================
cat("\n[6/9] Training final Bayesian-optimised XGBoost...\n")

if (n_stable < 1) stop("No stable genes found — check stability_cut parameter.")

X_fin      <- master_data %>% select(all_of(stable_genes))
preProc_fin <- caret::preProcess(X_fin, method = c("center", "scale"))
X_fin_sc   <- predict(preProc_fin, X_fin)
saveRDS(preProc_fin, file.path(dirs$models, "final_preProc_scaler.rds"))

y_fin_n <- as.numeric(master_data$Group == "Tumor")
pos_w   <- sum(y_fin_n == 0) / sum(y_fin_n == 1)
dtr_fin <- xgboost::xgb.DMatrix(data = as.matrix(X_fin_sc), label = y_fin_n)

scoring_fn <- function(eta, max_depth, subsample, colsample_bytree) {
  cv_r <- xgboost::xgb.cv(params=list(objective="binary:logistic", eval_metric="auc", eta=eta, max_depth=max_depth, subsample=subsample, colsample_bytree=colsample_bytree, scale_pos_weight=pos_w, nthread=1), data=dtr_fin, nfold=5, nrounds=200, early_stopping_rounds=15, verbose=0)
  best_nr <- if (is.null(cv_r$best_iteration)||is.na(cv_r$best_iteration[1])||cv_r$best_iteration[1]<1L) as.integer(nrow(cv_r$evaluation_log)) else as.integer(cv_r$best_iteration[1])
  list(Score=max(cv_r$evaluation_log$test_auc_mean), nrounds=best_nr)
}

opt_res <- tryCatch(ParBayesianOptimization::bayesOpt(FUN=scoring_fn, bounds=list(eta=c(0.01,0.3), max_depth=c(2L,6L), subsample=c(0.5,1.0), colsample_bytree=c(0.5,1.0)), initPoints=CONFIG$bayes_init, iters.n=CONFIG$bayes_iter, verbose=0), error=function(e) NULL)

if (!is.null(opt_res)) {
  bp <- ParBayesianOptimization::getBestPars(opt_res)
  final_params <- list(objective="binary:logistic", eval_metric="auc", eta=bp$eta, max_depth=bp$max_depth, subsample=bp$subsample, colsample_bytree=bp$colsample_bytree, scale_pos_weight=pos_w)
  final_rounds <- as.integer(opt_res$scoreSummary$nrounds[which.max(opt_res$scoreSummary$Score)][1])
  if (is.na(final_rounds) || final_rounds < 1L) final_rounds <- 150L
} else {
  final_params <- list(objective="binary:logistic", eval_metric="auc", eta=0.1, max_depth=3, scale_pos_weight=pos_w); final_rounds <- 150L
}

final_model <- xgboost::xgb.train(params=final_params, data=dtr_fin, nrounds=final_rounds, verbose=0)
saveRDS(final_model, file.path(dirs$models, "final_xgboost_model_Bayesian.rds"))

# Diagnostic ROC (5-fold CV)
cat("  Computing XGBoost Diagnostic ROC...\n")
set.seed(CONFIG$seed)
diag_folds <- caret::createFolds(master_data$Group, k=5, list=TRUE, returnTrain=FALSE)
diag_preds <- data.frame(Observed=integer(), Predicted=numeric())

for (dfi in seq_along(diag_folds)) {
  ti <- diag_folds[[dfi]]; ri <- setdiff(seq_len(nrow(master_data)), ti)
  if (length(unique(y_fin_n[ti])) < 2) next
  dm_tr <- xgboost::xgb.DMatrix(data=as.matrix(X_fin_sc[ri,]), label=y_fin_n[ri])
  dm_te <- xgboost::xgb.DMatrix(data=as.matrix(X_fin_sc[ti,]), label=y_fin_n[ti])
  set.seed(CONFIG$seed + dfi)
  dm <- xgboost::xgb.train(params=final_params, data=dm_tr, nrounds=final_rounds, verbose=0)
  diag_preds <- rbind(diag_preds, data.frame(Observed=y_fin_n[ti], Predicted=predict(dm, dm_te)))
}

roc_diag <- pROC::roc(diag_preds$Observed, diag_preds$Predicted, direction="<", quiet=TRUE)
ci_diag  <- pROC::ci.auc(roc_diag, conf.level=0.95)
diag_coords <- pROC::coords(roc_diag, "all", ret=c("specificity","sensitivity"), transpose=FALSE) %>% as.data.frame() %>% mutate(fpr=1-specificity)
youden <- pROC::coords(roc_diag, "best", best.method="youden", ret=c("specificity","sensitivity","threshold"), transpose=FALSE)[1,]

set.seed(CONFIG$seed)
diag_ci <- pROC::ci.se(roc_diag, specificities=seq(0,1,by=0.04), conf.level=0.95, method="bootstrap", boot.n=500, quiet=TRUE)
diag_ci_df <- data.frame(fpr=1-seq(0,1,by=0.04), lower=as.numeric(diag_ci[,1]), upper=as.numeric(diag_ci[,3]))

p_diag_roc <- ggplot() + geom_ribbon(data=diag_ci_df, aes(x=fpr, ymin=lower, ymax=upper), fill=COL$tumor, alpha=0.15) + geom_line(data=diag_coords, aes(x=fpr, y=sensitivity), colour=COL$tumor, linewidth=1.2) + geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") + geom_point(aes(x=1-youden$specificity, y=youden$sensitivity), colour=COL$tumor, size=3.5, shape=18) +  annotate("text", x=1-youden$specificity+0.03, y=youden$sensitivity-0.05, label=sprintf("Optimal\n(Se=%.2f, Sp=%.2f)", youden$sensitivity, youden$specificity), size=3.2, colour=COL$tumor, hjust=0) + scale_x_continuous(limits=c(0,1), expand=c(0.01,0.01)) + scale_y_continuous(limits=c(0,1), expand=c(0.01,0.01)) + labs(title="XGBoost Diagnostic ROC \u2014 Stable Gene Signature", subtitle=sprintf("AUC = %.4f  (95%% CI: %.4f \u2013 %.4f)", as.numeric(ci_diag[2]), as.numeric(ci_diag[1]), as.numeric(ci_diag[3])), x="1 - Specificity", y="Sensitivity") + theme_publication() + theme(aspect.ratio=1)
save_figure(p_diag_roc, "Fig03_XGB_Diagnostic_ROC", width=9, height=8.5)

# =============================================================================
# 7. BIOLOGICAL VALIDATION (DE ANALYSIS)
# =============================================================================
cat("\n[7/9] Biological validation (Wilcoxon + FDR)...\n")

de_res <- do.call(rbind, lapply(clean_stable_genes, function(g) {
  gc    <- gsub("-","_",g)
  n_exp <- master_data[[gc]][master_data$Group=="Normal"]
  t_exp <- master_data[[gc]][master_data$Group=="Tumor"]
  pv    <- tryCatch(wilcox.test(t_exp, n_exp)$p.value, error=function(e) NA)
  
  fc <- (mean(t_exp, na.rm=TRUE) + 0.001) / (mean(n_exp, na.rm=TRUE) + 0.001)
  data.frame(Gene=g, Mean_Normal=round(mean(n_exp,na.rm=TRUE),5), Mean_Tumor=round(mean(t_exp,na.rm=TRUE),5), Log2FC=round(log2(fc),5), P_Value=pv, stringsAsFactors=FALSE)
}))
de_res$FDR <- p.adjust(de_res$P_Value, method="BH")

# Best Practice: Biological Direction dictates the signature direction
de_res$DE_Direction <- ifelse(de_res$Log2FC > 0 & de_res$FDR < 0.05, "Up in Tumor", ifelse(de_res$Log2FC < 0 & de_res$FDR < 0.05, "Down in Tumor", "No Sig. Change"))
de_res$Dir_Symbol <- ifelse(de_res$DE_Direction == "Up in Tumor", "\u2191 Tumor", ifelse(de_res$DE_Direction == "Down in Tumor", "\u2193 Tumor", "-"))
de_res <- de_res %>% arrange(FDR)
write.csv(de_res, file.path(dirs$tables, "Biological_Validation_DE_Results.csv"), row.names=FALSE)
cat("  [SAVED] Biological_Validation_DE_Results.csv\n")

# ── 7B. Directional Score diagnostic ROC (Fig03B & 03C) ───────────────────────
cat("\n  Computing Clinical Directional Score diagnostic ROC...\n")

# Extract UP and DOWN genes directly from the biological truth
dir_up_genes   <- de_res$Gene[de_res$DE_Direction == "Up in Tumor"]
dir_down_genes <- de_res$Gene[de_res$DE_Direction == "Down in Tumor"]

# Score = Σ(Up genes) − Σ(Down genes) in raw log2(TPM+1) units
# We use X_fin (unscaled) to preserve original biological units
X_fin_hyphen <- X_fin
colnames(X_fin_hyphen) <- clean_stable_genes

dir_scores <- rep(0, nrow(X_fin_hyphen))
if(length(dir_up_genes) > 0) {
  dir_scores <- dir_scores + rowSums(X_fin_hyphen[, dir_up_genes, drop=FALSE])
}
if(length(dir_down_genes) > 0) {
  dir_scores <- dir_scores - rowSums(X_fin_hyphen[, dir_down_genes, drop=FALSE])
}

roc_dir  <- pROC::roc(y_fin_n, dir_scores, levels=c(0,1), direction="<", quiet=TRUE)
ci_dir   <- pROC::ci.auc(roc_dir, conf.level=0.95)

# Full XGBoost eval for DeLong test comparison
xgb_probs_full <- predict(final_model, xgboost::xgb.DMatrix(data=as.matrix(X_fin_sc)))
roc_xgb_full   <- pROC::roc(y_fin_n, xgb_probs_full, levels=c(0,1), direction="<", quiet=TRUE)
dl_test        <- pROC::roc.test(roc_dir, roc_xgb_full, method="delong")

dir_coords <- pROC::coords(roc_dir, "all", ret=c("specificity","sensitivity"), transpose=FALSE) %>% as.data.frame() %>% mutate(fpr=1-specificity)
youden_dir <- pROC::coords(roc_dir, "best", best.method="youden", ret=c("specificity","sensitivity","threshold"), transpose=FALSE)[1,]

set.seed(CONFIG$seed)
dir_ci_se <- pROC::ci.se(roc_dir, specificities=seq(0,1,by=0.04), conf.level=0.95, method="bootstrap", boot.n=500, quiet=TRUE)
dir_ci_df <- data.frame(fpr=1-seq(0,1,by=0.04), lower=as.numeric(dir_ci_se[,1]), upper=as.numeric(dir_ci_se[,3]))

p_dir_roc <- ggplot() +
  geom_ribbon(data=dir_ci_df, aes(x=fpr, ymin=lower, ymax=upper), fill="#1B7837", alpha=0.15) +
  geom_line(data=dir_coords, aes(x=fpr, y=sensitivity), colour="#1B7837", linewidth=1.2) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50", linewidth=0.6) +
  geom_point(aes(x=1-youden_dir$specificity, y=youden_dir$sensitivity), colour="#1B7837", size=3.5, shape=18) +
  annotate("text", x=1-youden_dir$specificity+0.03, y=youden_dir$sensitivity-0.05, label=sprintf("Optimal\n(Se=%.2f, Sp=%.2f)", youden_dir$sensitivity, youden_dir$specificity), size=3.2, colour="#1B7837", hjust=0) +
  scale_x_continuous(limits=c(0,1), expand=c(0.01,0.01)) + scale_y_continuous(limits=c(0,1), expand=c(0.01,0.01)) +
  labs(title="Diagnostic ROC \u2014 Clinical Directional Score", subtitle=sprintf("AUC = %.4f  (95%% CI: %.4f \u2013 %.4f)  |  Direction vs XGBoost DeLong: p%s", as.numeric(ci_dir[2]), as.numeric(ci_dir[1]), as.numeric(ci_dir[3]), ifelse(dl_test$p.value<0.001,"<0.001", sprintf("=%.4f",dl_test$p.value))), x="1 \u2212 Specificity  (False Positive Rate)", y="Sensitivity  (True Positive Rate)", caption="Score = \u03a3(Up) - \u03a3(Down). Evaluated on full cohort.") + theme_publication() + theme(aspect.ratio=1)
save_figure(p_dir_roc, "Fig03B_Directional_Diagnostic_ROC", width=9, height=8.5)

# Fig03C (Overlay)
ci_xgb_full <- pROC::ci.auc(roc_xgb_full, conf.level=0.95)
xgb_coords_full <- pROC::coords(roc_xgb_full, "all", ret=c("specificity","sensitivity"), transpose=FALSE) %>% as.data.frame() %>% mutate(fpr=1-specificity, Method="XGBoost")
dir_coords$Method <- "Directional Score"
roc_overlay <- rbind(xgb_coords_full, dir_coords)

p_roc_overlay <- ggplot(roc_overlay, aes(x=fpr, y=sensitivity, colour=Method, linetype=Method)) +
  geom_line(linewidth=1.1) + geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50", linewidth=0.5) +
  scale_colour_manual(values=c("XGBoost"=COL$tumor, "Directional Score"="#1B7837")) + scale_linetype_manual(values=c("XGBoost"="solid", "Directional Score"="solid")) +
  annotate("text", x=0.55, y=0.22, label=sprintf("XGBoost AUC = %.4f \nDirectional  AUC = %.4f\nDeLong p %s", as.numeric(ci_xgb_full[2]), as.numeric(ci_dir[2]), ifelse(dl_test$p.value<0.001,"<0.001",sprintf("=%.4f",dl_test$p.value))), size=3.5, hjust=0, fontface="bold", colour="grey20") +
  scale_x_continuous(limits=c(0,1), expand=c(0.01,0.01)) + scale_y_continuous(limits=c(0,1), expand=c(0.01,0.01)) +
  labs(title="XGBoost vs Directional Score", x="1 \u2212 Specificity  (False Positive Rate)", y="Sensitivity  (True Positive Rate)") + theme_publication() + theme(aspect.ratio=1, legend.position=c(0.97,0.05), legend.justification=c(1,0), legend.background=element_rect(fill="white",colour="grey80", linewidth=0.4))
save_figure(p_roc_overlay, "Fig03C_XGB_vs_Dir_ROC_Overlay", width=9, height=8.5)

# =============================================================================
# 8. SHAP SUITE — Publication Quality (DEG-Directed)
# =============================================================================
cat("\n[8/9] Generating publication-quality SHAP suite...\n")

X_fin_sc_plot           <- as.matrix(X_fin_sc)
colnames(X_fin_sc_plot) <- clean_stable_genes

shap_result <- SHAPforxgboost::shap.values(xgb_model=final_model, X_train=X_fin_sc_plot)
shap_matrix   <- shap_result$shap_score   
shap_abs_mean <- colMeans(abs(shap_matrix))

# Gene order: descending mean |SHAP| (Machine Learning Importance)
gene_order <- names(sort(shap_abs_mean, decreasing=TRUE))

shap_long <- SHAPforxgboost::shap.prep(shap_contrib=shap_matrix, X_train=X_fin_sc_plot) %>% mutate(variable=factor(variable, levels=rev(gene_order)))

# ── 8A. SHAP Beeswarm ────────────────────────────────────────────────────────
p_bee <- ggplot(shap_long, aes(x=value, y=variable, colour=rfvalue)) +
  ggbeeswarm::geom_quasirandom(groupOnX=FALSE, bandwidth=0.3, size=0.85, alpha=0.65, varwidth=FALSE) +
  geom_vline(xintercept=0, colour="black", linewidth=0.5) +
  scale_colour_gradientn(colours=c(COL$low_expr, COL$mid_expr, COL$high_expr), values=scales::rescale(c(0, 0.5, 1)), name="Feature value", breaks=c(0, 0.5, 1), labels=c("Low", "Mid", "High"), guide=guide_colorbar(title.position="top", title.hjust=0.5, barwidth=unit(6,"cm"))) +
  scale_y_discrete(labels=function(x) setNames(lapply(x, function(g) bquote(italic(.(g)))), x)) +
  labs(title="SHAP Summary: Gene Importance to Prediction", subtitle="Points coloured by relative expression. Ordered by Mean |SHAP|.", x="SHAP value", y=NULL) +
  theme_publication(base_size=13) + theme(legend.position="bottom", axis.text.y=element_text(face="italic", size=12))
save_figure(p_bee, "Fig04A_SHAP_Beeswarm", width=10, height=max(6, n_stable*0.65+2.5))

# ── 8B. SHAP Importance Bar (COLOURED BY DEG DIRECTION) ──────────────────────
aligned_deg_dir <- de_res$DE_Direction[match(gene_order, de_res$Gene)]
shap_imp_df <- data.frame(Gene=factor(gene_order, levels=gene_order), MeanAbs=shap_abs_mean[gene_order], Direction=aligned_deg_dir)

dir_col_map <- c("Up in Tumor"=COL$tumor, "Down in Tumor"=COL$normal, "No Sig. Change"="grey70")

p_bar <- ggplot(shap_imp_df, aes(x=reorder(Gene,MeanAbs), y=MeanAbs, fill=Direction)) +
  geom_col(colour="black", linewidth=0.35, width=0.72) + geom_text(aes(label=sprintf("%.4f", MeanAbs)), hjust=-0.1, size=3.5, colour="grey20") +
  coord_flip() + scale_fill_manual(values=dir_col_map, name="Biological Direction") +
  scale_y_continuous(expand=expansion(mult=c(0,0.20)), breaks=pretty_breaks(n=5)) +
  labs(title="Feature Importance: Mean |SHAP Value|", subtitle="Importance determined by XGBoost; Direction determined by DE Biology.", x=NULL, y="Mean |SHAP value|") +
  theme_publication(base_size=13) + theme(axis.text.y=element_text(face="italic", size=12), legend.position="right")
save_figure(p_bar, "Fig04B_SHAP_Importance_Bar", width=9, height=max(5, n_stable*0.5+2))

# ── 8C. SHAP Dependence Plots ────────────────────────────────────────────────
for (g in gene_order) {
  # Safe extraction using [[g]] from data.table
  dep_df <- data.frame(Expression = X_fin_sc_plot[, g], SHAP = shap_matrix[[g]], Group = master_data$Group)
  
  # Determine color from Biology (DEG)
  d_dir <- de_res$DE_Direction[de_res$Gene == g]
  tcol  <- ifelse(length(d_dir) > 0 && d_dir == "Up in Tumor", COL$tumor, ifelse(length(d_dir) > 0 && d_dir == "Down in Tumor", COL$normal, "grey50"))
  
  p_dep <- ggplot(dep_df, aes(x=Expression, y=SHAP)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    geom_point(aes(colour=Group), size=1.2, alpha=0.55) +
    geom_smooth(method="loess", span=0.75, se=TRUE, colour=tcol, fill=tcol, alpha=0.12, linewidth=1.0) +
    scale_colour_manual(values=c(Normal=COL$normal, Tumor=COL$tumor)) +
    labs(title=bquote(bold("SHAP Dependence: ") ~ italic(.(g))), subtitle=paste0("Biological Direction: ", d_dir), x=bquote("Expression of" ~ italic(.(g)) ~ "(scaled)"), y="SHAP value") + theme_publication()
  save_figure(p_dep, paste0("Fig04C_SHAP_Dependence_", gsub("-","_",g)), width=8, height=6)
}

# =============================================================================
# FINAL DIRECTIONAL SIGNATURE TABLE
# =============================================================================
tertiles <- quantile(shap_abs_mean, probs=c(0.333, 0.667))

# Using explicit if() return() block to prevent ifelse from stripping the names off the vector
shap_cat <- sapply(shap_abs_mean, function(x) {
  if (x >= tertiles[2]) return("High")
  if (x >= tertiles[1]) return("Medium")
  return("Low")
})

final_sig_df <- data.frame(
  Gene            = gene_order,
  Direction       = de_res$Dir_Symbol[match(gene_order, de_res$Gene)],
  log2FC          = de_res$Log2FC[match(gene_order, de_res$Gene)],
  SHAP_importance = shap_cat[gene_order]
)

write.csv(final_sig_df, file.path(dirs$tables, "Final_Directional_Signature.csv"), row.names=FALSE)
saveRDS(final_sig_df, file.path(dirs$tables, "Final_Directional_Signature.rds"))
cat("  [SAVED] Final_Directional_Signature.csv\n  [SAVED] Final_Directional_Signature.rds\n")
cat("\n  ── Final Validated Signature ──\n"); print(final_sig_df, row.names=FALSE)

# =============================================================================
# 9. SESSION INFO & WRAP UP
# =============================================================================
cat("\n[9/9] Saving session info and completing...\n")

sink(file.path(dirs$logs, "01_session_info.txt")); print(sessionInfo()); sink()
sink(log_path, append=TRUE, split=TRUE)

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 01_training_pipeline.R  COMPLETE\n")
cat(" End   :", as.character(t_end), "\n")
cat(" Runtime:", round(difftime(t_end, t_start, units="mins"), 2), "minutes\n")
cat("\n SAVED ARTEFACTS:\n")
cat("  Models  : final_xgboost_model_Bayesian.rds, final_preProc_scaler.rds\n")
cat("  Tables  : Final_Directional_Signature.rds/.csv\n")
cat("            Gene_Bootstrap_Stability.csv, Model_Performance_Summary.csv\n")
cat("            Feature_Selection_Summary.txt, Biological_Validation_DE_Results.csv\n")
cat("            01_master_data_snapshot.rds, 01_upset_gene_lists.rds\n")
cat("  Figures : Fig01A_UMAP, Fig01B_tSNE, Fig01C_Calibration, Fig01D_Feature_Stability,\n")
cat("            Fig01E_Nested_CV_ROC, Fig02_UpSet, Fig03_XGB_Diagnostic_ROC,\n")
cat("            Fig03B_Directional_Diagnostic_ROC, Fig03C_XGB_vs_Dir_ROC_Overlay,\n")
cat("            Fig04A_SHAP_Beeswarm, Fig04B_SHAP_Importance_Bar, Fig04C_SHAP_Dependence_[GENE]\n")
cat("=============================================================\n")
sink()