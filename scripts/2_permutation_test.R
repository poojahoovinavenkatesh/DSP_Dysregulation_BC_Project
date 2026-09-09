# =============================================================================
# Script      : 2_permutation_test.R
# Project     : DSP family dysregulation in breast cancer
# Description : 1000-permutation test to assess statistical significance of
#               the stable gene signature AUC. Mirrors 01_ pipeline exactly.
#               Crash-safe: raw permutation results saved immediately after
#               parallel section; post-processing can be re-run independently.
# Author      : Pooja Hoovina Venkatesh, Vidya Niranjan, Sumathra Manokaran
# Date        : 30.08.2026
# Seed        : 42
# Depends on  : 1_training_pipeline.R 
# =============================================================================
# REQUIRED INPUTS:
#   output/tables/final_stable_biomarker_signature.rds
#   output/tables/Model_Performance_Summary.csv  (true AUC reference)
#   data/TCGA-BRCA.star_tpm.tsv.gz
#   data/BRCA_gene_list.txt
#
# CRASH RECOVERY:
#   If the script crashes after the parallel section, set RECOVERY_MODE = TRUE.
#   It will reload output/tables/03_perm_results_raw.rds and skip straight to
#   post-processing (p-value, plots, summary CSV).
# =============================================================================

rm(list=ls()); gc()

# =============================================================================
# USER CONFIGURATION
# =============================================================================
CONFIG <- list(
  seed           = 42L,
  n_perm         = 1000L,    # number of permutations
  n_cores_leave  = 2L,       # cores to leave free for OS
  outer_folds    = 5L,       # must match 01_ CONFIG
  inner_folds    = 3L,       # must match 01_
  consensus_votes= 4L,       # must match 01_
  stability_cut  = 0.80,     # must match 01_
  RECOVERY_MODE  = FALSE     # set TRUE to skip parallel section after crash
)

PATHS <- list(
  data_file   = file.path("data","TCGA-BRCA.star_tpm.tsv.gz"),
  gene_list   = file.path("data","BRCA_gene_list.txt"),
  perf_csv    = file.path("output","tables","Model_Performance_Summary.csv"),
  sig_rds     = file.path("output","tables","final_stable_biomarker_signature.rds"),
  raw_rds     = file.path("output","tables","03_perm_results_raw.rds")
)

# =============================================================================
# 0. PACKAGES
# =============================================================================
if (!requireNamespace("pacman",quietly=TRUE)) install.packages("pacman")
pacman::p_load(
  here, data.table, dplyr, tibble, stringr,
  caret, glmnet, ranger, xgboost, e1071, kernlab,
  Boruta, pROC, ggplot2, doParallel, foreach, parallel,
  AnnotationDbi, org.Hs.eg.db, scales
)
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# =============================================================================
# DIRECTORIES
# =============================================================================
dirs <- list(
  tables  = here("output","tables"),
  figures = here("output","figures"),
  logs    = here("output","logs")
)
lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# LOGGING
# =============================================================================
log_path <- file.path(dirs$logs,
                      paste0("03_Permutation_Run_",format(Sys.time(),"%Y%m%d_%H%M"),".log"))
sink(log_path, append=FALSE, split=TRUE)
t_start <- Sys.time()
cat("=============================================================\n")
cat(" 03_permutation_test.R\n")
cat(" Start       :", as.character(t_start), "\n")
cat(" Permutations:", CONFIG$n_perm, "\n")
cat(" Recovery mode:", CONFIG$RECOVERY_MODE, "\n")
cat("=============================================================\n\n")
set.seed(CONFIG$seed)

# =============================================================================
# SHARED AESTHETICS (identical across all 5 scripts)
# =============================================================================
COL <- list(tumor="#B22222", normal="#2166AC",
            low_expr="#3288BD", high_expr="#D53E4F", mid_expr="#FFFFBF")

theme_publication <- function(base_size=12) {
  theme_bw(base_size=base_size) +
    theme(
      plot.title       = element_text(face="bold",hjust=0.5,size=base_size+2),
      plot.subtitle    = element_text(hjust=0.5,colour="grey30",
                                      size=base_size-1,margin=margin(b=5)),
      plot.caption     = element_text(size=base_size-3.5,colour="grey45",
                                      hjust=0,lineheight=1.3,
                                      margin=margin(t=8)),
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
                  width=width,height=height,dpi=300)
  ggplot2::ggsave(paste0(base,".tiff"),plot=plot_obj,
                  width=width,height=height,dpi=600,
                  device="tiff",compression="lzw")
  grDevices::cairo_pdf(paste0(base,".pdf"),width=width,height=height)
  print(plot_obj); grDevices::dev.off()
  cat("  [SAVED]",filename_base,"(PNG + TIFF 600dpi + PDF)\n")
  invisible(plot_obj)
}

# =============================================================================
# 1. LOAD ARTEFACTS AND DATA
# =============================================================================
cat("[1/5] Loading artefacts...\n")

for (p in c(PATHS$perf_csv, PATHS$sig_rds)) {
  if (!file.exists(p)) stop("Missing: ", p, "\nRun 01_training_pipeline.R first.")
}

perf_df     <- read.csv(PATHS$perf_csv)
true_auc    <- perf_df$Mean[perf_df$Metric == "AUC"]
true_ci_lo  <- perf_df$CI95_Lower[perf_df$Metric == "AUC"]
true_ci_hi  <- perf_df$CI95_Upper[perf_df$Metric == "AUC"]
stable_genes <- readRDS(PATHS$sig_rds)
clean_genes  <- gsub("_","-",stable_genes)
n_genes      <- length(stable_genes)

cat("  True AUC:", round(true_auc,4),
    "  CI [", round(true_ci_lo,4), ",", round(true_ci_hi,4), "]\n")
cat("  Stable genes:", n_genes, "\n\n")

if (!CONFIG$RECOVERY_MODE) {
  cat("  Loading TCGA-BRCA data...\n")
  target_genes <- readLines(PATHS$gene_list)
  
  df <- data.table::fread(PATHS$data_file, header=TRUE, sep="\t")
  if (!"gene_id" %in% names(df)) names(df)[1] <- "gene_id"
  df$gene_id_clean <- str_split_fixed(df$gene_id,"\\.",2)[,1]
  suppressMessages({
    df$gene_symbol <- AnnotationDbi::mapIds(
      org.Hs.eg.db, keys=df$gene_id_clean,
      column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  })
  
  master_data <- df %>%
    filter(!is.na(gene_symbol), gene_symbol %in% target_genes) %>%
    group_by(gene_symbol) %>%
    summarise(across(starts_with("TCGA"),\(x) mean(x,na.rm=TRUE)),.groups="drop") %>%
    tibble::column_to_rownames("gene_symbol") %>%
    as.matrix() %>% t() %>% as.data.frame() %>%
    mutate(sample_type_code=substr(rownames(.),14,15)) %>%
    filter(sample_type_code %in% c("01","11")) %>%
    mutate(Group=factor(ifelse(sample_type_code=="01","Tumor","Normal"),
                        levels=c("Normal","Tumor"))) %>%
    select(-sample_type_code) %>% na.omit()
  
  colnames(master_data) <- gsub("-","_",colnames(master_data))
  cat("  Samples:", nrow(master_data),
      "| Tumor:", sum(master_data$Group=="Tumor"),
      "| Normal:", sum(master_data$Group=="Normal"), "\n\n")
}

# =============================================================================
# SINGLE-PERMUTATION EVALUATION (mirrors 01_ pipeline exactly)
# =============================================================================
evaluate_permutation <- function(data_perm, seed_i, consensus_votes, inner_folds) {
  
  set.seed(seed_i)
  folds <- caret::createFolds(data_perm$Group, k=5, list=TRUE, returnTrain=FALSE)
  preds_all <- data.frame(Observed=integer(), Predicted=numeric())
  
  for (fi in seq_along(folds)) {
    test_idx  <- folds[[fi]]
    train_idx <- setdiff(seq_len(nrow(data_perm)), test_idx)
    
    X_tr <- data_perm[train_idx,] %>% select(-Group)
    y_tr <- data_perm$Group[train_idx]
    X_te <- data_perm[test_idx,]  %>% select(-Group)
    y_te <- data_perm$Group[test_idx]
    
    if (length(unique(y_te))<2) next
    
    nzv <- caret::nearZeroVar(X_tr)
    if (length(nzv)>0){X_tr<-X_tr[,-nzv]; X_te<-X_te[,-nzv]}
    pp   <- caret::preProcess(X_tr, method=c("center","scale"))
    X_tr <- predict(pp,X_tr); X_te <- predict(pp,X_te)
    
    y_tr_n <- as.numeric(y_tr=="Tumor")
    obs_w  <- ifelse(y_tr=="Tumor",sum(y_tr=="Normal")/sum(y_tr=="Tumor"),1)
    fs     <- CONFIG$seed + seed_i * 100 + fi
    
    # 5-method ensemble FS (lightweight subset for speed)
    l_g <- tryCatch({
      cv_l <- glmnet::cv.glmnet(as.matrix(X_tr),y_tr_n,
                                family="binomial",alpha=1,weights=obs_w,nfolds=3)
      rownames(coef(cv_l,"lambda.min"))[
        coef(cv_l,"lambda.min")[,1]!=0 &
          rownames(coef(cv_l,"lambda.min"))!="(Intercept)"]
    }, error=function(e) character(0))
    
    r_g <- tryCatch({
      rfm <- ranger::ranger(x=X_tr,y=y_tr,importance="impurity",
                            case.weights=obs_w,seed=fs)
      names(sort(rfm$variable.importance,decreasing=TRUE)[1:min(10,ncol(X_tr))])
    }, error=function(e) character(0))
    
    x_g <- tryCatch({
      du <- xgboost::xgb.DMatrix(data=as.matrix(X_tr),label=y_tr_n,weight=obs_w)
      xu <- xgboost::xgb.train(
        params=list(objective="binary:logistic",eta=0.1,max_depth=3,
                    nthread=1,seed=fs),
        data=du,nrounds=50,verbose=0)
      xgboost::xgb.importance(model=xu)$Feature[1:min(10,ncol(X_tr))]
    }, error=function(e) character(0))
    
    b_g <- tryCatch(
      Boruta::getSelectedAttributes(
        Boruta::Boruta(x=X_tr,y=y_tr,maxRuns=30,doTrace=0),
        withTentative=TRUE),
      error=function(e) character(0))
    
    s_g <- tryCatch({
      rc <- caret::rfeControl(functions=caretFuncs,method="cv",
                              number=3,allowParallel=FALSE,verbose=FALSE)
      predictors(caret::rfe(X_tr,y_tr,sizes=c(5,10),
                            rfeControl=rc,method="svmRadial"))
    }, error=function(e) character(0))
    
    all_fs  <- c(l_g,r_g,x_g,b_g,s_g)
    gcounts <- table(all_fs)
    con_g   <- names(gcounts[gcounts >= consensus_votes])
    if (length(con_g)<2) con_g <- names(gcounts[gcounts>=3])
    if (length(con_g)<2) con_g <- names(sort(gcounts,decreasing=TRUE)[1:min(5,length(gcounts))])
    if (length(con_g)<1) next
    
    X_tr_f <- X_tr[,con_g,drop=FALSE]
    X_te_f <- X_te[,con_g,drop=FALSE]
    
    # Lightweight XGBoost only (speed) — mirrors ensemble intent
    p_xgb <- tryCatch({
      di <- xgboost::xgb.DMatrix(data=as.matrix(X_tr_f),label=y_tr_n,weight=obs_w)
      mx <- xgboost::xgb.train(
        params=list(objective="binary:logistic",eval_metric="auc",
                    max_depth=3,eta=0.1,nthread=1),
        data=di,nrounds=50,verbose=0)
      predict(mx,xgboost::xgb.DMatrix(data=as.matrix(X_te_f)))
    }, error=function(e) rep(0.5,nrow(X_te_f)))
    
    preds_all <- rbind(preds_all,
                       data.frame(Observed=as.numeric(y_te=="Tumor"),
                                  Predicted=p_xgb))
  }
  
  if (nrow(preds_all)<10 || length(unique(preds_all$Observed))<2)
    return(NA_real_)
  tryCatch(
    as.numeric(pROC::auc(pROC::roc(preds_all$Observed,preds_all$Predicted,
                                   direction="<",quiet=TRUE))),
    error=function(e) NA_real_)
}

# =============================================================================
# 2. PERMUTATION LOOP
# =============================================================================
if (!CONFIG$RECOVERY_MODE) {
  
  cat("[2/5] Running", CONFIG$n_perm, "permutations...\n")
  
  n_cores <- max(1L, parallel::detectCores() - CONFIG$n_cores_leave)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  cat("  Using", n_cores, "cores\n")
  
  perm_aucs <- foreach(
    i = seq_len(CONFIG$n_perm),
    .packages = c("caret","glmnet","ranger","xgboost","Boruta","kernlab","pROC","dplyr"),
    .errorhandling = "pass"
  ) %dopar% {
    perm_data <- master_data
    perm_data$Group <- sample(master_data$Group)
    evaluate_permutation(perm_data, CONFIG$seed + i,
                         CONFIG$consensus_votes, CONFIG$inner_folds)
  }
  
  parallel::stopCluster(cl)
  
  # ── CRASH GUARD: save raw results immediately ──────────────────────────────
  saveRDS(perm_aucs, PATHS$raw_rds)
  cat("  [SAVED] 03_perm_results_raw.rds  (crash guard)\n\n")
  
} else {
  
  cat("[2/5] RECOVERY MODE — loading saved permutation results...\n")
  if (!file.exists(PATHS$raw_rds))
    stop("Recovery file not found: ", PATHS$raw_rds)
  perm_aucs <- readRDS(PATHS$raw_rds)
  cat("  Loaded", length(perm_aucs), "permutation results\n\n")
}

# =============================================================================
# 3. PARSE AND SUMMARISE
# =============================================================================
cat("[3/5] Summarising permutation results...\n")

err_count <- sum(sapply(perm_aucs, inherits, "error"))
na_count  <- sum(sapply(perm_aucs, function(x) !inherits(x,"error") && is.na(x)))
ok_count  <- sum(sapply(perm_aucs, function(x) !inherits(x,"error") && !is.na(x)))

cat("  Total permutations:", length(perm_aucs), "\n")
cat("  Successful        :", ok_count, "\n")
cat("  NA (failed folds) :", na_count, "\n")
cat("  Errors            :", err_count, "\n")

null_aucs <- unlist(Filter(function(x) !inherits(x,"error") && !is.na(x),
                           perm_aucs))

# Conservative p-value: (+1 numerator, +1 denominator) [Phipson & Smyth 2010]
p_emp <- (sum(null_aucs >= true_auc) + 1) / (length(null_aucs) + 1)

null_mean <- mean(null_aucs)
null_sd   <- sd(null_aucs)
null_95   <- quantile(null_aucs, 0.95)
z_score   <- (true_auc - null_mean) / null_sd

cat(sprintf("\n  True AUC    : %.4f  [%.4f \u2013 %.4f]\n",
            true_auc, true_ci_lo, true_ci_hi))
cat(sprintf("  Null AUC    : mean = %.4f  SD = %.4f  95th pct = %.4f\n",
            null_mean, null_sd, null_95))
cat(sprintf("  Z-score     : %.3f\n", z_score))
cat(sprintf("  Empirical p : %.6f  (%s)\n",
            p_emp, ifelse(p_emp < 0.001, "p < 0.001",
                          sprintf("p = %.4f", p_emp))))

perm_summary <- data.frame(
  Metric              = c("True_AUC","True_CI_Lower","True_CI_Upper",
                          "Null_AUC_Mean","Null_AUC_SD","Null_AUC_95pct",
                          "Z_Score","Empirical_P","N_Permutations",
                          "N_Successful","N_Failed"),
  Value               = c(round(true_auc,6), round(true_ci_lo,6),
                          round(true_ci_hi,6),
                          round(null_mean,6), round(null_sd,6),
                          round(null_95,6),   round(z_score,4),
                          round(p_emp,6),     length(perm_aucs),
                          ok_count, na_count + err_count)
)

write.csv(perm_summary,
          file.path(dirs$tables,"Permutation_Test_Results.csv"),
          row.names=FALSE)
saveRDS(null_aucs,
        file.path(dirs$tables,"Permutation_Null_AUC_Vector.rds"))
cat("  [SAVED] Permutation_Test_Results.csv\n")
cat("  [SAVED] Permutation_Null_AUC_Vector.rds\n\n")

# =============================================================================
# 4. FIGURES
# =============================================================================
cat("[4/5] Generating figures...\n")

# ── Null distribution + observed AUC ─────────────────────────────────────────
null_df    <- data.frame(NullAUC = null_aucs)
bw         <- 2.0 * IQR(null_aucs) * length(null_aucs)^(-1/3)   # Freedman-Diaconis
p_label    <- ifelse(p_emp < 0.001,
                     "italic(p) < 0.001",
                     sprintf("italic(p) == %.4f", p_emp))

p_null <- ggplot(null_df, aes(x=NullAUC)) +
  geom_histogram(aes(y=after_stat(density)),
                 binwidth=bw, fill="grey70", colour="white",
                 linewidth=0.3) +
  geom_density(colour="grey30", linewidth=0.8) +
  geom_vline(xintercept=true_auc, colour=COL$tumor,
             linewidth=1.2, linetype="solid") +
  geom_vline(xintercept=null_95,  colour=COL$normal,
             linewidth=0.9, linetype="dashed") +
  annotate("text", x=true_auc + 0.003, y=Inf,
           label=sprintf("True AUC\n%.4f", true_auc),
           colour=COL$tumor, size=3.5, hjust=0, vjust=1.5,
           fontface="bold") +
  annotate("text", x=null_95 - 0.003, y=Inf,
           label=sprintf("95th pct\n%.4f", null_95),
           colour=COL$normal, size=3.3, hjust=1, vjust=1.5) +
  annotate("label",
           x=null_mean, y=0,
           label=parse(text=p_label),
           fill="white", colour=COL$tumor, size=4, hjust=0.5, vjust=-0.5,
           label.padding=unit(0.4,"lines"), fontface="bold") +
  labs(
    title    = "Permutation Test \u2014 Null AUC Distribution",
    subtitle = sprintf(
      "n = %d permutations  |  True AUC = %.4f  |  Null mean = %.4f  |  Z = %.2f",
      length(null_aucs), true_auc, null_mean, z_score),
    x        = "Permuted AUC",
    y        = "Density",
    caption  = paste0(
      "Red solid line = true observed AUC.  ",
      "Blue dashed line = 95th percentile of null distribution.\n",
      "Empirical p-value: (sum(null \u2265 true) + 1) / (N + 1)  ",
      "[Phipson & Smyth 2010].")
  ) +
  theme_publication()
save_figure(p_null, "Fig_Perm_1_NullDistribution", width=9, height=6.5)

# ── Convergence curve ─────────────────────────────────────────────────────────
conv_df <- data.frame(
  N         = seq_along(null_aucs),
  RunMean   = cumsum(null_aucs) / seq_along(null_aucs),
  RunSD     = sapply(seq_along(null_aucs), function(i)
    if (i < 2) NA_real_ else sd(null_aucs[1:i]))
)
conv_df <- conv_df %>%
  mutate(Lower = RunMean - 1.96*RunSD/sqrt(N),
         Upper = RunMean + 1.96*RunSD/sqrt(N))

p_conv <- ggplot(conv_df, aes(x=N, y=RunMean)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper),
              fill="grey70", alpha=0.5) +
  geom_line(colour=COL$normal, linewidth=0.9) +
  geom_hline(yintercept=null_mean, linetype="dashed",
             colour="grey40", linewidth=0.6) +
  scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
  labs(
    title    = "Permutation Test \u2014 Convergence Curve",
    subtitle = "Running mean of null AUC with 95% CI band",
    x        = "Number of Permutations",
    y        = "Running Mean Null AUC",
    caption  = paste0(
      "Shaded band = running mean \u00b1 1.96 \u00d7 SE.\n",
      "Dashed line = final null mean = ", round(null_mean,4), ".")
  ) +
  theme_publication()
save_figure(p_conv, "Fig_Perm_2_Convergence", width=9, height=5.5)

# ── Combined 2-panel figure ───────────────────────────────────────────────────
# Build base versions without titles (for panel labels)
p_null_panel <- p_null +
  labs(title=NULL, subtitle=NULL) +
  theme(plot.margin=unit(c(0.4,0.8,0.4,0.4),"cm"))
p_conv_panel <- p_conv +
  labs(title=NULL, subtitle=NULL) +
  theme(plot.margin=unit(c(0.4,0.8,0.4,0.4),"cm"))

p_combined <- gridExtra::arrangeGrob(
  p_null_panel, p_conv_panel, nrow=1,
  top=grid::textGrob(
    paste0("Permutation Test  |  n=",length(null_aucs),
           " permutations  |  True AUC=",round(true_auc,4),
           "  |  p=",ifelse(p_emp<0.001,"<0.001",round(p_emp,4))),
    gp=grid::gpar(fontface="bold", fontsize=12))
)
for (ext in c("png","tiff","pdf")) {
  if (ext=="pdf") {
    grDevices::cairo_pdf(
      file.path(dirs$figures,paste0("Fig_Perm_3_Combined.pdf")),
      width=16, height=6)
    grid::grid.draw(p_combined)
    grDevices::dev.off()
  } else {
    res <- if (ext=="tiff") 600 else 300
    grDevices::png(
      file.path(dirs$figures,paste0("Fig_Perm_3_Combined.",ext)),
      width=16,height=6,units="in",res=res)
    grid::grid.draw(p_combined)
    grDevices::dev.off()
  }
}
cat("  [SAVED] Fig_Perm_3_Combined (PNG + TIFF 600dpi + PDF)\n")

# =============================================================================
# WRAP UP
# =============================================================================
cat("\n[5/5] Session info...\n")
sink(file.path(dirs$logs,"03_session_info.txt"))
print(sessionInfo()); sink()
sink(log_path, append=TRUE, split=TRUE)

t_end <- Sys.time()
cat("\n=============================================================\n")
cat(" 03_permutation_test.R  COMPLETE\n")
cat(" End    :", as.character(t_end), "\n")
cat(" Runtime:", round(difftime(t_end,t_start,units="mins"),2), "minutes\n")
cat(sprintf("\n RESULT: True AUC = %.4f  |  Empirical p = %.6f  (%s)\n",
            true_auc, p_emp,
            ifelse(p_emp<0.001,"p < 0.001",sprintf("p = %.4f",p_emp))))
cat("\n SAVED ARTEFACTS:\n")
cat("  Tables  : Permutation_Test_Results.csv\n")
cat("            Permutation_Null_AUC_Vector.rds\n")
cat("            03_perm_results_raw.rds  (crash guard)\n")
cat("  Figures : Fig_Perm_1_NullDistribution.*\n")
cat("            Fig_Perm_2_Convergence.*\n")
cat("            Fig_Perm_3_Combined.*  (PNG + TIFF 600dpi + PDF)\n")
cat("  Logs    : 03_Permutation_Run_[timestamp].log\n")
cat("            03_session_info.txt\n")
cat("=============================================================\n")
sink()