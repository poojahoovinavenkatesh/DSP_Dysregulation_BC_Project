# Changelog

## [v1.1] — Post-analysis corrections to `scripts/1_training_pipeline.R`

These corrections were identified during code review after the analysis
was completed and when appiled on other cancer types. **They were not applied to the submitted/published results in breast cancer manuscript.** The
figures, tables and model artefacts reported in the manuscript were
produced by the original script (`v1.0-submission`), recorded in
`logs/01_Pipeline_Run_20260312_1616.log`, with the frozen artefacts in
`models/`.

Every change is marked in the source with a `[FIX n]` comment.

### Corrections made to breast cancer manuscript figures before submission

| Figure | Change | Reason |
|--------|--------|--------|
| Figure 3A, left panel (XGBoost) | Label "Sp=0.96" corrected to "Sp=0.97" | FIX 8, verified against frozen model |
| Figure 3A, right panel (Directional Score) | DeLong clause removed from subtitle | In-sample comparator (see Known issues) |
| Figure 2 | UpSet panel (Fig02_UpSet) omitted | FIX 1; standalone analysis, not the CV feature selection |

No numerical result reported in the text changed as a consequence of these
corrections, apart from any quotation of XGBoost specificity at the optimal
operating point, which is 0.97.

---

### FIX 1 — Parallel backend not deregistered after cluster shutdown

`parallel::stopCluster(cl)` was not followed by `foreach::registerDoSEQ()`.
The stopped cluster remained registered as `foreach`'s backend with dead
worker sockets, so any later `%dopar%` call failed with
`"invalid connection"`.

**Affected:** the standalone UpSet visualisation block only. `caret::rfe`
invokes `foreach` internally; in the master process after
`stopCluster()`, SVM-RFE failed, returned an empty vector, and was
silently removed by `Filter(length > 0)`. `Feature_Selection_Summary.txt`
and `Fig02_UpSet` consequently showed four methods rather than five.

**Not affected:** the 100-fold nested cross-validation. It executes inside
`%dopar%` worker processes, which do not inherit the parent session's
`foreach` registration, so `caret::rfe` ran sequentially and completed
normally in every fold. All five feature selection methods contributed to
signature discovery. Re-executing the identical `caret::rfe` call after
`registerDoSEQ()` confirmed that SVM-RFE is functional.

### FIX 2 — Per-fold feature selection lists were computed but discarded

Each fold returned `FSLists` (the gene vector from each of the five
methods), but the object was never extracted or saved. It is now written
to `01_cv_fold_FSLists.rds`, and the log reports how many of the 100 folds
each method returned genes in, flagging any silent failures.

### FIX 3 — SVM-RFE errors were swallowed

Both `caret::rfe` calls used `error = function(e) character(0)`, discarding
the error message. Errors are now reported via `conditionMessage(e)`.

### FIX 4 — UpSet `rfeControl()` inconsistent with the CV call

The UpSet block's `rfeControl()` omitted `allowParallel = FALSE`, which the
cross-validation call sets explicitly. Now consistent.

### FIX 5 — Methods dropped from the UpSet list without warning

A warning is now emitted naming any method that returned zero genes and was
excluded from `upset_list`.

### FIX 6 — Hardcoded "5-METHOD" header

`Feature_Selection_Summary.txt` printed "5-METHOD FEATURE SELECTION
SUMMARY" regardless of how many methods were written. The header is now
generated from `length(upset_list)`, and a note states that the file
reflects the standalone single-pass analysis, not the cross-validation.

### FIX 7 — Ensemble models failed silently to a constant 0.5

Each of the four ensemble models (XGBoost, random forest, elastic net,
SVM) used `error = function(e) rep(0.5, nrow(X_te_f))`. A failing model
contributed a constant 0.5 to the ensemble average with nothing logged.
Failures are now reported. Given the observed nested-CV performance,
widespread silent failure in the published run is implausible, but no
record exists either way.

### FIX 8 — Wrong variable in the Fig03 annotation

The optimal-operating-point label in `Fig03_XGB_Diagnostic_ROC` referenced
`youden_dir` (the Directional Score's Youden point) instead of `youden`
(XGBoost's). The plotted diamond marker used the correct variable; only
the label's position and printed Se/Sp values were affected. Because
`youden_dir` is defined later in the script, this also caused a crash in
any clean-session run.

**Effect on the Submitted/published figure (manuscript Figure 3A, left panel).** In
the original run the XGBoost panel's label read "Se=0.98, Sp=0.96", which
are the Directional Score's values. The cross-validated XGBoost ROC was
rebuilt from the frozen artefacts in `models/` using the script's own fold
seeds; it reproduced the published AUC exactly (0.9966) and gave the
correct optimal operating point as Se = 0.983, Sp = 0.973. The label in the
manuscript figure was corrected to "Se=0.98, Sp=0.97" before submission.
The curve, confidence band, AUC and marker position were unaffected; the
correctly placed label would sit about 0.014 further left on the
false-positive-rate axis, which is not visible at figure scale. The
Directional Score panel's label (Se=0.98, Sp=0.96) was correct and is
unchanged.

### FIX 9 — No seed immediately before Bayesian optimisation

`ParBayesianOptimization::bayesOpt()` draws its initial design from
whatever random number state the preceding code leaves behind; no
`set.seed()` intervened between the UpSet block and the optimisation call.
Any change to random number consumption earlier in the script therefore
altered the initial points and, through them, the selected
hyperparameters. `set.seed(CONFIG$seed)` is now called immediately before
the `scoring_fn` definition, so Section 6 is reproducible independently of
upstream code.

### FIX 10 — In-sample XGBoost ROC used as the DeLong comparator

`roc_xgb_full` was computed from `predict(final_model, X_fin_sc)`, scoring
the model on the same cohort it was trained on. Compared against the
Directional Score, which fits no parameters beyond the gene direction
signs, this gives XGBoost an overfitting advantage; the resulting DeLong
test measured that advantage rather than a genuine difference in
discrimination.

The comparator is now the five-fold cross-validated XGBoost ROC
(`roc_diag`), which is out-of-fold and therefore comparable. A switch,
`XGB_COMPARATOR`, is provided at that point in the script and defaults to
`"cv"`; setting it to `"resubstitution"` restores the previous behaviour
and emits a warning. The annotation in `Fig03C_XGB_vs_Dir_ROC_Overlay` now
states which comparator was used.

Because the in-sample and cross-validated XGBoost AUCs differ
substantially (0.9998 versus 0.9966 in the submitted published run), this fix
changes the DeLong p-value reported by the script. Any comparison run on a
new dataset will use the cross-validated comparator by default.

---

### Reproducibility of submitted/published results

Re-running the corrected script will **not** reproduce the final XGBoost
model exactly. FIX 1 allows SVM-RFE to complete in the UpSet block, which
changes the amount of random number generation consumed before
`ParBayesianOptimization::bayesOpt()`. No `set.seed()` intervenes between
the two, so the Bayesian optimisation draws different initial points and
converges on different hyperparameters.

**Reproduce exactly from the corrected script:**
cohort composition, nested cross-validation performance, the stable gene
signature, differential expression results, and all Directional Score
results.

**Will differ after re-running:**
the Bayesian-optimised hyperparameters and `final_xgboost_model_Bayesian.rds`,
SHAP values and their tertiles, and any XGBoost-derived AUC in downstream
scripts that load the model.

To reproduce the submitted/published analyses exactly, use `v1.0-submission` with the
artefacts in `models/`, or run downstream scripts against the frozen
artefacts without re-running this script.

### Cross-platform numerical variation

Independently of these corrections, SHAP magnitudes vary by up to ~5%
between Linux and Windows builds, owing to differing BLAS/LAPACK
implementations affecting XGBoost tree construction and the optimisation
trajectory. Gene ranking and biological direction are unaffected.
Published values correspond to the Linux run documented in
`logs/01_session_info.txt`.

---

### Notes on applying this script to a new dataset

All ten fixes above are in force, so a run on a different dataset
reproduces its own Bayesian hyperparameters deterministically (FIX 9) and
reports a DeLong comparison against the cross-validated XGBoost ROC rather
than an in-sample one (FIX 10).

Two behaviours remain by design and are worth knowing before interpreting
output on new data.

**Method selectivity is uneven.** Boruta is run with
`withTentative = TRUE`, and `caret::rfe` evaluates the full-size model
alongside the requested sizes and keeps whichever performs best. On a
candidate panel of modest size, both can return nearly every gene,
contributing a vote to almost all of them. The per-fold consensus vote is
therefore not the principal selective step; the 0.80 stability threshold
across the 100 folds is. Report the signature as arising from stability
selection, not from the vote alone.

**The UpSet block is a separate analysis.** It runs once on the full
cohort with `upSample`-balanced classes and no elastic-net pre-filter, so
`Feature_Selection_Summary.txt` and `Fig02_UpSet` do not describe the
cross-validated feature selection. Use `01_cv_fold_FSLists.rds` (FIX 2)
for per-method contribution within the folds.

---

### Relationship to the submitted/published breast cancer manuscript results

`Fig03C_XGB_vs_Dir_ROC_Overlay` is not reported in the breast cancer manuscript. In
Figure 3A, the clause "Direction vs XGBoost DeLong: p=0.0041" was removed
from the subtitle of the Directional Score panel before submission,
because printed beside the cross-validated XGBoost panel it would imply a
comparison between the two curves shown, when it was computed against the
in-sample model. No DeLong p-value comparing the two scores is reported in
the breast cancer manuscript. FIX 10 removes the cause of that discrepancy for future
runs.
