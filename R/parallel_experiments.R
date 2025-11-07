#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  source("helpers.R"); source("gen_models.R"); source("metrics.R")
})

# ---- config (same as yours) ----
P <- 100
N <- 200
N_REPS <- 100
MODELS <- 1:6
METHODS <- c("glasso","rope","gelnet")
ALPHAS_DEFAULT <- c(glasso=1, rope=0, gelnet=0.5)
TARGETS <- c("None","Identity","vIdentity","Eigenvalue","MSC", "Regression", "TrueDiag")
PENALIZE_DIAG <- c(TRUE, FALSE)
lambda_grid_dense <- 10^seq(-3, 0, length.out=121)
alpha_grid <- c(0.2, 0.35, 0.5, 0.65, 0.8)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

# ---- array indexing: map SLURM_ARRAY_TASK_ID -> (model_id, rep_idx) ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # Allow local test: runner.R <task_id>
  task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
} else {
  task_id <- as.integer(args[[1]])
}
set.seed(42 + task_id)

# Create mapping once (same logic as in sbatch)
all_pairs <- expand.grid(
  model_id = MODELS,
  rep_idx  = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>% arrange(model_id, rep_idx)

stopifnot(task_id >= 1, task_id <= nrow(all_pairs))
model_id <- all_pairs$model_id[task_id]
rep_idx  <- all_pairs$rep_idx[task_id]

message(sprintf("Task %d -> model %d, rep %d", task_id, model_id, rep_idx))

# ---- generate data for this (model_id, rep_idx) ----
if (model_id == 7) {
  B_list <- make_ba_conditional_hamiltonian(P, m_attach = 2, add_budget = 0.02)
  sim <- sample_gsm(B_list, n = N, beta = 0.1, lambda0 = 0.5, a = 3, b = 3)
  Y <- sim$Y; Sigma <- sim$Sigma; Theta_true <- sim$Theta
} else {
  mm <- make_model(model_id, P)
  Sigma <- mm$Sigma; Theta_true <- mm$Theta

  Sigma <- sanitize_cov(Sigma)
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) {
    Sigma <- Sigma + (abs(min(ev)) + 1e-6) * diag(nrow(Sigma))
  }
  Y <- mvtnorm::rmvnorm(N, sigma = Sigma)
}

# ---- parameter grid inside task ----
param_grid <- tidyr::expand_grid(
  method = METHODS,
  penalize_diag = PENALIZE_DIAG,
  target = TARGETS,
  alpha = alpha_grid
) %>%
  mutate(target_eff = ifelse(!penalize_diag, "None", target))

# ---- (optional) light intra-task parallel over grid ----
# If you want this, set BLAS threads to 1 in sbatch (we do below).
use_parallel <- as.logical(Sys.getenv("R_INTRA_PAR", "FALSE"))

eval_one <- function(row, Y_local, Theta_true_local) {
  method <- row$method
  pen_diag <- row$penalize_diag
  target_eff <- row$target_eff
  alpha <- row$alpha

  lam <- cv_select_lambda(
    Y = Y_local, trueTheta = Theta_true_local, method = method, alpha = alpha,
    target_type = target_eff, penalize_diag = pen_diag,
    lambda_grid = lambda_grid_dense
  )

  fit <- fit_method(S = cor(Y_local), Y = Y_local, trueTheta = Theta_true_local, method = method, alpha = alpha,
                    target_type = target_eff, penalize_diag = pen_diag, lambda = lam)
  Theta_hat <- fit$Theta

  KL  <- kl_loss_safe(Sigma, Theta_hat)
  L2  <- l2_loss(Theta_true, Theta_hat)
  SP  <- sp_loss(Theta_true_local, Theta_hat)
  gs  <- graph_scores(Theta_true_local, Theta_hat, eps=1e-5)
  cur <- curve_points(Theta_true_local, Theta_hat)
  PRAUC <- pr_auc(cur); ROCAUC <- roc_auc(cur)

  tibble(
    model = model_id, rep = rep_idx,
    method = method, alpha = alpha,
    penalize_diag = pen_diag, target = target_eff, lambda = lam,
    KL = KL, L2 = L2, SP = SP,
    edges = count_edges(Theta_hat),
    F1 = gs["F1"], MCC = gs["MCC"],
    TP = gs["TP"], TN = gs["TN"], FP = gs["FP"], FN = gs["FN"],
    PRAUC = PRAUC, ROCAUC = ROCAUC
  )
}

if (use_parallel) {
  # Simple base-parallel (fork) on Linux
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
  suppressWarnings({
    res_list <- parallel::mclapply(
      split(param_grid, seq_len(nrow(param_grid))), 
      function(row) eval_one(row, Y_local = Y, Theta_true_local = Theta_true), mc.cores = max(1L, cores)
    )
  })
  err <- vapply(res_list, inherits, logical(1), "data.frame")  # check no errors
  if (any(!err)) stop("Error in eval_one:", paste(which(!err), collapse = ", "))
  df <- bind_rows(res_list)
} else {
  df <- param_grid %>% pmap_dfr(~eval_one(tibble(method=..1, penalize_diag=..2, target=..3, alpha=..4,
                                                 target_eff=ifelse(!..2, "None", ..3)),
                                          Y_local = Y, Theta_true_local = Theta_true))
}

# ---- write per-task CSV ----
outf <- sprintf("results/scores_model_%d_rep_%d.csv", model_id, rep_idx)
readr::write_csv(df, outf)
message("Wrote: ", outf)
