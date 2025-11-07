#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr); library(tidyr); library(purrr); library(readr)
  source("helpers.R")
  source("metrics.R")
  source("gen_models.R")
})

# ----- CLI options -----
option_list <- list(
  make_option(c("--model"),  type="integer", default=7),
  make_option(c("--edges"),  type="integer", default=200),
  make_option(c("--nreps"),  type="integer", default=100),
  make_option(c("--p"),      type="integer", default=100),
  make_option(c("--n"),      type="integer", default=200),
  make_option(c("--seed"),   type="integer", default=1),
  make_option(c("--boots"),  type="integer", default=200,  help="Bootstrap resamples B"),
  make_option(c("--workers"),type="integer", default=as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","4")),
              help="Parallel workers for bootstrap or grid"),
  make_option(c("--grid-par"),action="store_true", default=TRUE,
              help="Parallelize across method/target grid instead of bootstrap")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$grid_par)) opt$grid_par <- TRUE
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ----- Constants / grids -----
METHODS <- c("glasso","rope","gelnet")
ALPHAS  <- c(glasso=1, rope=0, gelnet=0.5)   # per-method default alpha
TARGETS <- c("None","Identity","vIdentity","Eigenvalue","MSC", "Regression", "TrueDiag")
PENALIZE_DIAG <- c(TRUE, FALSE)
lambda_grid_dense <- 10^seq(-3, 0, length.out=121)
alpha_grid <- c(0.2, 0.35, 0.5, 0.65, 0.8)   # if you want to sweep alpha too

# ----- Map SLURM_ARRAY_TASK_ID -> rep_idx so each array task does one rep -----
rep_idx <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(rep_idx)) rep_idx <- 1  # local test
stopifnot(rep_idx >= 1, rep_idx <= opt$nreps)

# Reproducible RNG streams per task
set.seed(opt$seed + rep_idx)
RNGkind("L'Ecuyer-CMRG")

# ----- Data generation for this rep -----
if (opt$model == 7) {
  B_list <- make_ba_conditional_hamiltonian(opt$p, m_attach = 2, add_budget = 0.02, seed = opt$seed + rep_idx)
  n_comms <- length(B_list) - 1
  message(sprintf("Rep %d: BA conditional Hamiltonian with %d communities", rep_idx, n_comms))
  sim <- sample_gsm(B_list, n = opt$n, beta = 0.1, lambda0 = 0.5, a = 3, b = 3, seed = opt$seed + rep_idx)
  Y <- sim$Y; Sigma <- sim$Sigma; Theta_true <- sim$Theta
} else {
  stop("This runner currently implemented for --model 7 only (BA conditional).")
}

Sfull <- cor(Y)

# ----- Stability selection helper (parallelizable) -----
# instead of cross-validating over lambda, we pick lambda to match edge budget
# this is important for comparing methods fairly and is the main point of the experiment
# cross-validation is very slow here and not needed for the experiment, because we want to control edges directly and not via likelihood ()
stability_one <- function(Y, fit_args, edges_target, B, workers, par_boot = TRUE) {
  Sfull <- cor(Y)

  # λ matching the edge budget on the full S
  fm <- lambda_for_edges(
    function(lam) {
      fit_method(
        S = Sfull, Y = Y, trueTheta = fit_args$trueTheta, lambda = lam,
        method = fit_args$method, alpha = fit_args$alpha,
        target_type = fit_args$target_type, penalize_diag = fit_args$penalize_diag
      )$Theta
    },
    target_edges = edges_target, grid = lambda_grid_dense
  )
  lam <- fm$lambda

  p <- ncol(Y)
  votes <- matrix(0, p, p)
  boot_eval <- function(b) {
    # independent bootstrap + RNG stream
    set.seed(opt$seed*100000 + rep_idx*1000 + b)
    idx <- sample.int(nrow(Y), replace = TRUE)
    Sb  <- cor(Y[idx, , drop = FALSE])
    Th  <- fit_method(
      S = Sb, Y = Y[idx, , drop = FALSE], trueTheta = fit_args$trueTheta, lambda = lam,
      method = fit_args$method, alpha = fit_args$alpha,
      target_type = fit_args$target_type, penalize_diag = fit_args$penalize_diag
    )$Theta
    A <- (abs(Th) > 1e-12) * 1; diag(A) <- 0
    A
  }

  if (par_boot) {
    # R-level parallel across bootstrap draws
    res_list <- parallel::mclapply(seq_len(B), boot_eval, mc.cores = max(1L, workers))
    for (A in res_list) votes <- votes + A
  } else {
    for (b in seq_len(B)) votes <- votes + boot_eval(b)
  }

  stab <- votes / B
  list(lambda = lam, stab = stab, Theta = fm$Theta, edges = fm$edges)
}

# ----- Parameter grid inside this rep -----
param_grid <- tidyr::expand_grid(
  method = METHODS,
  penalize_diag = PENALIZE_DIAG,
  target = TARGETS,
  alpha = alpha_grid
) %>%
  mutate(target_eff = ifelse(!penalize_diag, "None", target))

# ----- Evaluate grid (either parallel over grid OR parallel inside stability bootstrap) -----
eval_one <- function(row, Y_local, Theta_true_local) {
  method <- row$method
  pen_diag <- row$penalize_diag
  target_eff <- row$target_eff
  alpha <- row$alpha

  fm <- stability_one(
    Y_local,
    fit_args = list(method = method, alpha = alpha,
                    target_type = target_eff, penalize_diag = pen_diag, trueTheta = Theta_true_local),
    edges_target = opt$edges,
    B = opt$boots,
    workers = opt$workers,
    par_boot = !opt$grid_par
  )

  Theta_hat <- fm$Theta; lam <- fm$lambda; ecount <- fm$edges

  KL  <- kl_loss_safe(Sigma, Theta_hat)
  L2  <- l2_loss(Theta_true, Theta_hat)
  SP  <- sp_loss(Theta_true, Theta_hat)
  cur <- curve_points(Theta_true, Theta_hat)
  gs  <- graph_scores(Theta_true, Theta_hat, eps = 1e-5)
  PRAUC <- pr_auc(cur); ROCAUC <- roc_auc(cur)

  tibble(
    model = opt$model, rep = rep_idx,
    method = method, alpha = alpha,
    penalize_diag = pen_diag, target = target_eff,
    lambda = lam, edges = ecount,
    F1 = gs["F1"], MCC = gs["MCC"],
    TP = gs["TP"], TN = gs["TN"], FP = gs["FP"], FN = gs["FN"],
    KL = KL, L2 = L2, SP = SP, PRAUC = PRAUC, ROCAUC = ROCAUC
  )
}

if (opt$grid_par) {
  # Parallel across grid; set BLAS threads=1 in sbatch
  rows <- parallel::mclapply(split(param_grid, seq_len(nrow(param_grid))),
                             function(row) eval_one(row, Y_local, Theta_true_local), mc.cores = max(1L, opt$workers), mc.preschedule = TRUE)
  err <- vapply(rows, inherits, logical(1), "data.frame")  # check no errors
  if (any(!err)) stop("Error in eval_one:", paste(which(!err), collapse = ", "))
  df <- dplyr::bind_rows(rows)
} else {
  # Serial grid; stability_one parallelizes across bootstrap
  #df <- param_grid %>% pmap_dfr(eval_one)
  df <- param_grid %>% pmap_dfr(~eval_one(tibble(method=..1, penalize_diag=..2, target=..3, alpha=..4,
                                                 target_eff=ifelse(!..2, "None", ..3)),
                                          Y_local = Y, Theta_true_local = Theta_true))
}

outf <- sprintf("results/edge_matched_model_%d_rep_%03d.csv", opt$model, rep_idx)
readr::write_csv(df, outf)
cat("Saved:", outf, "\n")

# (Optional) Single-run aggregation is done in a separate script/job
