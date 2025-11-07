
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr); library(tidyr); library(purrr); library(readr)
  source("helpers.R")
  source("metrics.R")
  source("gen_models.R")
})
option_list <- list(
  make_option(c("--model"), type="integer", default=7),
  make_option(c("--edges"), type="integer", default=200),
  make_option(c("--nreps"), type="integer", default=1),
  make_option(c("--p"), type="integer", default=100),
  make_option(c("--n"), type="integer", default=200),
  make_option(c("--seed"), type="integer", default=1)
)

stability_one <- function(Y, fit_args, edges_target){
  Sfull <- cor(Y)
  # λ that matches edge budget on full S
  fm <- lambda_for_edges(function(lam){
    fit_method(S=Sfull, lambda=lam, method=fit_args$method, alpha=fit_args$alpha,
               target_type=fit_args$target_type, penalize_diag=fit_args$penalize_diag, trueTheta=fit_args$trueTheta)$Theta
  }, target_edges=edges_target, grid=lambda_grid_dense)
  lam <- fm$lambda

  p <- ncol(Y); up <- upper.tri(diag(p)); votes <- matrix(0,p,p)
  for (b in 1:B) {
    idx <- sample(1:nrow(Y), replace=TRUE)
    Sb  <- cor(Y[idx, , drop=FALSE])
    Th  <- fit_method(S=Sb, lambda=lam, method=fit_args$method, alpha=fit_args$alpha,
                      target_type=fit_args$target_type, penalize_diag=fit_args$penalize_diag, trueTheta=fit_args$trueTheta)$Theta
    A   <- (abs(Th) > 1e-12)*1; diag(A) <- 0
    votes <- votes + A
  }
  stab <- votes / B
  list(lambda=lam, stab=stab, Theta=fm$Theta, edges=fm$edges)
}

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)
dir.create("results", showWarnings = FALSE)
METHODS <- c("glasso","rope","gelnet")
ALPHAS  <- c(glasso=1, rope=0, gelnet=0.5)
TARGETS <- c("None","Identity","vIdentity","Eigenvalue","MSC")
PENALIZE_DIAG <- c(TRUE, FALSE)
lambda_grid_dense <- 10^seq(-3, 0, length.out=121)
alpha_grid <- c(0.2, 0.35, 0.5, 0.65, 0.8)
for (rep_idx in 1:opt$nreps) {
  if (opt$model == 7) {
    B_list <- make_ba_conditional_hamiltonian(opt$p, m_attach = 2, add_budget = 0.02, seed = opt$seed + rep_idx)
    B <- length(B_list) - 1
    print(sprintf("Generated BA conditional Hamiltonian with %d communities and edge budget %d", length(B_list)-1, B))
    sim <- sample_gsm(B_list, n = opt$n, beta = 0.1, lambda0 = 0.5, a = 3, b = 3, seed = opt$seed + rep_idx)
    Y <- sim$Y; Sigma <- sim$Sigma; Theta_true <- sim$Theta
  } else {
    break
  }
  Sfull <- cor(Y)
  rows <- list()
  for (method in METHODS) {
    alpha <- ALPHAS[[method]]
    for (pen_diag in PENALIZE_DIAG) {
      for (target in TARGETS) {
        target_eff <- if (!pen_diag) "None" else target
        #fit_fn <- function(lam) {
        #  fit <- fit_method(S=Sfull, method=method, alpha=alpha,
        #                    target_type=target_eff, penalize_diag=pen_diag, lambda=lam)
        #  fit$Theta
        #}
        #fm <- lambda_for_edges(fit_fn, target_edges = opt$edges, grid = lambda_grid_dense)
        fm <- stability_one(Y,
                              fit_args = list(method=method, alpha=alpha,
                                              target_type=target_eff, penalize_diag=pen_diag, trueTheta=Theta_true),
                              edges_target = opt$edges)
        Theta_hat <- fm$Theta; lam <- fm$lambda; ecount <- fm$edges
        KL  <- kl_loss_safe(Sigma, Theta_hat)
        L2  <- l2_loss(Theta_true, Theta_hat)
        SP  <- sp_loss(Theta_true, Theta_hat)
        cur <- curve_points(Theta_true, Theta_hat)
        graph <- graph_scores(Theta_true, Theta_hat)
        PRAUC <- pr_auc(cur); ROCAUC <- roc_auc(cur)
        rows[[length(rows)+1]] <- data.frame(
          model=opt$model, rep=rep_idx, method=method, alpha=alpha,
          penalize_diag=pen_diag, target=target_eff,
          lambda=lam, edges=ecount,
          F1=gs["F1"], MCC=gs["MCC"],
          TP=gs["TP"], TN=gs["TN"], FP=gs["FP"], FN=gs["FN"],
          KL=KL, L2=L2, SP=SP, PRAUC=PRAUC, ROCAUC=ROCAUC,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  df <- dplyr::bind_rows(rows)
  readr::write_csv(df, sprintf("results/edge_matched_model_%d_rep_%03d.csv", opt$model, rep_idx))
  cat("Saved results for rep", rep_idx, "\n")
}
all <- purrr::map_dfr(list.files("results", pattern=sprintf("^edge_matched_model_%d_rep_.*\\.csv$", opt$model), full.names = TRUE), readr::read_csv, show_col_types=FALSE)
readr::write_csv(all, sprintf("results/edge_matched_model_%d_all.csv", opt$model))
cat("Done edge-matched.\n")
