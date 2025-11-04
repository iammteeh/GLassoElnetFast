
# R/reproduce_simulations.R — main entry to reproduce Section 4 (simulations)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  source("helpers.R")
})

set.seed(42)

# Configuration
P <- 100
N <- 200          # sample size per replication (paper uses variety; pick a moderate N)
N_REPS <- 10     # set smaller for a quick smoke test
MODELS <- 1:6
METHODS <- c("glasso","rope","gelnet")
ALPHAS <- c(glasso=1, rope=0, gelnet=0.5)
TARGETS <- c("None","Identity","vIdentity","Eigenvalue","MSC")
PENALIZE_DIAG <- c(TRUE, FALSE)   # "F" in paper == FALSE
lambda_grid <- 0.9^(0:40)

dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

run_one <- function(model_id, rep_idx) {
  mm <- make_model(model_id, P)
  Sigma <- mm$Sigma; Theta_true <- mm$Theta
  Y <- mvtnorm::rmvnorm(N, sigma = Sigma)
  S <- cor(Y) # recommended scaling
  tibble()
}

results <- list()

for (model_id in MODELS) {
  cat("== Model", model_id, "==\n")
  mm <- make_model(model_id, P)
  Sigma <- mm$Sigma; Theta_true <- mm$Theta

  for (rep_idx in 1:N_REPS) {
    # wherever Sigma is created/loaded:
    Sigma <- sanitize_cov(Sigma)
    if (any(!is.finite(Sigma))) print(head(which(!is.finite(Sigma), arr.ind=TRUE), 10))
    ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    cat(sprintf("min eig = %.3e\n", min(ev)))
    Y <- mvtnorm::rmvnorm(N, sigma = Sigma)

    for (method in METHODS) {
      alpha <- ALPHAS[[method]]
      for (pen_diag in PENALIZE_DIAG) {
        for (target in TARGETS) {
          # TrueDiag target only valid in simulations
          if (target == "TrueDiag") {
            T_ok <- TRUE
          } else {
            T_ok <- TRUE
          }
          cat(sprintf("Model %d, Rep %d, Method %s, Alpha %.2f, PenDiag %s, Target %s\n",
                      model_id, rep_idx, method, alpha, pen_diag, target))
          # If diagonal not penalized and target != None, it's effectively no target (paper note)
          target_eff <- if (!pen_diag) "None" else target

          # Select lambda by 5-fold CV
          lam <- cv_select_lambda(Y=Y, method=method, alpha=alpha,
                                  target_type=target_eff, penalize_diag=pen_diag,
                                  lambda_grid=lambda_grid)

          # Fit on full data
          fit <- fit_method(S=cor(Y), method=method, alpha=alpha,
                            target_type=target_eff, penalize_diag=pen_diag, lambda=lam)
          Theta_hat <- fit$Theta

          # Metrics
          KL  <- kl_loss_safe(Sigma, Theta_hat)
          L2  <- l2_loss(Theta_true, Theta_hat)
          SP  <- sp_loss(Theta_true, Theta_hat)
          gs  <- graph_scores(Theta_true, Theta_hat, eps=1e-5)

          results[[length(results)+1]] <- data.frame(
            model=model_id, rep=rep_idx,
            method=method, alpha=alpha,
            penalize_diag=pen_diag, target=target_eff, lambda=lam,
            KL=KL, L2=L2, SP=SP,
            F1=gs["F1"], MCC=gs["MCC"],
            TP=gs["TP"], TN=gs["TN"], FP=gs["FP"], FN=gs["FN"],
            stringsAsFactors = FALSE
          )
          cat("  Done.\n")
          # show result row
          print(results[[length(results)]])
        }
      }
    }
  }

  df <- bind_rows(results)
  readr::write_csv(df, sprintf("results/scores_model_%d.csv", model_id))
  results <- list()
}

# Aggregate and plotting (lightweight summaries)
suppressPackageStartupMessages({ library(ggplot2) })

all_scores <- purrr::map_dfr(MODELS, ~readr::read_csv(sprintf("results/scores_model_%d.csv", .x),
                                                       show_col_types=FALSE))
summ <- all_scores %>%
  group_by(model, method, target, penalize_diag) %>%
  summarise(across(c(KL,L2,SP,F1,MCC), list(mean=mean, sd=sd), .names="{.col}_{.fn}"),
            .groups="drop")

readr::write_csv(summ, "results/summary.csv")

for (metric in c("KL","L2","SP","F1","MCC")) {
  p <- ggplot(summ, aes(x=interaction(method,target,penalize_diag),
                        y=.data[[paste0(metric,"_mean")]])) +
    geom_col() +
    geom_errorbar(aes(ymin=.data[[paste0(metric,'_mean')]] - .data[[paste0(metric,'_sd')]],
                      ymax=.data[[paste0(metric,'_mean')]] + .data[[paste0(metric,'_sd')]]),
                  width=0.2) +
    facet_wrap(~model, scales="free_y") +
    labs(x="Method–Target–PenDiag", y=metric, title=paste(metric, "summary by model")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(sprintf("plots/%s_summary.pdf", metric), p, width=12, height=8)
}

cat("Done. See results/ and plots/.\n")
