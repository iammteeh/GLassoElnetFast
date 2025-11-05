
# R/helpers.R — utilities for simulation, targets, metrics, and CV

suppressPackageStartupMessages({
  library(GLassoElnetFast)
  library(glasso)
  library(rags2ridges)
  library(gcdnet)
  library(igraph)
  library(Matrix)
  library(mvtnorm)
  library(dplyr); library(tidyr); library(purrr)
})

.symm <- function(M) (M + t(M)) / 2
.pd_nudge <- function(M, tol = 1e-12, add = 1e-8) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= tol) M <- M + (abs(min(ev)) + add) * diag(nrow(M))
  .symm(M)
}
.diag_inv <- function(D) diag(1 / diag(D))

sanitize_cov <- function(S) {
  cat("sanitize")
  if (any(!is.finite(S))) {
    bad <- which(!is.finite(S), arr.ind = TRUE)
    stop(sprintf("Sigma has non-finite entries; first at [%d,%d]", bad[1,1], bad[1,2]))
  }
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) {
    S <- S + (abs(min(ev)) + 1e-6) * diag(nrow(S))
  }
  
  return(S)
}

# ---------- Targets (diagonal) ----------
diag_target <- function(S, type = c("None","Identity","vIdentity","Eigenvalue","MSC","Reg","TrueDiag"),
                        Y = NULL, trueTheta = NULL, use_correlation = TRUE) {
  type <- match.arg(type)
  if (type == "None") return(NULL)
  cor_flag <- isTRUE(use_correlation)
  
  pass_Y <- !is.null(Y)
  switch(type,
         "Identity"    = if (pass_Y) target(Y = Y,        type = "Identity",   cor = cor_flag) else target(S = S, type = "Identity",   cor = cor_flag),
         "vIdentity"   = if (pass_Y) target(Y = Y,        type = "vI",         cor = cor_flag) else target(S = S, type = "vI",         cor = cor_flag),
         "Eigenvalue"  = if (pass_Y) target(Y = Y,        type = "Eigenvalue", cor = cor_flag) else target(S = S, type = "Eigenvalue", cor = cor_flag),
         "MSC"         = if (pass_Y) target(Y = Y,        type = "MSC",        cor = cor_flag) else target(S = S, type = "MSC",        cor = cor_flag),
         stop("Unknown target")
  )
}

# ---------- Simulation models (p=100 default) ----------
make_model <- function(model, p) {
  if (model == 1) {
    # Compound symmetry covariance (σij = 0.62 off-diag, 1 on diag); derive Θ = Σ^{-1}
    rho <- 0.62
    Sigma <- matrix(rho, p, p); diag(Sigma) <- 1
    Theta <- solve(Sigma)
    list(Sigma=Sigma, Theta=Theta)
  } else if (model %in% 2:4) {
    # Liu & Wang (2017): build adjacency A then Θ = D(A + (|λ_min(A)| + 0.2) I) D, D diag with 1 then 3
    if (model == 2) { # scale-free
      g <- sample_pa(p, power=1, m=1, directed=FALSE)
    } else if (model == 3) { # hub
      groups <- split(1:p, rep(1:(p/10), each=10))
      g <- make_empty_graph(p)
      for (grp in groups) {
        hub <- grp[1]
        for (v in grp[-1]) g <- add_edges(g, c(hub, v))
      }
    } else { # block graph model (model 4)
      stopifnot(p %% 10 == 0)
      block <- p %/% 10

      # Build one PD block: diag = 1, off-diag = 0.5
      B <- matrix(0.5, nrow = block, ncol = block)
      diag(B) <- 1

      # Assemble a dense block-diagonal Theta_tilde in base R (no Matrix::bdiag)
      Theta_tilde <- matrix(0, p, p)
      for (k in 0:9) {
          idx <- (k*block + 1):((k + 1)*block)
          Theta_tilde[idx, idx] <- B
      }

      # Random permutation preserves PD
      perm <- sample.int(p)
      Theta_tilde <- Theta_tilde[perm, perm, drop = FALSE]

      # Symmetrize + tiny PD nudge (shouldn't be needed, but safe)
      Theta_tilde <- (Theta_tilde + t(Theta_tilde)) / 2
      ev <- eigen(Theta_tilde, symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) <= 1e-12) {
          Theta_tilde <- Theta_tilde + (abs(min(ev)) + 1e-8) * diag(p)
      }

      # Σ = Θ^{-1}
      Sigma <- solve(Theta_tilde)

      # Heteroscedastic scaling: Σ <- D^{-1} Σ D^{-1}, but do it explicitly
      D  <- diag(c(rep(1, p/2), rep(1.5, p/2)))
      Di <- diag(1 / diag(D))
      Sigma <- Di %*% Sigma %*% Di

      # Clean up numerics
      Sigma <- (Sigma + t(Sigma)) / 2
      if (!all(is.finite(Sigma))) {
          bad <- which(!is.finite(Sigma), arr.ind = TRUE)
          stop(sprintf("Model 4: non-finite Sigma at [%d,%d]", bad[1,1], bad[1,2]))
      }
      evS <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
      if (min(evS) <= 1e-12) {
          Sigma <- Sigma + (abs(min(evS)) + 1e-8) * diag(p)
      }
      Theta <- solve(Sigma)

      return(list(Sigma=Sigma, Theta=Theta))
    }
    A <- as_adjacency_matrix(g, type="both", sparse=FALSE)
    A[upper.tri(A)] <- t(A)[upper.tri(A)]
    A[A!=0] <- 0.3; diag(A) <- 0 # threshold and zero diagonal
    lambda_min <- min(eigen(A, symmetric=TRUE, only.values=TRUE)$values)
    D <- diag(c(rep(1, p/2), rep(3, p/2)))
    Theta <- D %*% (A + (abs(lambda_min)+0.2)*diag(p)) %*% D
    Sigma <- solve(Theta)
    list(Sigma=Sigma, Theta=Theta)
  } else if (model == 5) {
    stopifnot(p >= 3)
    Theta_tilde <- matrix(0, p, p)
    for (i in 1:p) {
      Theta_tilde[i,i] <- 1
      if (i+1<=p) Theta_tilde[i,i+1] <- Theta_tilde[i+1,i] <- 0.6
      if (i+2<=p) Theta_tilde[i,i+2] <- Theta_tilde[i+2,i] <- 0.3
    }
    Theta_tilde <- .pd_nudge(Theta_tilde)
    D <- diag(runif(p, 1, 5))
    Sigma <- solve(Theta_tilde); Sigma <- .diag_inv(D) %*% Sigma %*% .diag_inv(D)
    Theta <- solve(Sigma)
    list(Sigma=Sigma, Theta=Theta)
  } else if (model == 6) {
    pmat <- matrix(0, p, p)
    for (i in 1:(p-1)) for (j in (i+1):p) {
      if (rbinom(1,1,0.05)==1) {
        u <- runif(1, 0.4, 0.8)
        pmat[i,j] <- pmat[j,i] <- u
      }
    }
    Theta_tilde <- .symm(pmat)
    ev <- eigen(Theta_tilde, symmetric=TRUE, only.values=TRUE)$values
    lambda_min <- min(ev)
    delta <- if (lambda_min <= 0) abs(lambda_min) + 0.05 else 0.05
    if (lambda_min <= 0) Theta_tilde <- Theta_tilde + delta * diag(p)
    Theta2 <- pmat
    lambda_min <- min(eigen(Theta2 + t(Theta2), symmetric=TRUE, only.values=TRUE)$values)
    Theta_tilde <- Theta2 + (abs(lambda_min)+0.05)*diag(p)
    D <- diag(runif(p, 1, 5))
    Sigma <- solve(Theta_tilde); Sigma <- .diag_inv(D) %*% Sigma %*% .diag_inv(D)
    Theta <- solve(Sigma)
    list(Sigma=Sigma, Theta=Theta)
  } else stop("Unknown model")
}
# ---------- Fitting wrappers ----------
fit_method <- function(S, Y=NULL, method=c("glasso","rope","gelnet"),
                       alpha=0.5, target_type="None", penalize_diag=TRUE, lambda) {
  method <- match.arg(method)
  # Use correlation input, as recommended
  if (!is.null(Y)) S <- cor(Y)
  Tmat <- NULL
  if (target_type!="None") Tmat <- diag_target(S, type=target_type, Y=Y, use_correlation=TRUE)
  if (method=="rope") alpha <- 0
  if (method=="glasso") alpha <- 1

  if (method=="rope") {
    # Closed form when alpha=0
    if (is.null(Tmat)) Tmat <- matrix(0, nrow(S), ncol(S))
    res <- rope(S=S, lambda=lambda, Target=Tmat)
    return(list(Theta=res))
  } else {
    res <- gelnet(S=S, lambda=lambda, alpha=alpha,
                  Target=Tmat, penalize.diag = penalize_diag)
    return(list(Theta=res$Theta))
  }
}

# ---------- Cross-validation over lambda ----------
cv_select_lambda <- function(Y, method, alpha, target_type, penalize_diag,
                             lambda_grid) {
  # 5-fold CV like the examples; use negative log-likelihood as score.
  n <- nrow(Y)
  folds <- sample(rep(1:5, length.out=n))
  scores <- numeric(length(lambda_grid))
  for (li in seq_along(lambda_grid)) {
    lam <- lambda_grid[li]
    loss_fold <- 0
    for (f in 1:5) {
      Ytr <- Y[folds!=f,,drop=FALSE]
      Yte <- Y[folds==f,,drop=FALSE]
      S_tr <- cor(Ytr)
      fit <- fit_method(S=S_tr, Y=NULL, method=method, alpha=alpha,
                        target_type=target_type, penalize_diag=penalize_diag, lambda=lam)
      Theta_hat <- fit$Theta
      # Evaluate test (Gaussian log-likelihood up to constant): tr(S_te Θ̂) - log det Θ̂
      S_te <- cor(Yte)
      loss_fold <- loss_fold + (sum(S_te * Theta_hat) - as.numeric(determinant(Theta_hat, logarithm=TRUE)$modulus))
    }
    scores[li] <- loss_fold/5
  }
  lambda_grid[ which.min(scores) ]
}

lambda_for_edges <- function(fit_fn, target_edges, grid) {
  last <- NULL
  for (lam in grid) {
    Th <- fit_fn(lam)
    e <- count_edges(Th)
    last <- list(lambda=lam, Theta=Th, edges=e)
    if (e <= target_edges) break
  }
  last
}
