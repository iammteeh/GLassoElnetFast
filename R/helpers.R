
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
    } else { # block graph
      block <- p/10
      Theta_tilde <- bdiag(replicate(10, {
        B <- matrix(0.5, block, block); diag(B) <- 1; B
      }, simplify=FALSE))
      perm <- sample(1:p); Theta_tilde <- as.matrix(Theta_tilde)[perm, perm]
      D <- diag(c(rep(1, p/2), rep(1.5, p/2)))
      Sigma <- solve(Theta_tilde); Sigma <- D^{-1} %*% Sigma %*% D^{-1}
      Theta <- solve(Sigma)
      return(list(Sigma=Sigma, Theta=Theta))
    }
    A <- as_adjacency_matrix(g, type="both", sparse=FALSE)
    A[upper.tri(A)] <- t(A)[upper.tri(A)]
    A[A!=0] <- 0.3; diag(A) <- 0
    lambda_min <- min(eigen(A, symmetric=TRUE, only.values=TRUE)$values)
    D <- diag(c(rep(1, p/2), rep(3, p/2)))
    Theta <- D %*% (A + (abs(lambda_min)+0.2)*diag(p)) %*% D
    Sigma <- solve(Theta)
    list(Sigma=Sigma, Theta=Theta)
  } else if (model == 5) {
    Theta_tilde <- matrix(0, p, p)
    for (i in 1:p) {
      Theta_tilde[i,i] <- 1
      if (i+1<=p) Theta_tilde[i,i+1] <- Theta_tilde[i+1,i] <- 0.6
      if (i+2<=p) Theta_tilde[i,i+2] <- Theta_tilde[i+2,i] <- 0.3
    }
    D <- diag(runif(p, 1, 5))
    Sigma <- solve(Theta_tilde); Sigma <- D^{-1} %*% Sigma %*% D^{-1}
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
    Theta2 <- pmat
    lambda_min <- min(eigen(Theta2 + t(Theta2), symmetric=TRUE, only.values=TRUE)$values)
    Theta_tilde <- Theta2 + (abs(lambda_min)+0.05)*diag(p)
    D <- diag(runif(p, 1, 5))
    Sigma <- solve(Theta_tilde); Sigma <- D^{-1} %*% Sigma %*% D^{-1}
    Theta <- solve(Sigma)
    list(Sigma=Sigma, Theta=Theta)
  } else stop("Unknown model")
}

# ---------- Metrics ----------
eps_adj <- function(M, eps=1e-8) {
  # Avoid numerical indefiniteness
  M <- (M + t(M))/2
  ev <- eigen(M, symmetric=TRUE)$values
  if (min(ev) <= eps) {
    M <- M + (abs(min(ev)) + eps) * diag(nrow(M))
  }
  M
}

kl_loss <- function(Sigma, Theta_hat) {
  p <- nrow(Sigma)
  val <- sum(Sigma * Theta_hat) - determinant(Sigma %*% solve(Sigma))$modulus # trick not needed; use ΣΘ̂
  val <- sum(Sigma * Theta_hat) - log(det(Sigma %*% solve(solve(Theta_hat)))) - p # but simpler:
  # safer:
  A <- Sigma %*% Theta_hat
  as.numeric(sum(diag(A)) - log(det(A)) - p)
}

kl_loss_safe <- function(Sigma, Theta_hat) {
  p <- nrow(Sigma)
  A <- Sigma %*% Theta_hat
  A <- eps_adj(A)
  as.numeric(sum(diag(A)) - log(det(A)) - p)
}

l2_loss <- function(Theta, Theta_hat) sqrt(sum((Theta - Theta_hat)^2))

sp_loss <- function(Theta, Theta_hat) {
  E <- Theta - Theta_hat
  sqrt(max(eigen(t(E)%*%E, symmetric=TRUE, only.values=TRUE)$values))
}

graph_scores <- function(Theta, Theta_hat, eps=1e-5) {
  A_true <- (abs(Theta) >= eps) * 1; diag(A_true) <- 0
  A_hat  <- (abs(Theta_hat) >= eps) * 1; diag(A_hat) <- 0
  upper <- upper.tri(A_true)
  TP <- sum((A_true==1 & A_hat==1)[upper])
  TN <- sum((A_true==0 & A_hat==0)[upper])
  FP <- sum((A_true==0 & A_hat==1)[upper])
  FN <- sum((A_true==1 & A_hat==0)[upper])
  F1 <- ifelse((2*TP+FN+FP)==0, 0, 2*TP/(2*TP+FN+FP))
  denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC <- ifelse(denom==0, 0, (TP*TN - FP*FN)/denom)
  c(F1=F1, MCC=MCC, TP=TP, TN=TN, FP=FP, FN=FN)
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
