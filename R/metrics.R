
suppressPackageStartupMessages({
  library(dplyr)
})

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


count_edges <- function(Theta, eps=1e-12) {
  A <- (abs(Theta) > eps) * 1
  diag(A) <- 0
  sum(A)/2
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

curve_points <- function(Theta_true, Theta_hat) {
  p <- nrow(Theta_true)
  A_true <- (abs(Theta_true) > 1e-8) * 1; diag(A_true) <- 0
  S <- abs(Theta_hat); diag(S) <- 0
  ths <- sort(unique(as.vector(S[upper.tri(S)])), decreasing = TRUE)
  if (length(ths) == 0) ths <- c(0)
  res <- vector("list", length(ths))
  up <- upper.tri(A_true)
  for (i in seq_along(ths)) {
    th <- ths[i]
    A_hat <- (S >= th) * 1
    TP <- sum((A_true==1 & A_hat==1)[up])
    TN <- sum((A_true==0 & A_hat==0)[up])
    FP <- sum((A_true==0 & A_hat==1)[up])
    FN <- sum((A_true==1 & A_hat==0)[up])
    prec <- ifelse((TP+FP)==0, 1, TP/(TP+FP))
    rec  <- ifelse((TP+FN)==0, 0, TP/(TP+FN))
    tpr  <- rec
    fpr  <- ifelse((FP+TN)==0, 0, FP/(FP+TN))
    res[[i]] <- data.frame(th=th, TP=TP, FP=FP, TN=TN, FN=FN,
                           precision=prec, recall=rec, tpr=tpr, fpr=fpr)
  }
  dplyr::bind_rows(res)
}
trapz_auc <- function(x, y) {
  if (length(x) < 2) return(0)
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  sum( (x[-1] - x[-length(x)]) * (y[-1] + y[-length(x)]) / 2 )
}
pr_auc <- function(df) { trapz_auc(df$recall, df$precision) }
roc_auc <- function(df) { trapz_auc(df$fpr, df$tpr) }
