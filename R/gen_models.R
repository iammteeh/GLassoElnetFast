make_ba_conditional_hamiltonian <- function(p = 100, m_attach = 2, seed = 42,
                                            add_budget = 0.02) {
  set.seed(seed)
  g0 <- igraph::sample_pa(p, power = 1, m = m_attach, directed = FALSE)
  # Communities (proxy for invariant cycle subgroups)
  comm <- igraph::cluster_louvain(g0)
  memb <- igraph::membership(comm)
  groups <- split(1:p, memb)

  # Representatives
  reps <- vapply(groups, function(idx) idx[which.max(igraph::degree(g0, idx))], 1)

  # Connector path between reps
  Eplus <- list()
  for (k in seq_len(length(reps)-1)) {
    Eplus[[k]] <- c(reps[k], reps[k+1])
  }
  g1 <- igraph::add_edges(g0, unlist(Eplus))

# Intra-community long paths (DFS heuristic)
add_intra <- list()
for (idx in groups) {
  # Global DFS restricted to this community
  ord <- igraph::dfs(g1, root = idx[1], mode = "all", unreachable = FALSE)$order
  ord <- ord[ord %in% idx]
  if (length(ord) >= 2) {
    path_edges <- as.vector(rbind(ord[-length(ord)], ord[-1]))
    add_intra[[length(add_intra)+1]] <- path_edges
  }
}
g2 <- igraph::add_edges(g1, unlist(add_intra))

  # Edge budget to boost algebraic connectivity
  B <- ceiling(add_budget * igraph::gsize(g2))
  if (B > 0) {
    for (b in 1:B) {
      uv <- sample(1:p, 2, replace = FALSE)
      if (!igraph::are_adjacent(g2, uv[1], uv[2])) {
        g2 <- igraph::add_edges(g2, uv)
      }
    }
  }

  # Build precision templates (normalized Laplacians)
  A <- as.matrix(igraph::as_adjacency_matrix(g2, sparse = FALSE))
  L0 <- diag(rowSums(A)) - A

  # Community path Laplacians
  B_list <- list(B0 = L0)
  for (gidx in seq_along(groups)) {
    idx <- groups[[gidx]]
    # Path within idx in g2
    subA <- A[idx, idx, drop = FALSE]
    # simple path Laplacian: chain over idx order
    P <- matrix(0, length(idx), length(idx))
    for (i in 1:(length(idx)-1)) {
      P[i,i] <- P[i,i] + 1; P[i+1,i+1] <- P[i+1,i+1] + 1
      P[i,i+1] <- P[i,i+1] - 1; P[i+1,i] <- P[i+1,i] - 1
    }
    Bi <- matrix(0, p, p); Bi[idx, idx] <- P
    B_list[[paste0("B", gidx)]] <- Bi
  }
  B_list
}

sample_gsm <- function(B_list, n = 200, beta = 0.1,
                       lambda0 = 0.5, a = 3, b = 3, seed = 1) {
  set.seed(seed)
  p <- nrow(B_list[[1]])
  # Gamma scales per community basis
  k <- length(B_list) - 1
  tau <- rgamma(k, shape = a, rate = b)
  Theta <- beta * diag(p) + lambda0 * B_list$B0
  for (i in 1:k) Theta <- Theta + tau[i] * B_list[[paste0("B", i)]]
  # Make SPD safe
  Theta <- (Theta + t(Theta))/2
  eig <- eigen(Theta, symmetric = TRUE)$values
  if (min(eig) <= 1e-6) Theta <- Theta + (abs(min(eig)) + 1e-6) * diag(p)
  Sigma <- solve(Theta)
  Y <- mvtnorm::rmvnorm(n, sigma = Sigma)
  list(Y = Y, Sigma = Sigma, Theta = Theta, tau = tau)
}