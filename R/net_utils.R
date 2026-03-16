prepare_net_graph <- function(Omega, sgn, add_intercept = FALSE) {
  sgn_net <- as.integer(sgn)

  if (add_intercept) {
    sgn_net <- c(1L, sgn_net)
    if (inherits(Omega, "sparseMatrix")) {
      Omega_net <- Matrix::bdiag(Matrix::Matrix(0, 1, 1, sparse = TRUE), Omega)
    } else {
      Omega_net <- rbind(0, cbind(0, Omega))
    }
  } else {
    Omega_net <- Omega
  }

  if (inherits(Omega_net, "sparseMatrix")) {
    Matrix::diag(Omega_net) <- 0
    if (length(Omega_net@x) > 0) {
      Omega_net@x <- abs(Omega_net@x)
    }
    W <- OmegaSC(Omega_net, sgn_net)
  } else {
    diag(Omega_net) <- 0
    Omega_net <- abs(Omega_net)
    W <- OmegaC(Omega_net, sgn_net)
  }

  W$loc <- W$loc + 1
  W
}
