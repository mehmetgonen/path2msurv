library(Rcplex)

solve_survival_svm <- function(K, y, delta, C, tube, epsilon) {
  N <- length(y)

  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200

  result <- Rcplex(cvec = c(tube - y, tube + y), 
                   Amat = matrix(c(rep(1, N), rep(-1, N)), nrow = 1, byrow = TRUE), 
                   bvec = 0, 
                   Qmat = rbind(cbind(+K, -K), cbind(-K, +K)), 
                   lb = c(rep(0, N), rep(0, N)),
                   ub = c(rep(C, N), C * (1 - delta)),
                   control = opts,
                   objsense = "min",
                   sense = "E")

  alpha <- result$xopt[1:N] - result$xopt[(N + 1):(2 * N)]
  alpha[alpha > 0 & alpha < +C * epsilon] <- 0
  alpha[alpha < 0 & alpha > -C * epsilon] <- 0
  alpha[alpha > 0 & alpha > +C * (1 - epsilon)] <- +C
  alpha[alpha < 0 & alpha < -C * (1 - epsilon)] <- -C
  objective <- t(y) %*% alpha - tube * sum(abs(alpha)) - 0.5 * (t(alpha) %*% K) %*% alpha
  objective <- objective * (objective >= 0)
  
  support_indices <- which(alpha != 0)
  active_indices <- which(alpha != 0 & abs(alpha) < C)
  if (length(active_indices) > 0) {
    b <- mean(y[active_indices] - tube * sign(alpha[active_indices])) - mean(K[active_indices, support_indices] %*% alpha[support_indices])
  } else {
    b <- 0
  }

  model <- list(alpha = alpha, b = b, objective = objective)
}
