group_lasso_multitask_multiple_kernel_survival_train <- function(Km, y, delta, parameters) {
  T <- length(Km)
  P <- dim(Km[[1]])[3]
  eta <- rep(1 / P, P)
  Keta <- vector("list", T)
  for (t in 1:T) {
    Keta[[t]] <- calculate_Keta(Km[[t]], eta)
  }

  models <- vector("list", T)
  for (t in 1:T) {
    models[[t]] <- solve_survival_svm(Keta[[t]], y[[t]], delta[[t]], parameters$C, parameters$tube, parameters$epsilon)
  }
  objectives <- sum(sapply(models, function(model) {model$objective}))
  print(c(length(objectives), sum(sapply(models, function(model) {model$objective}))))
  while(1) {
    start_objective <- sum(sapply(models, function(model) {model$objective}))
    start_eta <- eta

    eta_new <- rep(0, P)
    for (t in 1:T) {
      for (m in 1:P) {
        eta_new[m] <- eta_new[m] + eta[m] * sqrt(t(models[[t]]$alpha) %*% Km[[t]][,,m] %*% models[[t]]$alpha)
      }
    }
    eta <- eta_new
    eta <- eta / sum(eta)
    eta[eta < parameters$epsilon] <- 0
    eta <- eta / sum(eta)
    for (t in 1:T) {
      Keta[[t]] <- calculate_Keta(Km[[t]], eta)
    }

    for (t in 1:T) {
      models[[t]] <- solve_survival_svm(Keta[[t]], y[[t]], delta[[t]], parameters$C, parameters$tube, parameters$epsilon)
    }
    
    objectives <- c(objectives, sum(sapply(models, function(model) {model$objective})))
    print(c(length(objectives), sum(sapply(models, function(model) {model$objective}))))
    if (length(objectives) == parameters$iteration_count) {
      break
    }
  }

  state <- list(alpha = sapply(models, function(model) model$alpha, simplify = FALSE), b = sapply(models, function(model) model$b, simplify = FALSE), eta = eta, objectives = objectives, parameters = parameters)
}
