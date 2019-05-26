single_kernel_survival_train <- function(K, y, delta, parameters) {
  model <- solve_survival_svm(K, y, delta, parameters$C, parameters$tube, parameters$epsilon)
  
  state <- list(alpha = model$alpha, b = model$b, parameters = parameters)
}
