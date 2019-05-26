single_kernel_survival_test <- function(K, state) {
  y <- K %*% state$alpha + state$b
  
  prediction <- list(y = y)
}
