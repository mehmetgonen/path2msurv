group_lasso_multiple_kernel_survival_test <- function(Km, state) {
  Keta <- calculate_Keta(Km, state$eta)
  y <- Keta %*% state$alpha + state$b

  prediction <- list(y = y)
}
