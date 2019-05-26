source("survival_helper.R")
source("solve_survival_svm_cplex.R")
source("group_lasso_multiple_kernel_survival_train.R")
source("group_lasso_multiple_kernel_survival_test.R")

cohorts <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-COAD", "TCGA-ESCA",
             "TCGA-GBM",  "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", 
             "TCGA-LGG",  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-OV",
             "TCGA-PAAD", "TCGA-READ", "TCGA-SARC", "TCGA-STAD", "TCGA-UCEC")
  
data_path <- "./data"
result_path <- "./results"
pathway <- "hallmark" #replace with "pid" if you would like to use PID pathways

for (cohort in cohorts) {
  if (dir.exists(sprintf("%s/%s", result_path, cohort)) == FALSE) {
    dir.create(sprintf("%s/%s", result_path, cohort)) 
  }
  
  for (replication in 1:100) {
    if (file.exists(sprintf("%s/%s/glmkl_pathway_%s_measure_CI_replication_%d_result.RData", result_path, cohort, pathway, replication)) == FALSE) {
      load(sprintf("%s/%s.RData", data_path, cohort))
      
      patients_with_mrna <- rownames(TCGA$mrna)
      patients_with_survival <- rownames(TCGA$clinical)[which(is.na(TCGA$clinical$vital_status) == FALSE & ((TCGA$clinical$vital_status == "Dead" & is.na(TCGA$clinical$days_to_death) == FALSE & TCGA$clinical$days_to_death > 0) | (TCGA$clinical$vital_status == "Alive" & ((is.na(TCGA$clinical$days_to_last_followup) == FALSE & TCGA$clinical$days_to_last_followup > 0) | (is.na(TCGA$clinical$days_to_last_known_alive) == FALSE & TCGA$clinical$days_to_last_known_alive > 0)))))]
      common_patients <- intersect(patients_with_mrna, patients_with_survival)
      
      X <- log2(TCGA$mrna[common_patients,] + 1)
      Y <- TCGA$clinical[common_patients, c("vital_status", "days_to_death", "days_to_last_followup", "days_to_last_known_alive")]
      
      delta <- 1 * (Y$vital_status == "Alive")
      y <- matrix(NA, length(delta), 1)
      y[delta == 0] <- Y$days_to_death[delta == 0]
      y[delta == 1] <- sapply(which(delta == 1), function(patient) {max(Y$days_to_last_followup[patient], Y$days_to_last_known_alive[patient], na.rm = TRUE)})
      
      valid_features <- as.numeric(which(apply(X, 2, sd) != 0))
      X <- X[,valid_features]
      
      dead_indices <- which(delta == 0)
      alive_indices <- which(delta == 1)
      
      C_set <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000)
      epsilon <- 1e-3
      fold_count <- 4
      train_ratio <- 0.8
      iteration_count <- 200
      tube <- 0
      
      pathways <- read_pathways(pathway)
      gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
      X <- X[, which(colnames(X) %in% gene_names)]
      
      set.seed(1606 * replication)
      train_dead_indices <- sample(dead_indices, ceiling(train_ratio * length(dead_indices)))
      train_alive_indices <- sample(alive_indices, ceiling(train_ratio * length(alive_indices)))
      
      concordance_index_matrix <- matrix(NA, nrow = fold_count, ncol = length(C_set), dimnames = list(1:fold_count, sprintf("%g", C_set)))
      
      dead_allocation <- sample(rep(1:fold_count, ceiling(length(train_dead_indices) / fold_count)), length(train_dead_indices))
      alive_allocation <- sample(rep(1:fold_count, ceiling(length(train_alive_indices) / fold_count)), length(train_alive_indices))
      for (fold in 1:fold_count) {
        train_indices <- c(train_dead_indices[which(dead_allocation != fold)], train_alive_indices[which(alive_allocation != fold)])
        test_indices <- c(train_dead_indices[which(dead_allocation == fold)], train_alive_indices[which(alive_allocation == fold)])
        
        X_train <- X[train_indices,]
        X_test <- X[test_indices,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
        
        N_train <- nrow(X_train)
        N_test <- nrow(X_test)
        N_pathway <- length(pathways)
        K_train <- array(0, dim = c(N_train, N_train, N_pathway))
        K_test <- array(0, dim = c(N_test, N_train, N_pathway))
        for (m in 1:N_pathway) {
          feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
          D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
          D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
          sigma <- mean(D_train)
          K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
          K_test[,,m] <- exp(-D_test^2 / (2 * sigma^2))
        }
        
        y_train <- y[train_indices]
        delta_train <- delta[train_indices]
        y_test <- y[test_indices]
        delta_test <- delta[test_indices]
        
        for (C in C_set) {
          print(sprintf("running fold = %d, C = %g", fold, C))
          parameters <- list()
          parameters$C <- C
          parameters$epsilon <- epsilon
          parameters$tube <- tube
          parameters$iteration_count <- iteration_count
          
          state <- group_lasso_multiple_kernel_survival_train(K_train, y_train, delta_train, parameters)
          prediction <- group_lasso_multiple_kernel_survival_test(K_test, state)
          concordance_index_matrix[fold, sprintf("%g", C)] <- concordance_index(prediction$y, y_test, delta_test)
        }
      }
      
      C_star_CI <- C_set[max.col(t(colMeans(concordance_index_matrix)), ties.method = "last")]    
      
      train_indices <- c(train_dead_indices, train_alive_indices)
      test_indices <- setdiff(1:length(delta), train_indices)
      
      X_train <- X[train_indices,]
      X_test <- X[test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      
      N_train <- nrow(X_train)
      N_test <- nrow(X_test)
      N_pathway <- length(pathways)
      K_train <- array(0, dim = c(N_train, N_train, N_pathway))
      K_test <- array(0, dim = c(N_test, N_train, N_pathway))
      for (m in 1:N_pathway) {
        feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
        D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
        D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
        sigma <- mean(D_train)
        K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
        K_test[,,m] <- exp(-D_test^2 / (2 * sigma^2))
      }
      
      y_train <- y[train_indices]
      delta_train <- delta[train_indices]
      y_test <- y[test_indices]
      delta_test <- delta[test_indices]
      
      parameters <- list()
      parameters$C <- C_star_CI
      parameters$epsilon <- epsilon
      parameters$tube <- tube
      parameters$iteration_count <- iteration_count
      
      state <- group_lasso_multiple_kernel_survival_train(K_train, y_train, delta_train, parameters)
      prediction <- group_lasso_multiple_kernel_survival_test(K_test, state)
      result <- list()
      result$C <- C_star_CI
      result$y_test <- y_test
      result$delta_test <- delta_test
      result$y_predicted <- prediction$y
      result$CI <- concordance_index(prediction$y, y_test, delta_test)
      
      save("state", file = sprintf("%s/%s/glmkl_pathway_%s_measure_CI_replication_%d_state.RData", result_path, cohort, pathway, replication))
      save("result", file = sprintf("%s/%s/glmkl_pathway_%s_measure_CI_replication_%d_result.RData", result_path, cohort, pathway, replication))
    }
  }
}
