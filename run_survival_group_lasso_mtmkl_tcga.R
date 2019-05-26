source("survival_helper.R")
source("solve_survival_svm_cplex.R")
source("group_lasso_multitask_multiple_kernel_survival_train.R")
source("group_lasso_multitask_multiple_kernel_survival_test.R")

cohorts <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-COAD", "TCGA-ESCA",
             "TCGA-GBM",  "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", 
             "TCGA-LGG",  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-OV",
             "TCGA-PAAD", "TCGA-READ", "TCGA-SARC", "TCGA-STAD", "TCGA-UCEC")
  
data_path <- "./data"
result_path <- "./results"
pathway <- "hallmark" #replace with "pid" if you would like to use PID pathways

full_name <- paste0("TCGA-", paste0(gsub(cohorts, pattern = "TCGA-", replacement = ""), collapse = "-"))
T <- length(cohorts)

if (dir.exists(sprintf("%s/%s", result_path, full_name)) == FALSE) {
  dir.create(sprintf("%s/%s", result_path, full_name)) 
}

for (replication in 1:100) {
  if (file.exists(sprintf("%s/%s/glmtmkl_pathway_%s_measure_CI_replication_%d_result.RData", result_path, full_name, pathway, replication)) == FALSE) {
    X <- vector("list", T)
    y <- vector("list", T)
    delta <- vector("list", T)
    dead_indices <- vector("list", T)
    alive_indices <- vector("list", T)
    for (t in 1:T) {  
      load(sprintf("%s/%s.RData", data_path, cohorts[t]))
    
      patients_with_mrna <- rownames(TCGA$mrna)
      patients_with_survival <- rownames(TCGA$clinical)[which(is.na(TCGA$clinical$vital_status) == FALSE & ((TCGA$clinical$vital_status == "Dead" & is.na(TCGA$clinical$days_to_death) == FALSE & TCGA$clinical$days_to_death > 0) | (TCGA$clinical$vital_status == "Alive" & ((is.na(TCGA$clinical$days_to_last_followup) == FALSE & TCGA$clinical$days_to_last_followup > 0) | (is.na(TCGA$clinical$days_to_last_known_alive) == FALSE & TCGA$clinical$days_to_last_known_alive > 0)))))]
      common_patients <- intersect(patients_with_mrna, patients_with_survival)
    
      X[[t]] <- log2(TCGA$mrna[common_patients,] + 1)
      Y <- TCGA$clinical[common_patients, c("vital_status", "days_to_death", "days_to_last_followup", "days_to_last_known_alive")]
      
      delta[[t]] <- 1 * (Y$vital_status == "Alive")
      y[[t]] <- matrix(NA, length(delta[[t]]), 1)
      y[[t]][delta[[t]] == 0] <- Y$days_to_death[delta[[t]] == 0]
      y[[t]][delta[[t]] == 1] <- sapply(which(delta[[t]] == 1), function(patient) {max(Y$days_to_last_followup[patient], Y$days_to_last_known_alive[patient], na.rm = TRUE)})
      
      valid_features <- as.numeric(which(apply(X[[t]], 2, sd) != 0))
      X[[t]] <- X[[t]][,valid_features]
      
      dead_indices[[t]] <- which(delta[[t]] == 0)
      alive_indices[[t]] <- which(delta[[t]] == 1)
    }
  
    C_set <- c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000)
    epsilon <- 1e-3
    fold_count <- 4
    train_ratio <- 0.8
    iteration_count <- 200
    tube <- 0
    
    pathways <- read_pathways(pathway)
    gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
    for (t in 1:T) {
      X[[t]] <- X[[t]][, which(colnames(X[[t]]) %in% gene_names)] 
    }
  
    concordance_index_matrix <- array(NA, c(fold_count, length(C_set), T), dimnames = list(1:fold_count, sprintf("%g", C_set), 1:T))
    
    train_dead_indices <- vector("list", T)
    train_alive_indices <- vector("list", T)
    dead_allocation <- vector("list", T)
    alive_allocation <- vector("list", T)
    for (t in 1:T) {
      set.seed(1606 * replication)
      train_dead_indices[[t]] <- sample(dead_indices[[t]], ceiling(train_ratio * length(dead_indices[[t]])))
      train_alive_indices[[t]] <- sample(alive_indices[[t]], ceiling(train_ratio * length(alive_indices[[t]])))
      
      dead_allocation[[t]] <- sample(rep(1:fold_count, ceiling(length(train_dead_indices[[t]]) / fold_count)), length(train_dead_indices[[t]]))
      alive_allocation[[t]] <- sample(rep(1:fold_count, ceiling(length(train_alive_indices[[t]]) / fold_count)), length(train_alive_indices[[t]]))
    }
  
    K_train <- vector("list", T)
    K_test <- vector("list", T)
    y_train <- vector("list", T)
    y_test <- vector("list", T)
    delta_train <- vector("list", T)
    delta_test <- vector("list", T)
    for (fold in 1:fold_count) {
      for (t in 1:T) {
        train_indices <- c(train_dead_indices[[t]][which(dead_allocation[[t]] != fold)], train_alive_indices[[t]][which(alive_allocation[[t]] != fold)])
        test_indices <- c(train_dead_indices[[t]][which(dead_allocation[[t]] == fold)], train_alive_indices[[t]][which(alive_allocation[[t]] == fold)])
      
        X_train <- X[[t]][train_indices,]
        X_test <- X[[t]][test_indices,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      
        N_train <- nrow(X_train)
        N_test <- nrow(X_test)
        N_pathway <- length(pathways)
        K_train[[t]] <- array(0, dim = c(N_train, N_train, N_pathway))
        K_test[[t]] <- array(0, dim = c(N_test, N_train, N_pathway))
        for (m in 1:N_pathway) {
          feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
          D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
          D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
          sigma <- mean(D_train)
          K_train[[t]][,,m] <- exp(-D_train^2 / (2 * sigma^2))
          K_test[[t]][,,m] <- exp(-D_test^2 / (2 * sigma^2))
        }
      
        y_train[[t]] <- y[[t]][train_indices]
        delta_train[[t]] <- delta[[t]][train_indices]
        y_test[[t]] <- y[[t]][test_indices]
        delta_test[[t]] <- delta[[t]][test_indices]
      }
      
      for (C in C_set) {
        print(sprintf("running fold = %d, C = %g", fold, C))
        parameters <- list()
        parameters$C <- C
        parameters$epsilon <- epsilon
        parameters$tube <- tube
        parameters$iteration_count <- iteration_count
        
        state <- group_lasso_multitask_multiple_kernel_survival_train(K_train, y_train, delta_train, parameters)
        prediction <- group_lasso_multitask_multiple_kernel_survival_test(K_test, state)
        for (t in 1:T) {
          concordance_index_matrix[fold, sprintf("%g", C), t] <- concordance_index(prediction$y[[t]], y_test[[t]], delta_test[[t]])
        }
      }
    }
  
    C_star_CI <- C_set[max.col(t(colMeans(apply(concordance_index_matrix, c(1, 2), mean))), ties.method = "last")]    
    
    for (t in 1:T) {
      train_indices <- c(train_dead_indices[[t]], train_alive_indices[[t]])
      test_indices <- setdiff(1:length(delta[[t]]), train_indices)
    
      X_train <- X[[t]][train_indices,]
      X_test <- X[[t]][test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      
      N_train <- nrow(X_train)
      N_test <- nrow(X_test)
      N_pathway <- length(pathways)
      K_train[[t]] <- array(0, dim = c(N_train, N_train, N_pathway))
      K_test[[t]] <- array(0, dim = c(N_test, N_train, N_pathway))
      for (m in 1:N_pathway) {
        feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
        D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
        D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
        sigma <- mean(D_train)
        K_train[[t]][,,m] <- exp(-D_train^2 / (2 * sigma^2))
        K_test[[t]][,,m] <- exp(-D_test^2 / (2 * sigma^2))
      }
      
      y_train[[t]] <- y[[t]][train_indices]
      delta_train[[t]] <- delta[[t]][train_indices]
      y_test[[t]] <- y[[t]][test_indices]
      delta_test[[t]] <- delta[[t]][test_indices]
    }
  
    parameters <- list()
    parameters$C <- C_star_CI
    parameters$epsilon <- epsilon
    parameters$tube <- tube
    parameters$iteration_count <- iteration_count
  
    state <- group_lasso_multitask_multiple_kernel_survival_train(K_train, y_train, delta_train, parameters)
    prediction <- group_lasso_multitask_multiple_kernel_survival_test(K_test, state)
    result <- list()
    result$C <- C_star_CI
    result$y_test <- y_test
    result$delta_test <- delta_test
    result$y_predicted <- prediction$y
    result$CI <- rep(0, T)
    for (t in 1:T) {
      result$CI[t] <- concordance_index(prediction$y[[t]], y_test[[t]], delta_test[[t]])
    }
    
    save("state", file = sprintf("%s/%s/glmtmkl_pathway_%s_measure_CI_replication_%d_state.RData", result_path, full_name, pathway, replication))
    save("result", file = sprintf("%s/%s/glmtmkl_pathway_%s_measure_CI_replication_%d_result.RData", result_path, full_name, pathway, replication))
  }
}
