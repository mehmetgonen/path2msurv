library(randomForestSRC)
source("survival_helper.R")

cohorts <- c("TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-COAD", "TCGA-ESCA",
             "TCGA-GBM",  "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", 
             "TCGA-LGG",  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-OV",
             "TCGA-PAAD", "TCGA-READ", "TCGA-SARC", "TCGA-STAD", "TCGA-UCEC")

data_path <- "./data"
result_path <- "./results"
pathway <- "none" #replace with "hallmark" or "pid" if you would like to use genes in Hallmark gene sets or PID pathways only

for (cohort in cohorts) {
  if (dir.exists(sprintf("%s/%s", result_path, cohort)) == FALSE) {
    dir.create(sprintf("%s/%s", result_path, cohort)) 
  }
  
  for (replication in 1:100) {
    if (file.exists(sprintf("%s/%s/random_forest_pathway_%s_measure_CI_replication_%d_result.RData", result_path, cohort, pathway, replication)) == FALSE) {
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
      
      ntree_set <- c(1:5) * 500
      fold_count <- 4
      train_ratio <- 0.8
      
      if (pathway != "none") {
        pathways <- read_pathways(pathway)
        gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
        X <- X[, which(colnames(X) %in% gene_names)]
      }
      
      set.seed(1606 * replication)
      train_dead_indices <- sample(dead_indices, ceiling(train_ratio * length(dead_indices)))
      train_alive_indices <- sample(alive_indices, ceiling(train_ratio * length(alive_indices)))
      
      concordance_index_matrix <- matrix(NA, nrow = fold_count, ncol = length(ntree_set), dimnames = list(1:fold_count, sprintf("%g", ntree_set)))
      
      dead_allocation <- sample(rep(1:fold_count, ceiling(length(train_dead_indices) / fold_count)), length(train_dead_indices))
      alive_allocation <- sample(rep(1:fold_count, ceiling(length(train_alive_indices) / fold_count)), length(train_alive_indices))
      for (fold in 1:fold_count) {
        train_indices <- c(train_dead_indices[which(dead_allocation != fold)], train_alive_indices[which(alive_allocation != fold)])
        test_indices <- c(train_dead_indices[which(dead_allocation == fold)], train_alive_indices[which(alive_allocation == fold)])
        
        X_train <- X[train_indices,]
        X_test <- X[test_indices,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
        X_train[is.na(X_train)] <- 0
        X_test[is.na(X_test)] <- 0
        
        y_train <- y[train_indices]
        delta_train <- delta[train_indices]
        y_test <- y[test_indices]
        delta_test <- delta[test_indices]
        
        for (ntree in ntree_set) {
          state <- rfsrc(formula = as.formula(Surv(time, event) ~ .), data = as.data.frame(cbind("time" = y_train, "event" = 1 - delta_train, X_train)), ntree = ntree)
          prediction <- predict(state, newdata = as.data.frame(cbind("time" = y_test, "event" = 1 - delta_test, X_test)))
          concordance_index_matrix[fold, sprintf("%g", ntree)] <- 1 - prediction$err.rate[ntree]
        }
      }
      
      ntree_star_CI <- ntree_set[max.col(t(colMeans(concordance_index_matrix, na.rm = TRUE)), ties.method = "last")]
      
      train_indices <- c(train_dead_indices, train_alive_indices)
      test_indices <- setdiff(1:length(delta), train_indices)
      
      X_train <- X[train_indices,]
      X_test <- X[test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      X_train[is.na(X_train)] <- 0
      X_test[is.na(X_test)] <- 0
      
      y_train <- y[train_indices]
      delta_train <- delta[train_indices]
      y_test <- y[test_indices]
      delta_test <- delta[test_indices]
      
      state <- rfsrc(formula = as.formula(Surv(time, event) ~ .), data = as.data.frame(cbind("time" = y_train, "event" = 1 - delta_train, X_train)), ntree = ntree_star_CI)
      prediction <- predict(state, newdata = as.data.frame(cbind("time" = y_test, "event" = 1 - delta_test, X_test)))
      result <- list()
      result$ntree <- ntree_star_CI
      result$y_test <- y_test
      result$delta_test <- delta_test
      result$y_predicted <- prediction$predicted
      result$CI <- 1 - prediction$err.rate[ntree_star_CI]
      
      save("result", file = sprintf("%s/%s/random_forest_pathway_%s_measure_CI_replication_%d_result.RData", result_path, cohort, pathway, replication))
    }
  }
}
