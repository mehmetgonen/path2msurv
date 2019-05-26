concordance_index <- function(y_predicted, y, delta) {
  N <- length(y)
  Delta <- (matrix(delta == 0, N, N, byrow = FALSE) & matrix(delta == 0, N, N, byrow = TRUE)) | 
           (matrix(delta == 0, N, N, byrow = FALSE) & matrix(delta == 1, N, N, byrow = TRUE) & (matrix(y, N, N, byrow = FALSE) - matrix(y, N, N, byrow = TRUE) < 0))
  Diff <- (matrix(y_predicted, N, N, byrow = FALSE) - matrix(y_predicted, N, N, byrow = TRUE)) * (matrix(y, N, N, byrow = FALSE) - matrix(y, N, N, byrow = TRUE)) > 0
  CI <- (sum(Delta * Diff) - sum(diag(Delta * Diff))) / (sum(Delta) - sum(diag(Delta)))
  if (is.na(CI) == TRUE) {
    CI <- 0
  }
  return(CI)
}

pdist <- function(X1, X2) {
  if (identical(X1, X2) == TRUE) {
    D <- as.matrix(dist(X1))
  }
  else {
    D <- as.matrix(dist(rbind(X1, X2)))
    D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
  }
  return(D)
}

read_pathways <- function(name) {
  symbols_lines <- read.table(sprintf("msigdb/%s.gmt", name), header = FALSE, sep = ",", stringsAsFactor = FALSE)
  pathways <- vector("list", nrow(symbols_lines))
  for (line in 1:nrow(symbols_lines)) {
    symbols_entries <- strsplit(symbols_lines[line, 1], "\t")
    pathways[[line]]$name <- symbols_entries[[1]][1]
    pathways[[line]]$link <- symbols_entries[[1]][2]
    pathways[[line]]$symbols <- sort(symbols_entries[[1]][-2:-1])
  }
  return(pathways)
}

calculate_Keta <- function(Km, eta) {
  P <- dim(Km)[3]
  
  Keta <- eta[1] * Km[,,1]
  for (m in 2:P) {
    Keta <- Keta + eta[m] * Km[,,m]
  }
  
  return(Keta)
}
