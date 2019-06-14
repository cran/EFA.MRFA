eigenc <- function(S) {
  ## Compute the eigenvalue descomposition
  #
  # S <- simetric matrix

  n<-dim(S)[2]
  B<-eigen(S)$vectors
  D<-eigen(S)$values

  i <- 1

  while (i<=(n-1)){
    if (D[i] > D[i+1]){

      tmp <- D[i]
      D[i] <- D[i+1]
      D[i+1] <- tmp
      tmp <- B[,i]

      for (j in 1:n){
        B[j,i] <- B[j,i+1]
      }
      for (j in 1:n){
        B[j,i+1] <- tmp[j]
      }
      i<-0
    }
    i <- i+1
  }

  OUT<-list('eigenvector_m'=B,'eigenvalue_m'=D)
  return(OUT)
}
