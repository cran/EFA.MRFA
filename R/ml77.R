ml77 <- function(S,m,conv){
  ## Jennrich & Robinson Algorithm (Joreskog, 1977): simplified version

  n <- dim(S)[2]

  uni <- matrix(0,n,1)
  iter <- 1

  IS <- solve(S)

  Uo <- matrix(0,n,1)
  U <- matrix(0,n,1)

  T1 <- 4
  To <- T1 * 2

  while ((abs(To-T1)) > conv){

    if (iter > 1){

      U2o<-diag(c(Uo))
      t1<-matrix(0,n,1)
      for (i in 1:n){
        UoUo<-Uo[i]*Uo[i]
        if (UoUo> 10^(-5)){
          lUoUo <- log(UoUo)
        }
        else {
          lUoUo <- log(10^(-250))
        }
        t1[i]<-lUoUo
      }
    }
    else {
      t1 <- matrix(0,n,1)
      for (i in 1:n){
        if (IS[i,i]>0){
          t1[i] <- log((1-(m/(2*n))) / IS[i,i])
        }
        else {
          t1[i] <- log((1-(m/(2*n))) / 10^(-5))
        }
      }
      Uo <- matrix(0,n,1)
      for (i in 1:n){
        Uo[i] <- sqrt(exp(t1[i]))
      }
      U2o <- diag(c(Uo))
    }

    for (i in 1:n){
      if (is.infinite(U2o[i,i])){
        U2o[i,i]<-1
      }
      if (is.nan(U2o[i,i])){
        U2o[i,i]<-0
      }
    }

    #non-convergence
    if (max(t(max(is.nan(U2o * IS * U2o))))){
      T1 <- -1
      return(T1)
    }
    if (max(t(max(is.infinite(U2o * IS *U2o))))){
      T1 <- -1
      return(T1)
    }

    OUT <- eigenc(U2o %*% IS %*% U2o)
    B<-OUT$eigenvector_m
    L<-OUT$eigenvalue_m

    a <- matrix(0,n,1)

    for (i in 1:n){
      if (L[i]>0){
        a[i] <- 1 - (1/L[i])
      }
      else {
        a[i] <- 1 - (1/10^(-5))
      }
    }

    phi1 <- matrix(0,n,1)

    for (r in 1:n){
      for (i in (m+1):n){
        phi1[r] <- phi1[r] + (B[r,i]*B[r,i]) * a[i]
      }
    }

    PHI2 <- matrix(0,n,n)
    for (s in 1:n){
      for (r in 1:n){
        for (i in (m+1):n){
          PHI2[r,s] <- PHI2[r,s] + (B[r,i]*B[s,i])
        }
        PHI2[r,s] <- PHI2[r,s]^(2)
      }
    }

    nphi1 <-c()
    NPHI2 <- c()
    for (i in 1:n){
      if (uni[i] == 0){
        nphi1 <- rbind(nphi1,phi1[i])
        NPHI2 <- cbind(NPHI2, PHI2[,i])
      }
    }

    NNPHI2 <- c()

    for (i in 1:n){
      if (uni[i] == 0){
        NNPHI2 <- rbind(NNPHI2, NPHI2[i,])
      }
    }


    # If the matrix is unable to reach convergence, this function will crash
    d <- solve(NNPHI2) %*% nphi1

    nd <- matrix(0,n,1)
    j <- 1

    for (i in 1:n){
      if (uni[i] == 1){
        if (PHI2[i,i]>0){
          nd[i] <- phi1[i] / PHI2[i,i]
        }
        else {
          nd[i] <- phi1[i] / 10^(-5)
        }
      }
      else {
        nd[i] <- d[j]
        j <- j + 1
      }
    }

    t1 <- t1 - nd

    U <- matrix(0,n,1)

    for(i in 1:n){
      U[i] <- sqrt(exp(t1[i]))
    }

    Uo <- U

    for (i in 1:n){
      if (U[i] < .00000001){
        uni[i] <- 1
      }
      if (is.infinite(U[i])){
        U[i] <- 1
      }
      if (is.nan(U[i])){
        U[i] <- 0
      }
    }

    U <- diag(as.numeric(U))
    H <- (S - U*U)%*%solve(S)

    #non-convergence

    if (max(t(max(is.nan(H))))){
      T1 <- -1
      return(T1)
    }
    if (max(t(max(is.infinite(H))))){
      T1 <- -1
      return(T1)
    }

    OUT<-eigend(H)
    K<-OUT$eigenvector_m
    D<-diag(OUT$eigenvalue_m)

    d <- diag(as.matrix(D[1:m,1:m]))
    b <- K[,1:m]
    L <- c()

    if (m>0){
      if (m>1){
        for (i in 1:m){
          k <- t(b[,i]) %*% solve(S) %*% b[,i]
          u <- as.vector(sqrt(d[i]/k)) * b[,i]
          L <- cbind(L,u)
        }
      }
      else {
        k <- t(b) %*% solve(S) %*% b
        u <- as.vector(sqrt(d/k)) * b
        L <- cbind(L,u)
      }
    }

    if (m>0){
      Sr <- L%*%t(L) + (U*U)
    }
    else {
      Sr <- (U*U)
    }

    To <- T1
    T1 <- log(det(Sr)) + sum(diag(S%*%solve(Sr))) - log(det(S)) - n

    iter <- iter + 1

    if (iter ==20){
      To <- T1
    }


  } #fin WHILE

  OUT<-list('t'=T1,'A'=L)
  return(OUT)
}
