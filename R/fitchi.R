fitchi<-function(R,A,N,t){

  p<-dim(A)[1]
  k<-dim(A)[2]

  if (p==1){
    p<-size(R)[2]
    k<-0
    Re<-R
  }
  else {
    Re <- (A%*%t(A)) + diag(p) - diag(diag(A%*%t(A)))
  }

  F2 <- 0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      F2 <- F2 + R[i,j]^2
    }
  }

  # CFI

  Chio <- F2*(N-1)
  Do <- 1/2 * p * (p-1)
  Chik <- (N-1) * t
  Dk <- round(1/2 * ((p-k) * (p-k+1)) - p)

  CFI <- abs(((Chio - Do) - (Chik - Dk)) / (Chio - Do))

  #RMSEA

  t2 <- log(det(Re)) + sum(diag(R%*%solve(Re))) - log(det(R)) - p
  Chik2 <- (N-1) * t2
  Fo <- Chik2 / (N - 1)
  FFo <- Fo - (Dk/(N-1))
  RMSEA <- sqrt(FFo / Dk)
  if (is.nan(RMSEA)) {# the Dk value is so small that it is negative
    RMSEA<-0
  }


  OUT<-list('CFI'=CFI,'RMSEA'=RMSEA)

  return(OUT)

}
