fitchi_zero<-function(R,N,t){

  p<-dim(R)[2]
  k<-0

  Dk <- round((1/2 * ((p-k) * (p-k+1)) -p))
  Chik <- (N-1) * t
  Fo <- Chik / (N - 1)
  FFo <- Fo - (Dk/(N-1))
  if (FFo < 0){
    FFo <- 0
  }
  RMSEA <- sqrt(FFo / Dk)

  return(RMSEA)

}
