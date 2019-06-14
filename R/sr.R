sr<-function(R,N){

  delta<-.25/sqrt(N)

  #izfisher
  izdelta<-(exp(2*delta)-1)/(exp(2*delta)+1)

  siz<-size(izdelta)
  n<-siz[1]

  nR<-R

  while (min(eigen(nR)$values)<0){
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        r<-nR[i,j]
        nr<-r
        zr<-(exp(2*r)-1)/(exp(2*r)+1)

        if (r<-izdelta){
          nr<-(exp(2*(zr+delta))-1)/(exp(2*(zr+delta))+1)
        }
        if (abs(r)<=izdelta){
          nr<-0
        }
        if (r>izdelta){
          nr<-(exp(2*(zr-delta))-1)/(exp(2*(zr-delta))+1)
        }

        nR[i,j]=nr
        nR[j,i]=nr
      }
    }
  }
  return(nR)
}
