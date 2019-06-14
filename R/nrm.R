nrm<-function(A) {

  buff3<-dim(A)
  n1<-buff3[1]
  m1<-buff3[2]
  d<-apply(A^2,2,sum)
  w<-d>1e-30
  N<-optimbase::zeros(n1,m1)
  if(sum(w)>0){
    a1<-A[,w]
    a2<-optimbase::ones(n1,1)*d[w]^.5
    N[,w]<-a1/a2[col(a1)]
  }

}
