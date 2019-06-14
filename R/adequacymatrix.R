adequacymatrix<-function(R,n){

  d<-det(R)

  #Bartlett
  out<-psych::cortest.bartlett(R,n)
  chisq<-out$chisq
  p_value<-out$p.value
  df<-out$df

  #KMO
  kmo_index<-psych::KMO(R)$MSA

  adeq<-list('d'=d,'chisq'=chisq,'p_value'=p_value,'df'=df,'kmo_index'=kmo_index)
  return(adeq)
}
