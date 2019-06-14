LineCon<-function(f1,f2,f3,fp1,fp2,fp3){

  Ct <- ((fp2-fp1)*(f3-f1) )/ (fp3-fp1)
  ft <- f1+Ct
  if (ft>=f2){
    ret<-FALSE
  }
  else {
    ret<-TRUE
  }

  return(ret)

}
