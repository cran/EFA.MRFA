GreaterLowerBound<-function(C, conv = 0.000001, T, pwarnings = FALSE){

  m1<-size(C)[1]
  m<-size(C)[2]
  ######################################################################
  #  C : covariance/correlation matrix
  ######################################################################
  if (missing(C)){
    stop("The argument C is not optional, please provide a valid covariance/correlation matrix")
  }
  else if (m1!=m){
      stop("The covariance/correlation matrix has to be a square matrix")
  }

  ######################################################################
  #  conv argument: Convergence criterion for glb step
  ######################################################################

  #The user has provided the value, check if it is an OK value
  else if(conv>0.1){
    stop("conv (convergence value for computing glb) has to be lower or equal than 0.1")
  }
  else if(conv<=0){
    stop("conv (convergence value for computing glb) has to be greater than 0")
  }

  ######################################################################
  #  T argument: random matrix for start (can be omitted)
  ######################################################################
  if (missing(T)){
    #function nrm
    buff<-matrix(rnorm(m*m),m)
    T<-nrm(buff)
  }
  else if (m1!=m){
    stop("The T argument has to contain a square matrix")
  }

  ######################################################################
  #  pwarnings argument: determines if the warnings will be printed in the console
  ######################################################################


  if (pwarnings!=0 && pwarnings!=1){
    stop("pwarnings argument has to be logical (TRUE or FALSE, 0 or 1)")
  }

  ################################# Everything  OK #################################
  ################################# Begin Analysis #################################

  siz<-size(C);
  m<-siz[2]

  C0<-C-diag(diag(C))

  iter<-0
  #trace function
    #remove(buff)
    buff<-T%*%C%*%transpose(T)
    f2<-sum(diag(buff))

  fold<-f2+2*conv
  labda=-10000
  extra<-0
  iterextra<-0
  gamma1<-matrix(0,m,1)
  contador<-0

  while (((fold-f2)>conv)|(labda<(-.1*conv))){

    if (is.finite(f2)=='FALSE'){
      stop()
      geterrmessage()
    }
    fold<-f2
    iter<-iter+1
    iterextra<-iterextra+1
    for (i in 1:m){
      yi<-T%*%C0[,i]
      delta1<-(transpose(yi)%*%yi)^.5
      if (delta1==0){
        T[,i]=nrm(T[,i])
      }
      if ((delta1>0)&(delta1<C[i,i])){
        delta2<-matrix(0,m)
        delta2[,]<-delta1
        T[,i]=-yi/delta2
      }
      if (delta1>=C[i,i]){
        T[,i]=-yi/C[i,i]
      }
    }
      buff<-T%*%C%*%transpose(T)
      f2<-sum(diag(buff))
      if ((fold-f2)<=conv & (extra==0)){
        tt<-apply((T^2),2,sum)
        w<-tt<=1.0000000001
        gamma1<-diag(C)
        gg<-diag(C0%*%transpose(T)%*%T%*%C0)
        gg<-as.numeric(gg)
        gamma1[w]<-gg[w]^.5
        labda<-min(eigen(C0+diag(gamma1))$values)
        if(labda<(-.1*conv)){
          extra<-1
          iterextra<-0
        }
      }


      if((extra==1)&(round((iterextra+1)/10)==(iterextra+1)/10)){
        tt<-apply((T^2),2,sum)
        w<-tt<=1.0000000001
        gamma1<-diag(C)
        gg<-diag(C0%*%transpose(T)%*%T%*%C0)
        gg<-as.numeric(gg)
        gamma1[w]<-gg[w]^.5
        labda0<-labda
        labda<-min(eigen(C0+diag(gamma1))$values) #reevaluate labda every 10 cycles
        if(abs(labda-labda0)<.0001*abs(labda0)){
          buff<-matrix(rnorm(m*m),m)
          T<-nrm(buff)
          extra<-0
          extraiter<-0
          #trace function
          remove(buff)
          buff<-T%*%C%*%transpose(T)
          f2<-sum(diag(buff))
        }
      }

  }

  tt<-apply((T^2),2,sum)
  w<-tt<=1.0000000001
  gamma1<-diag(C)
  gg<-diag(C0%*%transpose(T)%*%T%*%C0)
  gg<-as.numeric(gg)
  gamma1[w]<-gg[w]^.5
  labda0<-labda
  labda<-min(eigen(C0+diag(gamma1))$values) #reevaluate labda one last time
  if(labda<(-.1*conv)){
    #stop()
    #geterrmessage()
    if (pwarnings==TRUE){
      print(sprintf('Too much negative eigenvalue (%.6f)',labda))
    }
  }

  gam<-as.numeric(gamma1)
  my_list<-list("gam"=gam)
  return(my_list)

}
