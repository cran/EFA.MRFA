mrfa<-function(SIGMA, dimensionality = 1, random = 10, conv1, conv2, display = TRUE, pwarnings = FALSE){

  m1<-size(SIGMA)[1]
  m<-size(SIGMA)[2]

  ######################################################################
  #  SIGMA : covariance/correlation matrix
  ######################################################################
  if (missing(SIGMA)){
    stop("The argument SIGMA is not optional, please provide a valid covariance/correlation matrix")
  }
  else if (m1!=m){
      stop("The covariance/correlation matrix has to be a square matrix")
  }

  ######################################################################
  #  dimensionality argument: common factors used to find communality estimates
  ######################################################################


  if ((dimensionality>=m)||(dimensionality<0)){
    stop("dimensionality argument has to be between 0 and the number of variables - 1")
  }

  ######################################################################
  #  random argument: number of random starts
  ######################################################################

  if(random<=0){
    stop("random argument has to be a positive value")
  }

  ######################################################################
  #  conv1 argument: Convergence criterion for MRFA (conv1)
  ######################################################################
  if (missing(conv1)){
    if(dimensionality==0){
      conv1 = 0.0000001
    }
    else{
      conv1 = 0.0001
    }
  }
  #The userhas provided the value, check if it is an OK value
  else if(conv1>0.1){
    stop("conv1 (convergence value for computing MRFA) has to be lower or equal than 0.1")
  }
  else if(conv1<=0){
    stop("conv1 (convergence value for computing MRFA) has to be greater than 0")
  }

  ######################################################################
  #  conv2 argument: Convergence criterion for glb step (conv2)
  ######################################################################
  if (missing(conv2)){
    if (dimensionality==0){
      conv2 = 0.000001
    }
    else{
      conv2 = 0.001
    }
  }
  #The user has provided the value, check if it is an OK value
  else if(conv2>0.1){
    stop("conv2 (convergence value for computing glb) has to be lower or equal than 0.1")
  }
  else if(conv2<=0){
    stop("conv2 (convergence value for computing glb) has to be greater than 0")
  }

  ######################################################################
  #  display argument: determines if the output will be displayed in the console
  ######################################################################


  if (display!=0 && display!=1){
    stop("display argument has to be logical (TRUE or FALSE, 0 or 1)")
  }

  ######################################################################
  #  pwarnings argument: determines if the warnings will be printed in the console
  ######################################################################

  if (pwarnings!=0 && pwarnings!=1){
    stop("pwarnings argument has to be logical (TRUE or FALSE, 0 or 1)")
  }

  ################################# Everything  OK #################################
  ################################# Begin Analysis #################################

  gam<--1
  iter<-0
  while(gam[1]==-1){
    if (iter<1000){
      gam<-mrfa_s(SIGMA,dimensionality,random,conv1,conv2,pwarnings)$optimal_communalities
      iter<-iter+1
    }
  }

  SIGMARED<-SIGMA-diag(diag(SIGMA))+diag(c(gam))
  out_eigen<-eigen(SIGMARED)
  r<-dimensionality
  VV<-out_eigen$vectors[,1:r]
  LL<-out_eigen$values[1:r]
  LL<-diag(LL,r)
  A<-VV%*%sqrt(LL)

  OUT<-list('A'=A,'Matrix'=SIGMARED,'gam'=gam)
  if (display==0){
    invisible(OUT)
  }
  else{
    cat('Minimum Rank Factor Analysis Output:\n\n')
    cat('Factor Structure Matrix:\n\n')
    prmatrix(SIGMARED)
    cat('\n\n')
    cat('Optimal Communalities:\n\n')
    buff1=''
    for (i in 1:m){
      cat(sprintf("X%3.0f  %5.4f \n",i,gam[i]))
    }
    cat('\n\n')

    #if (dimensionality==1){
      f2<-size(A)[1]
      f1<-size(A)[2]
    #}
    #else {
    #  f1<-size(A)[1]
    #  f2<-size(A)[2]
    #}
    cat('       ')

    for (k in 1:f1){
      cat(sprintf('   Factor %2.0f',k))
    }
    cat('\n')
    for (i in 1:f2){
      cat(sprintf("V% 3.0f",i))
      for (j in 1:f1){
        cat(sprintf('     % 5.4f',A[i,j]))
      }
      cat("\n")
    }
    cat('\n')
    invisible(OUT)

  }
}
