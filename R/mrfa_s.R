mrfa_s<-function(SIGMA,r,nstarts,conv1,conv2,pwarnings){

  siz<-size(SIGMA)
  m<-siz[2]

  func<-matrix(0,nstarts,1)
  GAM<-matrix(0,m,nstarts)

  ana<-1
  for(ana in 1:nstarts){
    gam<-transpose(runif(m))
    #function nrm
      buff<-matrix(rnorm(m*m),m)

      T<-nrm(buff)

    SIGMARED=SIGMA-diag(diag(SIGMA))+diag(c(gam));
    #function ed
      v<-eigen(SIGMARED)$values
      u<-eigen(SIGMARED)$vectors

      p<-length(v)
      v<-t(sort(v))
      i<-c(1:p)
      j<-c(p:1)
      v<-v[j]
      v<-diag(v)
      u<-u[,i[j]]
      #function nrm
      u<-nrm(u)
      K<-u
      L<-v
      l<-diag(L)

      f1<-sum(l[(r+1):m])
      fold<-f1+2*conv1*f1
      iter<-0

      while ((((fold-f1)>(f1*conv1)) && (f1>(.01*conv1))) || (iter<2)){
        fold<-f1
        iter<-iter+1
        if (r<(m-1)){
          buff1<-K[,(r+1):m]^2
          buff2<-(apply(buff1,1,sum))^.5
          w<-(transpose(buff2))^.5
        }
        if (r==(m-1)){
          w<-(K[,(r+1):m]^2)^.5
        }
        siz_w<-size(w)
        if (siz_w[2]==1){
          Sigma1<-(w%*%transpose(w))*SIGMA
        }
        else{
          Sigma1<-(transpose(w)%*%w)*SIGMA
        }

        resultat<-GreaterLowerBound(Sigma1,conv2,T)

        gam<-resultat$gam
        #T<-resultat$T
        gam<-gam/(w^2)
        SIGMARED=SIGMA-diag(diag(SIGMA))+diag(c(gam))
        #function ed
          l<-eigen(SIGMARED)$values
          K<-eigen(SIGMARED)$vectors

        tmp<-size(l)
        if(tmp[2]<m){
          gam<-(-1)
          f1<-(-1)
          stop()
        }

        #correction of gamma's because negative eigenvalues

        if (l[m]<0){
          gam<-gam-l[m]
          l>l-l[m]
          if (l[m]<(-.1)){
            if (pwarnings==TRUE){
              print(sprintf('correction factor=%12.8f',l[m]))
            }
          }

          f1<-sum(l[(r+1):m])
        }

      }
      prova<-ana
      func[ana]<-f1
      GAM[,ana]<-gam
  }

  f2<-min(func)
  mi<-which.min(func)
  gam<-GAM[,mi]

  OUT<-list('optimal_communalities'=gam,'function_value'=f2)
  return(OUT)

}
