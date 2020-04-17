testme<-function(example = FALSE){

  cat('\nThis is an example function designed for ilustrating the functions included in EFA.MRFA package\n')
  cat('For more info about the package and the available functions type ?"EFA.MRFA-package" \n\n')

  IDAQ<-IDAQ

  siz<-size(IDAQ)
  N<-siz[1]
  m<-siz[2]

  R<-cor(IDAQ)

  #check adequacy of the matrix (determinant, Bartlett & KMO)
  adeq<-adequacymatrix(R,N)

  #Parallel Anallysis
  if (example == FALSE){
    out_PA<-parallelMRFA(X=IDAQ, display = F, graph = T)
  }
  else {
    out_PA<-parallelMRFA(X=IDAQ,Ndatsets = 5, display = 0, graph = T)
  }

  #MRFA
  out_MRFA<-mrfa(SIGMA=R,dimensionality=3,random=10,display = 0)

  A<-out_MRFA$A

  resul<-PCovR::promin(A)
  P<-resul$loadings
  U<-resul$U

  PHI<-t(U)%*%U

  #Printing

  cat('Dataset used: IDAQ (100 obs., 23 variables)\n\n')

  cat('Adequacy of the Dispersion Matrix:\n\n')
  cat(sprintf('Determinant of the matrix     = %17.15f\n',adeq$d))
  cat(sprintf('Bartlett\'s statistic          = %7.1f (df = %5.0f; P = %7.6f)\n',adeq$chisq,adeq$df,adeq$p_value))
  cat(sprintf('Kaiser-Meyer-Olkin (KMO) test = %7.5f ',adeq$kmo_index))
  if     (adeq$kmo_index >= 0.9){ cat(sprintf('(very good)'))}
  else if (adeq$kmo_index >= 0.8){ cat(sprintf('(good)'))}
  else if (adeq$kmo_index >= 0.7){ cat(sprintf('(fair)'))}
  else if (adeq$kmo_index >= 0.6){ cat(sprintf('(mediocre)'))}
  else if (adeq$kmo_index >= 0.5){ cat(sprintf('(bad)'))}
  else                     { cat(sprintf('(inaceptable)'))}
  cat('\n\n')

  cat('Implementation details:\n\n')
  cat('  Correlation matrices analized:                Pearson correlation matrices\n')
  cat('  Number of random correlation matrices:        500\n')
  cat('  Method to obtain random correlation matrices: Permutation of the raw data\n\n')

  cat('Parallel Analysis (PA) based on Minimum Rank Factor Analysis\n\n')

  cat('Item      Real-data        Mean of random   95 percentile of random\n')
  cat('          % of variance    % of variance    % of variance\n\n')
  buff_realevals<-out_PA$Real_Data
  buff_means<-out_PA$Mean_random
  buff_percentiles<-out_PA$Percentile_random

  campana<-0
  nf_per<-0
  nf_mean<-0
  for (i in 1:m){
    buff<-sprintf(' %3.0f      ',i)
    if (campana==0){
      if (buff_realevals[i]>buff_percentiles[i]){
        buff2=sprintf('%5.2f**          %5.2f            %5.2f\n',buff_realevals[i], buff_means[i], buff_percentiles[i])
        nf_per=nf_per+1
        nf_mean=nf_mean+1
      }
      else if(buff_realevals[i]>buff_means[i]){
        buff2=sprintf('%5.2f*           %5.2f            %5.2f\n',buff_realevals[i], buff_means[i], buff_percentiles[i])
        nf_mean=nf_mean+1
      }
      else{
        buff2=sprintf('%5.2f            %5.2f            %5.2f\n',buff_realevals[i], buff_means[i], buff_percentiles[i])
        campana<-1
      }
    }
    else{
      buff2=sprintf('%5.2f            %5.2f            %5.2f\n',buff_realevals[i], buff_means[i], buff_percentiles[i])
    }
    cat(buff,buff2)
  }
  cat('\n')
  if (nf_mean>nf_per){
    cat(sprintf('**  Advised number of factors:   %.0f\n',nf_per))
    cat(sprintf('*   Advised number of factors:   %.0f\n\n',nf_mean))
  }
  else {
    cat(sprintf('**  Advised number of factors:   %.0f\n\n',nf_per))
  }

  hullEFA(IDAQ, extr = "ML", graph = FALSE)

  #FACTOR ANALYSIS
  cat('Rotated Loading Matrix using Promin (Lorenzo-Seva, 1999)\n\n')
  cat('       Factor 1   Factor 2   Factor 3\n')
  f1<-size(P)[1]
  for (i in 1:f1){
    cat(sprintf("V% 3.0f  % 5.4f    % 5.4f    % 5.4f\n",i,P[i,1],P[i,2],P[i,3]))
  }
  cat('\n')

  #Correlation matrix
  cat('Correlation between factors\n\n')
  cat('          Factor 1  Factor 2  Factor 3\n')
  f1<-size(PHI)[1]
  for (i in 1:f1){
    cat(sprintf('Factor %1.0f % 5.4f   % 5.4f   % 5.4f \n',i,PHI[1,i],PHI[2,i],PHI[3,i]))
  }
  cat('\n')

  #PLOT

  # min_val=min(buff_realevals,buff_means,buff_percentiles)
  # max_val=max(buff_realevals,buff_means,buff_percentiles)
  #
  # yrange=as.numeric(cbind(min_val,max_val))
  #
  # xrange=matrix(0,1,2)
  # xrange[1]=size(buff_realevals)[2]
  # xrange[2]=size(buff_realevals)[1]
  #
  # plot(xrange,yrange,type='n',xlab='Factors',ylab="% Explained Common Variance",xaxt='n',yaxt='n')
  #
  # if(max_val<10){
  #   axis_y<-1
  # }
  # else{
  #   axis_y<-round(max_val/10)
  # }
  # axis_y2<-seq(0,round(max_val+1),by=axis_y)
  # axis(2,axis_y2)
  #
  # colors<-rainbow(3)
  #
  # if(xrange[1]<10){
  #   axis_x<-1
  # }
  # else{
  #   axis_x<-round((xrange[1]+1)/10)
  # }
  # axis_x2<-seq(0,round(xrange[1]+1),by=axis_x)
  # axis(1,axis_x2)
  # axis(1,nf_mean,col="#000000")
  # grid(nx=NA,ny=NULL)
  # abline(v=nf_mean,col="#000000",lty=3)
  #
  # linetype<-c(1:3)
  # plotchar<-seq(18,18+3,1)
  #
  # data_to_plot=cbind(buff_realevals,buff_means,buff_percentiles)
  #
  # for (i in 1:3){
  #   lines(data_to_plot[,i],type='b',lwd=1.5,lty=linetype[i],col=colors[i],pch=plotchar[i])
  # }
  # title("Parallel Analysis")
  # buff=character(length=3)
  # buff[1]="Real-Data"
  # buff[2]="Mean of random"
  # buff[3]="Percentile of random"
  #
  # legend(xrange[1]-round(xrange[1]*0.4),yrange[2],buff,cex=0.8, col=colors,pch=plotchar,lty=linetype)

  cat('This is an example function designed for ilustrate the functions included in EFA.MRFA package\n')
  cat('For more info about the package and the available functions type ?"EFA.MRFA-package" \n\n')

  out<-c("testme output has been printed in to the console, it is not designed for returning any value")

}
