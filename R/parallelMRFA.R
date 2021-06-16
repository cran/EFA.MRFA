parallelMRFA<-function(X, Ndatsets = 500, percent = 95, corr = "Pearson", display = TRUE, graph = TRUE){

  #begin timer
  ptm <- proc.time()

  ######################################################################
  #  X : Raw sample scores
  ######################################################################

  if (missing(X)){
    stop("The argument X is not optional, please provide a valid raw sample scores")
  }

  ######################################################################
  #  Ndatsets: Number of random datasets used to compute the random distribution of eigenvalues
  ######################################################################

  if(Ndatsets<0){
    stop("Ndatsets argument has to be a positive value")
  }

  ######################################################################
  #  percent: Desired percentile of distribution of random eigenvalues (for example 95 for the 95th percentile) to be used as threshold
  ######################################################################

  if (percent>99.99 || percent<=0){
    stop("percent argument has to be between 1 and 99")
  }

  ######################################################################
  #  corr: Determine if Pearson or Polychoric matrix will be used (1: Computes Pearson correlation matrices 2: Computes Polychoric correlation matrices)
  ######################################################################

  if (corr!="Pearson" && corr!="Polychoric"){
      if (corr=="pearson"){
        corr = 1
      }
      else if (corr=="polychoric"){
        corr = 2
      }
    else {
      stop("corr argument has to be 'Pearson' for computing Pearson correlation or 'Polychoric' for computing Polychoric/Tetrachoric correlation)")
    }
  }
  else {
    if (corr=="Pearson"){
      corr = 1
    }
    else if (corr=="Polychoric"){
      corr = 2
    }
  }


  ######################################################################
  #  display argument: determines if the output will be displayed in the console
  ######################################################################

  if (display!=0 && display!=1){
    stop("display argument has to be logical (TRUE or FALSE, 0 or 1)")
  }

  ######################################################################
  #  graph argument: determines if the Scree Test plot will be printed
  ######################################################################

  if (graph!=0 && graph!=1){
      stop("graph argument has to be logical (TRUE or FALSE, 0 or 1)")
  }


  if (corr==1){
    corr_char = 'Pearson correlation matrices'
  }
  else if (corr==2){
    corr_char = 'Polychoric correlation matrices'
  }

  ################################# Everything  OK #################################
  ################################# Begin Analysis #################################

  siz<-size(X)
  N<-siz[1]
  m<-siz[2]

  if (corr==1){
    R<-cor(X)
  }
  else if (corr==2){
    #Polychoric matrix
    R<-(psych::polychoric(X, correct = FALSE))$rho
  }

  #check adequacy of the matrix (determinant, Bartlett & KMO)
  adeq<-adequacymatrix(R,N)

  check_adeq<-0
  if (adeq$kmo_index<0.7){
    check_adeq<-1
  }

  check_sr<-0
  D<-eigen(R)$values
  # Check if the matrix is positive-defined, if it is not, perform a smoothing
  if (min(D)<0){
    R<-sr(R,N)
    check_sr<-1
  }

  out<-mrfa(SIGMA=R,dimensionality=(m-1),random=10,display = 0)

  SIGMARED<-out$Matrix

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

  realevals<-diag(v)
  si<-sum(realevals)
  realevals_p<-(realevals/si)*100

  evals<-array(0,dim=c(0,0))
  evals_p<-a<-array(0,dim=c(0,0))

  ptm_one <- proc.time()

  for (a in 1:Ndatsets){

    # if (display==TRUE){
    #
    #   if (a==1){
    #     #calculate the time elapsed for computing the first one for estimating the total elapsed time
    #     ptm_one <- proc.time()
    #   }
    #
    #   if (a==2){
    #     #waitbar
    #     cat('Computing Parallel Analysis: Please wait \n')
    #     pb <- txtProgressBar(min = 0, max = Ndatsets-1, style = 3)
    #   }
    # }

    Xi=matrix(0,N,m)

    buff<-matrix(runif((N*m), min=0, max=1),N)
    buff2<-matrix(apply(buff,2,sort),N)
    J<-array(0,dim=c(0,0))
    buff_i<-array(0,dim=c(0,0))
    buff_n<-array(0,dim=c(0,0))

    for (ii in 1:m){
      remove(buff_i)
      remove(buff_n)
      remove(J)
      buff_i<-sample(buff2[,ii])
      J<-sort.list(buff_i,decreasing=TRUE)
      buff_n=as.numeric(unlist(X[J[],ii]))
      Xi[,ii]=buff_n
    }

    if (corr==1){
      Ri<-cor(Xi)
    }
    else if (corr==2){
      #polychoric matrix
      Ri<-(psych::polychoric(Xi,correct = FALSE))$rho
    }

    D<-eigen(Ri)$values

    SIGMARED<-mrfa(SIGMA = Ri, dimensionality = (m-1), random = 10,conv1 = 0.001, conv2 = 0.01, display = 0)$Matrix

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

    L<-diag(v)

    Lp<-L/sum(L)*100

    evals<-c(evals,L)
    evals_p<-c(evals_p,Lp)

    if (display==TRUE){

      compT <- proc.time() - ptm_one
      compT<-compT[3]
      compT<-compT*(Ndatsets-a)/a

      secondsInAMinute = 60
      secondsInAnHour = 60 * secondsInAMinute
      secondsInADay = 24 * secondsInAnHour

      days <- floor(compT / secondsInADay)

      hourSeconds <- compT %% secondsInADay
      hours <- floor(hourSeconds / secondsInAnHour)

      minuteSeconds <- hourSeconds %% secondsInAnHour
      minutes <- floor(minuteSeconds / secondsInAMinute)

      remainingSeconds <- minuteSeconds %% secondsInAMinute
      seconds <- ceiling(remainingSeconds)

      if (compT > 3600){
        if (days >= 1){ #Very very rare, but just to be sure
          cat("Computing PA. Time remaining: +24 hours                                                     \r")
          flush.console()
        }
        else {
          cat("Computing PA. Time remaining: ", hours,"hours, ",minutes, "minutes and ",seconds, "seconds \r")
          flush.console()
        }
      }
      else{
        if (compT >= 60){
          cat("Computing PA. Time remaining: ", minutes, "minutes and ",seconds,"seconds \r")
          flush.console()
        }
        if (compT < 60) {
          cat("Computing PA. Time remaining ",seconds,"seconds                                                                  \r")
          flush.console()
        }
      }

      #if (a==1){
      #   time.taken.mrfa <- proc.time() - ptm_one
      #
      #   time.taken.mrfa<-time.taken.mrfa[3]
      #
      #   et.seconds<-time.taken.mrfa*Ndatsets
      #
      #
      #   if (corr==1){
      #     et.minutes<-(et.seconds/60)*1.8
      #   }
      #   else {
      #     et.minutes<-et.seconds/60
      #   }
      #
      #   if (et.minutes<=1){
      #     cat('Estimated time for the analysis: less than a minute')
      #   }
      #   if (et.minutes>1 && et.minutes<=1.5) {
      #     cat(sprintf('Estimated time for the analysis: %3.0f minute',round(et.minutes)))
      #   }
      #   if (et.minutes>1.5) {
      #     cat(sprintf('Estimated time for the analysis: %3.0f minutes',round(et.minutes)))
      #   }
      #
      #   cat('\n\n')
      # }
      # else{
      #   Sys.sleep(0.1)
      #   # update progress bar
      #   setTxtProgressBar(pb, a-1)
      # }
      }

  }
  if (display==TRUE){
    cat("\r","                                                                                                  ","\r")

    # close(pb)
  }

  evals<-matrix(evals,m)
  evals<-apply(evals,2,sort,decreasing=T)
  evals_buff<-apply(evals,1,sort,decreasing=T)
  evals_decreasing<-apply(evals_buff,1,sort,decreasing=T)
  means<-apply(evals_decreasing,1,mean)

  evals_p<-matrix(evals_p,m)
  evals_p<-apply(evals_p,2,sort,decreasing=T)
  evals_buff_p<-apply(evals_p,1,sort,decreasing=T)
  evals_decreasing_p<-apply(evals_buff_p,1,sort,decreasing=T)
  means_p<-apply(evals_decreasing_p,1,mean)

  buff_percent=(Ndatsets+1)-(round((percent*Ndatsets)/100)) #revert because evals is decreasing

  percentiles<-evals_decreasing[,buff_percent]

  percentiles_p<-evals_decreasing_p[,buff_percent]

  #realevals_P contains the real-data % of variance
  #means_p contains the mean of random % of variance
  #percentiles_p contains the percentile of random % of variance

  campana<-0
  nf_per<-0
  nf_mean<-0
  for (i in 1:m){
    if (campana==0){
      if (realevals_p[i]>percentiles_p[i]){
        nf_per=nf_per+1
        nf_mean=nf_mean+1
      }
      else if(realevals_p[i]>means_p[i]){
        nf_mean=nf_mean+1
      }
      else{
        campana<-1
      }
    }
  }

  time.taken <- proc.time() - ptm

  seconds<-time.taken[3]
  hours<-floor(seconds/3600)
  minutes<-floor(seconds/60)
  seconds<-floor(seconds-(minutes*60))
  total_time<-sprintf('%02.0f:%02.0f:%02.0f',hours,minutes,seconds)

  OUT<-list('Real_Data'=realevals_p,'Real_Data_eigenv'=realevals, 'Mean_random'=means_p,'Percentile_random'=percentiles_p,'N_factors_mean'=nf_mean,'N_factors_percentiles'=nf_per)

  if (graph==1){

    buff_realevals<-realevals_p
    buff_means<-means_p
    buff_percentiles<-percentiles_p

    min_val=min(buff_realevals,buff_means,buff_percentiles)
    max_val=max(buff_realevals,buff_means,buff_percentiles)

    ##### old plot

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

    data_to_plot=cbind(buff_realevals,buff_means,buff_percentiles)

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

    ################# new plot with ggplot2 #########################

    max_y_val=max(buff_realevals,buff_means,buff_percentiles)

    Factors=Variance=Source=NULL #for preventing a NOTE

    colnames(data_to_plot)<-c('Real-Data','Mean of random','Percentile of random')
    #rownames(data_to_plot)<-
    new_data_to_plot<-reshape2::melt(data_to_plot)
    colnames(new_data_to_plot)[1]<-'Factors'
    colnames(new_data_to_plot)[2]<-'Source'
    colnames(new_data_to_plot)[3]<-'Variance'
    p<-ggplot2::ggplot(new_data_to_plot, ggplot2::aes(x=Factors, y=Variance, col=Source, shape=Source)) + ggplot2::ggtitle("Parallel Analysis") + ggplot2::geom_vline(xintercept = nf_mean, linetype="dashed") + ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::scale_y_continuous(name ="% Explained Common Variance", breaks=seq(0,max_y_val,round(max_y_val*.1))) + ggplot2::scale_x_continuous(breaks=seq(0,m,round(m*.1))) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), strip.text.x = ggplot2::element_blank(), strip.background = ggplot2::element_rect(colour="white", fill="white"),legend.position=c(.85,.75))
    print(p)
  }
  if (display==0){
    invisible(OUT)
  }
  else if (display==1){
    cat('\n')
    cat('Parallel Analysis (PA) based on Minimum Rank Factor Analysis\n\n')

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
    cat(sprintf('  Correlation matrices analized:                %s\n',corr_char))
    cat(sprintf('  Number of random correlation matrices:        %.f\n',Ndatsets))
    cat('  Method to obtain random correlation matrices: Permutation of the raw data\n\n')

    cat(sprintf('Item      Real-data        Mean of random   %.f percentile of random\n',percent))
    cat('          % of variance    % of variance    % of variance\n\n')
    buff_realevals<-realevals_p
    buff_means<-means_p
    buff_percentiles<-percentiles_p

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
        buff2=sprintf('%5.2f            %5.2f            %5.2f\n',abs(realevals_p[i]), abs(means_p[i]), abs(percentiles_p[i]))
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

    #cat('References \n\n')
    #cat('Buja, A., & Eyuboglu, N. (1992). Remarks on parallel analysis. Multivariate Behavioral Research, 27(4), 509-540.\n\n')
    #cat('Timmerman, M. E., & Lorenzo-Seva, U. (2011). Dimensionality Assessment of Ordered Polytomous Items with Parallel Analysis. Psychological Methods, 16, 209-220.\n')

    if (check_sr==1){
      cat('WARNING: The matrix was not positive-defined, a smoothing procedure has been applied (Devlin, Gnanadesikan & Kettenring, 1981)\n\n')
    }
    if (check_adeq==1){
      cat('WARNING: The matrix is not suitable for performing factor analysis (KMO <0.7)\n\n')
    }

    #cat(sprintf('Computing time: %s \n\n',total_time))

    invisible(OUT)
  }

}
