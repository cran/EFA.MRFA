hullEFA <- function(X, maxQ, extr = "ULS", index_hull = "CAF",  display = TRUE, graph = TRUE, details = TRUE){

  ######################################################################
  #  X : Raw sample scores
  ######################################################################

  if (missing(X)){
    stop("The argument X is not optional, please provide a valid raw sample scores")
  }

  ######################################################################
  #  extr: Determine the extraction method: ML (Maximum Likehood) or ULS (Unweighted Least Squares)
  ######################################################################

  if (extr!="ML" && extr!="ULS"){
    stop("extr argument has to be ML or ULS")
  }

  ######################################################################
  #  index_hull argument: determines which index will be used in the hull method
  ######################################################################


  if (index_hull!="CAF" && index_hull!="CFI" && index_hull!="RMSEA"){
    stop("index_hull argument has to be one of the availables indices (see documentation)")
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
    corr_char = 'Pearson correlation matrices'


  ################################# Everything  OK #################################
  ################################# Begin Analysis #################################

  N<-dim(X)[1]
  p<-dim(X)[2]

  x<-as.matrix(X)
  SIGMA<-cor(X)


  if (missing(maxQ)){

#    if (method_maxQ=="PA+1"){
      # determine the number of maximum number of dimensions suggested by PA
      realevals<-eigen(SIGMA)$values
      evals<-matrix(0,p,500)

      for (a in 1:500){
        evals[,a]<- svd(cor(matrix(x[round(matrix(runif(N*p),N,p)*(N-1)+1)],N,p)))$d
      }

      means <- apply(evals,1,mean)

      for (root in 1:p){
        if (realevals[root] < means[root]){
          maxQ <- root - 1
          break
        }
      }

#    }
#    else if (method_maxQ=="Ledermann"){
#      maxQ <-(2*p+1-sqrt(8*p+1))/2
#    }
#    else {
#      stop("The method_maxQ has to be 'PA+1' or 'Ledermann'")
#    }
  }
  else {
    if (is.numeric(maxQ)){
    #the user has defined a maximum number of factors to be retained, check if it is plausible
      if (maxQ>p){
        stop('The maxQ value has to  be equal or lower than the number of variables')
      }
      maxQ = maxQ - 1
    }
    else {
      stop("The maxQ argument has to be a numeric one, indicating the maximum number of factors to be retained.")
    }

  }

  q<-0

  f_CAF <- 1-psych::KMO(SIGMA)$MSA
  f_CFI <- 0
  f_NNFI <- 0
  f_GFI <- 0
  f_AGFI <- 0

  t<-0
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      t <- t + SIGMA[i,j]^2
    }
  }

  if (Re(t) > (-1)){
    RMSEA <- fitchi_zero(SIGMA,N,t)
    f_RMSEA <- 1-Re(RMSEA)
  }
  else {
    r_CAF<-(-1)
    r_CFI<-(-1)
    r_RMSEA<-(-1)
    stop("Convergence was not possible, the analysis was stopped.")
  }

  #g <- p*q + p-q * (q-1)/2 No rippe anymore
  g <- (1 / 2 * ((p - q) * (p - q + 1)) - p)

  again<-1
  q<-1

  while (again==1){

    if (extr=='ML'){

    #try to reach convergence
    out <- tryCatch(
      {
        out<-ml77(SIGMA,q,0.001)
      },
      error=function(cond) {
        # Non-convergence
        out<-NA
        cat(sprintf('The maximum number of factors available for extraction were %.0f \n\n',q-1))
        return(out)
      }
    )

    # if (is.na(out)){
    #   cat(sprintf('The maximum number of factors available for extraction were %.0f \n\n',q-1))
    #   break()
    # }
    t<-out$t
    A<-out$A

    SIGMAREP <- A%*%t(A)
    #Re2<- (A%*%t(A)) + diag(1,p) - diag(diag(A%*%t(A)))

    if (Re(t) > 0){
      SIGMA_ERROR <- SIGMA - SIGMAREP
      CAF <- 1 - psych::KMO(SIGMA_ERROR - diag(diag(SIGMA_ERROR)) + diag(p))$MSA
      f_CAF <- rbind(f_CAF,CAF)

      OUT<-fitchi(SIGMA,A,N,t)

      CFI<-OUT$CFI
      RMSEA <- 1 - Re(OUT$RMSEA)
    }
    else {
      r_CAF<-(-1)
      r_CFI<-(-1)
      r_RMSEA<-(-1)
      stop("Convergence was not possible, the analysis was stopped.")
    }

    }
    else { # ULS
      out_eigen<-eigen(SIGMA)
      VV<-out_eigen$vectors[,1:q]
      LL<-diag(out_eigen$values[1:q])
      if (q==1){
        VV<-transpose(VV)
        LL<-out_eigen$values[1]
      }

      A<-VV%*%sqrt(LL)

      A <-psych::fa(X, fm="minres",nfactors = q, rotate = "none", min.err = 0.000001)$loadings

      SIGMAREP <- A%*%t(A)
      RES <- SIGMA - SIGMAREP
      t<-0
      for (i in 1:(p-1)){
        for (j in (i+1):p){
          t <- t + RES[i,j]^2
        }
      }
      #Re2<- (A%*%t(A)) + diag(1,p) - diag(diag(A%*%t(A)))

      if (Re(t) > 0){
        SIGMA_ERROR <- SIGMA - SIGMAREP
        CAF <- 1 - psych::KMO(SIGMA_ERROR - diag(diag(SIGMA_ERROR)) + diag(p))$MSA
        f_CAF <- rbind(f_CAF,CAF)

        OUT<-fitchi(SIGMA,A,N,t)

        CFI<-OUT$CFI
        RMSEA <- 1 - Re(OUT$RMSEA)
      }
      else {
        r_CAF<-(-1)
        r_CFI<-(-1)
        r_RMSEA<-(-1)
        stop("Convergence was not possible, the analysis was stopped.")
      }
    }

    f_CFI <- rbind(f_CFI, CFI)
    f_RMSEA <- rbind(f_RMSEA, RMSEA)


    # not rippe
    #g <- rbind(g, p*q + p - q * (q-1)/2)
    g<- rbind(g,g <- (1 / 2 * ((p - q) * (p - q + 1)) - p))
    q<-q+1

    if ((q>maxQ+1)){ # && (CFI<=max(f_CFI)) for better adjustment
      again<-0
    }
    if (q>p){
      again<-0
    }

  }

  maxQ<-q-1

  ### CAF ###

  if (index_hull == "CAF"){
    f_ <- f_CAF
  }
  if (index_hull == "CFI"){
    f_ <- f_CFI
  }
  if (index_hull == "RMSEA"){
    f_ <- f_RMSEA
  }

    out<-cbind(c(0:maxQ),f_,g)
    out_best<-out[1,]
    for (i in 2:(maxQ+1)){
      if (i==2){
        if (max(out_best[2]) < out[i,2]){
          out_best <- rbind(out_best, out[i,])
        }
      }
      else {
        if( size(out_best)[1] == 1){ #fix 17/4/2020
          if ((out_best[2]) < out[i,2]){
            out_best <- rbind(out_best, out[i,])
          }
        }
        else {
          if (max(out_best[,2]) < out[i,2]){
            out_best <- rbind(out_best, out[i,])
          }
        }
      }
    }
    out_c <- out #comlpete output
    out <- out_best

    # Testing triplets of adjacent solutions
    tmp4<-dim(out)[1]

    i<-2

    while (i < tmp4-1){
      f1 <- out[i-1,2]
      f2 <- out[i,2]
      f3 <- out[i+1,2]
      fp1 <- out[i-1,3]
      fp2 <- out[i,3]
      fp3 <- out[i+1,3]

      if (LineCon(f1,f2,f3,fp1,fp2,fp3)==FALSE){
        out <- rbind(out[1:(i-1),],out[(i+1):tmp4,])
        tmp4 <- dim(out)[1]
        i<-1
      }
      i<-i+1
    }

    # st calculation
    tmp4<-dim(out)[1]
    st<-c(matrix(0,tmp4,1))

    for (i in 2:(tmp4-1)){
      fi <- out[i,2]
      fi_p <- out[i-1,2]
      fi_n <- out[i+1,2]
      pi <- out[i,3]
      pi_p <- out[i-1,3]
      pi_n <- out[i+1,3]
      st[i] <- ((fi - fi_p) / (pi - pi_p)) / ((fi_n - fi) / (pi_n - pi))
    }
    out <- cbind(out,st)

    #########

    j<-which(out[,4]==max(out[,4]))

    if (index_hull == "CAF"){
      r_CAF<-out[j,1]
    }
    if (index_hull == "CFI"){
      r_CFI<-out[j,1]
    }
    if (index_hull == "RMSEA"){
      r_RMSEA<-out[j,1]
    }

    r<-out[j,1]


    dif<-(dim(out_c)[1])-(dim(out)[1])
    out_red <- matrix(0,dif,3)
    j<-2
    h<-1
    k2<-dim(out)[1]
    for (i in 2:maxQ){
      if (j<=k2){
        if (out_c[i,1]==out[j,1]){
          #the n factor is good, contains st
          j=j+1
        }
        else {
          # the i factor is under the line
          out_red[h,] <- out_c[i,]
          h<-h+1
        }
      }
      else {
        if (out_c[i,1]==out[j-1,1]){
          #the n factor is good, contains st
        }
        else {
          # the i factor is under the line
          out_red[h,] <- out_c[i,]
          h<-h+1
        }
      }
    }
    if (dif>1 && (out_red[dif,1]==0)){ #sometimes, if the last value is outside the hull, it does not get included in out_red
      out_red[dif,] <- out_c[(dim(out_c)[1]),]
    }



    if (index_hull == "CAF"){
      t_CAF_red <- out_red
      t_CAF<-out
      colnames(t_CAF)<-c('q','f','g','st')
      rownames(t_CAF)<-NULL
      colnames(t_CAF_red)<-c('q','f','g')
      rownames(t_CAF_red)<-NULL
      OUT<-list('Matrix'=t_CAF,'n_factors'=r_CAF)
    }
    if (index_hull == "CFI"){
      t_CFI_red <- out_red
      t_CFI<-out
      colnames(t_CFI)<-c('q','f','g','st')
      rownames(t_CFI)<-NULL
      colnames(t_CFI_red)<-c('q','f','g')
      rownames(t_CFI_red)<-NULL
      OUT<-list('Matrix'=t_CFI,'n_factors'=r_CFI)
    }
    if (index_hull == "RMSEA"){
      t_RMSEA_red <- out_red
      t_RMSEA<-out
      colnames(t_RMSEA)<-c('q','f','g','st')
      rownames(t_RMSEA)<-NULL
      colnames(t_RMSEA_red)<-c('q','f','g')
      rownames(t_RMSEA_red)<-NULL
      OUT<-list('Matrix'=t_RMSEA,'n_factors'=r_RMSEA)
    }


  # print output

  if (display==T){

    if (index_hull=="CAF"){

      cat('HULL METHOD - CAF INDEX\n')
      cat('\n')
      cat('        q      f          g       st\n')
      if (details==FALSE){
        f1<-dim(t_CAF)[1]
        for (i in 1:f1){
          cat(sprintf('       %2.0f      %.4f  %4.0f       %6.4f \n',t_CAF[i,1],t_CAF[i,2],t_CAF[i,3],t_CAF[i,4]))
        }
      }
      else {
        f1<-dim(out_c)[1]
        g1<-dim(t_CAF)[1]
        j <- 1
        h <- 0
        for (i in 1:f1){
          if ((i-h)<=g1){
            if (t_CAF[i-h,1] == (i-j+h)){ #starts in 0
                cat(sprintf('       %2.0f      %.4f  %4.0f       %6.4f \n',t_CAF[i-h,1],t_CAF[i-h,2],t_CAF[i-h,3],t_CAF[i-h,4]))
            }
            else { # outside convex hull
              cat(sprintf('       %2.0f*     %.4f  %4.0f       \n',t_CAF_red[j,1],t_CAF_red[j,2],t_CAF_red[j,3]))
              j <- j+1
              h <- h+1
            }
          }
          else { # outside convex hull
            cat(sprintf('       %2.0f*     %.4f  %4.0f       \n',t_CAF_red[j,1],t_CAF_red[j,2],t_CAF_red[j,3]))
            j <- j+1
            h <- h+1
          }
        }

      }
      cat('\n')
      cat(sprintf('Number of advised dimensions: %.0f \n',r_CAF))
      if (details==TRUE){
        if (g1!=f1){
          cat(sprintf('* Value outside the convex Hull \n'))
        }
      }
      cat('\n')
      cat('-----------------------------------------------\n')
      cat('\n')

    }

    if (index_hull == "CFI"){

      cat('HULL METHOD - CFI INDEX\n')
      cat('\n')
      cat('        q      f          g       st\n')
      if (details==FALSE){
        f1<-dim(t_CFI)[1]
        for (i in 1:f1){
          cat(sprintf('       %2.0f      %.4f  %4.0f       %6.4f \n',t_CFI[i,1],t_CFI[i,2],t_CFI[i,3],t_CFI[i,4]))
        }
      }
      else {
        f1<-dim(out_c)[1]
        g1<-dim(t_CFI)[1]
        j <- 1
        h <- 0
        for (i in 1:f1){
          if ((i-h)<=g1){
            if (t_CFI[i-h,1] == (i-j+h)){ #starts in 0
              cat(sprintf('       %2.0f      %.4f  %4.0f       %6.4f \n',t_CFI[i-h,1],t_CFI[i-h,2],t_CFI[i-h,3],t_CFI[i-h,4]))
            }
            else { # outside convex hull
              cat(sprintf('       %2.0f*      %.4f  %4.0f       \n',t_CAF_red[j,1],t_CAF_red[j,2],t_CAF_red[j,3]))
              j <- j+1
              h <- h+1
            }
          }
          else { # outside convex hull
            cat(sprintf('       %2.0f*     %.4f  %4.0f       \n',t_CFI_red[j,1],t_CFI_red[j,2],t_CFI_red[j,3]))
            j <- j+1
            h <- h+1
          }
        }

      }
      cat('\n')
      cat(sprintf('Number of advised dimensions: %.0f \n',r_CFI))
      if (details==TRUE){
        if (g1!=f1){
          cat(sprintf('* Value outside the convex Hull \n'))
        }
      }

      cat('\n')
      cat('-----------------------------------------------\n')
      cat('\n')

    }

    if (index_hull == "RMSEA"){

      cat('HULL METHOD - RMSEA INDEX\n')
      cat('\n')
      cat('        q      f          g       st\n')
      if (details==FALSE){
        f1<-dim(t_RMSEA)[1]
        for (i in 1:f1){
          cat(sprintf('       %2.0f      %.4f  %4.0f       %6.4f \n',t_RMSEA[i,1],t_RMSEA[i,2],t_RMSEA[i,3],t_RMSEA[i,4]))
        }
      }
      else {
        f1<-dim(out_c)[1]
        g1<-dim(t_RMSEA)[1]
        j <- 1
        h <- 0
        for (i in 1:f1){
          if ((i-h)<=g1){
            if (t_RMSEA[i-h,1] == (i-j+h)){ #starts in 0
              cat(sprintf('       %2.0f      %.4f  %4.0f       %6.4f \n',t_RMSEA[i-h,1],t_RMSEA[i-h,2],t_RMSEA[i-h,3],t_RMSEA[i-h,4]))
            }
            else { # outside convex hull
              cat(sprintf('       %2.0f*     %.4f  %4.0f       \n',t_RMSEA_red[j,1],t_RMSEA_red[j,2],t_RMSEA_red[j,3]))
              j <- j+1
              h <- h+1
            }
          }
          else { # outside convex hull
            cat(sprintf('       %2.0f*     %.4f  %4.0f       \n',t_RMSEA_red[j,1],t_RMSEA_red[j,2],t_RMSEA_red[j,3]))
            j <- j+1
            h <- h+1
          }
        }

      }
      cat('\n')
      cat(sprintf('Number of advised dimensions: %.0f \n',r_RMSEA))
      if (details==TRUE){
        if (g1!=f1){
          cat(sprintf('* Value outside the convex Hull \n'))
        }
      }
      cat('\n')
      cat('-----------------------------------------------\n')
      cat('\n')

    }

    invisible(OUT)
  }

  # graph

  if (graph==T){

    #OLD GRAPH

    ## COLZES
    #names1<-c('t_CAF','t_CFI','t_NNFI','t_GFI','t_AGFI','t_RMSEA','t_RMSR')
    #names2<-c('r_CAF','r_CFI','r_NNFI','r_GFI','r_AGFI','r_RMSEA','r_RMSR')
    #names3<-c('CAF','CFI','NNFI','GFI','AGFI','RMSEA','RMSR')

    #data_to_plot<-cbind(r_CAF,r_CFI,r_NNFI,r_GFI,r_AGFI,r_RMSEA,r_RMSR)

    #buff_new <- data.frame()
    #for (i in 1:7){
    #  buff0<-as.data.frame(get(names1[i])[,1:2])
    #  Index<-names3[i]
    #  colnames(buff0)[1]<-'Factors'
    #  buff0<-cbind(buff0,Index)
    #  buff_new<-rbind(buff_new,buff0)
    #}

    #p<-ggplot2::ggplot(data=buff_new, ggplot2::aes(x=Factors, y=f, group=Index, colour=Index)) +
    #        ggplot2::geom_line() +
    #        ggplot2::geom_point() +
    #        ggplot2::ggtitle("Hull Method") +
    #        ggplot2::geom_vline(xintercept = unique(c(data_to_plot)), linetype="dashed")
    #print(p)

    title<-sprintf('Hull Method: %s Index',index_hull)

    f1<-dim(out_c)[1]
    g1<-dim(out)[1]

    out_c<-as.data.frame(out_c)
    buff<-character(length = f1)
    buff[]<-"f1"
    out_c2<-cbind(out_c,buff)
    colnames(out_c)<-c('Factors','f','g')

    out<-as.data.frame(out)
    buff<-character(length = g1)
    buff[]<-"g1"
    out2<-cbind(out[,1:3],buff)
    colnames(out)<-c('Factors','f','g','st')

    out_final<-rbind(out_c2,out2)
    colnames(out_final)<-c('Factors','f','g','buff')

    ## Colour final graphic

    Factors=f=NULL #for preventing a NOTE

    p<-ggplot2::ggplot(data=out_final, ggplot2::aes(x=Factors, y=f, colour=buff)) +
      ggplot2::geom_line(data=out) +
      ggplot2::geom_point() +
      ggplot2::ggtitle(title) +
      ggplot2::geom_vline(xintercept = r, linetype="dashed") +
      ggplot2::scale_x_continuous(breaks=scales::pretty_breaks(n=f1)) +
      ggplot2::theme(legend.position="none")
    print(p)


    # Plot containing all the proposed factors to be retained regarding all the indices
    # data_to_plot<-cbind(r_CAF,r_CFI,r_RMSEA)
    #
    # if (length(unique(c(data_to_plot)))>1){
    #   colnames(data_to_plot)<-c('CAF','CFI','RMSEA')
    #   new_data_to_plot<-reshape2::melt(data_to_plot)
    #   new_data_to_plot<-new_data_to_plot[,2:3]
    #   colnames(new_data_to_plot)[1]<-'Index'
    #   colnames(new_data_to_plot)[2]<-'Factors'
    #
    #   p<-ggplot2::ggplot(data=new_data_to_plot, ggplot2::aes(x=Index, y=Factors, fill=Index)) +
    #     ggplot2::geom_bar(colour="black", stat="identity")+
    #     ggplot2::guides(fill=FALSE)
    #   print(p)
    # }
    # else {
    #   cat(sprintf('WARNING: Plot was not generated because all the indices suggest the same number of factors (%.0f)\n\n',r_CAF))
    # }
  }

  if (display==F){
    invisible(OUT)
  }

}
