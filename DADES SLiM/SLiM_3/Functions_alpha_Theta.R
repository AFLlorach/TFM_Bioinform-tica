source("./asymptoticMK_local.R")

#############
# Functions #
#############

#calculation of alpha from standard equation
calc.alpha <- function(Ps,Ds,Pn,Dn) {
  alpha <- 1 - (Ds/Dn)*(Pn/Ps)
  return(as.numeric(alpha))
}

#calculation of different theta estimators using Achaz (2009) approach
CalcThetaUnfolded <- function(sfs,w) {
  th <- 0
  for(i in 1:length(sfs)) {
    th <- th + w[i] * i * sfs[i]
  }
  if(sum(w)) th <- th/(sum(w))
  return(th)
}

#weights depending on defined estimators of variability: each estimator focus on different frequencies. 
#Taj=Tajima theta estimator. Watt=Watterson theta estimator. FW=Fay and Wu theta estimator. FL=Fu and Li theta estimator. vhf="very high frequency" theta estimator.
weight.stats.unfolded <-function(nsam,estimator="Taj") {
  w <- array(0,dim=c(floor(nsam-1)))
  if(estimator=="FL") {
    for(i in 1:length(w)) {
      if(i==1) 
        w[i] <- i
      else 
        w[i] <- 0
    }
  }
  if(estimator=="Watt")
    for(i in 1:length(w)) {
      w[i] <- 1/i
    }
  if(estimator=="Taj")
    for(i in 1:length(w)) {
      w[i] <- nsam-i
    }
  if(estimator=="FW")
    for(i in 1:length(w)) {
      w[i] <- i
    }
  if(estimator=="vhf")
    for(i in 1:length(w)) {
      w[i] <- i/(nsam-i)
    }
  return(w)
}

#weights to estimate theta given a frequency interval
#be careful with intervals at small sample sizes. Intervals will be unequal because of discrete ranges
weight.finterval.unfolded <-function(nsam,finit,fend) {
  nf <- nsam-1
  nf.init <- floor(nsam*finit)
  nf.end  <- floor(nsam*fend)
  w <- array(0,dim=c(nf))
  for(i in 1:nf) {
    if(i>=nf.init && i<nf.end)
      w[i] <- 1
    else 
      w[i] <- 0
  }
  if(sum(w)) w <-w/sum(w)
  return(w)
}

#MKT asymptotic calculation considering the unfolded SFS for syn and nsyn in five Theta estimators (Fu&Li, Wattersom, Tajima, Fay&Wu, vhf(very high frequency))
mkt.five.stats.sfs <- function(nsam,sfs.s,sfs.n,div.s,div.n,true.alpha=F,output="quartz",prefix.path="", ...) {
  
  #output="table",output="pdf",output="quartz"
  if((nsam != length(sfs.s)+1) || (nsam != length(sfs.n)+1))
    stop("sfs.s and sfs.n must contain the absolute number of variants for each frequency, from 1 to (nsam-1).")
  if(length(div.s)!=1 || length(div.n)!=1)
    stop("div.s and div.n must contain a single numerical value indicating the divergence in absolute numbers.")
  
  n.nsam <- nsam
  sfs.s <- as.numeric(sfs.s)
  sfs.n <- as.numeric(sfs.n)
  div.syn <- as.numeric(div.s)
  div.nsyn <- as.numeric(div.n)
  
  #Definitions of vectors for theta estimators
  name.stats <- c("FL","Watt","Taj","FW","vhf")
  n.stats <- length(name.stats)
  Theta.stat.Nsyn <- array(0,dim=c(n.stats,1))
  rownames(Theta.stat.Nsyn)  <- name.stats
  colnames(Theta.stat.Nsyn)  <- "theta.n"
  Theta.stat.Syn  <- array(0,dim=c(n.stats,1))
  rownames(Theta.stat.Syn)  <- name.stats
  colnames(Theta.stat.Syn)  <- "theta.s"
  alpha.mkt.stat  <- array(0,dim=c(n.stats,1))
  rownames(alpha.mkt.stat)  <- name.stats
  colnames(alpha.mkt.stat)  <- "alpha"
  daf.stat.theta <- array(0,dim=c(n.stats,3))
  colnames(daf.stat.theta)  <- c("daf","Pi","P0")
  
  #calculate stats
  for(int in 1:n.stats) {
    #weights
    w.stat <- weight.stats.unfolded(n.nsam,estimator=name.stats[int])
    #estimates of theta for syn and nsyn
    Theta.stat.Syn[int,1]   <- CalcThetaUnfolded(sfs=sfs.s,w=w.stat)
    Theta.stat.Nsyn[int,1]  <- CalcThetaUnfolded(sfs=sfs.n,w=w.stat)
    #estimate mean frequency values for each estimate (on x-axis).
    if(sum(w.stat)) daf.stat.theta[int,1] <- round(sum(w.stat/sum(w.stat)*c(1:(n.nsam-1)))/n.nsam,4)
    daf.stat.theta[int,3] <- round(Theta.stat.Syn[int],0)
    daf.stat.theta[int,2] <- round(Theta.stat.Nsyn[int],0)
    #calculation of alpha
    alpha.mkt.stat[int,1]   <- calc.alpha(Ps=Theta.stat.Syn[int],Ds=div.syn,Pn=Theta.stat.Nsyn[int],Dn=div.nsyn)
  }
  #show th alpha values for each stats
  show(alpha.mkt.stat)
  
  #calculate the aMKT using the stats
  aa <- NULL
  tryCatch(
    {
      aa <- asymptoticMK.path(d0=div.syn, d=div.nsyn, xlow=0, xhigh=1, df=as.data.frame(daf.stat.theta), true_alpha=true.alpha, output=output,prefix.path=prefix.path) 
    },
    error = function(e) {
      message(sprintf("Error calculating MKTa for observed data in file %s",f))
    }
  )

  #Result
  cat("Result MKT using theta stats estimators\n")
  show(aa)
  
  myres <- list()
  myres[[1]] <- Theta.stat.Syn
  myres[[2]] <- Theta.stat.Nsyn
  rownames(daf.stat.theta)  <- name.stats
  myres[[3]] <- daf.stat.theta
  myres[[4]] <- alpha.mkt.stat
  myres[[5]] <- aa
  return(myres)
}

#MKT asymptotic calculation considering the unfolded SFS for syn and nsyn in interval Theta estimators (weight of 1 inside the interval and 0 outside)
mkt.intervals.sfs <- function(n.int=5,nsam,sfs.s,sfs.n,div.s,div.n,true.alpha=F,output="quartz",prefix.path="", ...) {
  
  #output="table",output="pdf",output="quartz"
  if((nsam != length(sfs.s)+1) || (nsam != length(sfs.n)+1))
    stop("sfs.s and sfs.n must contain the absolute number of variants for each frequency, from 1 to (nsam-1).")
  if(length(div.s)!=1 || length(div.n)!=1)
    stop("div.s and div.n must contain a single numerical value indicating the divergence in absolute numbers.")
  if(length(n.int) != 1 || !is.numeric(n.int) || n.int >= nsam)
    stop("The number of intervals n.int must be a single numrical value smaller than nsam.")
  
  n.nsam <- nsam
  sfs.s <- as.numeric(sfs.s)
  sfs.n <- as.numeric(sfs.n)
  div.syn <- as.numeric(div.s)
  div.nsyn <- as.numeric(div.n)
  
  #Definitions of vectors for interval theta estimators
  Theta.names <- NULL
  for(i in seq(floor(n.nsam/n.int),n.nsam,floor(n.nsam/n.int))) {
    Theta.names <- c(Theta.names, sprintf("Theta.i%s.Syn",i), sprintf("Theta.i%s.NSyn",i))
  }
  Theta.names <- c(Theta.names,"Div.Syn","Div.Nsyn")
  Theta.int.results <- array(0,dim=c(1+n.int*2+2))
  Theta.int.results <- as.data.frame(Theta.int.results)
  Alpha.theta.int.results <- array(0,dim=c(1,5))
  colnames(Alpha.theta.int.results) <- c("alpha.real","alpha.mkt","alpha.asym","alphaAs.CIL","alphaAs.CIH")
  Alpha.theta.int.results <- as.data.frame(Alpha.theta.int.results)
  
  #Definitions for theta intervals:    
  name.int <- sprintf("%.0f/%.0f",1:n.int,n.int)
  w.int  <- array(0,dim=c(n.int,n.nsam))
  Theta.int.Nsyn <- array(0,dim=c(n.int,1))
  rownames(Theta.int.Nsyn)  <- name.int
  colnames(Theta.int.Nsyn)  <- "theta.n"
  Theta.int.Syn  <- array(0,dim=c(n.int,1))
  rownames(Theta.int.Syn)  <- name.int
  colnames(Theta.int.Syn)  <- "theta.s"
  alpha.mkt.int  <- array(0,dim=c(n.int,1))
  rownames(alpha.mkt.int)  <- name.int
  colnames(alpha.mkt.int)  <- "alpha"
  daf.int.theta <- array(0,dim=c(n.int,3))
  colnames(daf.int.theta)  <- c("daf","Pi","P0")
  
  #calculate stats
  for(int in 1:n.int) {
    f.init=(int-1)/(n.int)
    f.end =(int  )/(n.int)
    w.int <- weight.finterval.unfolded(n.nsam,f.init,f.end)
    #estimates of theta for syn and nsyn
    Theta.int.Syn[int,1]  <- CalcThetaUnfolded(sfs=sfs.s,w=w.int)
    Theta.int.Nsyn[int,1]  <- CalcThetaUnfolded(sfs=sfs.n,w=w.int)
    #estimate mean frequency values for each estimate (on x-axis).
    if(sum(w.int)) daf.int.theta[int,1] <- round(sum(w.int/sum(w.int)*c(1:(n.nsam-1)))/n.nsam,4)
    daf.int.theta[,3] <- sapply(c(Theta.int.Syn),round,0)
    daf.int.theta[,2] <- sapply(c(Theta.int.Nsyn),round,0)
    daf.int.theta <- as.data.frame(daf.int.theta)
    #calculation of alpha
    alpha.mkt.int[int,1]   <- calc.alpha(Ps=Theta.int.Syn[int],Ds=div.syn,Pn=Theta.int.Nsyn[int],Dn=div.nsyn)
  }
  show(alpha.mkt.int)
  
  #calculate the aMKT using the intervals
  aa <- NULL
  tryCatch(
    {
      aa <- asymptoticMK.path(d0=div.syn, d=div.nsyn, xlow=0, xhigh=1, df=daf.int.theta, true_alpha=true.alpha, output=output,prefix.path=prefix.path)
    },
    error = function(e) {
      message(sprintf("Error calculating MKTa for observed data in file %s",f))
    }
  )
  
  #Result
  cat("Result MKT using theta interval estimators\n")
  show(aa)
  
  myres <- list()
  myres[[1]] <- data.frame(Theta.int.Syn=Theta.int.Syn)
  myres[[2]] <- data.frame(Theta.int.Nsyn=Theta.int.Nsyn)
  rownames(daf.int.theta)  <- name.int
  myres[[3]] <- daf.int.theta
  myres[[4]] <- alpha.mkt.int
  myres[[5]] <- aa
  return(myres)
}

#MKT asymptotic calculation considering the unfolded frequency and sample size for each position for syn and nsyn in five Theta estimators (Fu&Li, Wattersom, Tajima, Fay&Wu, vhf(very high frequency))
mkt.five.stats.table <- function(nsam,tab.s,tab.n,div.s,div.n,true.alpha=F,output="quartz",prefix.path="", ...) {

  #output="table",output="pdf",output="quartz"
  if((nsam < max(tab.s[,2])) || (nsam < max(tab.n[,2])))
    stop("tab.s and tab.n must contain a maximum number of samples of nsam.")
  if(length(div.s)!=1 || length(div.n)!=1)
    stop("div.s and div.n must contain a single numerical value indicating the divergence in absolute numbers.")
  
  n.nsam <- nsam
  tab.s <- as.matrix(tab.s)
  tab.n <- as.matrix(tab.n)
  div.syn <- as.numeric(div.s)
  div.nsyn <- as.numeric(div.n)
  
  #Definitions for theta estimators
  name.stats <- c("FL","Watt","Taj","FW","vhf")
  n.stats <- length(name.stats)
  Theta.stat.Nsyn <- array(0,dim=c(n.stats,1))
  rownames(Theta.stat.Nsyn)  <- name.stats
  colnames(Theta.stat.Nsyn)  <- "theta.n"
  Theta.stat.Syn  <- array(0,dim=c(n.stats,1))
  rownames(Theta.stat.Syn)  <- name.stats
  colnames(Theta.stat.Syn)  <- "theta.s"
  alpha.mkt.stat  <- array(0,dim=c(n.stats,1))
  rownames(alpha.mkt.stat)  <- name.stats
  colnames(alpha.mkt.stat)  <- "alpha"
  daf.stat.theta <- array(0,dim=c(n.stats,3))
  colnames(daf.stat.theta)  <- c("daf","Pi","P0")
  
  #calculation of theta and alpha for each theta statistic
  for(int in 1:n.stats) {
    ws <- 0
    #theta for syn:
    for(l.s in 1:dim(tab.s)[1]) {
      n.nsam <- tab.s[l.s,2]
      n.freq <- tab.s[l.s,1]
      if(n.freq>0 && n.freq<n.nsam) {
        w.stat <- weight.stats.unfolded(n.nsam,estimator=name.stats[int]); ws <- ws + 1
        sfs <- array(0,dim=c(n.nsam-1)); sfs[n.freq] <- 1;
        Theta.stat.Syn[int]  <- Theta.stat.Syn[int] + CalcThetaUnfolded(sfs=sfs,w=w.stat)
        #estimate mean frequency values for each estimate (on x-axis).
        if(sum(w.stat)) daf.stat.theta[int,1] <- daf.stat.theta[int,1] + sum(w.stat/sum(w.stat)*c(1:(n.nsam-1))) * 1/n.nsam
      }
    }
    #calculate the average frequency for each statistic using syn data, considering the weights
    if(ws) daf.stat.theta[int,1] <- round(daf.stat.theta[int,1] / ws,4)
    daf.stat.theta[int,3] <- round(Theta.stat.Syn[int],0)
    
    #theta for nsyn:
    for(l.s in 1:dim(tab.n)[1]) {
      n.nsam <- tab.n[l.s,2]
      n.freq <- tab.n[l.s,1]
      if(n.freq>0 && n.freq<n.nsam) {
        w.stat <- weight.stats.unfolded(n.nsam,estimator=name.stats[int]);
        sfs <- array(0,dim=c(n.nsam-1)); sfs[n.freq] <- 1;
        Theta.stat.Nsyn[int]  <- Theta.stat.Nsyn[int] + CalcThetaUnfolded(sfs=sfs,w=w.stat)
      }
    }
    daf.stat.theta[int,2] <- round(Theta.stat.Nsyn[int],0)
    
    #calculation of theta
    alpha.mkt.stat[int]   <- calc.alpha(Ps=Theta.stat.Syn[int],Ds=div.syn,Pn=Theta.stat.Nsyn[int],Dn=div.nsyn)
  }
  show(alpha.mkt.stat)
  
  ##############################################################
  #calculate the aMKT using the stats
  aa <- NULL
  tryCatch(
    {
      aa <- asymptoticMK.path(d0=div.syn, d=div.nsyn, xlow=0, xhigh=1, df=as.data.frame(daf.stat.theta), true_alpha=true.alpha, output=output,prefix.path=prefix.path)
    },
    error = function(e) {
      message(sprintf("Error calculating MKTa for observed data in file %s",f))
    }
  )
  
  #Result
  cat("Result MKT using theta stats estimators\n")
  show(aa)
  
  myres <- list()
  myres[[1]] <- Theta.stat.Syn
  myres[[2]] <- Theta.stat.Nsyn
  rownames(daf.stat.theta)  <- name.stats
  myres[[3]] <- daf.stat.theta
  myres[[4]] <- alpha.mkt.stat
  myres[[5]] <- aa
  return(myres)
}
  
#MKT asymptotic calculation considering the unfolded frequency and sample size for each position for syn and nsyn in interval Theta estimators (weight of 1 inside the interval and 0 outside)
mkt.five.int.table <- function(n.int=5,nsam,tab.s,tab.n,div.s,div.n,true.alpha=F,output="quartz",prefix.path="", ...) {
  
  #output="table",output="pdf",output="quartz"
  if((nsam < max(tab.s[,2])) || (nsam < max(tab.n[,2])))
    stop("tab.s and tab.n must contain a maximum number of samples of nsam.")
  if(length(div.s)!=1 || length(div.n)!=1)
    stop("div.s and div.n must contain a single numerical value indicating the divergence in absolute numbers.")
  if(length(n.int) != 1 || !is.numeric(n.int) || n.int >= nsam)
    stop("The number of intervals n.int must be a single numrical value smaller than nsam.")
  
  n.nsam <- nsam
  tab.s <- as.matrix(tab.s)
  tab.n <- as.matrix(tab.n)
  div.syn <- as.numeric(div.s)
  div.nsyn <- as.numeric(div.n)
  
  #Definitions of vectors for interval theta estimators
  Theta.names <- NULL
  for(i in seq(floor(n.nsam/n.int),n.nsam,floor(n.nsam/n.int))) {
    Theta.names <- c(Theta.names, sprintf("Theta.i%s.Syn",i), sprintf("Theta.i%s.NSyn",i))
  }
  Theta.names <- c(Theta.names,"Div.Syn","Div.Nsyn")
  Theta.int.results <- array(0,dim=c(1+n.int*2+2))
  Theta.int.results <- as.data.frame(Theta.int.results)
  Alpha.theta.int.results <- array(0,dim=c(1,5))
  colnames(Alpha.theta.int.results) <- c("alpha.real","alpha.mkt","alpha.asym","alphaAs.CIL","alphaAs.CIH")
  Alpha.theta.int.results <- as.data.frame(Alpha.theta.int.results)
  
  #Definitions for theta intervals:    
  name.int <- sprintf("%.0f/%.0f",1:n.int,n.int)
  w.int  <- array(0,dim=c(n.int,n.nsam))
  Theta.int.Nsyn <- array(0,dim=c(n.int,1))
  rownames(Theta.int.Nsyn)  <- name.int
  colnames(Theta.int.Nsyn)  <- "theta.n"
  Theta.int.Syn  <- array(0,dim=c(n.int,1))
  rownames(Theta.int.Syn)  <- name.int
  colnames(Theta.int.Syn)  <- "theta.s"
  alpha.mkt.int  <- array(0,dim=c(n.int,1))
  rownames(alpha.mkt.int)  <- name.int
  colnames(alpha.mkt.int)  <- "alpha"
  daf.int.theta <- array(0,dim=c(n.int,3))
  colnames(daf.int.theta)  <- c("daf","Pi","P0")
  
  for(int in 1:n.int) {
    ws <- 0
    #theta for syn:
    for(l.s in 1:dim(tab.s)[1]) {
      n.nsam <- tab.s[l.s,2]
      n.freq <- tab.s[l.s,1]
      if(n.freq>0 && n.freq<n.nsam && n.nsam>n.int) {
        f.init=(int-1)/(n.int)
        f.end =(int  )/(n.int)
        w.int <- weight.finterval.unfolded(n.nsam,f.init,f.end); ws <- ws + 1
        sfs <- array(0,dim=c(n.nsam-1)); sfs[n.freq] <- 1;
        Theta.int.Syn[int]  <- Theta.int.Syn[int] + CalcThetaUnfolded(sfs=sfs,w=w.int)
        #estimate mean frequency values for each estimate (on x-axis).
        if(sum(w.int)) daf.int.theta[int,1] <- daf.int.theta[int,1] + 1/sum(w.int) * sum(w.int*c(1:(n.nsam-1))) * 1/n.nsam
      }
    }
    if(ws) daf.int.theta[int,1] <- round(daf.int.theta[int,1]/ws, 4)
    
    #theta for nsyn:
    for(l.s in 1:dim(tab.n)[1]) {
      n.nsam <- tab.n[l.s,2]
      n.freq <- tab.n[l.s,1]
      if(n.freq>0 && n.freq<n.nsam && n.nsam>n.int) {
        f.init=(int-1)/(n.int)
        f.end =(int  )/(n.int)
        w.int <- weight.finterval.unfolded(n.nsam,f.init,f.end);
        sfs <- array(0,dim=c(n.nsam-1)); sfs[n.freq] <- 1;
        Theta.int.Nsyn[int]  <- Theta.int.Nsyn [int] + CalcThetaUnfolded(sfs=sfs,w=w.int)
      }
    }
    #calculation of alpha
    alpha.mkt.int[int]   <- calc.alpha(Ps=Theta.int.Syn[int],Ds=div.syn,Pn=Theta.int.Nsyn[int],Dn=div.nsyn)
  }
  show(alpha.mkt.int)
  
  #Round the variability estimates
  daf.int.theta[,2] <- sapply(c(Theta.int.Nsyn),round,0)
  daf.int.theta[,3] <- sapply(c(Theta.int.Syn),round,0)
  daf.int.theta <- as.data.frame(daf.int.theta)
  
  ##############################################################
  #calculate the aMKT using the stats
  aa <- NULL
  tryCatch(
    {
      aa <- asymptoticMK.path(d0=div.syn, d=div.nsyn, xlow=0, xhigh=1, df=daf.int.theta, true_alpha=true.alpha, output=output,prefix.path=prefix.path)
    },
    error = function(e) {
      message(sprintf("Error calculating MKTa for observed data in file %s",f))
    }
  )
  
  #Result
  cat("Result MKT using theta interval estimators\n")
  show(aa)
  
  myres <- list()
  myres[[1]] <- data.frame(Theta.int.Syn=Theta.int.Syn)
  myres[[2]] <- data.frame(Theta.int.Nsyn=Theta.int.Nsyn)
  rownames(daf.int.theta)  <- name.int
  myres[[3]] <- daf.int.theta
  myres[[4]] <- alpha.mkt.int
  myres[[5]] <- aa
  return(myres)
}

#Function from Messer and Haller (2017) that allows to include a pre-name path in the pdf files 
#Asymptotic MK analysis including plots and tables. Here a path for plots is included from the original version
createOutputDevice.path <- function(output, plotNum, width=3.5, height=3.5, bg="white", prefix.path="", ...)
{
  titleString <- paste0(prefix.path,"Figure", plotNum)
  filenameString <- paste0(titleString, ".pdf")
  
  if (output == "quartz")
  {
    quartz(title=titleString, width=width, height=height, bg=bg, ...)
  }
  else if (output == "pdf")
  {
    pdf(file=filenameString, width=width, height=height, bg=bg, ...)
  }
  else
  {
    stop(paste0("Unsupported output device name: ", output))
  }
}
asymptoticMK.path <- function(d0, d, xlow, xhigh, df, true_alpha=NA, output="table",force_linear=F,force_exp=F,prefix.path="", ...)
{
  require(nls2)
  require(MASS)
  
  if (is.na(d0) || is.null(d0))
    stop("Malformed d0 (must be numeric).")
  if (is.na(d) || is.null(d))
    stop("Malformed d (must be numeric).")
  if (is.na(xlow) || is.null(xlow))
    stop("Malformed xlow (must be numeric).")
  if (is.na(xhigh) || is.null(xhigh))
    stop("Malformed xhigh (must be numeric).")
  
  #	Bounds-check response variables
  #
  if (d0 <= 0)
    stop("d0 must greater than zero.")
  if (d <= 0)
    stop("d must greater than zero.")
  if ((xlow < 0.0) || (xlow > 1.0))
    stop("xlow must be in the interval [0,1].")
  if ((xhigh < 0.0) || (xhigh > 1.0))
    stop("xhigh must be in the interval [0,1].")
  if (xlow >= xhigh)
    stop("xlow must be less than xhigh.")
  
  #	Read in the file and check its format
  #
  if (NCOL(df) != 3)
    stop("Dataframe df does not contain exactly three tab-separated columns.")
  if (NROW(df) <= 0)
    stop("Dataframe df contains no data rows.")
  
  cols <- names(df)
  
  suppressWarnings(	# the goal is to generate NAs here, so we don't want to see the warnings...
    if (!is.na(as.numeric(cols[1])) || !is.na(as.numeric(cols[2])) || !is.na(as.numeric(cols[3])))
      stop("Dataframe df has a numeric column name; probably the required header row is missing.")
  )
  
  f <- df[[1]]
  p <- df[[2]]
  p0 <- df[[3]]
  
  if (!is.numeric(f))
    stop("The first column of the dataframe df, frequency, is not numeric.")
  if (!is.numeric(p))
    stop("The second column of the dataframe df, p, is not numeric.")
  if (!is.numeric(p0))
    stop("The third column of the dataframe df, p0, is not numeric.")
  if (any(is.na(f)))
    stop("The first column of the dataframe df, frequency, contains NA values (not allowed).")
  if (any(is.na(p)))
    stop("The second column of the dataframe df, p, contains NA values (not allowed).")
  if (any(is.na(p0)))
    stop("The third column of the dataframe df, p0, contains NA values (not allowed).")
  if (any(f < 0.0) || any(f > 1.0))
    stop("The first column of the dataframe df, frequency, contains values out of the required range [0,1].")
  if (any(p < 0))		# note that zero is allowed, although not recommended
    stop("The second column of the dataframe df, p, contains values < 0 (not allowed).")
  if (all(p == 0))		# not all can be zero, however
    stop("The second column of the dataframe df, p, contains all values == 0 (not allowed).")
  if (any(p0 <= 0))
    stop("The third column of the dataframe df, p0, contains values <= 0 (not allowed).")
  
  if (NROW(df) < 3)
    stop("At least three data rows are required, to constrain the fit.")
  
  #	Compute alpha values and trim
  #
  alpha <- 1 - (d0/d) * (p/p0)
  cutoff_f1 <- xlow
  cutoff_f2 <- xhigh
  
  trim <- ((f >= cutoff_f1) & (f <= cutoff_f2))
  
  if (sum(trim) < 3)
    stop("At least three data points are required after trimming the frequency range, to constrain the fit.")
  
  f_trimmed <- f[trim]
  alpha_trimmed <- alpha[trim]
  
  #	Compute the original McDonald-Kreitman alpha; we decided to use the trimmed data for this.
  #
  alpha_nonasymp <- 1 - (d0/d) * (sum(p[trim])/sum(p0[trim]))			# from trimmed data
  #alpha_nonasymp <- 1 - (d0/d) * (sum(p)/sum(p0))						# from untrimmed data
  
  #	Fit models
  #
  mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 10)
  
  if (length(mod1) == 0)
  {
    # try a deeper scan for a decent fit
    mod1 <- fitMKmodel(alpha_trimmed, f_trimmed, 20)
  }
  
  tryCatch({
    mod2 <- lm(alpha_trimmed ~ f_trimmed)
  },
  error=function(cond) {})
  
  linear_better <- FALSE
  
  if ((length(mod1) == 0) || (AIC(mod2) < AIC(mod1)))
    linear_better <- TRUE
  
  if(force_linear) linear_better <- TRUE
  if(force_exp) linear_better <- FALSE
  
  if (!linear_better)
  {
    # if we're leaning toward the exponential model, check for ridiculously wide confidence intervals; sometimes
    # we should reject the exponential model for that reason, because it is basically just a linear model with
    # a "cheat" of a swing up or down to fit one additional data point perfectly, which is lame :->
    ci_pred <- predictNLS(mod1, newdata=data.frame(f_trimmed=1.0))
    alpha_1_low <- ci_pred[6]
    alpha_1_high <- ci_pred[7]
    
    if ((alpha_1_low < -100) || (alpha_1_high > 100))
      linear_better <- TRUE
  }
  
  # Prepare for output and plotting
  full_seq <- seq(from=min(f), to=max(f), by=0.001)
  trimmed_seq <- seq(from=min(f_trimmed), to=max(f_trimmed), by=0.001)
  
  if (linear_better)
  {
    alpha_1_est <- predict(mod2, newdata=data.frame(f_trimmed=1.0))
    ci_pred <- predict(mod2, newdata=data.frame(f_trimmed=1.0), interval="confidence")	# we want confidence, not prediction
    alpha_1_low <- ci_pred[2]
    alpha_1_high <- ci_pred[3]
    const_a <- coef(mod2)["(Intercept)"]
    const_b <- coef(mod2)["f_trimmed"]
    const_c <- NA
    
    full_predicts <- predict(mod2, newdata=data.frame(f_trimmed=full_seq))
    trimmed_predicts <- predict(mod2, newdata=data.frame(f_trimmed=trimmed_seq))
    fit_color <- "red"
  }
  else
  {
    alpha_1_est <- predict(mod1, newdata=data.frame(f_trimmed=1.0))
    const_a <- coef(mod1)["const_a"]
    const_b <- coef(mod1)["const_b"]
    const_c <- coef(mod1)["const_c"]
    
    full_predicts <- predict(mod1, newdata=data.frame(f_trimmed=full_seq))
    trimmed_predicts <- predict(mod1, newdata=data.frame(f_trimmed=trimmed_seq))
    fit_color <- "red"
  }
  
  
  #	BEGIN OUTPUT
  #
  
  result_df <- data.frame(model=(if ((length(mod1) == 0) || linear_better) "linear" else "exponential"), a=const_a, b=const_b, c=const_c, alpha_asymptotic=alpha_1_est, CI_low=alpha_1_low, CI_high=alpha_1_high, alpha_original=alpha_nonasymp, row.names=NULL)
  
  if (output == "table")
  {
    # Just return the results dataframe; no plots or other output
    return(result_df)
  }
  else
  {
    cat("\nAsymptotic McDonald-Kreitman Evaluator: Results\n")
    cat("\n")
    cat("   from: http://benhaller.com/messerlab/asymptoticMK.html\n")
    cat("\n")
    cat("Analysis dataset:\n")
    cat("\n")
    cat(paste0("   d0 = ", format(d0, scientific=FALSE)), "\n")
    cat(paste0("   d = ", format(d, scientific=FALSE)), "\n")
    cat(paste0("   x interval = [", format(xlow, digits=3, nsmall=3, scientific=FALSE), ", ", format(xhigh, digits=3, nsmall=3, scientific=FALSE), "]"), "\n")
    cat(paste0("   f = (", paste0(format(f, scientific=FALSE), collapse=", "), ")\n"))
    cat(paste0("   p0 = (", paste0(format(p0, scientific=FALSE), collapse=", "), ")\n"))
    cat(paste0("   p = (", paste0(format(p, scientific=FALSE), collapse=", "), ")\n"))
    cat("\n")
    
    if ((length(mod1) == 0) || linear_better)
      cat("Fitted model: linear, alpha(x) = a + bx\n")
    else
      cat("Fitted model: exponential, alpha(x) = a + b * exp(−cx)\n")
    if (length(mod1) == 0)
      cat("   (The exponential fit failed to converge, usually because the data are not exponential in shape.)\n")
    
    cat("\n")
    
    cat(paste0("a\t", format(const_a, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("b\t", format(const_b, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("c\t", format(const_c, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("alpha_asymptotic\t", format(alpha_1_est, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("95% CI(lower)\t", format(alpha_1_low, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("95% CI(upper)\t", format(alpha_1_high, digits=5, nsmall=5, scientific=FALSE)), "\n")
    cat(paste0("alpha_original\t", format(alpha_nonasymp, digits=5, nsmall=5, scientific=FALSE)), "\n")
    
    cat("\nIf you use this service, please cite our paper:\n\n   B.C. Haller, P.W. Messer. (2017). asymptoticMK: A web-based tool\n      for the asymptotic McDonald-Kreitman test. G3: Genes, Genomes,\n      Genetics 7(5), 1569-1575. doi:10.1534/g3.117.039693\n\nPlease let us know of any issues with this service at philipp {dot} messer <at> gmail [dot] com.  Thanks!\n\n")
    
    
    #	Output plots
    #
    
    #	PLOT 1: Frequency spectrum: p and p0 versus x
    #
    createOutputDevice.path(output=output,prefix.path=prefix.path, plotNum=1)
    
    par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
    plot(x=c(0,1), y=range(c(p,p0)), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab="polymorphism counts")
    points(x=f, y=p0, col="black", pch=19, cex=0.7)
    points(x=f, y=p, col="red", pch=19, cex=0.7)
    legend(x="topright", legend=c("p0 : neutral region", "p : test region"), col=c("black", "red"), pch=19, cex=0.9, pt.cex=0.7)
    
    closeOutputDevice(output=output, plotNum=1)
    
    
    #	PLOT 2: Frequency spectrum: p and p0 versus x
    #
    createOutputDevice.path(output=output,prefix.path=prefix.path, plotNum=2)
    
    normalized_p0 <- p0 / sum(p0)
    normalized_p <- p / sum(p)
    
    par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
    plot(x=c(0,1), y=range(c(normalized_p,normalized_p0)), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab="normalized SFS")
    points(x=f, y=normalized_p0, col="black", pch=19, cex=0.7)
    points(x=f, y=normalized_p, col="red", pch=19, cex=0.7)
    legend(x="topright", legend=c("p0 : neutral region", "p : test region"), col=c("black", "red"), pch=19, cex=0.9, pt.cex=0.7)
    
    closeOutputDevice(output=output, plotNum=2)
    
    
    #	PLOT 3: alpha(x) ~ x
    #
    createOutputDevice.path(output=output,prefix.path=prefix.path, plotNum=3)
    
    par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
    plot(x=c(0,1), y=range(alpha), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab=expression(paste("MK ", alpha, "(", x, ")")))
    points(x=f, y=alpha, col="black", pch=19, cex=0.7)
    
    closeOutputDevice(output=output, plotNum=3)
    
    
    #	PLOT 4: alpha(x) ~ x plus fit to that data
    #
    createOutputDevice.path(output=output,prefix.path=prefix.path, plotNum=4)
    
    yr <- na.omit(c(alpha, alpha_1_est, max(0.0, alpha_1_low), min(1.0, alpha_1_high), alpha_nonasymp))
    if (!is.na(true_alpha))
      yr <- c(yr, true_alpha)
    
    par(mar=c(3.1, 3.1, 2, 2), tcl=-0.3, mgp=c(1.9, 0.4, 0), family="serif")
    plot(x=c(0,1), y=range(yr), cex.axis=0.8, cex.lab=1.0, type="n", xlab=expression(paste("derived allele frequency, ", italic(x))), ylab=expression(paste("MK ", alpha, "(", x, ")")))
    polygon(x=c(-1, -1, 2.0, 2.0), y=c(alpha_1_low, alpha_1_high, alpha_1_high, alpha_1_low), col="#DDDDDD", border=NA)
    
    if (!is.na(alpha_nonasymp))
      abline(h=alpha_nonasymp, lty=3, col="#777777")
    
    points(x=f, y=alpha, col=ifelse(trim, "black", "#999999"), pch=19, cex=0.7)
    abline(v=cutoff_f1, col="#5555FF")
    abline(v=cutoff_f2, col="#5555FF")
    
    if (!is.na(true_alpha))
      abline(h=true_alpha, col="#00AA00")
    
    lines(x=full_seq, full_predicts, col="#333333", lwd=1)
    lines(x=trimmed_seq, trimmed_predicts, col=fit_color, lwd=2)
    abline(h=alpha_1_est, lty=2, col=fit_color)
    
    box()
    closeOutputDevice(output=output, plotNum=4)
    
    # return the same table as for output="table", but using invisible()
    return(invisible(result_df))
  }
}

################
#    data      #
################
data.sfs <- function() { 
  #data
  n.nsam <- 50
  sfs.n <- c(790,382,229,144,138,89,70,77,47,62,52,46,37,49,35,38,32,29,27,23,30,24,32,23,22,11,14,24,8,19,8,7,16,24,15,9,14,6,12,7,12,8,7,11,15,12,10,5,7)	
  sfs.s <- c(564,318,190,132,128,93,75,93,63,65,64,44,33,47,31,28,39,33,28,28,28,28,22,26,15,25,17,27,26,12,13,15,19,20,12,15,14,21,9,16,8,16,12,13,11,13,11,5,11)
  pos.pol.s <- c(166666)
  pos.pol.n <- c(333333)
  dat.sfs <- rbind(c(sfs.n,pos.pol.n),c(sfs.s,pos.pol.s))
  colnames(dat.sfs) <- c(sprintf("fr%.0f",1:(n.nsam-1)),"PosP")
  rownames(dat.sfs) <- c("sfs.n","sfs.s")
  
  divergence.n <- c(2675)
  divergence.s <- c(2916)
  pos.div.n <- c(333333)
  pos.div.s <- c(166666)
  div <- data.frame(mi=333333,Di=2675,m0=166666,D0=2916)
  
  fix.ben <- data.frame(FixBen=146)
  true.alpha <- fix.ben$FixBen/divergence.n
  
  myres <- list()
  myres[[1]] <- dat.sfs
  myres[[2]] <- div
  myres[[3]] <- fix.ben
  myres[[4]] <- data.frame(true.alpha=true.alpha)
  return(myres)
}  
data.table <- function() {
  
  nsam <- 50
  sfs.n <- c(790,382,229,144,138,89,70,77,47,62,52,46,37,49,35,38,32,29,27,23,30,24,32,23,22,11,14,24,8,19,8,7,16,24,15,9,14,6,12,7,12,8,7,11,15,12,10,5,7,333333,2675,333333,146)	
  sfs.s <- c(564,318,190,132,128,93,75,93,63,65,64,44,33,47,31,28,39,33,28,28,28,28,22,26,15,25,17,27,26,12,13,15,19,20,12,15,14,21,9,16,8,16,12,13,11,13,11,5,11,166666,2916,166666,0)
  dat.sfs <- rbind(sfs.n,sfs.s) 
  #keep divergence in "divergence"
  divergence <- data.frame(mi=dat.sfs[1,nsam+2],Di=dat.sfs[1,nsam+1],m0=dat.sfs[2,nsam+2],D0=dat.sfs[2,nsam+1])
  
  set.seed(123)
  miss.val <- 0.5 
  
  #sample a value of nsam per position for each missing mean value, assuming a binomial distribution of missing. keep nsam
  new.nsam.syn <- rbinom(n=sum(dat.sfs[2,c(1:(nsam-1),nsam+1)]),size=nsam,prob=c(1-miss.val))
  new.nsam.nsyn <- rbinom(n=sum(dat.sfs[1,c(1:(nsam-1),nsam+1)]),size=nsam,prob=c(1-miss.val))
  
  #for each frequency and variant, sample without replacement the new number of samples: keep freq
  new.fr.syn <- NULL
  new.fr.nsyn <- NULL
  i.s <- 1; i.n <- 1; frec<-1;
  for(frec in 1:(nsam-1)) {
    new.fr.syn <- c(new.fr.syn,rhyper(nn=dat.sfs[2,frec],m=frec,n=nsam-frec,k=new.nsam.syn[(i.s):(i.s+dat.sfs[2,frec]-1)])); i.s <- i.s + dat.sfs[2,frec]
    new.fr.nsyn <- c(new.fr.nsyn,rhyper(nn=dat.sfs[1,frec],m=frec,n=nsam-frec,k=new.nsam.nsyn[(i.n):(i.n+dat.sfs[1,frec]-1)])); i.n <- i.n + dat.sfs[1,frec] 
  }
  new.fr.syn <- c(new.fr.syn,rhyper(nn=dat.sfs[2,nsam+1],m=nsam,n=0,k=new.nsam.syn[(i.s):(i.s+dat.sfs[2,nsam+1]-1)])); i.s <- i.s + dat.sfs[2,nsam+1]-1
  new.fr.nsyn <- c(new.fr.nsyn,rhyper(nn=dat.sfs[1,nsam+1],m=nsam,n=0,k=new.nsam.nsyn[(i.n):(i.n+dat.sfs[1,nsam+1]-1)])); i.n <- i.n + dat.sfs[1,nsam+1]-1
  
  tab.syn <- cbind(new.fr.syn,new.nsam.syn)
  colnames(tab.syn) <- c("freq.s","nsam.s")
  tab.nsyn <- cbind(new.fr.nsyn,new.nsam.nsyn)
  colnames(tab.nsyn) <- c("freq.n","nsam.n")
  
  dd.s <- which(tab.syn[,1]==tab.syn[,2])
  tab.syn <- tab.syn[-dd.s,]
  dd.n <- which(tab.nsyn[,1]==tab.nsyn[,2])
  tab.nsyn <- tab.nsyn[-dd.n,]
  
  #Estimate the length of the remaining sequence (needs at least two samples)
  len.syn <- data.frame(len.s=sum(rbinom(n=dat.sfs[2,nsam]-sum(tab.syn[,2]<2),size=nsam,prob=c(1-miss.val)) > 1))
  len.nsyn <- data.frame(len.n=sum(rbinom(n=dat.sfs[1,nsam]-sum(tab.nsyn[,2]<2),size=nsam,prob=c(1-miss.val)) > 1))
  #Estimate the number of fixed positions of the remaining sequence (needs at least two samples)
  div.syn  <- data.frame(div.s=length(dd.s) + sum(rbinom(n=dat.sfs[2,nsam+1],size=1,prob=len.syn$len.s/dat.sfs[2,nsam+2])))
  div.nsyn <- data.frame(div.n=length(dd.n) + sum(rbinom(n=dat.sfs[1,nsam+1],size=1,prob=len.nsyn$len.n/dat.sfs[1,nsam+2])))
  
  myres <- list()
  myres[[1]] <- tab.syn
  myres[[2]] <- tab.nsyn
  myres[[3]] <- div.syn
  myres[[4]] <- div.nsyn
  myres[[5]] <- len.syn
  myres[[6]] <- len.nsyn
  myres[[7]] <- data.frame(true.alpha=146/2675)
  return(myres)
}

#include another data with VCF and annotation file to obtain syn and nsyn (using another libraries)

################
#   Examples   #
################
#Show examples: keep positions for syn and for nsyn and make MKT for all three examples in a vignette.

#Example 1: no missing data
dat <- data.sfs()
sfs <- dat[[1]]
nsam <- length(sfs[1,])
sfs.s <- sfs[2,1:(nsam-1)]
sfs.n <- sfs[1,1:(nsam-1)]
div <- dat[[2]]
div.s <- div$D0
div.n <- div$Di
daf <- data.frame(daf=c(1:(nsam-1))/nsam,Pi=sfs.n,P0=sfs.s)
true.alpha <- dat[[4]]$true.alpha

#plot SFS:
#par(mfrow=c(1,1))
#plot(x=1:(nsam-1),y=sfs[2,1:(nsam-1)]/sfs[2,nsam],type="l",xlab="SFS",ylab="freq",ylim=c(0,max(sfs[1,1:(nsam-1)]/sfs[1,nsam],sfs[2,1:(nsam-1)]/sfs[2,nsam])),main="SFS / positions")
#lines(x=1:(nsam-1),y=sfs[1,1:(nsam-1)]/sfs[1,nsam],type="l",col="red")
#legend("topright",legend=c("Syn","Nsyn"),col=c("black","red"),lty=c(1,1))  s

#Asymptotic SFS (Messer and Petrov 2013, Haller and Messer 2017)
aa <- asymptoticMK.path(d0=div.s, d=div.n, xlow=0, xhigh=1, df=daf, true_alpha=true.alpha, output="pdf",prefix.path="aMKT.")
mkt.sfs <- aa

#aMKT calculated with Theta stats
mkt5s.sfs <- mkt.five.stats.sfs(nsam,sfs.s=sfs.s,sfs.n=sfs.n,div.s=div.s,div.n=div.n,true.alpha,output="pdf",prefix.path="aMKT5s.")
mkt5s.sfs
#aMKT calculated with Theta intervals
mkt5i.sfs <- mkt.intervals.sfs(n.int=5,nsam,sfs.s,sfs.n,div.s,div.n,true.alpha,output="pdf",prefix.path="aMKT5i.")
mkt5i.sfs

#Example 2: unequal sample sizes (missing data)
dat.t <- data.table()
div.s <- dat.t[[3]]$div.s
div.n <- dat.t[[4]]$div.n
true.alpha <- dat.t[[7]]$true.alpha
tab.s <- dat.t[[1]]
tab.n <- dat.t[[2]]

#aMKT calculated with Theta stats in case of unequal sample size
mkt5s.table <- mkt.five.stats.table(nsam,tab.s,tab.n,div.s,div.n,true.alpha,output="pdf",prefix.path="aMKTtabs.")
mkt5s.table
#aMKT calculated with Theta intervals in case of unequal sample size
mkt5i.table <- mkt.five.int.table(n.int=5,nsam,tab.s,tab.n,div.s,div.n,true.alpha,output="pdf",prefix.path="aMKTtabi.")
mkt5i.table

