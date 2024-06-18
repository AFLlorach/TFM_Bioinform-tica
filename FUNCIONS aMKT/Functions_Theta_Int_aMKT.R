setwd("./")
#getwd()

source("./asymptoticMK_local.R")

############################################
# FUNCTIONS
############################################

#calculation of alpha from standard equation
calc.alpha <- function(Ps,Ds,Pn,Dn) {
  alpha <- 1 - (Ds/Dn)*(Pn/Ps)
  return(alpha)
} 
#calculation of different theta estimators
CalcThetaUnfolded <- function(sfs,w) {
  th <- 0
  for(i in 1:length(sfs)) {
    th <- th + w[i] * i * sfs[i]
  }
  if(sum(w)) th <- th/(sum(w))
  return(th)
}
#weights to estimate theta given a frequency interval
weight.finterval.unfolded <-function(nsam,finit,fend) {
  nf <- nsam-1
  nf.init <- (nsam*finit)
  nf.end  <- (nsam*fend)
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

############################################
# DEFINITIONS
############################################

#define sample and size of sfs (to avoid zeros)
n.int <- 5 #nombre d'intervals
nsam <- c(25*2) #sample size
miss.p <- c(0,0.1,0.25,0.5,0.75,0.9) #proportion of missing data
#read files from the directory: looking for slim outputs
slim_files <- system("ls SLiM_*_sum*.txt",intern=T)

#Define data frame to keep results from variability estimates using INTERVALS: theta + alpha
Theta.names <- NULL
for(i in seq(floor(100/n.int),100,floor(100/n.int))) {
  Theta.names <- c(Theta.names, sprintf("Theta.i%s.Syn",i), sprintf("Theta.i%s.NSyn",i))
}
Theta.names <- c(Theta.names,"Div.Syn","Div.Nsyn")
Theta.int.results <- array(0,dim=c(length(slim_files)*length(miss.p),2+n.int+n.int*2+2))
colnames(Theta.int.results) <- c("scenario","missing",sprintf("freq.st%.0f",1:n.int),Theta.names)
Theta.int.results <- as.data.frame(Theta.int.results)

Alpha.theta.int.results <- array(0,dim=c(length(slim_files)*length(miss.p),7))
colnames(Alpha.theta.int.results) <- c("scenario","missing","alpha.real","alpha.mkt","alpha.asym","alphaAs.CIL","alphaAs.CIH")
Alpha.theta.int.results <- as.data.frame(Alpha.theta.int.results)

############################
#RUN: for each slim file calculate the variability and the alpha, considering missing proportions
############################
i <- 1
f <- slim_files[1] #in case you want to check without loop
for(f in slim_files) {
  #Each SLiM file contains 1 SFS for syn and 1 SFS for nonsyn + divergence and finally the true alpha value
  ##########################
  #keep SFS in "dat.sfs"
  dat.sfs <- read.table(file=f,header=T,row.names=1) #<- read a file with SFS, divergence and true beneficial from slim
  #keep divergence in "divergence"
  divergence <- data.frame(mi=dat.sfs[1,nsam+2],Di=dat.sfs[1,nsam+1],m0=dat.sfs[2,nsam+2],D0=dat.sfs[2,nsam+1])
  
  ##########################
  #TO DO:  Re-calculate nsam SFSs assuming that there is a proportion of missing values!!
  #############################################
  miss.val <- miss.p[4] #in case you want to check without loop
  for(miss.val in miss.p) {

    #sample a value of nsam per position for each missing mean value, assuming a binomial distribution of missing. keep nsam
    new.nsam.syn <- rbinom(n=sum(dat.sfs[2,c(1:(nsam-1),nsam+1)]),size=nsam,prob=c(1-miss.val))
    new.nsam.nsyn <- rbinom(n=sum(dat.sfs[1,c(1:(nsam-1),nsam+1)]),size=nsam,prob=c(1-miss.val))
    
    #for each frequency and variant, sample without replacement the new number of samples: keep freq
    new.fr.syn <- NULL
    new.fr.nsyn <- NULL
    i.s <- 1; i.n <- 1; frec <- 1
    for(frec in 1:(nsam-1)) {
      new.fr.syn <- c(new.fr.syn,rhyper(nn=dat.sfs[2,frec],m=frec,n=nsam-frec,k=new.nsam.syn[(i.s):(i.s+dat.sfs[2,frec]-1)])); i.s <- i.s + dat.sfs[2,frec]
      new.fr.nsyn <- c(new.fr.nsyn,rhyper(nn=dat.sfs[1,frec],m=frec,n=nsam-frec,k=new.nsam.nsyn[(i.n):(i.n+dat.sfs[1,frec]-1)])); i.n <- i.n + dat.sfs[1,frec] 
    }
    new.fr.syn <- c(new.fr.syn,rhyper(nn=dat.sfs[2,nsam+1],m=nsam,n=0,k=new.nsam.syn[(i.s):(i.s+dat.sfs[2,nsam+1]-1)])); i.s <- i.s + dat.sfs[2,nsam+1]-1
    new.fr.nsyn <- c(new.fr.nsyn,rhyper(nn=dat.sfs[1,nsam+1],m=nsam,n=0,k=new.nsam.nsyn[(i.n):(i.n+dat.sfs[1,nsam+1]-1)])); i.n <- i.n + dat.sfs[1,nsam+1]-1
    
    tab.syn <- cbind(new.fr.syn,new.nsam.syn)
    colnames(tab.syn) <- c("freq","nsam")
    tab.nsyn <- cbind(new.fr.nsyn,new.nsam.nsyn)
    colnames(tab.nsyn) <- c("freq","nsam")
    
    dd.s <- which(tab.syn[,1]==tab.syn[,2])
    tab.syn <- tab.syn[-dd.s,]
    dd.n <- which(tab.nsyn[,1]==tab.nsyn[,2])
    tab.nsyn <- tab.nsyn[-dd.n,]
    
    #Now, estimate the length of the remaining sequence (needs at least two samples)
    len.syn <- sum(rbinom(n=dat.sfs[2,nsam]-sum(tab.syn[,2]<2),size=nsam,prob=c(1-miss.val)) > 1)
    len.nsyn <- sum(rbinom(n=dat.sfs[1,nsam]-sum(tab.nsyn[,2]<2),size=nsam,prob=c(1-miss.val)) > 1)
    #Estimate the number of fixed positions of the remaining sequence (needs at least two samples)
    div.syn  <- length(dd.s)# + sum(rbinom(n=dat.sfs[2,nsam+1],size=1,prob=len.syn/dat.sfs[2,nsam+2]))
    div.nsyn <- length(dd.n)# + sum(rbinom(n=dat.sfs[1,nsam+1],size=1,prob=len.nsyn/dat.sfs[1,nsam+2]))
    
    #Now We have a new dataset that includes missing (different sample sizes and lost variants and positions)
    
    #Estimation of alpha from theta estimators
    #Theta estimators are LIKE summary of sections of the SFS: 
    #############################################

    #Definitions for theta intervals:    
    w.int  <- array(0,dim=c(n.int,nsam))
    Theta.int.Nsyn <- array(0,dim=c(n.int))
    Theta.int.Syn  <- array(0,dim=c(n.int))
    alpha.mkt.int  <- array(0,dim=c(n.int))
    daf.int.theta <- array(0,dim=c(n.int,3))
    colnames(daf.int.theta)  <- c("daf","Pi","P0")
    
    #calculation of theta and alpha for each theta interval
    int <- 5
    for(int in 1:n.int) {
      ws <- 0
      #theta for syn:
      for(l.s in 1:dim(tab.syn)[1]) {
        n.nsam <- tab.syn[l.s,2]
        n.freq <- tab.syn[l.s,1]
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
      for(l.s in 1:dim(tab.nsyn)[1]) {
        n.nsam <- tab.nsyn[l.s,2]
        n.freq <- tab.nsyn[l.s,1]
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
    #Round the variability estimates
    daf.int.theta[,2] <- sapply(c(Theta.int.Nsyn),round,0)
    daf.int.theta[,3] <- sapply(c(Theta.int.Syn),round,0)
    daf.int.theta <- as.data.frame(daf.int.theta)
    
    ##############################################################
    #calculate the aMKT using the stats
    aa <- NULL
    tryCatch(
      {
        aa <- asymptoticMK(d0=div.syn, d=div.nsyn, xlow=0, xhigh=1, df=daf.int.theta, true_alpha=NA, output="table")
      },
      error = function(e) {
        message(sprintf("Error calculating MKTa for observed data in file %s",f))
      }
    )
    #aa
    #aa$alpha_asymptotic
    #aa$alpha_original
    
    #true.alpha
    #sum beneficial / total.nsyn
    true.alpha <- dat.sfs[1,nsam+3] / dat.sfs[1,nsam+1]

    #TO DO: 
    #Keep all results:
    Theta.int.results$scenario[i] <- f
    Theta.int.results$missing[i]  <- miss.val
    Theta.int.results[i,3:(3+n.int-1)]  <- daf.int.theta[,1]
    Theta.int.results[i,seq(3+n.int,3+n.int-1+n.int*2,2)]  <- Theta.int.Syn
    Theta.int.results[i,seq(3+n.int+1,3+n.int+n.int*2,2)]  <- Theta.int.Nsyn
    Theta.int.results$Div.Nsyn[i]           <- div.nsyn
    Theta.int.results$Div.Syn[i]            <- div.syn
    
    Alpha.theta.int.results$scenario[i]    <- f
    Alpha.theta.int.results$missing[i]     <- miss.val
    Alpha.theta.int.results$alpha.real[i]  <- true.alpha
    Alpha.theta.int.results$alpha.mkt[i]   <- aa$alpha_original
    Alpha.theta.int.results$alpha.asym[i]  <- aa$alpha_asymptotic
    Alpha.theta.int.results$alphaAs.CIL[i] <- aa$CI_low
    Alpha.theta.int.results$alphaAs.CIH[i] <- aa$CI_high
    
    i <- i + 1
  }
}
write.table(x=Theta.int.results,file=sprintf("Theta.int.results.txt"),row.names=F,quote=F)
write.table(x=Alpha.theta.int.results,file=sprintf("Alpha.theta.int.results.txt"),row.names=F,quote=F)
