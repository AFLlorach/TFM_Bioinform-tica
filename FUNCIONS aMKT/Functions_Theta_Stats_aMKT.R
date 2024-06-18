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
weight.stats.unfolded <-function(nsam,estimator="Taj") {
  w <- array(0,dim=c(floor(nsam-1)))
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
  if(estimator=="FL") {
    for(i in 1:length(w)) {
      if(i==1) 
        w[i] <- i
      else 
        w[i] <- 0
    }
  }
  return(w)
}

############################################
# DEFINITIONS
############################################

#define sample and size of sfs (to avoid zeros)
nsam <- c(25*2) #sample size
miss.p <- c(0,0.1,0.25,0.5,0.75,0.9) #proportion of missing data
#read files from the directory: looking for slim outputs
slim_files <- system("ls SLiM_*_sum*.txt",intern=T)

#Define data frame to keep results from variability estimates using FuLi, Watt, Taji, FW, vhf: theta + alpha
Theta.stat.results <- array(0,dim=c(length(slim_files)*length(miss.p),14+5))
colnames(Theta.stat.results) <- c("scenario","missing",sprintf("freq.st%.0f",1:5),"Theta.FuLi.Syn","Theta.FuLi.Nsyn",
                             "Theta.Watt.Syn","Theta.Watt.Nsyn",
                             "Theta.Taji.Syn","Theta.Taji.Nsyn",
                             "Theta.FayWu.Syn","Theta.FayWu.Nsyn",
                             "Theta.vhf.Syn","Theta.vhf.Nsyn",
                             "Div.Syn","Div.Nsyn")
Theta.stat.results <- as.data.frame(Theta.stat.results)

Alpha.theta.stat.results <- array(0,dim=c(length(slim_files)*length(miss.p),7))
colnames(Alpha.theta.stat.results) <- c("scenario","missing","alpha.real","alpha.mkt","alpha.asym","alphaAs.CIL","alphaAs.CIH")
Alpha.theta.stat.results <- as.data.frame(Alpha.theta.stat.results)

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
    i.s <- 1; i.n <- 1; frec<-1;
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
    
    #Estimate the length of the remaining sequence (needs at least two samples)
    len.syn <- sum(rbinom(n=dat.sfs[2,nsam]-sum(tab.syn[,2]<2),size=nsam,prob=c(1-miss.val)) > 1)
    len.nsyn <- sum(rbinom(n=dat.sfs[1,nsam]-sum(tab.nsyn[,2]<2),size=nsam,prob=c(1-miss.val)) > 1)
    #Estimate the number of fixed positions of the remaining sequence (needs at least two samples)
    div.syn  <- length(dd.s)# + sum(rbinom(n=dat.sfs[2,nsam+1],size=1,prob=len.syn/dat.sfs[2,nsam+2]))
    div.nsyn <- length(dd.n)# + sum(rbinom(n=dat.sfs[1,nsam+1],size=1,prob=len.nsyn/dat.sfs[1,nsam+2]))
    
    #Now We have a new dataset that includes missing (different sample sizes and lost variants and positions)
    
    #Estimation of alpha from theta estimators
    #Theta estimators are LIKE summary of sections of the SFS: 
    #############################################

    #Definitions for theta estimators
    name.stats <- c("FL","Watt","Taj","FW","vhf")
    n.stats <- length(name.stats)
    Theta.stat.Nsyn <- array(0,dim=n.stats)
    Theta.stat.Syn  <- array(0,dim=n.stats)
    alpha.mkt.stat  <- array(0,dim=n.stats)
    daf.stat.theta <- array(0,dim=c(n.stats,3))
    colnames(daf.stat.theta)  <- c("daf","Pi","P0")
    
    #calculation of theta and alpha for each theta statistic
    int <- 5
    for(int in 1:n.stats) {
      ws <- 0
      #theta for syn:
      for(l.s in 1:dim(tab.syn)[1]) {
        n.nsam <- tab.syn[l.s,2]
        n.freq <- tab.syn[l.s,1]
        if(n.freq>0 && n.freq<n.nsam) {
          w.stat <- weight.stats.unfolded(n.nsam,estimator=name.stats[int]); ws <- ws + 1
          sfs <- array(0,dim=c(n.nsam-1)); sfs[n.freq] <- 1;
          Theta.stat.Syn[int]  <- Theta.stat.Syn[int] + CalcThetaUnfolded(sfs=sfs,w=w.stat)
          #estimate mean frequency values for each estimate (on x-axis).
          if(sum(w.stat)) daf.stat.theta[int,1] <- daf.stat.theta[int,1] + 1/sum(w.stat) * sum(w.stat*c(1:(n.nsam-1))) * 1/n.nsam
        }
      }
      #calculate the average frequency for each statistic using syn data, considering the weights
      if(ws) daf.stat.theta[int,1] <- round(daf.stat.theta[int,1] / ws,4)
      daf.stat.theta[int,3] <- round(Theta.stat.Syn[int],0)
      
      #theta for nsyn:
      for(l.s in 1:dim(tab.nsyn)[1]) {
        n.nsam <- tab.nsyn[l.s,2]
        n.freq <- tab.nsyn[l.s,1]
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
    
    ##############################################################
    #calculate the aMKT using the stats
    aa <- NULL
    tryCatch(
      {
        aa <- asymptoticMK(d0=div.syn, d=div.nsyn, xlow=0, xhigh=1, df=as.data.frame(daf.stat.theta), true_alpha=NA, output="table")
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

    #Keep all results:
    Theta.stat.results$scenario[i] <- f
    Theta.stat.results$missing[i]  <- miss.val
    Theta.stat.results[i,3:(3+5-1)]  <- daf.stat.theta[,1]
    Theta.stat.results$Theta.FuLi.Nsyn[i]  <- Theta.stat.Nsyn[1]
    Theta.stat.results$Theta.FuLi.Syn[i]   <- Theta.stat.Syn[1]
    Theta.stat.results$Theta.Watt.Nsyn[i]  <- Theta.stat.Nsyn[2]
    Theta.stat.results$Theta.Watt.Syn[i]   <- Theta.stat.Syn[2]
    Theta.stat.results$Theta.Taji.Nsyn[i]  <- Theta.stat.Nsyn[3]
    Theta.stat.results$Theta.Taji.Syn[i]   <- Theta.stat.Syn[3]
    Theta.stat.results$Theta.FayWu.Nsyn[i] <- Theta.stat.Nsyn[4]
    Theta.stat.results$Theta.FayWu.Syn[i]  <- Theta.stat.Syn[4]
    Theta.stat.results$Theta.vhf.Nsyn[i]   <- Theta.stat.Nsyn[5]
    Theta.stat.results$Theta.vhf.Syn[i]    <- Theta.stat.Syn[5]
    Theta.stat.results$Div.Nsyn[i]         <- div.nsyn
    Theta.stat.results$Div.Syn[i]          <- div.syn
    
    Alpha.theta.stat.results$scenario[i]    <- f
    Alpha.theta.stat.results$missing[i]     <- miss.val
    Alpha.theta.stat.results$alpha.real[i]  <- true.alpha
    Alpha.theta.stat.results$alpha.mkt[i]   <- aa$alpha_original
    Alpha.theta.stat.results$alpha.asym[i]  <- aa$alpha_asymptotic
    Alpha.theta.stat.results$alphaAs.CIL[i] <- aa$CI_low
    Alpha.theta.stat.results$alphaAs.CIH[i] <- aa$CI_high
    
    i <- i + 1
  }
}
write.table(x=Theta.stat.results,file=sprintf("Theta.stats.results.txt"),row.names=F,quote=F)
write.table(x=Alpha.theta.stat.results,file=sprintf("Alpha.theta.stats.results.txt"),row.names=F,quote=F)
