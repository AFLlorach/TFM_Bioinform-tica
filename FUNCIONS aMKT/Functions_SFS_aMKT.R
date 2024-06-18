setwd("./")
#getwd()

source("./asymptoticMK_local.R")

#define sample and size of sfs (to avoid zeros)
nsam <- c(25*2) #sample size
#read files from the directory: looking for slim outputs
slim_files <- system("ls SLiM_*_sum*.txt",intern=T)

#Define data frame for the asymptoticMKT results: 
daf.results <- array(0,dim=c(length(slim_files),nsam*2+1))
colnames(daf.results) <- c("scenario",sprintf("Pi.Syn%.0f",1:(nsam-1)),sprintf("Pi.NSyn%.0f",1:(nsam-1)),"Div.Syn","Div.Nsyn")
daf.results <- as.data.frame(daf.results)

Alpha.daf.results <- array(0,dim=c(length(slim_files),6))
colnames(Alpha.daf.results) <- c("scenario","alpha.real","alpha.mkt","alpha.asym","alphaAs.CIL","alphaAs.CIH")
Alpha.daf.results <- as.data.frame(Alpha.daf.results)

############################
#RUN: for each slim file calculate the variability and the alpha, considering missing proportions
############################
i <- 1
f <- slim_files[1] #in case you want to check without loop
for(f in slim_files) {
  #Each SLiM file contains 1 SFS for syn and 1 SFS for nonsyn + divergence and finally the true alpha value
  ##########################
  #keep SFS in "dat.sfs"
  dat.sfs <- as.matrix(read.table(file=f,header=T,row.names=1)) #<- read a file with SFS, divergence and true beneficial from slim
  #keep divergence in "divergence"
  divergence <- data.frame(mi=dat.sfs[1,nsam+2],Di=dat.sfs[1,nsam+1],m0=dat.sfs[2,nsam+2],D0=dat.sfs[2,nsam+1])
  
  #make the data frame daf (f,Pi,Po)
  daf <- data.frame(f=c(1:(nsam-1)/nsam),p=dat.sfs[1,c(1:(nsam-1))],p0=dat.sfs[2,c(1:(nsam-1))])
    
  ##############################################################
  #calculate the aMKT using the stats
  aa <- NULL
  tryCatch(
    {
      aa <- asymptoticMK(d0=as.numeric(divergence[4]), d=as.numeric(divergence[2]), xlow=0, xhigh=1, df=daf, true_alpha=NA, output="table")
    },
    error = function(e) {
      message(sprintf("Error calculating MKTa for observed data in file %s",f))
    }
  )
  #aa
  #aa$alpha_asymptotic
  #aa$alpha_original
  
  #TO DO: 
  #true.alpha
  #sum beneficial / total.nsyn
  true.alpha <- dat.sfs[1,nsam+3] / dat.sfs[1,nsam+1]
  
  #TO DO: 
  #Keep all results:
  daf.results$scenario[i] <- f
  daf.results[i,2:(nsam)] <- daf$p0
  daf.results[i,(nsam+1):(2*nsam-1)] <- daf$p
  daf.results[i,(2*nsam)] <- as.numeric(divergence[4])
  daf.results[i,(2*nsam+1)] <- as.numeric(divergence[2])
  
  Alpha.daf.results$scenario[i]    <- f
  Alpha.daf.results$alpha.real[i]  <- true.alpha
  Alpha.daf.results$alpha.mkt[i]   <- aa$alpha_original
  Alpha.daf.results$alpha.asym[i]  <- aa$alpha_asymptotic
  Alpha.daf.results$alphaAs.CIL[i] <- aa$CI_low
  Alpha.daf.results$alphaAs.CIH[i] <- aa$CI_high
  
  i <- i + 1
  
}
write.table(x=daf.results,file=sprintf("daf.results.txt"),row.names=F,quote=F)
write.table(x=Alpha.daf.results,file=sprintf("Alpha.daf.results.txt"),row.names=F,quote=F)
