setwd("C:/Users/albaf/Documents/ALBA/ESTUDIS/BIOINFORMÀTICA/TFM/SLiM/SLiM_2/Bootstrap")
#getwd()

#define sample and size of sfs (to avoid zeros)
nsam <- c(25*2) #sample size
#read files from the directory: looking for slim outputs
slim_files <- system("ls *SLiM_02_sum*.txt",intern=T)
#number of bootstrap samples
nboot <- 2

f <- slim_files[1] #in case you want to check without loop
for(f in slim_files) {
  
  #Each SLiM file contains 1 SFS for syn and 1 SFS for nonsyn + divergence and finally the true alpha value
  ##########################
  #keep SFS in "dat.sfs"
  dat.sfs <- read.table(file=f,header=T,row.names=1) #<- read a file with SFS, divergence and true beneficial from slim
  #copy data frame 
  new.dat.sfs <- dat.sfs
  new.dat.sfs[,c(1:(nsam-1),nsam+1,nsam+3)] <- 0 # all values to zero except length
  #new.dat.sfs[1,nsam+3]<-NA
  
  ##########
  #bootstrap
  ##########
  #for each frequency and variant, create a vector of freq values with all positions
  pos.fr.syn <- NULL
  pos.fr.nsyn <- NULL
  for(freq in 1:(nsam-1)) {
    pos.fr.syn <- c(pos.fr.syn,rep(freq,dat.sfs[2,freq]))
    pos.fr.nsyn <- c(pos.fr.nsyn,rep(freq,dat.sfs[1,freq]))
  }
  #for fixed
  pos.fr.syn <- c(pos.fr.syn,rep(nsam,dat.sfs[2,nsam+1]))
  pos.fr.nsyn <- c(pos.fr.nsyn,rep(nsam,dat.sfs[1,nsam+1]))
  #for non-variable positions
  pos.fr.syn <- c(pos.fr.syn,rep(0,dat.sfs[2,nsam+2]-sum(dat.sfs[2,c(1:(nsam-1),nsam+1)])))
  pos.fr.nsyn <- c(pos.fr.nsyn,rep(0,dat.sfs[1,nsam+2]-sum(dat.sfs[1,c(1:(nsam-1),nsam+1)])))
  
  #Sample (bootstrap)
  for(iter in 1:nboot) {
    new.pos.syn  <- sample(1:dat.sfs[2,nsam+2],replace=TRUE)
    new.pos.nsyn <- sample(1:dat.sfs[1,nsam+2],replace=TRUE)
    #recalculate SFS
    for(freq in 1:(nsam-1)) {
      new.dat.sfs[2,freq] <- sum(pos.fr.syn[new.pos.syn]==freq)
      new.dat.sfs[1,freq] <- sum(pos.fr.nsyn[new.pos.nsyn]==freq)
    }
    new.dat.sfs[2,nsam+1] <- sum(pos.fr.syn[new.pos.syn]==nsam)
    new.dat.sfs[1,nsam+1] <- sum(pos.fr.nsyn[new.pos.nsyn]==nsam)
    
    #keep results in a file with name *_slim_SFS*.txt_boot_*.txt
    #BE CAREFUL. If you run again the script and forget delete bootstrap files, all these may be used as sources!!!
    write.table(x=new.dat.sfs,file=sprintf("%s_boot_%03.0f.txt",f,iter),quote=F,sep="\t")
    #check
    #r.new.dat.sfs <- read.table(file=sprintf("%s_boot_%03.0f.txt",f,iter),header=T,row.names=1)
  }
}

