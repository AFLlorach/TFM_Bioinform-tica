#calculation of alpha from standard equation
calc.alpha <- function(Ps,Ds,Pn,Dn) {
  ratiok <- Ds/Dn
  ratiop <- Ps/Pn
  alpha <- 1 - ratiok/ratiop
  return(alpha)
} 

setwd("./")
setwd("/Users/sramos/Desktop/TFM-UOC-Proposals-2023/NEW_simulations_aMKT/DADES SLiM/")

pdf("Plots_MKT_Bootstrap.pdf")
nsam <- 50
nmodels <- c(0:7)
miss.val <- c(0,0.1,0.25,0.5,0.75,0.9)
for(m in nmodels) {
  #path.res <- sprintf("./SLiM_%.0f/Boostrap_%02.0f/Results_SLiM_%02.0f",m,m,m)
  path.res <- sprintf("./SLiM_%.0f/Boostrap_%02.0f",m,m);#/Results_SLiM_%02.0f",m,m,m)
  
  # alpha.daf <- read.table(file=sprintf("%s/Alpha.daf.results_%02.0f.txt",path.res,m),header=T)
  # alpha.int <- read.table(file=sprintf("%s/Alpha.theta.int.results_%02.0f.txt",path.res,m),header=T)
  # alpha.stt <- read.table(file=sprintf("%s/Alpha.theta.stats.results_%02.0f.txt",path.res,m),header=T)
  # 
  # daf.data <- read.table(file=sprintf("%s/daf.results_%02.0f.txt",path.res,m),header=T)
  # int.data <- read.table(file=sprintf("%s/Theta.int.results_%02.0f.txt",path.res,m),header=T)
  # stt.data <- read.table(file=sprintf("%s/Theta.stats.results_%02.0f.txt",path.res,m),header=T)
  
  alpha.daf <- read.table(file=sprintf("%s/Alpha.daf.results.txt",path.res),header=T)
  alpha.int <- read.table(file=sprintf("%s/Alpha.theta.int.results.txt",path.res),header=T)
  alpha.stt <- read.table(file=sprintf("%s/Alpha.theta.stats.results.txt",path.res),header=T)
  
  daf.data <- read.table(file=sprintf("%s/daf.results.txt",path.res),header=T)
  int.data <- read.table(file=sprintf("%s/Theta.int.results.txt",path.res),header=T)
  stt.data <- read.table(file=sprintf("%s/Theta.stats.results.txt",path.res),header=T)
  
  #calculate alpha from SFS
  cur.alpha <- calc.alpha(Ps=daf.data[,c(2:(nsam))],Ds=daf.data[,2*nsam],Pn=daf.data[,c((nsam+1):(2*nsam-1))],Dn=daf.data[,2*nsam+1])
  colnames(cur.alpha) <- sprintf("alpha_%.0f",c(1:(nsam-1)))
  cur.int.alpha <- cbind(int.data[,1:7],calc.alpha(Ps=int.data[,seq(8,17,2)],Ds=int.data[,18],Pn=int.data[,seq(9,17,2)],Dn=int.data[,19]))
  colnames(cur.int.alpha) <- c(colnames(int.data)[1:7],sprintf("alpha_%.0f",c(1:5)))
  cur.stt.alpha <- cbind(stt.data[,1:7],calc.alpha(Ps=stt.data[,seq(8,17,2)],Ds=stt.data[,18],Pn=stt.data[,seq(9,17,2)],Dn=stt.data[,19]))
  colnames(cur.stt.alpha) <- c(colnames(stt.data)[1:7],sprintf("alpha_%.0f",c(1:5)))
  
  #plot first (sfs and sfs*freq)
  par(mfrow=c(2,2))
  plot(x=c(1:(nsam-1)),y=daf.data[1,c(2:(nsam))],xlab="freq",ylab="SFS.Syn",type="l",ylim=c(0,max(daf.data[1,c(2:(nsam))])),main=sprintf("Model_%.0f",m))
  plot(x=c(1:(nsam-1)),y=daf.data[1,c((nsam+1):(2*nsam-1))],xlab="freq",ylab="SFS.NSyn",type="l",ylim=c(0,max(daf.data[1,c((nsam+1):(2*nsam-1))])),main=sprintf("Model_%.0f",m))
  plot(x=c(1:(nsam-1)),y=daf.data[1,c(2:(nsam))]*c(1:(nsam-1)),xlab="freq",ylab="i*SFS.Syn",type="l",ylim=c(0,max(daf.data[1,c(2:(nsam))]*c(1:(nsam-1)))),main=sprintf("Model_%.0f",m))
  plot(x=c(1:(nsam-1)),y=daf.data[1,c((nsam+1):(2*nsam-1))]*c(1:(nsam-1)),xlab="freq",ylab="i*SFS.NSyn",type="l",ylim=c(0,max(daf.data[1,c((nsam+1):(2*nsam-1))]*c(1:(nsam-1)))),main=sprintf("Model_%.0f",m))
  
  #Plot results: alpha from daf
  par(mfrow=c(2,1))
  
  plot(x=c(1:(nsam-1))/nsam,y=cur.alpha[1,],ylim=c(-1,1),type="b",pch=20,main=sprintf("Model_%.0f: alpha (SFS)",m),xlab="freq",ylab="alpha")
  ci <- apply(cur.alpha[-1,],2,quantile,probs=c(0.025,0.975)) #confidence interval and plot
  segments(x0=as.numeric(c(1:(nsam-1))/nsam),x1=as.numeric(c(1:(nsam-1))/nsam),y0=as.numeric(cur.alpha[1,]),y1=as.numeric(ci[1,]))
  segments(x0=as.numeric(c(1:(nsam-1))/nsam),x1=as.numeric(c(1:(nsam-1))/nsam),y0=as.numeric(cur.alpha[1,]),y1=as.numeric(ci[2,]))
  abline(h=alpha.daf$alpha.real[1],col="red",lwd=1)
  abline(h=alpha.daf$alpha.asym[1],col="blue",lwd=1)
  lines(x=1,y=alpha.daf$alpha.asym[1],type="p",pch=20,col="blue")
  ci.asy <-quantile(alpha.daf$alpha.asym[-1],probs=c(0.025,0.975))
  segments(x0=1,x1=1,y0=alpha.daf$alpha.asym[1],y1=ci.asy[1],col="blue")
  segments(x0=1,x1=1,y0=alpha.daf$alpha.asym[1],y1=ci.asy[2],col="blue")
  #legend("topleft",legend=c("True","Asymptotic"),lwd=c(1,1),col=c("red","blue"))
  
  da <- density(alpha.daf$alpha.asym[-1])
  plot(da,main=sprintf("Model_%.0f (SFS)\n Asymptotic alpha",m),xlim=c(-1,1))
  abline(v=alpha.daf$alpha.real[1],col="red",lwd=1)
  abline(v=alpha.daf$alpha.asym[1],col="blue",lwd=1)
  legend("topleft",legend=c("True","Asymptotic"),lwd=c(1,1),col=c("red","blue"))
  
  #plot with zoom
  plot(x=c(1:(nsam-1))/nsam,y=cur.alpha[1,],type="b",pch=20,main=sprintf("Model_%.0f: alpha (SFS)",m),
       ylim=c(min(ci,alpha.daf$alpha.real[1],alpha.daf$alpha.asym[1]),max(ci,alpha.daf$alpha.real[1],alpha.daf$alpha.asym[1])),
       xlab="freq",ylab="alpha",)
  segments(x0=as.numeric(c(1:(nsam-1))/nsam),x1=as.numeric(c(1:(nsam-1))/nsam),y0=as.numeric(cur.alpha[1,]),y1=as.numeric(ci[1,]))
  segments(x0=as.numeric(c(1:(nsam-1))/nsam),x1=as.numeric(c(1:(nsam-1))/nsam),y0=as.numeric(cur.alpha[1,]),y1=as.numeric(ci[2,]))
  abline(h=alpha.daf$alpha.real[1],col="red",lwd=1)
  abline(h=alpha.daf$alpha.asym[1],col="blue",lwd=1)
  lines(x=1,y=alpha.daf$alpha.asym[1],type="p",pch=20,col="blue")
  segments(x0=1,x1=1,y0=alpha.daf$alpha.asym[1],y1=ci.asy[1],col="blue")
  segments(x0=1,x1=1,y0=alpha.daf$alpha.asym[1],y1=ci.asy[2],col="blue")

  plot(da,main=sprintf("Model_%.0f (SFS)\n Asymptotic alpha",m),xlim=c(min(da$x,alpha.daf$alpha.real[1],alpha.daf$alpha.asym[1]),max(da$x,alpha.daf$alpha.real[1],alpha.daf$alpha.asym[1])))
  abline(v=alpha.daf$alpha.real[1],col="red",lwd=1)
  abline(v=alpha.daf$alpha.asym[1],col="blue",lwd=1)
  legend("topleft",legend=c("True","Asymptotic"),lwd=c(1,1),col=c("red","blue"))
  
  
  for(miss in miss.val) {
    #Plot results: alpha interval theta stats
    row.m <- which(cur.int.alpha[,2]==miss)
    ci <- apply(cur.int.alpha[row.m[-1],8:12],2,quantile,probs=c(0.025,0.975))
    plot(x=as.numeric(cur.int.alpha[row.m[1],3:7]),y=as.numeric(cur.int.alpha[row.m[1],8:12]),
         ylim=c(-1,1),xlim=c(0,1),type="b",pch=20,
         main=sprintf("Model_%.0f: Missing %.2f \n alpha (Theta Intervals)",m,miss),
         xlab="freq",ylab="alpha")
    segments(x0=as.numeric(apply(cur.int.alpha[row.m,3:7],2,mean)),x1=as.numeric(apply(cur.int.alpha[row.m,3:7],2,mean)),y0=as.numeric(cur.int.alpha[row.m[1],8:12]),y1=as.numeric(ci[1,]))
    segments(x0=as.numeric(apply(cur.int.alpha[row.m,3:7],2,mean)),x1=as.numeric(apply(cur.int.alpha[row.m,3:7],2,mean)),y0=as.numeric(cur.int.alpha[row.m[1],8:12]),y1=as.numeric(ci[2,]))
    abline(h=alpha.daf$alpha.real[1],col="red",lwd=1)
    abline(h=alpha.int$alpha.asym[row.m[1]],col="blue",lwd=1)
    lines(1,alpha.int$alpha.asym[row.m[1]],pch=20,type="p",col="blue")
    ci.asy <-quantile(alpha.int$alpha.asym[row.m[-1]],probs=c(0.025,0.975))
    segments(x0=1,x1=1,y0=alpha.int$alpha.asym[row.m[1]],y1=ci.asy[1],col="blue")
    segments(x0=1,x1=1,y0=alpha.int$alpha.asym[row.m[1]],y1=ci.asy[2],col="blue")

    #plot with zoom 
    da <- density(alpha.int$alpha.asym[row.m[-1]])
    plot(da,main=sprintf("Model_%.0f (Theta Intervals)\n Asymptotic alpha",m),xlim=c(-1,1))
    abline(v=alpha.daf$alpha.real[1],col="red",lwd=1)
    abline(v=alpha.int$alpha.asym[row.m[1]],col="blue",lwd=1)
    legend("topleft",legend=c("True","Asymptotic"),lwd=c(1,1),col=c("red","blue"))

    plot(x=as.numeric(cur.int.alpha[row.m[1],3:7]),y=as.numeric(cur.int.alpha[row.m[1],8:12]),
         ylim=c(min(ci,alpha.daf$alpha.real[1],alpha.int$alpha.asym[row.m[1]]),max(ci,alpha.daf$alpha.real[1],alpha.int$alpha.asym[row.m[1]])),
         xlim=c(0,1),type="b",pch=20,
         main=sprintf("Model_%.0f: Missing %.2f \n alpha (Theta Intervals)",m,miss),
         xlab="freq",ylab="alpha")
    segments(x0=as.numeric(apply(cur.int.alpha[row.m,3:7],2,mean)),x1=as.numeric(apply(cur.int.alpha[row.m,3:7],2,mean)),y0=as.numeric(cur.int.alpha[row.m[1],8:12]),y1=as.numeric(ci[1,]))
    segments(x0=as.numeric(apply(cur.int.alpha[row.m,3:7],2,mean)),x1=as.numeric(apply(cur.int.alpha[row.m,3:7],2,mean)),y0=as.numeric(cur.int.alpha[row.m[1],8:12]),y1=as.numeric(ci[2,]))
    abline(h=alpha.daf$alpha.real[1],col="red",lwd=1)
    abline(h=alpha.int$alpha.asym[row.m[1]],col="blue",lwd=1)
    lines(1,alpha.int$alpha.asym[row.m[1]],pch=20,type="p",col="blue")
    segments(x0=1,x1=1,y0=alpha.int$alpha.asym[row.m[1]],y1=ci.asy[1],col="blue")
    segments(x0=1,x1=1,y0=alpha.int$alpha.asym[row.m[1]],y1=ci.asy[2],col="blue")
    
    plot(da,main=sprintf("Model_%.0f (Theta Intervals)\n Asymptotic alpha",m),xlim=c(min(da$x,alpha.daf$alpha.real[1],alpha.int$alpha.asym[row.m[1]]),max(da$x,alpha.daf$alpha.real[1],alpha.int$alpha.asym[row.m[1]])))
    abline(v=alpha.daf$alpha.real[1],col="red",lwd=1)
    abline(v=alpha.int$alpha.asym[row.m[1]],col="blue",lwd=1)
    legend("topleft",legend=c("True","Asymptotic"),lwd=c(1,1),col=c("red","blue"))
    
  }
  for(miss in miss.val) {
    #Plot results: alpha from classical theta stats
    row.m <- which(cur.stt.alpha[,2]==miss)
    plot(x=as.numeric(cur.stt.alpha[row.m[1],3:7]),y=as.numeric(cur.stt.alpha[row.m[1],8:12]),
         ylim=c(-1,1),xlim=c(0,1),type="b",pch=20,
         main=sprintf("Model_%.0f: Missing %.2f \n alpha (Theta Stats)",m,miss),
         xlab="freq",ylab="alpha")
    ci <- apply(cur.stt.alpha[row.m[-1],8:12],2,quantile,probs=c(0.025,0.975))
    segments(x0=apply(cur.stt.alpha[row.m,3:7],2,mean),x1=apply(cur.stt.alpha[row.m,3:7],2,mean),y0=as.numeric(cur.stt.alpha[row.m[1],8:12]),y1=ci[1,])
    segments(x0=apply(cur.stt.alpha[row.m,3:7],2,mean),x1=apply(cur.stt.alpha[row.m,3:7],2,mean),y0=as.numeric(cur.stt.alpha[row.m[1],8:12]),y1=ci[2,])
    abline(h=alpha.daf$alpha.real[1],col="red",lwd=1)
    abline(h=alpha.stt$alpha.asym[row.m[1]],col="blue",lwd=1)
    lines(1,alpha.stt$alpha.asym[row.m[1]],pch=20,type="p",col="blue")
    ci.asy <-quantile(alpha.stt$alpha.asym[row.m[-1]],probs=c(0.025,0.975))
    segments(x0=1,x1=1,y0=alpha.stt$alpha.asym[row.m[1]],y1=ci.asy[1],col="blue")
    segments(x0=1,x1=1,y0=alpha.stt$alpha.asym[row.m[1]],y1=ci.asy[2],col="blue")

    da <- density(alpha.stt$alpha.asym[row.m[-1]])
    plot(da,main=sprintf("Model_%.0f (Theta Stats)\n Asymptotic alpha",m),xlim=c(-1,1))
    abline(v=alpha.daf$alpha.real[1],col="red",lwd=1)
    abline(v=alpha.stt$alpha.asym[row.m[1]],col="blue",lwd=1)
    legend("topleft",legend=c("True","Asymptotic"),lwd=c(1,1),col=c("red","blue"))
    
    #plot with zoom 
    plot(x=as.numeric(cur.stt.alpha[row.m[1],3:7]),y=as.numeric(cur.stt.alpha[row.m[1],8:12]),
         ylim=c(min(ci,alpha.daf$alpha.real[1],alpha.stt$alpha.asym[row.m[1]]),max(ci,alpha.daf$alpha.real[1],alpha.stt$alpha.asym[row.m[1]])),
         xlim=c(0,1),type="b",pch=20,
         main=sprintf("Model_%.0f: Missing %.2f \n alpha (Theta Stats)",m,miss),
         xlab="freq",ylab="alpha")
    segments(x0=apply(cur.stt.alpha[row.m,3:7],2,mean),x1=apply(cur.stt.alpha[row.m,3:7],2,mean),y0=as.numeric(cur.stt.alpha[row.m[1],8:12]),y1=ci[1,])
    segments(x0=apply(cur.stt.alpha[row.m,3:7],2,mean),x1=apply(cur.stt.alpha[row.m,3:7],2,mean),y0=as.numeric(cur.stt.alpha[row.m[1],8:12]),y1=ci[2,])
    abline(h=alpha.daf$alpha.real[1],col="red",lwd=1)
    abline(h=alpha.stt$alpha.asym[row.m[1]],col="blue",lwd=1)
    lines(1,alpha.stt$alpha.asym[row.m[1]],pch=20,type="p",col="blue")
    segments(x0=1,x1=1,y0=alpha.stt$alpha.asym[row.m[1]],y1=ci.asy[1],col="blue")
    segments(x0=1,x1=1,y0=alpha.stt$alpha.asym[row.m[1]],y1=ci.asy[2],col="blue")
    
    plot(da,main=sprintf("Model_%.0f (Theta Stats)\n Asymptotic alpha",m),xlim=c(min(da$x,alpha.daf$alpha.real[1],alpha.stt$alpha.asym[row.m[1]]),max(da$x,alpha.daf$alpha.real[1],alpha.stt$alpha.asym[row.m[1]])))
    abline(v=alpha.daf$alpha.real[1],col="red",lwd=1)
    abline(v=alpha.stt$alpha.asym[row.m[1]],col="blue",lwd=1)
    legend("topleft",legend=c("True","Asymptotic"),lwd=c(1,1),col=c("red","blue"))
  }
}
dev.off()

