numero<-c(1:53)
syn<-rep(c(0), each=53)
nsyn<-rep(c(0), each=53)

sumatori<-data.frame(numero,nsyn,syn)
View(sumatori)

setwd("C:/Users/albaf/Documents/SLiM/SLiM_7")

files <- list.files(path=".", pattern="*.txt", full.names=TRUE, recursive=FALSE)
# sum
for (i in files) {
  fil <- read.table(i, sep="", header = FALSE, row.names = NULL) # load file
  gir <- data.frame(t(fil), row.names = NULL) #change rows to columns
  gir= gir[-1, ]# delete the first row
  row.names(gir)<-NULL
  colnames(gir)[1] <- 'fr' 
  colnames(gir)[2] <- 'nsyn'
  colnames(gir)[3] <- 'syn'
  gir$nsyn<-as.numeric(gir$nsyn)
  gir$syn<-as.numeric(gir$syn)
  for (i in 1:53) {
    sumatori$nsyn[i]<-sum(sumatori$nsyn[i],gir$nsyn[i])
    sumatori$syn[i]<-sum(sumatori$syn[i],gir$syn[i])
  }
}
str(gir)
View(sumatori)

Sm <- data.frame(t(sumatori), row.names = NULL)
View(Sm)

