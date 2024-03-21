

install.packages('powerEQTL')
library(powerEQTL)
library(tidyverse)
library(data.table)


##################################
### POWER ANALYSIS 

# sample size
N <- c(100,200,300,400,500)
nn <- length(N)

# MAF
MAF <- seq(0.5,50,0.1)/100
nq <- length(MAF)

nTest=6139877*736 #MAF5% cutoff

# get average SD of module-PC
files <- list.files(pattern="*txtsigpc_sd.gz", full.names=FALSE)
ldf <- list()
for (i in 1:length(files)){
	ldf[i] <- read.table(files[i],header=T)
}
names(ldf) <- files
master <- bind_rows(ldf)
master <- as.data.frame(t(master))
# mean module-PC SD in my own data is 4.850


# obtain power at various values for true slope

power_slr_1 <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
    for (j in 1:nq){
        result <- powerEQTL.SLR(MAF=MAF[j], 
        						slope=1,
                                  nTests=nTest, 
                                  n=N[i], 
                                  power=NULL,
                                  sigma.y=4.850,
                                  FWER=0.05)
        power_slr_1[i,j] <-result;
    }
}

power_slr_2 <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
    for (j in 1:nq){
        result <- powerEQTL.SLR(MAF=MAF[j], 
        						slope=2,
                                  nTests=nTest, 
                                  n=N[i], 
                                  power=NULL,
                                  sigma.y=4.850,
                                  FWER=0.05)
        power_slr_2[i,j] <-result;
    }
}


power_slr_3 <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
    for (j in 1:nq){
        result <- powerEQTL.SLR(MAF=MAF[j], 
                    slope=3,
                                  nTests=nTest, 
                                  n=N[i], 
                                  power=NULL,
                                  sigma.y=4.850,
                                  FWER=0.05)
        power_slr_3[i,j] <-result;
    }
}


power_slr_4 <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
    for (j in 1:nq){
        result <- powerEQTL.SLR(MAF=MAF[j], 
        						slope=4,
                                  nTests=nTest, 
                                  n=N[i], 
                                  power=NULL,
                                  sigma.y=4.850,
                                  FWER=0.05)
        power_slr_4[i,j] <-result;
    }
}

power_slr_5 <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
    for (j in 1:nq){
        result <- powerEQTL.SLR(MAF=MAF[j], 
        						slope=5,
                                  nTests=nTest, 
                                  n=N[i], 
                                  power=NULL,
                                  sigma.y=4.850,
                                  FWER=0.05)
        power_slr_5[i,j] <-result;
    }
}

power_slr_6 <- array(numeric(nn*nq), dim=c(nn,nq))
for (i in 1:nn){
    for (j in 1:nq){
        result <- powerEQTL.SLR(MAF=MAF[j], 
        						slope=6,
                                  nTests=nTest, 
                                  n=N[i], 
                                  power=NULL,
                                  sigma.y=4.850,
                                  FWER=0.05)
        power_slr_6[i,j] <-result;
    }
}



# set up graph
xrange <- range(MAF*100)
yrange <- c(0:1)
colors <- rainbow(length(N))

pdf('~/www/test.pdf',width=7,height=5)

plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main="Sig=0.05, nSNP=6,139,877, Effect Size=1") + abline(v=0, h=seq(0,1,.1), lty=2, col="grey89") + abline(h=0, v=c(1:10), lty=2,col="grey89") + lines(MAF*100, power_slr_1[1,], type="l", lwd=4, col=colors[1])+ lines(MAF*100, power_slr_1[2,], type="l", lwd=4, col=colors[2])+ lines(MAF*100, power_slr_1[3,], type="l", lwd=4, col=colors[3])+ lines(MAF*100, power_slr_1[4,], type="l", lwd=4, col=colors[4])+ lines(MAF*100, power_slr_1[5,], type="l", lwd=4, col=colors[5]) #+ legend("topleft", title="Sample size (n)", as.character(N),fill=colors, bty='n')

plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main="Sig=0.05, nSNP=6,139,877, Effect Size=2") + abline(v=0, h=seq(0,1,.1), lty=2, col="grey89") + abline(h=0, v=c(1:10), lty=2,col="grey89") + lines(MAF*100, power_slr_2[1,], type="l", lwd=4, col=colors[1])+ lines(MAF*100, power_slr_2[2,], type="l", lwd=4, col=colors[2])+ lines(MAF*100, power_slr_2[3,], type="l", lwd=4, col=colors[3])+ lines(MAF*100, power_slr_2[4,], type="l", lwd=4, col=colors[4])+ lines(MAF*100, power_slr_2[5,], type="l", lwd=4, col=colors[5]) #+ le

plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main="Sig=0.05, nSNP=6,139,877, Effect Size=3") + abline(v=0, h=seq(0,1,.1), lty=2, col="grey89") + abline(h=0, v=c(1:10), lty=2,col="grey89") + lines(MAF*100, power_slr_3[1,], type="l", lwd=4, col=colors[1])+ lines(MAF*100, power_slr_3[2,], type="l", lwd=4, col=colors[2])+ lines(MAF*100, power_slr_3[3,], type="l", lwd=4, col=colors[3])+ lines(MAF*100, power_slr_3[4,], type="l", lwd=4, col=colors[4])+ lines(MAF*100, power_slr_3[5,], type="l", lwd=4, col=colors[5]) #+ le

plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main="Sig=0.05, nSNP=6,139,877, Effect Size=4") + abline(v=0, h=seq(0,1,.1), lty=2, col="grey89") + abline(h=0, v=c(1:10), lty=2,col="grey89") + lines(MAF*100, power_slr_4[1,], type="l", lwd=4, col=colors[1])+ lines(MAF*100, power_slr_4[2,], type="l", lwd=4, col=colors[2])+ lines(MAF*100, power_slr_4[3,], type="l", lwd=4, col=colors[3])+ lines(MAF*100, power_slr_4[4,], type="l", lwd=4, col=colors[4])+ lines(MAF*100, power_slr_4[5,], type="l", lwd=4, col=colors[5]) #+ legend("topleft", title="Sample size (n)", as.character(N),fill=colors, bty='n')

plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main="Sig=0.05, nSNP=6,139,877, Effect Size=5") + abline(v=0, h=seq(0,1,.1), lty=2, col="grey89") + abline(h=0, v=c(1:10), lty=2,col="grey89") + lines(MAF*100, power_slr_5[1,], type="l", lwd=4, col=colors[1])+ lines(MAF*100, power_slr_5[2,], type="l", lwd=4, col=colors[2])+ lines(MAF*100, power_slr_5[3,], type="l", lwd=4, col=colors[3])+ lines(MAF*100, power_slr_5[4,], type="l", lwd=4, col=colors[4])+ lines(MAF*100, power_slr_5[5,], type="l", lwd=4, col=colors[5]) #+ legend("topleft", title="Sample size (n)", as.character(N),fill=colors, bty='n')

plot(xrange, yrange, log='x', type="n",
     xlab="MAF (%)",
     ylab="Power",
     main="Sig=0.05, nSNP=6,139,877, Effect Size=6") + abline(v=0, h=seq(0,1,.1), lty=2, col="grey89") + abline(h=0, v=c(1:10), lty=2,col="grey89") + lines(MAF*100, power_slr_6[1,], type="l", lwd=4, col=colors[1])+ lines(MAF*100, power_slr_6[2,], type="l", lwd=4, col=colors[2])+ lines(MAF*100, power_slr_6[3,], type="l", lwd=4, col=colors[3])+ lines(MAF*100, power_slr_6[4,], type="l", lwd=4, col=colors[4])+ lines(MAF*100, power_slr_6[5,], type="l", lwd=4, col=colors[5]) #+ legend("topleft", title="Sample size (n)", as.character(N),fill=colors, bty='n')


dev.off()





