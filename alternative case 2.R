library(deconvolveR)
data(surg)
surg
## KS two-sample tes
nSIM<-100
T_hat<-rep(0,nSIM)
tau = seq(0,1,length = 100)
m <- length(tau)
g<-matrix(0,nrow=m,ncol=nSIM)
f<-matrix(0,nrow=m,ncol=nSIM)
g_pool<-matrix(0,nrow=m,ncol=nSIM)
B <- 99
w <- seq(0,0.09,0.01)
power<-rep(0,length(w))
P_value<-matrix(0,nrow=nSIM,ncol=length(w))
library(gplots)
set.seed(123456)

ptm <- proc.time()
for ( k in 1: length(w)){
  for (i in 1:nSIM){
    index = sample(c(rep(1,length(surg$n)/2), rep(0,length(surg$n)/2)),length(surg$n))
    G1 = surg[which(index==1),]
    G2 = surg[which(index==0),]
    dataX <- G1
    resultX = deconv(tau = tau, X = dataX, family = "Binomial", c0 = 1)
    g[,i]<-resultX$stats[, "g"]  ## estimated density function
    G_est<-resultX$stats[, "G"] ## estimated cumulative density function
    #unstim
    dataY <- G2
    resultY = deconv(tau = tau, X = dataY, family = "Binomial",  c0 = 1)
    f[,i]<-resultY$stats[, "g"] ## estimated density function
    F_est<-resultY$stats[, "G"] ## estimated cumulative density function
    NX = length(dataX[,1])
    NY = length(dataY[,1])
    ## alternative
    Theta_est <- sample(tau, size = NX, replace = TRUE, prob = f[,i])
    Theta_new <- (1-w[k])*Theta_est+w[k]*ifelse(Theta_est>0,1,0)
    count_alt <- rbinom(NY,size = G2[,1], prob = Theta_new)
    dataY <- data.frame(cbind(G2[,1],count_alt))
    colnames(dataY)<-c("n","s")
    resultY = deconv(tau = tau, X = dataY, family = "Binomial",  c0 = 1)
    f[,i]<-resultY$stats[, "g"] ## estimated density function
    F_est<-resultY$stats[, "G"] ## estimated cumulative density function
    ## Under null hypothesis: pull X and Y together to estimate g(theta)
    sample = rbind(dataX, dataY) 
    results = deconv(tau = tau, X = sample, family = "Binomial",  c0 = 1)
    g_pool[,i]<- results$stats[, "g"]
    Difference_hat = abs(G_est - F_est)
    T_hat = max(Difference_hat)
    ### bootstrap
    
    genBootSampleX <- function() {
      thetaXStar = sample(tau, size = NX, replace = TRUE, prob = g_pool[,i])
      xStar <- sapply(seq_len(NX),
                      function(j){
                        rbinom(n = 1, prob = thetaXStar[j], size = dataX[j,1])
                      })
      data.frame(cbind(dataX[,1], xStar))
    }
    
    genBootSampleY <- function() {
      taoStar <- sample(tau, size = NY, replace = TRUE, prob = g_pool[,i])
      yStar <- sapply(seq_len(NY),
                      function(j){
                        rbinom(n = 1, prob = taoStar[j], size = dataY[j,1])
                      })
      data.frame(cbind(dataY[,1], yStar))
    }
    
    bootResults = lapply(seq_len(B),
                         function(j) {
                           XBoot = genBootSampleX()
                           resultXboot = deconv(tau = tau, X = XBoot, family = "Binomial",  c0 = 1)
                           gBoot = resultXboot$stats[, "g"]
                           GBoot = resultXboot$stats[, "G"]
                           YBoot <- genBootSampleY()
                           resultYboot = deconv(tau = tau, X = YBoot, family = "Binomial",  c0 = 1)
                           fBoot = resultYboot$stats[, "g"]
                           FBoot <- resultYboot$stats[, "G"]
                           list(GBoot=GBoot, FBoot=FBoot)
                         })
    GBoot <- t(data.frame(sapply(bootResults, function(x) x[1])))
    FBoot <- t(data.frame(sapply(bootResults, function(x) x[2])))
    Difference = abs(GBoot-FBoot)
    T_bsp = apply(Difference, 1, max)
    P_value[i,k]<-(length(which(T_bsp>=T_hat))+1)/(B+1)
  }
  power[k]<-length(which(P_value[,k]<0.05))/nSIM
}
proc.time() - ptm

## power plot

power = round(power, digits=2)
pdf("Alternative case 2 100 simulation 99 bootstrap power plot w=seq(0,0.08,0.01).pdf")
#power[1] = 0.05 
plot(x = w, y = power, type = "l",
     xlab = "Change in combination rate", ylab = "Power",
     main = "Power plot")
text(w, power, power, adj = c(0.3, -0.5))
dev.off()