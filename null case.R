library(deconvolveR)
data(surg)
surg
## KS two-sample tes
P_value<-rep(0,1000)
T_hat<-rep(0,1000)
tau = seq(0,1,length = 100)
m <- length(tau)
g<-matrix(0,nrow=m,ncol=1000)
f<-matrix(0,nrow=m,ncol=1000)
g_pool<-matrix(0,nrow=m,ncol=1000)
B <- 99
library(gplots)
set.seed(123456)

ptm <- proc.time()
for (i in 1:1000){
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
  ## Under null hypothesis: pull X and Y together to estimate g(theta)
  sample = rbind(dataX, dataY) 
  results = deconv(tau = tau, X = sample, family = "Binomial",  c0 = 1)
  g_pool[,i]<- results$stats[, "g"]
  NX = length(dataX[,1])
  NY = length(dataY[,1])
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
                       function(k) {
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
  P_value[i]<-(length(which(T_bsp>=T_hat))+1)/(B+1)
  }
proc.time() - ptm

hist(P_value)
library(qqtest)
qqtest(P_value,dist="uniform")


