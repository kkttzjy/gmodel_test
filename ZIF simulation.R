set.seed(238923)
N <- 300           ##  set # of data points
nSIM <- 100      ##  set # of simulations
tau <- seq(1, 32)   ##  sample space for unobserved parameter theta
P_value=rep(0,nSIM) ##  initialize P-value matrix
Tstat=rep(0,nSIM) ##  initialize LRT statistics
piX<-0.5
piY<-0.5
B <- 99
m <- length(tau)
g<-matrix(0,nrow=m,ncol=nSIM)
f<-matrix(0,nrow=m,ncol=nSIM)
g_pool<-matrix(0,nrow=m,ncol=nSIM)
pi_X<-rep(0,nSIM)
pi_Y<-rep(0,nSIM)
pi_pool<-rep(0,nSIM)
g_b<-matrix(0,nrow=m,ncol=B)
f_b<-matrix(0,nrow=m,ncol=B)
g_poolb<-matrix(0,nrow=m,ncol=B)
pi_Xb<-rep(0,B)
pi_Yb<-rep(0,B)
pi_poolb<-rep(0,B)
library("VGAM")

T_hat = rep(0,nSIM)
T_asy = rep(0,nSIM)
P_value_ASY=rep(0,nSIM)
P_value_KS=rep(0,nSIM)
P_value1=rep(0,nSIM)
P_value2=rep(0,nSIM)
N_asy = B


#library(deconvolveR)
library(MASS)

for ( i in 1:nSIM){
  #X
  Theta <- rchisq(n=N,  df = 10)
  ## generate d
  dX <- runif(n=N,  min=0.5, max=1)
  ## generate X
  dataX <- rzipois(N, lambda = Theta*dX, pstr0 = piX)
  ## estimate the dist of theta (g(theta)) & need r function deconv_integrate
  resultX <- deconv_new(tau=tau,X=dataX, d=dX, c0=1)
  pi_X[i]<-resultX$pi
  g[,i]<-resultX$stats[, "g"]  ## estimated density function
  G_est<-resultX$stats[, "G"] ## estimated cumulative density function
  X.bias = resultX$bias.all
  X.Cov= resultX$Cov.all

  #Y
  Tao <-rchisq(n=N,  df = 10)
  ## generate d
  dY<-runif(n=N,  min=0.5, max=1)
  ## sample data Y from Poisson dist with mean=tao*dY
  dataY <- rzipois(N, lambda = Tao*dY, pstr0 = piY)
  ## estimate the dist of theta (g(theta)) & need r function deconv_integrate
  resultY<-deconv_new(tau=tau,X=dataY, d=dY, c0=1)
  pi_Y[i]<-resultY$pi
  f[,i]<-resultY$stats[, "g"] ## estimated density function
  F_est<-resultY$stats[, "G"] ## estimated cumulative density function
  Y.bias = resultY$bias.all
  Y.Cov= resultY$Cov.all

  ## Under null hypothesis: pull X and Y together to estimate g(theta)
  sample <- c(dataX, dataY)
  d<-c(dX, dY)
  results <- deconv_new(tau=tau,X=sample,d=d,c0=1)
  pi_pool[i]<-results$pi
  g_pool[,i]<- results$stats[, "g"]

  ### Accelerated bootstrap
  alpha_pool <- results$mle
  resultsX<-deconv_new(tau=tau,X=dataX, d=dX, c0=1, mle = alpha_pool)
  resultsY<-deconv_new(tau=tau,X=dataY, d=dY, c0=1, mle = alpha_pool)

  Pool.bias.X = resultsX$bias.all
  Pool.bias.Y = resultsY$bias.all

  Pool.Cov.X= resultsX$Cov.all
  Pool.Cov.Y= resultsY$Cov.all


  pi_diff_hat = pi_X[i] - pi_Y[i]
  bias.pi_diff = X.bias[1] - Y.bias[1]
  E_hat = G_est - F_est
  bias.E = X.bias[-1] - Y.bias[-1]
  Difference_hat <- abs(G_est-F_est)
  Difference_pi <- abs(pi_X[i]-pi_Y[i])
  T_hat[i]=max(max(Difference_hat),Difference_pi)
  bias.dff = Pool.bias.X-Pool.bias.Y
  Cov.dff = Pool.Cov.X+Pool.Cov.Y
  All_asy = abs(mvrnorm(N_asy, bias.dff, Cov.dff, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
  pi_diff_asy = All_asy[,1]
  T_pre = cbind(pi_diff_asy, apply(All_asy[,-1], 1, max))
  T_asy = apply(T_pre, 1, max)
  P_value_ASY[i]<-(length(which(T_asy>=T_hat[i]))+1)/(N_asy+1)

  ### Simple bootstrap (KS)
  genBootSampleX <- function() {
    thetaStar <- sample(tau, size = N, replace = TRUE, prob = g_pool[,i])
    xStar <- sapply(seq_len(N),
                    function(j){
                      dX_boot<-dX[j]
                      dataX_boot<-rzipois(n = 1, lambda=thetaStar[j]*dX_boot, pstr0 = pi_pool[i])
                    })
  }

  genBootSampleY <- function() {
    taoStar <- sample(tau, size = N, replace = TRUE, prob = g_pool[,i])
    yStar <- sapply(seq_len(N),
                    function(j){
                      dY_boot<-dY[j]
                      dataY_boot<-rzipois(n = 1, lambda=taoStar[j]*dY_boot, pstr0 = pi_pool[i])})
  }
  bootResults <-lapply(seq_len(B),
                       function(k) {
                         XBoot <- genBootSampleX()
                         resultXboot <- deconv_new(tau = tau, X = XBoot, d=dX)
                         GBoot <- resultXboot$stats[, "G"]
                         g_boot <-resultXboot$stats[, "g"]
                         pi_Xboot<-resultXboot$pi
                         YBoot <- genBootSampleY()
                         resultYboot <- deconv_new(tau = tau, X = YBoot, d=dY)
                         FBoot <- resultYboot$stats[, "G"]
                         f_boot <-resultYboot$stats[, "g"]
                         pi_Yboot<-resultYboot$pi
                         sampleBoot <- c(XBoot, YBoot)
                         d_boot <- c(dX, dY)
                         resultsboot <- deconv_new(tau = tau, X = sampleBoot, d=d_boot)
                         GpBoot <- resultsboot$stats[, "G"]
                         g_poolboot<- resultsboot$stats[, "g"]
                         pi_poolboot<-resultsboot$pi
                         list(GBoot=GBoot, FBoot=FBoot, GpBoot=GpBoot,
                              g_boot=g_boot, pi_Xboot=pi_Xboot, f_boot=f_boot, pi_Yboot=pi_Yboot,
                              g_poolboot=g_poolboot, pi_poolboot=pi_poolboot)
                       })

  GBoot <- t(data.frame(sapply(bootResults, function(x) x[1])))
  FBoot <- t(data.frame(sapply(bootResults, function(x) x[2])))
  GpBoot <- t(data.frame(sapply(bootResults, function(x) x[3])))
  g_b <- t(data.frame(sapply(bootResults, function(x) x[4])))
  pi_Xb <- t(data.frame(sapply(bootResults, function(x) x[5])))
  f_b <- t(data.frame(sapply(bootResults, function(x) x[6])))
  pi_Yb <- t(data.frame(sapply(bootResults, function(x) x[7])))
  g_poolb <- t(data.frame(sapply(bootResults, function(x) x[8])))
  pi_poolb <- t(data.frame(sapply(bootResults, function(x) x[9])))

  Difference_cdf <- abs(GBoot-FBoot)
  Difference_max <-apply(Difference_cdf, 1, max)
  Difference_pi <- data.frame(abs(pi_Xb-pi_Yb))
  T1_bsp = Difference_max
  T2_bsp = Difference_pi
  #Difference <- data.frame(cbind(Difference_max, Difference_pi))
  #T_bsp=apply(Difference,1,max)
  P_value1[i]<-(length(which(T1_bsp>=T1_hat))+1)/(B+1)
  P_value2[i]<-(length(which(T2_bsp>=T2_hat))+1)/(B+1)
  P_value_KS[i]<-min(P_value1[i],P_value2[i])

}

