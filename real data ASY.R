## load whole data
genedata <- read.table("counts.txt")
E.day = colnames(genedata)
E.day = strsplit(E.day, split='.', fixed=TRUE)
E.day = do.call(rbind, E.day)[,1]
index_E3 = which(E.day == "E3")
index_E4 = which(E.day == "E4")

## E3 and E4 data seperate
E3.data = genedata[,index_E3]
E4.data = genedata[,index_E4]
data = cbind(E3.data, E4.data)
## calculate d
d.E3<-colSums(E3.data)
d.E4<-colSums(E4.data)
d<-colSums(data)

## select 1000 genes with highest mean(sd) (mean for gene across cells)
genemean<-as.data.frame(rowMeans(data))
library(matrixStats)
genedata1<-as.matrix(data)
genemean$sd<-rowSds(genedata1)
colnames(genemean)<-c("gemean","genesd")
genemean$geneid<-seq(1,length(data[,1]),1)
genemean<-as.data.frame(genemean)
genemean2<-genemean[order(-genemean$gemean, -genemean$genesd),]
finalgene<-genemean2[1:2000,]


select.E3 = E3.data[finalgene$geneid,]
select.E4 = E4.data[finalgene$geneid,]
select_gene = cbind(select.E3, select.E4)

prop<-matrix(0,2000,length(d))
for (j in 1:length(d)){
  for (i in 1:2000){
    prop[i,j]<-(select_gene[i,j]+1)/d[j]
  }
}
max<-apply(prop,1,max)
min<-apply(prop,1,min)
range<-cbind(min,max)

## KS two-sample test
P_value_asy<-rep(0,2000)
T_hat<-rep(0,2000)
T_index<-rep(0,2000)
m<-100
g<-matrix(0,nrow=m,ncol=2000)
f<-matrix(0,nrow=m,ncol=2000)
g_pool<-matrix(0,nrow=m,ncol=2000)
pi_X<-rep(0,2000)
pi_Y<-rep(0,2000)
pi_pool<-rep(0,2000)
N_asy <- 999
library(gplots)
library("VGAM")
library(MASS)
set.seed(238923)

for (i in 498:2000){
  tau <- exp(seq(log(min[i]), log(max[i]), length.out = 100))
  m <- length(tau)
  distance <- tau-c(0,tau[-100])
  NX <- length(d.E3)
  NY <- length(d.E4)
  dataX <- t(select.E3[i,])
  dX <- d.E3
  # remove outliers
  quant_X <- quantile(dataX/dX, probs=c(.25, .75), na.rm = FALSE)
  iqr_X <- IQR(dataX/dX)
  up_X <-  quant_X[2]+1.5*iqr_X # Upper Range
  low_X<- quant_X[1]-1.5*iqr_X # Lower Range
  X.index <- which(dataX/dX>= (quant_X[1] - 1.5*iqr_X) & dataX/dX <= (quant_X[2]+1.5*iqr_X))
  dataX <- dataX[X.index]
  dX <- dX[X.index]
  resultX<-deconv_new(tau=tau,X=dataX, d=dX, c0=1)
  pi_X[i]<-resultX$pi
  g[,i]<-resultX$stats[, "g"]  ## estimated density function
  G_est<-resultX$stats[, "G"] ## estimated cumulative density function
  X.bias = resultX$bias.all
  X.Cov= resultX$Cov.all
  #unstim
  dataY <- t(select.E4[i,])
  dY <- d.E4
  quant_Y <- quantile(dataY/dY, probs=c(.25, .75), na.rm = FALSE)
  iqr_Y <- IQR(dataY/dY)
  up_Y<-  quant_Y[2]+1.5*iqr_Y # Upper Range
  low_Y<- quant_Y[1]-1.5*iqr_Y # Lower Range
  Y.index <- which(dataY/dY>= (quant_Y[1] - 1.5*iqr_Y) & dataY/dY <= (quant_Y[2]+1.5*iqr_Y))
  dataY <- dataY[Y.index]
  dY <- dY[Y.index]
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

  while( det(Cov.dff) <= 0){
    diag(Cov.dff) = diag(Cov.dff)+1e-5
  }

  All_asy = abs(mvrnorm(N_asy, bias.dff, Cov.dff, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
  pi_diff_asy = All_asy[,1]
  T_pre = cbind(pi_diff_asy, apply(All_asy[,-1], 1, max))
  T_asy = apply(T_pre, 1, max)
  P_value_asy[i]<-(length(which(T_asy>=T_hat[i]))+1)/(N_asy+1)

}






P_value2<-P_value_asy
P_value_adj<-p.adjust(P_value2, method="fdr")
P_value3<-rep(0,2000)
P_value3<-P_value_adj
P_value3<-as.data.frame(cbind(P_value3,pi_X, pi_Y, T_index))
P_value3$geneid<-rownames(finalgene)

write.csv(P_value3, file = "G model 2000 genes results with pi asy p-value.csv")

