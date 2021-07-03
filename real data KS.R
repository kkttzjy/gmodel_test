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

## select 2000 genes with highest mean(sd) (mean for gene across cells)
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

## KS two-sample tes
P_value<-rep(0,2000)
T_hat<-rep(0,2000)
T_index<-rep(0,2000)
m<-100
g<-matrix(0,nrow=m,ncol=2000)
f<-matrix(0,nrow=m,ncol=2000)
g_pool<-matrix(0,nrow=m,ncol=2000)
pi_X<-rep(0,2000)
pi_Y<-rep(0,2000)
pi_pool<-rep(0,2000)
B <- 999
library(gplots)
library(sm)
library("VGAM")
set.seed(238923)

pdf("plot check 2000 test new KS test with early stopping rule.pdf")
for (i in 1:2000){
  tau <- exp(seq(log(min[i]), log(max[i]), length.out = 100))
  m <- length(tau)
  distance <- tau-c(0,tau[-100])
  NX <- length(d.E3)
  NY <- length(d.E4)
  dataX <- t(select.E3[i,])
  dX <- d.E3
  resultX<-deconv_new(tau=tau,X=dataX, d=dX, c0=1)
  pi_X[i]<-resultX$pi
  g[,i]<-resultX$stats[, "g"]  ## estimated density function
  G_est<-resultX$stats[, "G"] ## estimated cumulative density function
  #unstim
  dataY <- t(select.E4[i,])
  dY <- d.E4
  resultY<-deconv_new(tau=tau,X=dataY, d=dY, c0=1)
  pi_Y[i]<-resultY$pi
  f[,i]<-resultY$stats[, "g"] ## estimated density function
  F_est<-resultY$stats[, "G"] ## estimated cumulative density function
  ## Under null hypothesis: pull X and Y together to estimate g(theta)
  sample <- c(dataX, dataY)
  d<-c(dX, dY)
  results <- deconv_new(tau=tau,X=sample,d=d,c0=1)
  pi_pool[i]<-results$pi
  g_pool[,i]<- results$stats[, "g"]
  Difference_hat <- abs(G_est-F_est)
  Difference_pi <- abs(pi_X[i]-pi_Y[i])
  T_hat=max(max(Difference_hat),Difference_pi)
  T1_hat = max(Difference_hat)
  T2_hat = Difference_pi
  T_index[i] = ifelse(T1_hat>T2_hat,1,0)
  ### Plots
  E3_dataest<-as.data.frame(dataX/dX)
  names(E3_dataest)<-"theta"
  E4_dataest<-as.data.frame(dataY/dY)
  names(E4_dataest)<-"theta"
  E3_dataest$trt<-"E3"
  E4_dataest$trt<-"E4"
  data_est<-rbind(E3_dataest, E4_dataest)
  data_est$trt<-as.factor(data_est$trt)
  par(mfrow=c(2,2))
  plot(y=g[,i], x=tau, type="l",
       xlim=c(min[i], max[i]),
       xlab=expression(theta/tau),
       ylab="density",
       main="Density ")
  lines(y=f[,i], x=tau, col=2, lty=2, lwd=2)
  legend('topright',cex = 0.6, pt.cex(2), c("E3", "E4"), col=c(1,2), lty=c(1,2))

  plot(y=G_est, x=tau, type="l",
       xlim=c(min[i], max[i]),
       xlab=expression(theta/tau),
       ylab="CDF",
       main="CDFs ")
  lines(y=F_est, x=tau, col=2, lty=2, lwd=2)
  legend("bottomright", cex = 0.6, c("E3", "E4"), col=c(1,2), lty=c(1,2))

  trt.f <- factor(data_est$trt, levels= c(1,2),
                  labels = c("E3","E4"))

  sm.density.compare(data_est$theta, data_est$trt,xlab="theta/tau",col=c(1,2),xlim=c(min[i], max[i]))
  legend("topright", cex = 0.6, levels(data_est$trt), col =c(1,2), lty=c(1,2) )



  plot(x=dX,y=dataX)
  points(x=dY,y=dataY, col='red')
  legend('topright', cex = 0.6, c("E3", "E4"), col=c(1,2),pch=c(1,1))

  ### bootstrap

  genBootSampleX <- function() {
    thetaStar <- sample(tau, size = NX, replace = TRUE, prob = g_pool[,i])
    xStar <- sapply(seq_len(NX),
                    function(j){
                      dX_boot<-dX[j]
                      dataX_boot<-rzipois(n = 1, lambda=thetaStar[j]*dX_boot, pstr0 = pi_pool[i])
                    })
  }

  genBootSampleY <- function() {
    taoStar <- sample(tau, size = NY, replace = TRUE, prob = g_pool[,i])
    yStar <- sapply(seq_len(NY),
                    function(j){
                      dY_boot<-dY[j]
                      dataY_boot<-rzipois(n = 1, lambda=taoStar[j]*dY_boot, pstr0 = pi_pool[i])})
  }

  I<-0
  delta<-0.4
  p_0<-0.01
  c<-(1+delta)*p_0/(1-p_0)
  a<-5
  P_value_B<-0
  for (j in 1:B){
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
    Difference_hat <- abs(GBoot-FBoot)
    Difference_pi <- abs(pi_Xboot-pi_Yboot)
    T_bsp=max(max(Difference_hat),Difference_pi)
    T1_bsp = max(Difference_hat)
    T2_bsp = Difference_pi
    if ( T_bsp >= T_hat){
      I<-I+1
    }
    P_value_B<-I/j
    if (P_value_B > ((a/j+c)/(1+c)) ){
      P_value[i]<-P_value_B
      break
    }
    if (j==B){
      P_value[i]<-(I+1)/(B+1)
    }
  }
  mtext(text = paste(rownames(select_gene)[i],P_value[i], pi_X[i], pi_Y[i]), outer=TRUE,  cex=1, line=-1)
}
dev.off()






P_value2<-P_value
P_value_adj<-p.adjust(P_value2, method="fdr")

P_value3<-as.data.frame(cbind(P_value, P_value_adj, pi_X, pi_Y, T_index))
P_value3$geneid<-rownames(finalgene)

write.csv(P_value3, file = "G model 2000 genes results with pi new codes.csv")



