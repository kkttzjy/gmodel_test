deconv_new <- function(tau, X, Q, P, d,   
                       c0 = 1,
                       scale = TRUE,
                       pDegree = 5,
                       aStart = 0.5, ...) {
  
  m <- length(tau)
  y=1
  theta_opt<-X/d
  Q <- cbind(sqrt((m - 1) / (m)), scale(splines::ns(tau, pDegree)))
  Q <- apply(Q, 2, function(w) w / sqrt(sum(w * w)))
  P_pre<-matrix(0,nrow=length(X),ncol=m)
  P<-matrix(0,nrow=length(X),ncol=m)
  p <- ncol(Q)+1
  pGiven <- length(aStart)
  if (pGiven == 1) {
    aStart <- rep(aStart, p)
  } 
  P_pre <-t(sapply(seq(1,length(X),1),function(i) dpois(X[i],lambda=tau*d[i])))
  
  ## precalculate part of P matrix
  for ( i in 1:length(X)){
    if (theta_opt[i] <= tau[1] ){
      P_pre[i,1] <- dpois(X[i],lambda=theta_opt[i]*d[i])
    }
    for ( j in 2:m){
      if ( tau[j-1] <= theta_opt[i] & tau[j] >= theta_opt[i] ){
        if (abs(theta_opt[i]-tau[j-1]) < abs(theta_opt[i]-tau[j])){
          P_pre[i,j-1] <- dpois(X[i],lambda=theta_opt[i]*d[i])
        } else{
          P_pre[i,j] <- dpois(X[i],lambda=theta_opt[i]*d[i])
        }
      }
    }
  }
  
  V<-matrix(0,nrow=length(X),ncol=m)
  V <-t(sapply(seq(1,length(X),1),function(i) 
  {as.numeric(I(X[i]==0))*(1-dpois(0,lambda=tau*d[i]))-(1-as.numeric(I(X[i]==0)))*dpois(X[i],lambda=tau*d[i])
  }))
  
  
  distance <- tau-c(0,tau[-m])
  statsFunction <- function(a) {
    # for ( i in 1:length(X)){
    #   for (j in 1:m){
    #     P[i,j]= ifelse (tau[j]==0,
    #                     P_pre[i,j],
    #                     ifelse(
    #                       X[i]==0,
    #                       a[1]+(1-a[1])*P_pre[i,j],
    #                       (1-a[1])*P_pre[i,j]
    #                     )
    #     )}}
    P = (1-a[1])*P_pre
    P[which(X==0),] = a[1] + P[which(X==0),] 
    pi <- a[1]
    g <- as.vector(exp(Q %*% a[-1]))
    g <- g / sum(g)
    G <- cumsum(g)
    f <- as.vector(P %*% g)
    l<--sum(y * log(f)) + c0 * sum(a[-1]^2)^.5
    yHat <- if (length(y) == 1 && y == 1) y else sum(y) * f
    Pt <- P / f
    W <- g * (t(Pt) - 1 )
    qw <- t(Q) %*% W
    ywq <- (yHat * t(W)) %*% Q
    I1 <- qw %*% ywq   ## Fisher Information Matrix I(\alpha)
 
    ### pi part
    dPpi <- matrix(0,nrow=length(X),ncol=m)
    dpg <-rep(0,length(X))
    
    # for ( i in 1:length(X)){
    #   dPpi[i,] <- ifelse(X[i]==0, 1, 0)*(1-exp(-tau*d[i])) - ifelse(X[i]==0, 0, 1)*dpois(X[i],lambda=tau*d[i])
    #   dpg[i] <- t(dPpi[i,]) %*% g
    # }
    # 
    for ( i in 1:length(X)){
      dPpi[i,] <- ifelse(X[i]==0, 1, 0)*(1-P_pre[i,]) - ifelse(X[i]==0, 0, 1)*P_pre[i,]
      dpg[i] <- t(dPpi[i,]) %*% g
    }
    Ipi <- (1/(f^2)) %*% (dpg^2)
    
    Ipia <- rep(0, p-1)
    Ipia_pre1 <- matrix(0, nrow = m, ncol = p-1)
    Ipia_pre2 <- matrix(0, nrow = length(X), ncol = p-1)
    
    for ( i in 1:length(X)){
      for ( j in 1:m){
        Ipia_pre1[j,] <- t(Q[j,]) * (g[j]/(f[i]^2)*(P[i,j]*dpg[i]-dPpi[i,j]*f[i])) 
      }
      Ipia_pre2[i,] <- colSums(Ipia_pre1)
    }
    Ipia <- colSums(Ipia_pre2)
    I1 <- rbind(t(c(Ipi, Ipia)), cbind(Ipia, I1))
    
    aa <- sqrt(sum(a^2))
    sDot <- c0 * a / aa
    sDotDot <- (c0 / aa) * ( diag(length(a)) - outer(a, a) /aa^2 )
 
    ## The R value
    R <- sum(diag(sDotDot)) / sum(diag(I1))
    
    I2 <- solve(I1 + sDotDot)
    bias <- as.vector(-I2 %*% sDot)
    Cov <- I2 %*% (I1 %*% t(I2))
    ## Se <- diag(Cov)^.5
    Dq <- (diag(g) - outer(g, g)) %*% Q
   
    
    Tran_der <- rbind(t(c(1, rep(0, p-1))), cbind(rep(0, m), Dq))
    
    bias.T <- Tran_der %*% bias
    Cov.T <- Tran_der %*% Cov %*% t(Tran_der)
    
    
    
    Tran_G <- diag(length(tau))
    Tran_G[lower.tri(Tran_G)] <- 1
    Tran2 <- rbind(t(c(1, rep(0, length(tau)))), cbind(rep(0, length(tau)), Tran_G))
    bias.T2 <- Tran2 %*% bias.T
    Cov.T2 <- Tran2 %*% Cov.T %*% t(Tran2)
    

    
    mat <- cbind(tau, g, G)
    colnames(mat) = c("theta", "g", "G")
    list(S = R, bias = bias, cov = Cov, 
         mat = mat, 
         pi = pi,
         bias.all = bias.T2,
         Cov.all = Cov.T2
        )
  }
  
  loglik <- function(a) {
    P = (1-a[1])*P_pre
    P[which(X==0),] = a[1] + P[which(X==0),] 
    
    g <- exp(Q %*% a[-1])
    g <- as.vector(g / sum(g))
    f <- as.vector(P %*% g)
    value <- -sum(y * log(f)) + c0 * sum(a^2)^.5
    # Pt <- P / f
    # W <- g * (t(Pt) - 1 )
    # qw <- t(Q) %*% W
    # aa <- sqrt(sum(a^2))
    # sDot <- c0 * a / aa
    # ldot.pi <- (V/f) %*% g
    # ldot <- rbind(t(ldot.pi), qw)
    # attr(value, "gradient") <- -rowSums(ldot) + sDot
    value
  }
  grr <- function(a) {
    P = (1-a[1])*P_pre
    P[which(X==0),] = a[1] + P[which(X==0),] 
    
    g <- exp(Q %*% a[-1])
    g <- as.vector(g / sum(g))
    f <- as.vector(P %*% g)
    value <- -sum(y * log(f)) + c0 * sum(a^2)^.5
    Pt <- P / f
    W <- g * (t(Pt) - 1 )
    qw <- t(Q) %*% W
    aa <- sqrt(sum(a^2))
    sDot <- c0 * a / aa
    ldot.pi <- (V/f) %*% g
    ldot <- rbind(t(ldot.pi), qw)
    -rowSums(ldot) + sDot
  }
  
  result <- stats::nlminb(aStart, loglik, gradient = grr, 
                          lower=c(0, rep(-Inf,6)), upper=c(0.9999, rep(Inf,6)))
  
  mle <- result$par
  stats <- statsFunction(mle)
  list(mle = mle,
       Q = Q,
       S = stats$S,
       cov = stats$cov,
       stats = stats$mat,
       loglik = loglik,
       pi=stats$pi,
       bias.all = stats$bias.all,
       Cov.all = stats$Cov.all,
       l=stats$l)
}

