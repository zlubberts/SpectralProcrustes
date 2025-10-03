#Experiment 3: Rank-1 perturbation, theta=pi/4 (perturbation partly in subspace of X), 
#             measure of strength of HA is given by sigma = sv of rank 1

verybeginning <- Sys.time()
library("JuliaCall")
#julia_install_package_if_needed("Manopt")
#julia_install_package_if_needed("Manifolds")
#julia_install_package_if_needed("LinearAlgebra")
#julia_install_package_if_needed("RecursiveArrayTools")
path_to_julia = "Add path to Julia here"

library(irlba)
require(doParallel)

procrustesFrob <- function(A,B) {
  #A, B are d x n
  #return W = argmin ||WB-A||_F
  tmp <- svd(A %*% t(B))
  return(tmp$u %*% t(tmp$v))
}

ase <- function(A,r){
  tmp <- irlba(A,r)
  return(tmp$u %*% diag(sqrt(tmp$d)))
}

getX <- function(n,r){
  blocksizes <- rmultinom(1,n-r,rep(1/r,r))+rep(1,r)
  X <- matrix(0,nrow=n,ncol=r)
  if(r>2){
    X[,1] <- c(rep(1,blocksizes[1]),rep(0,n-blocksizes[1]))
    for(i in 2:(r-1)){
      X[,i] <- c(rep(0,sum(blocksizes[1:(i-1)])),rep(1,blocksizes[i]),rep(0,n-sum(blocksizes[1:i])))
    }
    X[,r] <- c(rep(0,n-blocksizes[r]),rep(1,blocksizes[r]))
  }else{
    X <- cbind(c(rep(1,blocksizes[1]),rep(0,blocksizes[2])),c(rep(0,blocksizes[1]),rep(1,blocksizes[2])))
  }
  X <- X %*% diag(as.numeric(blocksizes^(-1/2)))
  B <- matrix(rep(0,r^2),nrow=r)
  B[upper.tri(B)] <- runif((r-1)*r/2)/r
  B <- B+t(B)
  diag(B) <- runif(r)*(r-1)/r
  B <- diag(as.numeric(sqrt(blocksizes))) %*% B %*% diag(as.numeric(sqrt(blocksizes)))
  tmp <- svd(B/sqrt(n))
  return(X %*% tmp$u %*% diag(as.numeric(sqrt(tmp$d))))
}

perturb <- function(X,n,r){
  qrX <- qr(X)
  Q <- qr.Q(qrX)
  a <- rnorm(r)
  a <- a/sqrt(sum(a^2))
  a <- Q%*%a
  b <- rnorm(n)
  y <- t(Q)%*%b
  b <- b-Q%*%y
  b <- b/sqrt(sum(b^2))
  u <- a+b
  #tested, this actually does have length sqrt(2) (to mach prec)
  u <- u/sqrt(2)
  v <- rnorm(r)
  v <- v/sqrt(sum(v^2))
  return(X+(u %*% t(v))/2)
}

drawAdjacency <- function(X,n){
  A <- matrix(0,nrow=n,ncol=n)
  P <- X %*% t(X)
  P[P<0] <- 0
  P[P>1] <- 1
  A[upper.tri(A)]<-(runif(n*(n-1)/2)<=P[upper.tri(P)])*1
  A <- A+t(A)
  return(A)
}

n <- 1000
r <- 3

nmc <- 100
bs <- 1000

mult <- 1.96/sqrt(nmc)

tvec <- (0:20)/20
num_t <- length(tvec)

allpowercurves <- array(0,dim=c(num_t,5,nmc),dimnames=list(t=tvec,statistic=1:5,trials=1:nmc))



cluster <- makeCluster(4)
registerDoParallel(cluster)
clusterExport(cluster,"path_to_julia")
clusterEvalQ(cluster, {
  library("JuliaCall")
  library("irlba")
  julia_setup(JULIA_HOME=path_to_julia)
  julia_source("procrustes_functions.jl")
})
cat("Setup time: ",difftime(Sys.time(),verybeginning,units="secs")," seconds.\n")
loopstart <- Sys.time()
for(i in 1:nmc){
  #generate problem instance
  X <- getX(n,r)
  Y <- perturb(X,n,r)
  #get critical value
  bootresults <- matrix(0,nrow=bs,ncol=5)
  bootresults <- foreach (j = 1:bs,.combine=rbind) %dopar% {
    A <- drawAdjacency(X,n)
    B <- drawAdjacency(X,n) #this is under the null!
    Xhat <- ase(A,r)
    Yhat <- ase(B,r)
    p0 <- procrustesFrob(t(Xhat),t(Yhat))
    T1 <- norm(p0%*%t(Yhat)-t(Xhat),type="F")
    T2 <- norm(p0%*%t(Yhat)-t(Xhat),type="2")
    T4 <- sum(apply(p0%*%t(Yhat)-t(Xhat),2,function(x) sqrt(sum(x^2))))
    pS <- julia_call("spectral_procrustes", t(Xhat), t(Yhat), p0)
    T3 <- norm(pS%*%t(Yhat)-t(Xhat),type="2")
    pR <- julia_call("robust_procrustes",t(Xhat),t(Yhat),p0)
    T5 <- sum(apply(pR%*% t(Yhat)-t(Xhat),2,function(x) sqrt(sum(x^2))))
    c(T1,T2,T3,T4,T5)
  }
  cat("i:",i,"cv")
  critvals <- apply(bootresults,2,quantile,probs=0.95)
  #generate power curve
  allpowercurves[1,,i] <- rep(.05,5)
  for(tidx in 2:num_t){
    t <- tvec[tidx]
    Yt <- (1-t)*X+t*Y
    bootresults <- matrix(0,nrow=bs,ncol=5)
    bootresults <- foreach (j = 1:bs,.combine=rbind) %dopar% {
      A <- drawAdjacency(X,n)
      B <- drawAdjacency(Yt,n)
      Xhat <- ase(A,r)
      Yhat <- ase(B,r)
      p0 <- procrustesFrob(t(Xhat),t(Yhat))
      T1 <- norm(p0%*%t(Yhat)-t(Xhat),type="F")
      T2 <- norm(p0%*%t(Yhat)-t(Xhat),type="2")
      T4 <- sum(apply(p0%*%t(Yhat)-t(Xhat),2,function(x) sqrt(sum(x^2))))
      pS <- julia_call("spectral_procrustes", t(Xhat), t(Yhat), p0)
      T3 <- norm(pS%*%t(Yhat)-t(Xhat),type="2")
      pR <- julia_call("robust_procrustes",t(Xhat),t(Yhat),p0)
      T5 <- sum(apply(pR%*% t(Yhat)-t(Xhat),2,function(x) sqrt(sum(x^2))))
      ((c(T1,T2,T3,T4,T5)-critvals)>0)*1
    }
    allpowercurves[tidx,,i] <- apply(bootresults,2,mean)
    cat(tidx)
  }
  cat("done.\n")
  cat("Estimated time per curve: ",difftime(Sys.time(),loopstart,units="secs")[[1]]/i," seconds.\n")
}
stopCluster(cl=cluster)
powercurves <- array(0,dim=c(num_t,3,5),dimnames=list(t=tvec,value=c("lower","mean","upper"),statistic=1:5))
powercurves[,2,] <- apply(allpowercurves,c(1,2),mean)
errorbars <- apply(allpowercurves,c(1,2),sd)*mult
powercurves[,1,] <- powercurves[,2,]-errorbars
powercurves[,3,] <- powercurves[,2,]+errorbars
save(powercurves,file="powercurves3.RData")

