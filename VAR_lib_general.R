library(flare)
library(igraph)
library(lpSolve)
library(mnormt)
library(matrixcalc)
library(pracma)
library(MASS)
library(doParallel)
library(foreach)
library(parallel)
library(abind)
library(compiler)
enableJIT(3)

# Summary of the library
# This library includes several functions: 
# 1.weight_epanechnikovker(t, m, n, b)
# 2.tvvar_fit(it,X,tau,b)
# 3.stvar_fit(X,tau)
# 4.lsvar_fit(it,X,lam,b)
# 5.mlevar_fit(it,X,b)
# 6.tvvar_gen(T,d,nSim,tol,graph)#data generation
# 7.tvvar_gen2()#data generation and cal threshold

weight_epanechnikovker <- function(t, m, n, b)
{
  v = (t-m/n)/b
  if(abs(v) > 1)
  {
    return(0)
  }
  else
  {
    numerator = 0.75 * (1 - v^2)
    vi=sapply(1:n,function(x){(t-x/n)/b})
    indx=abs(vi)<=1
    denumerator=sum(0.75*(1-vi[indx]^2))
    return(numerator/denumerator)
  }
}
vc_weight_epanechnikovker<-Vectorize(weight_epanechnikovker,vectorize.args='m')#vectorization of a scalar function
# weight_epanechnikovker=cmpfun(weight_epanechnikovker)
vc_weight_epanechnikovker=cmpfun(vc_weight_epanechnikovker)


# -----------------------------------------------------------------------------
# fit a TV_VAR model
# tvvar_fit is a function who fit a TVVAR model.
# there are four input: 
# it: estimate the transition matrix A_t at time t (X_t=X_{t-1}A_t), t can be larger than T
# X is a T*d array which stores the data of the time series. T is the length and d is the dimension.
# tau is the tuning parameters
# b is the bandwidth

tvvar_fit<-function(it,X,tau,b){
  T=dim(X)[1]
  d=dim(X)[2]
  #compute the kernel weights for t_i
  ker_weight=vc_weight_epanechnikovker(it/T,1:T,T,b)
  #compute the kernel weights for t_{i-1}
  ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
  #sigma_{i-1,0}
  s10=crossprod(X*drop(ker_weight_1),X)
  #sigma_{i,0}
  s00=crossprod(X*drop(ker_weight),X)
  #sigma_{i-1,1}
  s11=crossprod(X[1:(T-1),]*drop(ker_weight_1[1:(T-1)]),X[2:T,])
  #sigma_{i,1}
  s01=crossprod(X[1:(T-1),]*drop(ker_weight[1:(T-1)]),X[2:T,])
  #sigma_{i-1,-1}
  s1_1=crossprod(X[2:T,]*drop(ker_weight_1[2:T]),X[1:(T-1),])
  #sigma_{i,-1}
  s0_1=crossprod(X[2:T,]*drop(ker_weight[2:T]),X[1:(T-1),])
  #solve linear programming
  f.obj=rep(1,2*d)
  # f.con1=cbind(-s10,s10)
  # f.con2=cbind(s00,-s10)
  # f.con=rbind(f.con1,f.con2)
  f.con=matrix(rep(0,4*d^2),nrow=2*d)
  f.con[1:d,1:d]=-s10
  f.con[1:d,(d+1):(2*d)]=s10
  f.con[(d+1):(2*d),1:d]=s00
  f.con[(d+1):(2*d),(d+1):(2*d)]=-s10
  f.dir=rep("<=",2*d)
  f.rhs=matrix(0,2*d,d)
  A_hat_tmp=sapply(1:d,FUN=function(x){
    # f.rhs[,x]=matrix(c(tau-pmax(s11[,x],t(s1_1)[,x]),
    #                    tau+pmin(s01[,x],t(s0_1)[,x])),2*d,1)
    f.rhs[,x]=matrix(c(tau-pmax(s11[,x],s1_1[x,]),
                       tau+pmin(s01[,x],s0_1[x,])),2*d,1)
    lp.out=lp('min',f.obj,f.con,f.dir,f.rhs[,x])$solution
    return(matrix(lp.out[1:d]-lp.out[(1+d):(d+d)],d,1))
  })
  A_hat=t(A_hat_tmp)
  return(A_hat)
}


# -----------------------------------------------------------------------------
# fit a STAT_VAR model based on Fang Han's paper
# stvar_fit is a function who fit a STAT_VAR model.
# it has two input: 
# X is a T*d array which stores the data of the time series. T is the length and d is the dimension.
# tau is the tuning parameters

stvar_fit<-function(X,tau){
  T=dim(X)[1]
  d=dim(X)[2]
  #marginal covariance and lag 1 covariance
  S0=cov(X)
  S1=scale(crossprod(X[1:(T-1),],X[2:T,]),center=FALSE,scale=rep((T-1),d))
  # S1=crossprod(X[1:(T-1),],X[2:T,])/(T-1)
  # f.con.1=cbind(-S0,S0)
  # f.con.2=cbind(S0,-S0)
  # f.con=rbind(f.con.1,f.con.2)
  f.con=matrix(rep(0,4*d^2),nrow=2*d)
  f.con[1:d,1:d]=-S0
  f.con[1:d,(d+1):(2*d)]=S0
  f.con[(d+1):(2*d),1:d]=S0
  f.con[(d+1):(2*d),(d+1):(2*d)]=-S0
  f.rhs=matrix(0,2*d,d)
  f.obj=rep(1,2*d)
  f.dir=rep("<=",2*d)
  A_hat_stat_tmp=sapply(1:d,FUN=function(x){
    f.rhs[,x]=matrix(c(tau-S1[,x],tau+S1[,x]),2*d,1)
    lp.out=lp('min',f.obj,f.con,f.dir,f.rhs[,x])$solution
    # lp.out=solveLP( f.obj, f.rhs[,x], f.con, maximum = FALSE,
    #                 const.dir = rep( "<=", length( f.rhs[,x] ) ),
    #                 maxiter = 1000, zero = 1e-9, tol = 1e-6, dualtol = tol,
    #                 lpSolve = TRUE, solve.dual = FALSE, verbose = 0 )$solution
    out=matrix(lp.out[1:d]-lp.out[(1+d):(d+d)],d,1)
    return(out)
  })
  A_hat_stat=t(A_hat_stat_tmp)
  return(A_hat_stat)
}

# -----------------------------------------------------------------------------
# fit a Ridge_VAR model based on LS with l2 penalty
# lsvar_fit is a function who fit a Ridge_VAR model.
# there are four input: 
# it: estimate the transition matrix A_t at time t (X_t=X_{t-1}A_t), t can be larger than T
# X is a T*d array which stores the data of the time series. T is the length and d is the dimension.
# lam is the shrinkage parameter
# b is the bandwidth

lsvar_fit<-function(it,X,lam,b){
  T=dim(X)[1]
  d=dim(X)[2]
  #compute the kernel weights for t_i
  ker_weight=vc_weight_epanechnikovker(it/T,1:T,T,b)
  #compute the kernel weights for t_{i-1}
  ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
  W1=crossprod(X[2:T,]*drop(ker_weight[2:T]),X[1:(T-1),])
  W2=crossprod(X[1:T,]*drop(ker_weight_1[1:T]),X[1:T,])
  # A_hat_LS=solve(W2+lam*diag(d))%*%W1
  # A_hat_LS=crossprod(chol2inv(chol(W2+lam*diag(d))),W1)#W2 is symmetric
  A_hat_LS=tcrossprod(W1,chol2inv(chol(W2+lam*diag(d))))#W2 is symmetric
  return(A_hat_LS)
}


# -----------------------------------------------------------------------------
# fit a MLE_VAR model based on MLE
# mlevar_fit is a function who fit a MLE_VAR model.
# there are three input: 
# it: estimate the transition matrix A_t at time t (X_t=X_{t-1}A_t), t can be larger than T
# X is a T*d array which stores the data of the time series. T is the length and d is the dimension.
# b is the bandwidth
mlevar_fit<-function(it,X,b){
  T=dim(X)[1]
  d=dim(X)[2]
  #compute the kernel weights for t_i
  ker_weight=vc_weight_epanechnikovker(it/T,1:T,T,b)
  #compute the kernel weights for t_{i-1}
  ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
  W1=crossprod(X[2:T,]*drop(ker_weight[2:T]),X[1:(T-1),])
  W2=crossprod(X[1:T,]*drop(ker_weight_1[1:T]),X[1:T,])
  A_hat_MLE=tcrossprod(W1,chol2inv(chol(W2)))
  # A_hat_MLE=tcrossprod(W1,solve(W2))
  # A_hat_MLE=crossprod(chol2inv(chol(W2)),W1)
  # A_hat_MLE=W1%*%chol2inv(chol(W2))
  return(A_hat_MLE)
}

# -----------------------------------------------------------------------------
# fit a Lasso_VAR model based on LS with l1 penalty
# lavar_fit is a function who fit a Lasso_VAR model.
# there are four input: 
# it: estimate the transition matrix A_t at time t (X_t=X_{t-1}A_t), t can be larger than T
# X is a T*d array which stores the data of the time series. T is the length and d is the dimension.
# lam is the shrinkage parameter
# b is the bandwidth
# alpha is the stepsize
# A_old is the initial matrix of A_k and A_oldold is A_{k-1}

abso<-function(x){
  x[x<0]=0
  return(x)
}
lavar_fit<-function(it,X,lam,b,alpha,A_old,A_oldold){
  # alpha=1/max(svd(crossprod(X,X))$d)
  T=dim(X)[1]
  d=dim(X)[2]
  #compute the kernel weights for t_i
  ker_weight=vc_weight_epanechnikovker(it/T,1:T,T,b)
  #compute the kernel weights for t_{i-1}
  ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
  W1=crossprod(X[2:T,]*drop(ker_weight[2:T]),X[1:(T-1),])
  # W2=crossprod(X[1:T,]*drop(ker_weight_1[1:T]),X[1:T,])
  W2=crossprod(X[1:T,]*drop(ker_weight_1[1:T]),X[1:T,])
  dis=1
  step=2
  alpha=1/max(svd(2*W2)$d)
  # alpha=1/max(svd(crossprod(X,X))$d)
  while (dis>10^(-5)){
    print(c(step,dis))
    y=A_old+(step-1)/(step+2)*(A_old-A_oldold)
    tmp=y+2*alpha*(W1-tcrossprod(y,W2))
    # tmp=y+2*alpha*(W1-crossprod(W2,y))
    A_new=abso(abs(tmp)-alpha*lam)*sign(tmp)
    # dis=sum((A_new-A_old)^2)
    dis=norm(A_new-A_old,type="F")
    # dis=norm(A_new-y,type="F")
    # isallzero=(sum(A_new-A_old)==d*d)
    A_oldold=A_old
    A_old=A_new
    step=step+1
  }
  return(A_new)
  # return(A_old)
}











# -----------------------------------------------------------------------------
# generate time-varying vector AR(1) time series
# graph= "hub","cluster","random","band","scale-free"
tvvar_gen<-function(T,d,nSim,v=v,u=u,prob=prob,g=g_gen,graph1=graph1,graph2=graph2
                    ,rand_seed1,rand_seed2,rho_A0=rho_A0,rho_A1=rho_A1,t1,t2){
  #store time series
  TV_X=array(0,c(T,d,nSim))
  TV_A=array(0,c(d,d,T,nSim))
  # Generate TV-VAR(k) processes
  # Generate the initial transition matrix A:
  # this serves as the baseline of the tv coefficent matrices
  set.seed(rand_seed1)
  A0 = sugm.generator(d = d, graph = graph1,v=v,u=u,prob=prob, g=g, vis = FALSE, seed=rand_seed1)$sigma
  A1 = sugm.generator(d = d, graph = graph2,v=v,u=u,prob=prob, g=g, vis = FALSE, seed=rand_seed1+1)$sigma
  # Normalize A's spectral norm
  A0_2norm = max(svd(A0)$d)
  A0 = rho_A0 * A0 / A0_2norm
  A1_2norm = max(svd(A1)$d)
  A1 = rho_A1 * A1 / A1_2norm
  # Generate the marginal covariance matrix
  Sigma = diag(d)
  # Solve for the innovation covariance matrix
  Psi = Sigma - A0 %*% Sigma %*% t(A0)
  
  for(nn in 1:nSim)
  {
    set.seed(rand_seed2[nn])
    # Generate the iid random innovations: mean zero and normalized covariance matrix
    e = mvrnorm(n = T, mu = rep(0,d), Sigma = Psi)
    # e = rmt(n = T, mean = rep(0, d), sqrt(3/5) * Psi, df = 5)
    # e = rmt(n = T, mean = rep(0, d), sqrt(1/3) * Psi, df = 3)
    A = array(0, c(d, d, T))
    # Generate the TV-VAR(1) time series X, with sparse transition matrix A_i
    X = matrix(0,T,d)
    X[1,] = mvrnorm(n=1,mu=rep(0,d),Sigma = Psi)
    # X[1,] =rmt(n = 1, mean = rep(0, d), sqrt(3/5) * Psi, df = 5)
    # X[1,] =rmt(n = 1, mean = rep(0, d), sqrt(1/3) * Psi, df = 3)
    for (i in 1:T){
      A_tmp = matrix(0, d, d)
      A_tmp=(i/T)^t1*A1 +(1-(i/T))^t2* A0
      # A_tmp=A0+0.5* sin(2*pi*i/T)
      A[,,i] = A_tmp
      if (i>=2){
        X[i,]=A_tmp%*%X[i-1,]+matrix(e[i,], d, 1)
      }
    }
    TV_X[,,nn]=X
    TV_A[,,,nn]=A
  }
  return(list(TV_X=TV_X,TV_A=TV_A))
}

tvvar_gen2<-function(T,d,nSim,v=v,u=u,prob=prob,g=g_gen,graph1=graph1,graph2=graph2
                     ,rand_seed1,rand_seed2,rho_A0=rho_A0,rho_A1=rho_A1,t1,t2){
  #store time series
  TV_X=array(0,c(T,d,nSim))
  TV_A=array(0,c(d,d,T,nSim))
  # Generate TV-VAR(k) processes
  # Generate the initial transition matrix A:
  set.seed(rand_seed1)
  A0 = sugm.generator(d = d, graph = graph1,v=v,u=u,prob=prob, g=g, vis = FALSE, seed=rand_seed1)$sigma
  A1 = sugm.generator(d = d, graph = graph2,v=v,u=u,prob=prob, g=g, vis = FALSE, seed=rand_seed1+1)$sigma
  # Normalize A's spectral norm
  A0_2norm = max(svd(A0)$d)
  A0 = rho_A0 * A0 / A0_2norm
  A1_2norm = max(svd(A1)$d)
  A1 = rho_A1 * A1 / A1_2norm
  Sigma = diag(d)
  # Solve for the innovation covariance matrix
  Psi = Sigma - A0 %*% Sigma %*% t(A0)
  l1_thr=array(0,c(T-1,nSim))#store the l1 norm to calculate thresholds
  for(nn in 1:nSim)
  {
    set.seed(rand_seed2[nn])
    # Generate the iid random innovations: mean zero and normalized covariance matrix
    e = mvrnorm(n = T, mu = rep(0,d), Sigma = Psi)
    # tdf=7
    # e = rmt(n = T, mean = rep(0, d), sqrt((tdf-2)/tdf) * Psi, df = tdf)
    A = array(0, c(d, d, T))
    # Generate the TV-VAR(1) time series X, with sparse transition matrix A_i
    X = matrix(0,T,d)
    X[1,] = mvrnorm(n=1,mu=rep(0,d),Sigma = Psi)
    # X[1,] =rmt(n = 1, mean = rep(0, d), sqrt((tdf-2)/tdf) * Psi, df = tdf)
    for (i in 1:T){
      A_tmp = matrix(0, d, d)
      A_tmp=(i/T)^t1*A1 +(1-(i/T))^t2* A0
      A[,,i] = A_tmp
      if (i>=2){
        X[i,]=A_tmp%*%X[i-1,]+matrix(e[i,], d, 1)
      }
    }
    TV_X[,,nn]=X
    TV_A[,,,nn]=A
    
    #calculate the true threshold for pattern recovery
    l1_thr[1,nn]=norm(solve(Psi),type="1")
    l1_thr[2:(T-1),nn]=sapply(3:T,FUN=function(i){
      #A_{i-1}*A_{i-2}*...*A_3*A_2;
      #A_{i-1}*A_{i-2}*...*A_3;
      #...;
      #A_{i-1}*A_{i-2}
      #A_{i-1}
      acc_At=Reduce("%*%",lapply((i-1):2,FUN=function(x){return(A[,,x])}),accumulate=TRUE)
      tmp_mat=lapply(1:(i-2),FUN=function(x){return(tcrossprod((acc_At[[x]]%*%Psi),acc_At[[x]]))})
      Sigma_all=Reduce("+",tmp_mat,accumulate=FALSE)+Psi
      return(norm(solve(Sigma_all),type="1"))
    })
  }
  return(list(TV_X=TV_X,TV_A=TV_A,l1_thr=l1_thr))
}



