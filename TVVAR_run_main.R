rm(list=ls())
library(flare)
library(igraph)
library(lpSolve)
library(mnormt)
library(matrixcalc)
library(pracma)
library(MASS)
library(doParallel)
library(foreach)
# library(parallel)
library(abind)
library(compiler)
library(mvtnorm)



setwd()
source("VAR_lib_general.R")
# -------------------------------------------
# Necessary function
# -------------------------------------------
norm_vec <- function(x) sqrt(sum(x^2))/length(x)
acomb <- function(...) abind(..., along=3)#function for combination of array
acomb2 <- function(...) abind(..., along=4)#function for combination of array
enableJIT(3)
# -------------------------------------------
# Simulation setup
# -------------------------------------------
ncores=8
rand_seed2=1001:1100
# rand_seed2=1:10
nSim=length(rand_seed2)
d =20
T = 100
#bandwidth fitting
b_tvvar=0.8*T^(-1/5)
b_las=0.8*T^(-1/5)
b_Rid=0.8*T^(-1/5)
b_mle=0.8
#bandwidth training
b_train_tvvar=0.35
b_train_las=0.35
b_train_Rid=0.35
#####################
#data generation
graph1="random"
graph2="random"
rho_A0=0.2
rho_A1=1
t1=2
t2=4
g_gen=NULL
prob=0.01
v=0.01
u=10
filename=paste("Results_d",d,"_",graph1,"_g",g_gen,"_prob",prob,"_v",v,"_u",u,".RData",sep="")
####################
L_tau=30#number of tau
L_lam=30#number of lambda (shrinkage parameter for ridge method)
tau_grid=seq(0.001,1,length.out=L_tau)#all tau
lam_grid=seq(0.001,2,length.out=L_lam)#all lambda
# #Theoretically guarranteed region [b,1-b]
# t_lower=round(T*b_tvvar)+1
# t_upper=round(T*(1-b_tvvar))-1

#---------------------------------------------
# Data generation
#---------------------------------------------
gen=tvvar_gen(T,d,nSim,v=v,u=u,prob=prob,g=g_gen,graph1=graph1,graph2=graph2,
              rand_seed1=200,rand_seed2=rand_seed2,rho_A0=rho_A0,rho_A1=rho_A1,t1,t2)
TV_X=gen$TV_X
TV_A=gen$TV_A


cl <- makeCluster(ncores)
registerDoParallel(cl)
ptm <- proc.time()
#---------------------------------------------
# Cross-validation to select tau and lamdba
#---------------------------------------------
# seperate the data to 2 part evenly and do one-step-ahead predition and compare the errors
# initialization
X_train=matrix(rep(0,T*d/2),nrow=T/2)
A_hat_tvvar=A_hat_statvar=A_hat_Rid=A_hat_las=matrix(rep(0,d*d),nrow=d)
X_test=X_pred_tvvar=X_pred_statvar=X_pred_Rid=X_pred_las=matrix(rep(0,d),nrow=d)
# b_train=0.8*(T/2)^(-1/5)
# b_train=0.8
error_pred_tvvar_statvar=array(0,c(T/2-1,L_tau,2,nSim))
error_pred_Rid_las=array(0,c(T/2-1,L_lam,2,nSim))
selected_tau_tvvar=selected_tau_statvar=selected_lam_Rid=selected_lam_las=matrix(rep(0,nSim),nrow=nSim)

###########
#TV_VAR AND stat_VAR
###########
for (nn in 1:nSim){
  cat(paste("TV&STAT",nn))
  X=TV_X[,,nn]
  error_pred_tmp<-foreach(itau=1:L_tau, .combine='acomb',.multicombine=TRUE,.packages=c('lpSolve'))%:%
    foreach(it=1:(T/2-1), .combine='cbind',.packages=c('lpSolve'))%dopar%{
      tau=tau_grid[itau]
      X_train=X[it:(T/2+it-1),]
      X_test=X[T/2+it,]
      #TVVAR estimation
      A_hat_tvvar=tvvar_fit((T/2+1),X_train,tau,b_train_tvvar)
      X_pred_tvvar=A_hat_tvvar%*%X_train[T/2,]
      #statVAR estimation
      A_hat_statvar=stvar_fit(X_train,tau)
      X_pred_statvar=A_hat_statvar%*%X_train[T/2,]
      error_tmp=c(norm_vec(X_pred_tvvar-X_test),norm_vec(X_pred_statvar-X_test))
      # error_tmp=c(max(abs((X_pred_tvvar-X_test))),max(abs(X_pred_statvar-X_test)))
    }
  error_pred_tvvar_statvar[,,,nn]=aperm(error_pred_tmp, c(2,3,1))
  selected_tau_tvvar[nn]=tau_grid[which.min(colMeans(error_pred_tvvar_statvar[,,1,nn]))]
  selected_tau_statvar[nn]=tau_grid[which.min(colMeans(error_pred_tvvar_statvar[,,2,nn]))]
}

###########
#LS
###########
for (nn in 1:nSim){
  cat(paste("LS",nn))
  X=TV_X[,,nn]
  # LS
  error_pred_tmp<-foreach(ilam=1:L_lam, .combine='acomb',.multicombine=TRUE)%:%
    foreach(it=1:(T/2-1), .combine='cbind')%dopar%{
      lam=lam_grid[ilam]
      X_train=X[it:(T/2+it-1),]
      X_test=X[T/2+it,]
      #Ridge estimation
      A_hat_Rid=lsvar_fit((T/2+1),X_train,lam,b_train_las)
      X_pred_Rid=A_hat_Rid%*%X_train[T/2,]
      #Lasso estimation
      A_hat_las=lavar_fit((T/2+1),X_train,lam,b_train_Rid,alpha=0.0001,A_old=array(0.5,c(d,d)),A_oldold=array(0.5,c(d,d)))
      X_pred_las=A_hat_las%*%X_train[T/2,]
      # norm_vec(X_pred_LS-X_test)
      error_tmp=c(norm_vec(X_pred_Rid-X_test),norm_vec(X_pred_las-X_test))
    }
  # selected_lam_LS[nn]=lam_grid[which.min(colMeans(error_pred_tmp))]
  error_pred_Rid_las[,,,nn]=aperm(error_pred_tmp, c(2,3,1))
  selected_lam_Rid[nn]=lam_grid[which.min(colMeans(error_pred_Rid_las[,,1,nn]))]
  selected_lam_las[nn]=lam_grid[which.min(colMeans(error_pred_Rid_las[,,2,nn]))]
}

#---------------------------------------------
# Fit the model using selected parameters
#---------------------------------------------
# initialization
A_hat_tvvar=A_hat_Rid=A_hat_las=A_hat_MLE=array(0,c(d,d,T,nSim))
A_hat_statvar=array(0,c(d,d,nSim))
mean_err_tvvar_all=mean_err_stat_all=mean_err_Rid_all=mean_err_las_all=mean_err_MLE_all=array(0,c(4,nSim))#4 kinds of matrix norms
tmp=tmp5=tmp2=tmp3=tmp4=array(0,c(4,T,nSim))


###########
#TV_VAR
###########
#fit model
print("fit_TVVAR")
A_hat_tvvar<-foreach(nn=1:nSim, .combine='acomb2',.multicombine=TRUE,.packages=c('lpSolve'))%:%
  foreach(it=1:T, .combine='acomb',.packages=c('lpSolve'))%dopar%{
    tau1=selected_tau_tvvar[nn]
    output=tvvar_fit(it,TV_X[,,nn],tau1,b_tvvar)
  }

#calculate error
tmp<-foreach(nn=1:nSim, .combine='acomb',.multicombine=TRUE)%:%
  foreach(it=1:T, .combine='cbind')%dopar%{
    At=TV_A[,,it,nn]
    A_hat_tmp=A_hat_tvvar[,,it,nn]
    #     St=abs(At)>tol
    #     S_hat_tmp=abs(A_hat_tmp)>tol
    #Error
    error_Inf_tmp=Matrix::norm(A_hat_tmp-At,"I")
    error_L1_tmp=Matrix::norm(A_hat_tmp-At,"1")
    error_spec_tmp=max(svd(A_hat_tmp-At)$d)
    error_F_tmp=Matrix::norm(A_hat_tmp-At,"F")/sqrt(d)
    matrix(c(error_Inf_tmp,error_L1_tmp,error_spec_tmp,error_F_tmp),
           nrow=4,ncol=1)
  }
#Theoretically guarranteed region [b,1-b]
t_lower=round(T*b_tvvar)+1
t_upper=round(T*(1-b_tvvar))-1
mean_err_tvvar_all=sapply(1:nSim,FUN=function(nn){
  rowMeans(tmp[,t_lower:t_upper,nn])
})

###########
#stat_VAR
###########
#fit model
print("fit stat_VAR")
A_hat_statvar<-foreach(nn=1:nSim, .combine='acomb',.multicombine=TRUE,.packages=c('lpSolve'))%dopar%{
  tau2=selected_tau_statvar[nn]
  output=stvar_fit(TV_X[,,nn],tau2)
}
#calculate error
tmp2=array(0,c(4,T,nSim))
tmp2<-foreach(nn=1:nSim, .combine='acomb',.multicombine=TRUE)%:%
  foreach(it=1:T, .combine='cbind')%dopar%{
    At=TV_A[,,it,nn]
    A_hat_tmp=A_hat_statvar[,,nn]
    #     St=abs(At)>tol
    #     S_hat_tmp=abs(A_hat_tmp)>tol
    #Error
    error_Inf_tmp=Matrix::norm(A_hat_tmp-At,"I")
    error_L1_tmp=Matrix::norm(A_hat_tmp-At,"1")
    error_spec_tmp=max(svd(A_hat_tmp-At)$d)
    error_F_tmp=Matrix::norm(A_hat_tmp-At,"F")/sqrt(d)
    matrix(c(error_Inf_tmp,error_L1_tmp,error_spec_tmp,error_F_tmp),
           nrow=4,ncol=1)
  }
#Theoretically guarranteed region [b,1-b]
t_lower=round(T*b_tvvar)+1
t_upper=round(T*(1-b_tvvar))-1
mean_err_stat_all=sapply(1:nSim,FUN=function(nn){
  rowMeans(tmp2[,t_lower:t_upper,nn])
})

###########
#Lasso
##########
#fit model
print("fit Las")
A_hat_las<-foreach(nn=1:nSim, .combine='acomb2',.multicombine=TRUE)%:%
  foreach(it=1:T, .combine='acomb')%dopar%{
    lam=selected_lam_las[nn]
    output=lavar_fit(it,TV_X[,,nn],lam,b_las,alpha=0.0001,A_old=array(0.5,c(d,d)),A_oldold=array(0.5,c(d,d)))
  }
#calculate error
tmp3<-foreach(nn=1:nSim, .combine='acomb',.multicombine=TRUE,.packages=c('lpSolve'))%:%
  foreach(it=1:T, .combine='cbind',.packages=c('lpSolve'))%dopar%{
    At=TV_A[,,it,nn]
    A_hat_tmp=A_hat_las[,,it,nn]
    #Error
    error_Inf_tmp=Matrix::norm(A_hat_tmp-At,"I")
    error_L1_tmp=Matrix::norm(A_hat_tmp-At,"1")
    error_spec_tmp=max(svd(A_hat_tmp-At)$d)
    error_F_tmp=Matrix::norm(A_hat_tmp-At,"F")/sqrt(d)
    matrix(c(error_Inf_tmp,error_L1_tmp,error_spec_tmp,error_F_tmp),
           nrow=4,ncol=1)
  }
#Theoretically guarranteed region [b,1-b]
t_lower=round(T*b_las)+1
t_upper=round(T*(1-b_las))-1
mean_err_las_all=sapply(1:nSim,FUN=function(nn){
  rowMeans(tmp3[,t_lower:t_upper,nn])
})

###########
#Ridge
###########
#fit model
print("fit Rid")
A_hat_Rid<-foreach(nn=1:nSim, .combine='acomb2',.multicombine=TRUE)%:%
  foreach(it=1:T, .combine='acomb')%dopar%{
    lam=selected_lam_Rid[nn]
    output=lsvar_fit(it,TV_X[,,nn],lam,b_Rid)
  }
#calculate error
tmp4<-foreach(nn=1:nSim, .combine='acomb',.multicombine=TRUE)%:%
  foreach(it=1:T, .combine='cbind')%dopar%{
    At=TV_A[,,it,nn]
    A_hat_tmp=A_hat_Rid[,,it,nn]
    #Error
    error_Inf_tmp=Matrix::norm(A_hat_tmp-At,"I")
    error_L1_tmp=Matrix::norm(A_hat_tmp-At,"1")
    error_spec_tmp=max(svd(A_hat_tmp-At)$d)
    error_F_tmp=Matrix::norm(A_hat_tmp-At,"F")/sqrt(d)
    matrix(c(error_Inf_tmp,error_L1_tmp,error_spec_tmp,error_F_tmp),
           nrow=4,ncol=1)
  }
#Theoretically guarranteed region [b,1-b]
t_lower=round(T*b_Rid)+1
t_upper=round(T*(1-b_Rid))-1
mean_err_Rid_all=sapply(1:nSim,FUN=function(nn){
  rowMeans(tmp4[,t_lower:t_upper,nn])
})

###########
#MLE
###########
#fit model
print("fit MLE")
A_hat_MLE<-foreach(nn=1:nSim, .combine='acomb2',.multicombine=TRUE)%:%
  foreach(it=1:T, .combine='acomb')%dopar%{
    output=mlevar_fit(it,TV_X[,,nn],b_mle)
  }
#calculate error
tmp5<-foreach(nn=1:nSim, .combine='acomb',.multicombine=TRUE)%:%
  foreach(it=1:T, .combine='cbind')%dopar%{
    At=TV_A[,,it,nn]
    A_hat_tmp=A_hat_MLE[,,it,nn]
    #     St=abs(At)>tol
    #     S_hat_tmp=abs(A_hat_tmp)>tol
    #Error
    error_Inf_tmp=Matrix::norm(A_hat_tmp-At,"I")
    error_L1_tmp=Matrix::norm(A_hat_tmp-At,"1")
    error_spec_tmp=max(svd(A_hat_tmp-At)$d)
    error_F_tmp=Matrix::norm(A_hat_tmp-At,"F")/sqrt(d)
    matrix(c(error_Inf_tmp,error_L1_tmp,error_spec_tmp,error_F_tmp),
           nrow=4,ncol=1)
  }
#Theoretically guarranteed region [b,1-b]
t_lower=round(T*b_mle)+1
t_upper=round(T*(1-b_mle))-1
mean_err_MLE_all=sapply(1:nSim,FUN=function(nn){
  rowMeans(tmp5[,t_lower:t_upper,nn])
})

t=proc.time() - ptm
stopCluster(cl)

avg_err=cbind(rowMeans(mean_err_tvvar_all),rowMeans(mean_err_stat_all),rowMeans(mean_err_las_all),
              rowMeans(mean_err_Rid_all),rowMeans(mean_err_MLE_all))
avg_err
sd_err=cbind(apply(mean_err_tvvar_all,1,sd),apply(mean_err_stat_all,1,sd),apply(mean_err_las_all,1,sd),
             apply(mean_err_Rid_all,1,sd),apply(mean_err_MLE_all,1,sd))
save.image(filename)