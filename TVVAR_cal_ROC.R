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
library(parallel)
library(abind)
library(compiler)
library(mvtnorm)


setwd("D:\\OneDrive\\Time-varying VAR\\code\\simulation\\R code")
# setwd("C:\\Users\\DX\\OneDrive\\Time-varying VAR\\code\\simulation\\R code")
# setwd("C:\\Users\\Xin\\OneDrive\\Time-varying VAR\\code\\simulation\\R code")
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
rand_seed2=c(1:9,11:26,101:125)
# rand_seed2=21:34
nSim=length(rand_seed2)
d = 20
T = 100
b=0.8*T^(-1/5)
#####################
#data generation
graph1="hub"
graph2="hub"
rho_A0=0.2
rho_A1=1
t1=2
t2=4
g_gen=8
prob=NULL
v=0.001
u=10
if (graph1=="band"){
  tol_true=10^-3
}else{
  tol_true=0
}
filename=paste("ROC_d",d,"_",graph1,"_g",g_gen,"_prob",prob,"_v",v,"_u",u,".Rdata",sep="")
####################
L_tau=30#number of tau
# L_lam=30#number of lambda (shrinkage parameter for ridge method)
tau_grid=seq(0.001,0.4,length.out=L_tau)#all tau
# lam_grid=seq(0.001,2,length.out=L_lam)#all lambda
#Theoretically guarranteed region [b,1-b]
t_lower=round(T*b)+1
t_upper=round(T*(1-b))-1

#---------------------------------------------
# Data generation
#---------------------------------------------
gen=tvvar_gen2(T,d,nSim,v=v,u=u,prob=prob,g=g_gen,graph1=graph1,graph2=graph2,
              rand_seed1=50,rand_seed2=rand_seed2,rho_A0=rho_A0,rho_A1=rho_A1,t1,t2)
TV_X=gen$TV_X
TV_A=gen$TV_A
l1_thr=gen$l1_thr


FPR_all=TPR_all=matrix(0,L_tau,nSim)
# A_prod=array(0,c())

cl <- makeCluster(ncores)
registerDoParallel(cl)

for (nn in 1:nSim)
{
  print(nn)
  A=TV_A[,,,nn]
  X=TV_X[,,nn]
  
  results_tmp<-foreach(itau=1:L_tau, .combine='acomb',.multicombine=TRUE,
                       .packages=c('lpSolve'))%:%
    foreach(it=2:T, .combine='cbind',.packages=c('lpSolve'))%dopar%{
      tau=tau_grid[itau]
      ker_weight_1=vc_weight_epanechnikovker((it-1)/T,1:T,T,b)
      s10=crossprod(X*drop(ker_weight_1),X)
      s10[s10<10^-3]=0
      tol=2*tau*base::norm(solve(s10),type="1")
      # tol=2*tau*l1_thr[it-1,nn]
      At=TV_A[,,it,nn]
      A_hat_tmp=tvvar_fit(it,TV_X[,,nn],tau,b)
      St=abs(At)>tol_true
      S_hat_tmp=abs(A_hat_tmp)>tol
      #St==0
      if (sum(St)==0){
        TPR_tmp=1
      }else{
        TPR_tmp=1-sum(St*(1-S_hat_tmp))/sum(St)
      }
      #St==1
      if (sum(St)==d*d){
        FPR_tmp=0
      }else{
        FPR_tmp=sum((1-St)*S_hat_tmp)/sum(1-St)
      }
      output=matrix(c(FPR_tmp,TPR_tmp),nrow=2,ncol=1)
    }
  FPR_all[,nn]=colMeans(results_tmp[1,t_lower:t_upper,])
  TPR_all[,nn]=colMeans(results_tmp[2,t_lower:t_upper,])
}
stopCluster(cl)
# TPR_all[is.na(TPR_all)]=0
FPR_all_mean=rowMeans(FPR_all)
TPR_all_mean=rowMeans(TPR_all)
plot(FPR_all_mean[2:30],TPR_all_mean[2:30])
plot(FPR_all[,1],TPR_all[,1])