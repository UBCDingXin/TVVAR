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
library(abind)
library(compiler)
library(mvtnorm)
library(pracma)

setwd("C:\\Users\\Xin\\OneDrive\\Time-varying VAR\\code\\simulation\\R code")
source("VAR_lib_general.R")
load('C:\\Users\\Xin\\OneDrive\\Time-varying VAR\\code\\Real data\\exchange rate\\exchange_rate_all.RData')

###############################################
# set parameters
##############################################
ncores=8
T_test=50
# tv-var and stat-var parameters
L_tau = 10
tau_tvvar_grid = seq(0.01, 0.1, length.out=L_tau)
tau_statvar_grid=seq(0.01, 0.8, length.out=L_tau)
# ls parameters
L_lam = 10
lam_grid_rid = seq(0.0003, 0.001, length.out=L_lam)
lam_grid_las = seq(0.0003, 0.002, length.out=L_lam)
# -------------------------------------------
# Necessary function
# -------------------------------------------
norm_vec <- function(x) sqrt(sum(x^2))/length(x)
norm_vec2 <- function(x) sqrt(sum(x^2))
# norm_vec2 <- function(x) max(abs(x))
acomb <- function(...) abind(..., along=3)#function for combination of array
acomb2 <- function(...) abind(..., along=4)#function for combination of array
cal_dis_vec<-function(x){
  leng_x=length(x)
  x2=c(0,x[1:(leng_x-1)])
  diff=abs(x-x2)
  return(diff[2:leng_x])
}
enableJIT(3)
##########################################
#data processing
#########################################
exchange_rate[,1]=NULL

# indx = c(1,2,7,4,5,6,13,9,10,21);
# indx = c(3,4,6,7,8,9,11,20);
# indx = c(3,4,5,6,7,8,9,11,12,13,14,15,17,20,21,22);
indx= c(3,4,5,6,7,8,9,11,12,13,14,15,17,20,21);
stock50.data = exchange_rate[1:100,indx]

stock50.data = as.matrix(stock50.data )

# Centering data
stock50.mean = colMeans(stock50.data)
stock50.sd=apply(stock50.data, 2, sd)
T = nrow(stock50.data)
for(i in 1:T)
{
  stock50.data[i,] = (stock50.data[i,] - stock50.mean)/stock50.sd
  # stock50.data[i,] = (stock50.data[i,] - stock50.mean)
}

X=stock50.data;
rawX=X
T_train=nrow(rawX)-T_test
# b = 0.8* T_train^(-1/5)
b=0.3
b_mle=0.4
# b = 2* T_train^(-1/5)
# d=ncol(rawX)
# ############################################
# detrend
# ############################################
for(i in 1:length(X[1,]))
{
  # rawX[,i]=detrend(X[,i], tt = 'linear',bp=round(1158*0.19))
  rawX[,i]=detrend(X[,i], tt = 'linear',bp=round(20))
}

# p_val=c()
# for (i in 1:length(indx)){
#   fit=ar(rawX[,i], aic = FALSE, order.max = 1, method = "yw")
#   BOXtest=Box.test(fit$resid,lag=1,type="Ljung-Box")
#   p_val[i]=BOXtest$p.value
# }
# p_val
#####################################################
#fit data
#####################################################
d=ncol(rawX)
error_pred_tvvar_statvar=array(0,c(T_test,L_tau,2))
error_pred_Rid_las=array(0,c(T_test,L_lam,2))
error_pred_MLE=array(0,T_test)
A_hat_tvvar=A_hat_statvar=array(0,c(d,d,T_test,L_tau))
A_hat_rid=A_hat_las=array(0,c(d,d,T_test,L_lam))
A_hat_mle=array(0,c(d,d,T_test))
X_tvvar_pred=X_statvar_pred=X_rid_pred=X_las_pred=X_mle_pred=array(0,c(T_test,d))


cl <- makeCluster(ncores)
registerDoParallel(cl)
ptm <- proc.time()

###########################################
#VAR
###########################################
########
#TVVAR
A_hat_tvvar<-foreach(itau=1:L_tau,.combine='acomb2',.multicombine=TRUE,.packages=c("lpSolve"))%:%
  foreach(it=1:T_test, .combine='acomb',.packages=c('lpSolve'))%dopar%{
    tau=tau_tvvar_grid[itau]
    #traning data
    X_train=rawX[it:(T_train+it-1),]
    tvvar_fit((T_train+1),X_train,tau,b)
  }
########
#STAT-VAR
A_hat_statvar<-foreach(itau=1:L_tau,.combine='acomb2',.multicombine=TRUE,.packages=c("lpSolve"))%:%
  foreach(it=1:T_test, .combine='acomb',.packages=c('lpSolve'))%dopar%{
    tau=tau_statvar_grid[itau]
    #traning data
    X_train=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
    stvar_fit(X_train,tau)
  }
##################
#error for tvvar and statvar
error_pred_tmp1<-foreach(itau=1:L_tau, .combine='acomb',.multicombine=TRUE,.packages=c('lpSolve'))%:%
  foreach(it=1:T_test, .combine='cbind',.packages=c('lpSolve'))%dopar%{
    #traning data for tvvar
    X_train1=rawX[it:(T_train+it-1),]
    #traning data for statvar
    X_train2=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
    X_test=rawX[T_train+it,]
    #TVVAR estimation
    A_hat_tvvar_tmp=A_hat_tvvar[,,it,itau]
    X_pred_tvvar=A_hat_tvvar_tmp%*%X_train1[T_train,]
    #statVAR estimation
    A_hat_statvar_tmp=A_hat_statvar[,,it,itau]
    X_pred_statvar=A_hat_statvar_tmp%*%X_train2[nrow(X_train2),]
    # error_tmp=c(norm_vec(X_pred_tvvar-X_test),norm_vec(X_pred_statvar-X_test))
    error_tmp=c(norm_vec2(X_pred_tvvar-X_test),norm_vec2(X_pred_statvar-X_test))
  }
error_pred_tvvar_statvar=aperm(error_pred_tmp1, c(2,3,1))
tau_tvvar_indx=which.min(colMeans(error_pred_tvvar_statvar[,,1]))
tau_statvar_indx=which.min(colMeans(error_pred_tvvar_statvar[,,2]))
######
#tvvar prediction
X_tvvar_pred<-foreach(it=1:T_test, .combine='rbind',.packages=c('lpSolve'))%dopar%{
  X_train=rawX[it:(T_train+it-1),]
  matrix(A_hat_tvvar[,,it,tau_tvvar_indx]%*%X_train[T_train,],nrow=1)
}
######
#statvar prediction 
X_statvar_pred<-foreach(it=1:T_test, .combine='rbind',.packages=c('lpSolve'))%dopar%{
  X_train=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
  matrix(A_hat_statvar[,,it,tau_statvar_indx]%*%X_train[nrow(X_train),],nrow=1)
}

##############################
#LS
##############################
########
#Rid
A_hat_rid<-foreach(ilam=1:L_lam, .combine='acomb2',.multicombine=TRUE)%:%
  foreach(it=1:T_test, .combine='acomb')%dopar%{
    lam=lam_grid_rid[ilam]
    X_train=rawX[it:(T_train+it-1),]
    lsvar_fit((T_train+1),X_train,lam,b)
  }
########
#Las
A_hat_las<-foreach(ilam=1:L_lam, .combine='acomb2',.multicombine=TRUE)%:%
  foreach(it=1:T_test, .combine='acomb')%dopar%{
    lam=lam_grid_las[ilam]
    X_train=rawX[it:(T_train+it-1),]
    lavar_fit((T_train+1),X_train,lam,b,alpha=0.0001,A_old=array(1,c(d,d)),A_oldold=array(0.5,c(d,d)))
  }
error_pred_tmp2<-foreach(ilam=1:L_lam, .combine='acomb',.multicombine=TRUE)%:%
  foreach(it=1:T_test, .combine='cbind')%dopar%{
    # lam=lam_grid[ilam]
    X_train=rawX[it:(T_train+it-1),]
    X_test=rawX[T_train+it,]
    #Ridge estimation
    A_hat_Rid_tmp=A_hat_rid[,,it,ilam]
    X_pred_Rid=A_hat_Rid_tmp%*%X_train[T_train,]
    #Lasso estimation
    A_hat_las_tmp=A_hat_las[,,it,ilam]
    X_pred_las=A_hat_las_tmp%*%X_train[T_train,]
    # error_tmp=c(norm_vec(X_pred_Rid-X_test),norm_vec(X_pred_las-X_test))
    error_tmp=c(norm_vec2(X_pred_Rid-X_test),norm_vec2(X_pred_las-X_test))
  }
error_pred_Rid_las=aperm(error_pred_tmp2, c(2,3,1))
lam_rid_indx=which.min(colMeans(error_pred_Rid_las[,,1]))
lam_las_indx=which.min(colMeans(error_pred_Rid_las[,,2]))
######
#ridge prediction
X_rid_pred<-foreach(it=1:T_test, .combine='rbind',.packages=c('lpSolve'))%dopar%{
  X_train=rawX[it:(T_train+it-1),]
  matrix(A_hat_rid[,,it,lam_rid_indx]%*%X_train[T_train,],nrow=1)
}
######
#lasso prediction 
X_las_pred<-foreach(it=1:T_test, .combine='rbind',.packages=c('lpSolve'))%dopar%{
  X_train=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
  matrix(A_hat_las[,,it,lam_las_indx]%*%X_train[nrow(X_train),],nrow=1)
}


#################################################
#MLE
###############################################
A_hat_mle<-foreach(it=1:T_test, .combine='acomb')%dopar%{
  X_train=rawX[it:(T_train+it-1),]
  X_test=rawX[T_train+it,]
  #mle
  A_hat_MLE=mlevar_fit((T_train+1),X_train,b_mle)
}
error_pred_MLE<-foreach(it=1:T_test, .combine='cbind')%dopar%{
  X_train=rawX[it:(T_train+it-1),]
  X_test=rawX[T_train+it,]
  #mle
  X_pred_MLE=A_hat_mle[,,it]%*%X_train[T_train,]
  # norm_vec(X_pred_MLE-X_test)
  norm_vec2(X_pred_MLE-X_test)
}
X_mle_pred<-foreach(it=1:T_test, .combine='rbind')%dopar%{
  X_train=rawX[it:(T_train+it-1),]
  matrix(A_hat_mle[,,it]%*%X_train[T_train,],nrow=1)
}

######################
#naive
#######################
error_pred_naive<-foreach(it=1:T_test, .combine='cbind')%dopar%{
  X_train=rawX[it:(T_train+it-1),]
  X_test=rawX[T_train+it,]
  #mle
  norm_vec2(X_train[T_train,]-X_test)
}


t=proc.time() - ptm
stopCluster(cl)

err=c(min(colMeans(error_pred_tvvar_statvar[,,1])),
      min(colMeans(error_pred_tvvar_statvar[,,2])),
      min(colMeans(error_pred_Rid_las[,,1])),
      min(colMeans(error_pred_Rid_las[,,2])),
      mean(error_pred_MLE),
      mean(error_pred_naive))
round(err,digits=4)
err_sd=c(sd(error_pred_tvvar_statvar[,which.min(colMeans(error_pred_tvvar_statvar[,,1])),1]),
         sd(error_pred_tvvar_statvar[,which.min(colMeans(error_pred_tvvar_statvar[,,2])),2]),
         sd(error_pred_Rid_las[,which.min(colMeans(error_pred_Rid_las[,,1])),1]),
         sd(error_pred_Rid_las[,which.min(colMeans(error_pred_Rid_las[,,2])),2]),
         sd(error_pred_MLE))
round(err_sd,digits=4)


save.image("exchange_data_results.Rdata")