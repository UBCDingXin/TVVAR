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
# library(tseries)
library(pracma)

setwd("D:\\OneDrive\\Time-varying VAR\\code\\simulation\\R code")
source("VAR_lib_general.R")
load('D:\\OneDrive\\Time-varying VAR\\code\\Real data\\stock data\\stockdata.RData')
###############################################
# set parameters
##############################################
ncores=8
T_test=100
# tv-var and stat-var parameters
L_tau = 50
tau_tvvar_grid = seq(0.0001, 0.01, length.out=L_tau)
tau_statvar_grid=seq(0.0001, 0.5, length.out=L_tau)
# ls parameters
L_lam = 50
lam_grid = seq(0.0001, 0.05, length.out=L_lam)


# -------------------------------------------
# Necessary function
# -------------------------------------------
norm_vec <- function(x) sqrt(sum(x^2))/length(x)
norm_vec2 <- function(x) sqrt(sum(x^2))
acomb <- function(...) abind(..., along=3)#function for combination of array
acomb2 <- function(...) abind(..., along=4)#function for combination of array
enableJIT(3)
##########################################
#data processing
#########################################
stock.data = stock$data
# normalize data
stock.mean = colMeans(stock$data)
stock.sd=apply(stock$data, 2, sd)
T_raw = nrow(stock$data)
for(i in 1:T_raw)
{
  stock.data[i,] = (stock$data[i,] - stock.mean)/stock.sd
  # stock.data[i,] = (stock$data[i,] - stock.mean)
}

# indx=c(231 ,233 ,339  ,95 ,344 ,141 ,251 ,324 ,194 ,156)
indx=c(231 , 395, 63  ,95 ,344 ,141 ,251 ,324 ,194 ,156)
old_rawX=stock.data[,indx]
rawX=old_rawX
T_train=nrow(rawX)-T_test
b = 0.8* T_train^(-1/5)
d=ncol(rawX)
# ############################################
# detrend
# ############################################
for(i in 1:length(old_rawX[1,]))
{
  # rawX[,i]=detrend(old_rawX[,i], tt = 'linear',bp=round((nrow(old_rawX)-T_test)*b))
  rawX[,i]=detrend(old_rawX[,i], tt = 'linear',bp=round(50))
}


#####################################################
#fit data
#####################################################
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
    lam=lam_grid[ilam]
    X_train=rawX[it:(T_train+it-1),]
    lsvar_fit((T_train+1),X_train,lam,b)
  }
########
#Las
A_hat_las<-foreach(ilam=1:L_lam, .combine='acomb2',.multicombine=TRUE)%:%
  foreach(it=1:T_test, .combine='acomb')%dopar%{
    lam=lam_grid[ilam]
    X_train=rawX[it:(T_train+it-1),]
    lavar_fit((T_train+1),X_train,lam,b,alpha=0.0001,A_old=array(0.5,c(d,d)),A_oldold=array(0.5,c(d,d)))
  }
error_pred_tmp2<-foreach(ilam=1:L_lam, .combine='acomb',.multicombine=TRUE)%:%
  foreach(it=1:T_test, .combine='cbind')%dopar%{
    lam=lam_grid[ilam]
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
  A_hat_MLE=mlevar_fit((T_train+1),X_train,b)
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

t=proc.time() - ptm
stopCluster(cl)

err=c(min(colMeans(error_pred_tvvar_statvar[,,1])),
      min(colMeans(error_pred_tvvar_statvar[,,2])),
      min(colMeans(error_pred_Rid_las[,,1])),
      min(colMeans(error_pred_Rid_las[,,2])),
      mean(error_pred_MLE))
round(err,digits=3)
err_sd=c(sd(error_pred_tvvar_statvar[,which.min(colMeans(error_pred_tvvar_statvar[,,1])),1]),
         sd(error_pred_tvvar_statvar[,which.min(colMeans(error_pred_tvvar_statvar[,,2])),2]),
         sd(error_pred_Rid_las[,which.min(colMeans(error_pred_Rid_las[,,1])),1]),
         sd(error_pred_Rid_las[,which.min(colMeans(error_pred_Rid_las[,,2])),2]),
         sd(error_pred_MLE))
round(err_sd,digits=3)


save.image("stock_data_results.Rdata")
