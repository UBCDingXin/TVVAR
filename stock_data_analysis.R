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
# library(RSEIS)
library(pracma)



setwd("C:\\Users\\Xin\\OneDrive\\Time-varying VAR\\code\\simulation\\R code")
source("VAR_lib_general.R")
load('C:\\Users\\Xin\\OneDrive\\Time-varying VAR\\code\\Real data\\stock data\\stockdata.RData')
###############################################
# set parameters
##############################################
ncores=8
T_test=100
# tv-var and stat-var parameters
L_tau = 15
tau_tvvar_grid = seq(0.005, 0.1, length.out=L_tau)
tau_statvar_grid=seq(0.001, 0.1, length.out=L_tau)
# ls parameters
L_lam = 15
lam_grid = seq(0.2, 1, length.out=L_lam)


# -------------------------------------------
# Necessary function
# -------------------------------------------
norm_vec <- function(x) sqrt(sum(x^2))/length(x)
norm_vec2 <- function(x) sqrt(sum(x^2))
# norm_vec2 <- function(x) max(abs(x))
# norm_vec2 <- function(x) sum(abs(x))
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
stock.data = (stock$data)[1:150,]
# stock.data = (stock$data)
# min_val=apply(stock.data,MARGIN=2,FUN=min)
# min_mat=matrix(min_val,byrow=TRUE,nrow=dim(stock.data)[1],ncol=dim(stock.data)[2])
# stock.data=log(stock.data+1+abs(min_mat))
# normalize data
stock.mean = colMeans(stock$data)
stock.sd=apply(stock$data, 2, sd)
T_raw = nrow(stock.data)
for(i in 1:T_raw)
{
  stock.data[i,] = (stock$data[i,] - stock.mean)/stock.sd
}

# indx=c(231 ,233 ,339  ,95 ,344 ,141 ,251 ,324 ,194 ,156) #0.197 1.121 0.221 0.195 0.216
# indx=c(231 , 395, 63  ,95 ,344 ,141 ,251 ,324 ,194 ,156)
# indx=c(231, 395 , 63,233 ,339  ,95 ,344 ,141 ,251 ,324 ,194 ,156)
# indx=c(2,3,7,13,16,26,44,47,60,63,64,79,95,128,136,141,147,156,194,205,207,231,233,251,324,344,395)
# indx=sort(stock.sd,decreasing = TRUE,index.return=TRUE)$ix[2:51]

# old_rawX=stock.data[,indx]
old_rawX=stock.data
rawX=old_rawX
T_train=nrow(rawX)-T_test
# b = 0.8* T_train^(-1/5)
b = 0.3
b_mle=0.5
# b_mle=0.8* T_train^(-1/5)
# b=b_mle=1
# d=ncol(rawX)
# min_val=apply(rawX,MARGIN=2,FUN=min)
# min_mat=matrix(min_val,byrow=TRUE,nrow=dim(rawX)[1],ncol=dim(rawX)[2])
# rawX=log(rawX+1+abs(min_mat))
############################################
#detrend
############################################
for(i in 1:length(old_rawX[1,]))
{
  rawX[,i]=pracma::detrend(old_rawX[,i], tt = 'linear',bp=50)
}
# select stocks
# p_val=c()
# for (i in 1:dim(rawX)[2]){
#   fit=ar(rawX[,i], aic = FALSE, order.max = 1, method = "yw")
#   BOXtest=Box.test(fit$resid,lag=1,type="Ljung-Box")
#   p_val[i]=BOXtest$p.value
# }
# indx_tmp=which(p_val<0.05)
# all_max=c()
# for (i in 1:length(indx_tmp)){
#   all_max[i]=max(cal_dis_vec(rawX[,indx_tmp[i]]))
# }
# indx=indx_tmp[which(all_max<0.7)]#0.3

all_max=c()
for (i in 1:dim(rawX)[2]){
  all_max[i]=max(cal_dis_vec(rawX[,i]))
}
indx_tmp=which(all_max<0.4)
p_val=c()
for (i in 1:length(indx_tmp)){
  fit=ar(rawX[,indx_tmp[i]], aic = FALSE, order.max = 1, method = "yule-walker")
  BOXtest=Box.test(fit$resid,lag=1,type="Ljung-Box")
  p_val[i]=BOXtest$p.value
}
indx=indx_tmp[which(p_val<0.05)]
indx=unique(c(indx,c(231 , 395, 63  ,95 ,344 ,141 ,251 ,324 ,194 ,156)))


rawX=rawX[,indx]

#####################################################
#fit data
#####################################################
d=ncol(rawX)
error_pred_tvvar_statvar=array(0,c(T_test,L_tau,2))
error_pred_Rid_las=array(0,c(T_test,L_lam,2))
error_pred_stat_Rid_las=array(0,c(T_test,L_lam,2))
error_pred_MLE=array(0,T_test)
A_hat_tvvar=A_hat_statvar=array(0,c(d,d,T_test,L_tau))
A_hat_rid=A_hat_las=A_hat_statrid=A_hat_statlas=array(0,c(d,d,T_test,L_lam))
A_hat_mle=array(0,c(d,d,T_test))
X_tvvar_pred=X_statvar_pred=X_rid_pred=X_las_pred=X_mle_pred=X_statlas_pred=X_statrid_pred=array(0,c(T_test,d))


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
  X_train=rawX[it:(T_train+it-1),]
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

# ##############################
# #Stationary LS
# ##############################
# ########
# #stat Rid
# A_hat_statrid<-foreach(ilam=1:L_lam, .combine='acomb2',.multicombine=TRUE)%:%
#   foreach(it=1:T_test, .combine='acomb')%dopar%{
#     lam=lam_grid[ilam]
#     X_train=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
#     statlsvar_fit(X_train,lam)
#   }
# ########
# #stat Las
# A_hat_statlas<-foreach(ilam=1:L_lam, .combine='acomb2',.multicombine=TRUE)%:%
#   foreach(it=1:T_test, .combine='acomb')%dopar%{
#     lam=lam_grid[ilam]
#     X_train=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
#     statlavar_fit(X_train,lam,alpha=0.0001,A_old=array(0.5,c(d,d)),A_oldold=array(0.5,c(d,d)))
#   }
# error_pred_tmp3<-foreach(ilam=1:L_lam, .combine='acomb',.multicombine=TRUE)%:%
#   foreach(it=1:T_test, .combine='cbind')%dopar%{
#     lam=lam_grid[ilam]
#     X_train=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
#     X_test=rawX[T_train+it,]
#     #Ridge estimation
#     A_hat_statRid_tmp=A_hat_statrid[,,it,ilam]
#     X_pred_statRid=A_hat_statRid_tmp%*%X_train[nrow(X_train),]
#     #Lasso estimation
#     A_hat_statlas_tmp=A_hat_statlas[,,it,ilam]
#     X_pred_statlas=A_hat_statlas_tmp%*%X_train[nrow(X_train),]
#     error_tmp=c(norm_vec2(X_pred_statRid-X_test),norm_vec2(X_pred_statlas-X_test))
#   }
# error_pred_stat_Rid_las=aperm(error_pred_tmp3, c(2,3,1))
# lam_stat_rid_indx=which.min(colMeans(error_pred_stat_Rid_las[,,1]))
# lam_stat_las_indx=which.min(colMeans(error_pred_stat_Rid_las[,,2]))
# ######
# #stat ridge prediction
# X_statrid_pred<-foreach(it=1:T_test, .combine='rbind',.packages=c('lpSolve'))%dopar%{
#   X_train=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
#   matrix(A_hat_statrid[,,it,lam_stat_rid_indx]%*%X_train[nrow(X_train),],nrow=1)
# }
# ######
# #stat lasso prediction 
# X_statlas_pred<-foreach(it=1:T_test, .combine='rbind',.packages=c('lpSolve'))%dopar%{
#   X_train=rawX[(T_train+it-round(T_train*b-1)):(T_train+it-1),]
#   matrix(A_hat_statlas[,,it,lam_stat_las_indx]%*%X_train[nrow(X_train),],nrow=1)
# }


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
      min(colMeans(error_pred_stat_Rid_las[,,1])),
      min(colMeans(error_pred_stat_Rid_las[,,2])),
      mean(error_pred_naive))
round(err,digits=4)
err_sd=c(sd(error_pred_tvvar_statvar[,which.min(colMeans(error_pred_tvvar_statvar[,,1])),1]),
         sd(error_pred_tvvar_statvar[,which.min(colMeans(error_pred_tvvar_statvar[,,2])),2]),
         sd(error_pred_Rid_las[,which.min(colMeans(error_pred_Rid_las[,,1])),1]),
         sd(error_pred_Rid_las[,which.min(colMeans(error_pred_Rid_las[,,2])),2]),
         sd(error_pred_MLE),
         sd(error_pred_naive))
round(err_sd,digits=4)


save.image("stock_data_results.Rdata")
