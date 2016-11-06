#prediction plot
n_comp=27
lwdn=1.7
name=stock$info[,3][indx[n_comp]]
testX=rawX[(T_train+1):T_raw,]
plot(testX[,n_comp],type='o',xlab='Time',ylab='Detrended Close Price',main=name,lty=1,col=2,lwd=2)
lines(X_tvvar_pred[,n_comp],lty=1,col=4,lwd=1.85);points(X_tvvar_pred[,n_comp],pch=3)
lines(X_las_pred[,n_comp],lty=2,col=16,lwd=1.85)
lines(X_rid_pred[,n_comp],lty=6,col=36,lwd=lwdn)
lines(X_mle_pred[,n_comp],lty=4,col=6,lwd=lwdn)
lines(X_statvar_pred[,n_comp],lty=5,col=10,lwd=lwdn)

legend(-2, -0.75, c("True value","TV-VAR","Lasso", "Ridge",'MLE','Stat. VAR'),
       col = c(2, 4,16, 36,6,10),
       text.col = "black", lty = c(1, 1,2, 6,4,5),lwd=c(2,1.85,1.85,1.7,1.7,1.7),
       pch=c(1,3,NA,NA,NA,NA),
       merge = TRUE, bg = "gray90")



#prediction plot
n_comp=3
lwdn=1.7
currname=c('Euro/U.S.','U.K./U.S.','Australia/U.S.','New Zealand/U.S.','Canada/U.S.',
           'Singapore/U.S','Brazil/U.S.','Norway/U.S.')
name=paste('Exchange Rate:',currname[n_comp])
testX=rawX[1204:1303,]
plot(testX[,n_comp],type='o',xlab='Time',ylab='Detrended Exchange Rate',main=name,lty=1,col=2,lwd=2)
lines(X_tvvar_pred[,n_comp],lty=1,col=4,lwd=1.85);points(X_tvvar_pred[,n_comp],pch=3)
lines(X_las_pred[,n_comp],lty=2,col=16,lwd=1.85)
lines(X_rid_pred[,n_comp],lty=6,col=36,lwd=lwdn)
lines(X_mle_pred[,n_comp],lty=4,col=6,lwd=lwdn)
lines(X_statvar_pred[,n_comp],lty=5,col=10,lwd=lwdn)

legend(0, -0.3, c("True value","TV-VAR","Lasso", "Ridge",'MLE','Stat. VAR'),
       col = c(2, 4,16, 36,6,10),
       text.col = "black", lty = c(1, 1,2, 6,4,5),lwd=c(2,1.85,1.85,1.7,1.7,1.7),
       pch=c(1,3,NA,NA,NA,NA),
       merge = TRUE, bg = "gray90")




mean1=apply(error1,MARGIN=2, FUN=mean)
mean2=apply(error2,MARGIN=2, FUN=mean)
mean3=apply(error3,MARGIN=2, FUN=mean)
mean4=apply(error4,MARGIN=2, FUN=mean)
mean_all=cbind(mean1[1],mean2[1],mean3,mean4[1])
sd1=apply(error1,MARGIN=2, FUN=sd)
sd2=apply(error2,MARGIN=2, FUN=sd)
sd3=apply(error3,MARGIN=2, FUN=sd)
sd4=apply(error4,MARGIN=2, FUN=sd)
sd_all=cbind(sd1[1],sd2[1],sd3,sd4[1])
mean_all
sd_all

#finance
#ajacency matrix
library(igraph)

# adjmatrix1= A_hat_pred[,,25,1]
# adjmatrix2= A_hat_pred[,,50,1]
# adjmatrix3= A_hat_pred[,,75,1]
# adjmatrix4= A_hat_pred[,,100,1]
no_lam=2#which lambda we will use
adjmatrix1= A_hat_tvvar[,,25,no_lam]
adjmatrix2= A_hat_tvvar[,,65,no_lam]
adjmatrix3= A_hat_tvvar[,,85,no_lam]
adjmatrix4= A_hat_tvvar[,,100,no_lam]
sector.info=stock$info[,2]
indx1=which(stock$info[,2]=="Consumer Staples")
indx2=which(stock$info[,2]=="Consumer Discretionary")
indx3=which(stock$info[,2]=="Industrials")
indx4=which(stock$info[,2]=="Financials")
indx5=which(stock$info[,2]=="Utilities")
indx6=which(stock$info[,2]=="Energy")
sector.info[indx1]="CS"
sector.info[indx2]="CD"
sector.info[indx3]="IND"
sector.info[indx4]="FIN"
sector.info[indx5]="UTI"
sector.info[indx6]="ENE"
coln=c()
for (i in 1:10){
  coln[i]=paste(stock$info[,3][indx[i]],"(",sector.info[indx[i]],")",sep="")
}
# coln=c(stock$info[,3][indx])
# coln=c('Euro','U.K.','Australia','New Zealand','Canada','Singapore','Brazil','Norway')
roln=coln
colnames(adjmatrix1)=coln
colnames(adjmatrix2)=coln
colnames(adjmatrix3)=coln
colnames(adjmatrix4)=coln
z1=graph.adjacency(adjmatrix1, mode=c("undirected"), weighted=TRUE, diag=FALSE,
                   add.colnames=NULL, add.rownames=NA)
z2=graph.adjacency(adjmatrix2, mode=c("undirected"), weighted=TRUE, diag=FALSE,
                   add.colnames=NULL, add.rownames=NA)
z3=graph.adjacency(adjmatrix3, mode=c("undirected"), weighted=TRUE, diag=FALSE,
                   add.colnames=NULL, add.rownames=NA)
z4=graph.adjacency(adjmatrix4, mode=c("undirected"), weighted=TRUE, diag=FALSE,
                   add.colnames=NULL, add.rownames=NA)
z1=set_vertex_attr(z1, 'sector', 
                   value=c(1,2,3,4,4,5,4,1,4,6))
z2=set_vertex_attr(z2, 'sector', 
                   value=c(1,2,3,4,4,5,4,1,4,6))
z3=set_vertex_attr(z3, 'sector', 
                   value=c(1,2,3,4,4,5,4,1,4,6))
z4=set_vertex_attr(z4, 'sector', 
                   value=c(1,2,3,4,4,5,4,1,4,6))

V(z1)[V(z1)$sector==1]$color <- 1;V(z1)[V(z1)$sector==2]$color <- 2
V(z1)[V(z1)$sector==3]$color <- 3;V(z1)[V(z1)$sector==4]$color <- 4
V(z1)[V(z1)$sector==5]$color <- 5;V(z1)[V(z1)$sector==6]$color <- 6

V(z2)[V(z2)$sector==1]$color <- 1;V(z2)[V(z2)$sector==2]$color <- 2
V(z2)[V(z2)$sector==3]$color <- 3;V(z2)[V(z2)$sector==4]$color <- 4
V(z2)[V(z2)$sector==5]$color <- 5;V(z2)[V(z2)$sector==6]$color <- 6

V(z3)[V(z3)$sector==1]$color <- 1;V(z3)[V(z3)$sector==2]$color <- 2
V(z3)[V(z3)$sector==3]$color <- 3;V(z3)[V(z3)$sector==4]$color <- 4
V(z3)[V(z3)$sector==5]$color <- 5;V(z3)[V(z3)$sector==6]$color <- 6

V(z4)[V(z4)$sector==1]$color <- 1;V(z4)[V(z4)$sector==2]$color <- 2
V(z4)[V(z4)$sector==3]$color <- 3;V(z4)[V(z4)$sector==4]$color <- 4
V(z4)[V(z4)$sector==5]$color <- 5;V(z4)[V(z4)$sector==6]$color <- 6

# fin_adja_t25
plot.igraph(z1,vertex.label=V(z1)$name,
            layout=layout.fruchterman.reingold, 
            edge.color="black",edge.width=E(z1)$weight,
            vertex.color=V(z1)$color)

plot.igraph(z2,vertex.label=V(z2)$name,
            layout=layout.fruchterman.reingold, 
            edge.color="black",edge.width=E(z2)$weight,
            vertex.color=V(z2)$color)
plot.igraph(z3,vertex.label=V(z3)$name,
            layout=layout.fruchterman.reingold, 
            edge.color="black",edge.width=E(z3)$weight,
            vertex.color=V(z3)$color)
plot.igraph(z4,vertex.label=V(z4)$name,
            layout=layout.fruchterman.reingold, 
            edge.color="black",edge.width=E(z4)$weight,
            vertex.color=V(z4)$color)



