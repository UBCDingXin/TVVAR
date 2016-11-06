rm(list=ls())
require(MESS)
pattern="hub"
nfolder=7
folder=c("\\original1","\\uv1","\\large_b1","\\uv2","\\uv3","\\uv4","\\uv5","\\uv6","\\uv7")
#uv1:b=0.31; one side;0.001:0.3
#uv2:b=0.31; two side;0.001:0.3
#uv3:b=0.31; two side;0.0001:0.2
#uv4:b=0.31; two side;0.001:0.25
#uv5:b=0.31; two side:0,001:0.8; threshold=0.001
#uv6:b=0.31; two side:0,001:0.45; threshold=0.001
#uv7:b=0.31; two side:0.0001:0.05;threshold=rule out small entry in s10 and use the formula

# path=paste("D:\\OneDrive\\Time-varying VAR\\code\\result\\final results\\ROC\\tmp2\\",pattern,sep="")
# path=paste("D:\\OneDrive\\Time-varying VAR\\code\\result\\final results\\ROC\\",pattern,folder[nfolder],sep="")
path=paste("C:\\Users\\Xin\\OneDrive\\Time-varying VAR\\code\\result\\final results\\ROC\\",pattern,folder[nfolder],sep="")
setwd(path)
all_files=list.files()
FPR_plot=array(0,c(29,4))
TPR_plot=array(0,c(29,4))

for (i in 1:4){
  load(all_files[i])
  FPR_all_mean=rowMeans(FPR_all)
  TPR_all_mean=rowMeans(TPR_all)
  FPR_plot[,i]=FPR_all_mean[2:30]
  TPR_plot[,i]=TPR_all_mean[2:30]
}

plot(NULL,xlab="False Positive Rate (1-Specificity)", ylab='True Positive Rate (Sensitivity)',main=pattern,
     ylim=c(0,0.9),xlim=c(0,0.9))
lines(TPR_plot[,1]~FPR_plot[,1],lty=1,lwd=2.5,col='red');points(TPR_plot[,1]~FPR_plot[,1],pch=1,col='red')
lines(TPR_plot[,2]~FPR_plot[,2],lty=1,lwd=1.5,col='blue');points(TPR_plot[,2]~FPR_plot[,2],pch=2,col='blue')
lines(TPR_plot[,3]~FPR_plot[,3],lty=1,lwd=1.5,col='lightsalmon');points(TPR_plot[,3]~FPR_plot[,3],pch=3,col='lightsalmon')
lines(TPR_plot[,4]~FPR_plot[,4],lty=1,lwd=1.5,col='gray8');points(TPR_plot[,4]~FPR_plot[,4],pch=4,col='gray8')
abline(a = 0, b = 1,lty=2,col="black")
legend(0.71, 0.23, c("d=20","d=30", "d=40","d=50"),
       col = c('red', 'blue', "lightsalmon","gray8"),
       text.col = "black", lty = c(1, 1, 1, 1),lwd=c(2.5,1.5,1.5,1.5),pch=c(1,2,3,4),
       merge = TRUE, bg = "gray90",cex=1.1)



