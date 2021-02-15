
DgetAUC<-function(est.alpha.1,Z.i,Z.a,D,M,N, Y.ia){



  est.P.ia=matrix(NA,N,M)




  for(i in 1:N){
    for(a in 1:M){
      est.P.ia[i,a]=inv.logit(est.alpha.1-sum((Z.i[i,] - Z.a[a,]) ^ 2))
    }
  }



  c=roc(c(Y.ia),c(est.P.ia),auc.polygon=FALSE, grid=FALSE,
        ylab="true positive rate",xlab="false positive rate",xlim=c(1,0),plot=TRUE,auc=TRUE,main=NULL, font.axis = 2,bty="n", cex.axis =1.2,ann=FALSE,  col="black",legacy.axes = TRUE, lwd=5)
  text(x=0.45, y=.1,cex=1.6, labels=paste("AUC = ", as.character(round(as.numeric(c$auc),digits = 4)),sep = ""), font= 2, col="black")

  return(c$auc)

}

