
DVlsmh<-function(Niter,Y.i,Y.ia,M,N,D,Zs){

  mxzmatKL.new=matrix(NA,nrow=Niter,ncol=1)
  mxzmatKL.new[1,]=10000000
  
  
  Lambda.0=diag(D)*.0025
  Lambda.1=diag(D)*.0025
  
  
  alpha.0=as.numeric(glm(c(Y.i)~c(as.matrix(dist(Zs$Z.i)^2)))$coeff[1])
  alpha.1=as.numeric(glm(c(Y.ia)~c(-as.matrix(Zs$Z.i%*% t(Zs$Z.a))))$coeff[1])
  Z.i=Zs$Z.i
  Z.a=Zs$Z.a
  
  dis.new=.5
  i=2
  while( i<Niter ){

    rbZ.a=Z.a
    
    for(each in 2:nrow(Z.i)){rbZ.a=rbind(Z.a,rbZ.a)}
    
    rbZ.i=Z.i[rep(1:nrow(Z.i), times = rep(nrow(Z.a),nrow(Z.i))), ]
    dnj.ia <- t(sqrtm(Lambda.1)%*%t(rbZ.i)+sqrtm(Lambda.0)%*%t(rbZ.a))
    B <- solve(diag(nrow(Lambda.1))- 2* sqrtm(Lambda.0) %*% sqrtm(Lambda.1))
    C=solve(diag(D)+4*Lambda.0)
    
    
    
    cont.i <- dist(Z.i %*% chol(C))^2
    
    alpha.0=DVupalpha.0(alpha.0, C, Y.i, cont.i)
    alpha.1=DVupalpha.1(alpha.1, B, Y.ia, dnj.ia, Z.i, Z.a)

    Lambda.0=DVupLambda.0(alpha.0, alpha.1, C,B,  Y.i, cont.i, Z.i, Z.a, Lambda.0, Lambda.1 ,p.lambda.0=1)
    Lambda.0[upper.tri(Lambda.0)] <- t(Lambda.0)[upper.tri(t(Lambda.0))]
    
    if(!is.positive.definite(Lambda.0)){
      Lambda.0= nearPD(Lambda.0)$mat
    }
    
    B <- solve(diag(nrow(Lambda.1))- 2* sqrtm(Lambda.0) %*% sqrtm(Lambda.1))
    
    Lambda.1=DVupLambda.1.ia(alpha.1, B, D, M, Z.i,Z.a, Lambda.0, Lambda.1 ,p.lambda.1=1)

    Lambda.1[upper.tri(Lambda.1)] <- t(Lambda.1)[upper.tri(t(Lambda.1))]
    
    if(!is.positive.definite(Lambda.1)){
      Lambda.1= nearPD(Lambda.1)$mat
    }

    Z.i=DVgetZ.i.t(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i, Z.a,alpha.0,p.lambda.0=1,alpha.1)
    
    Z.a=DVgetZ.a.t.i.ia(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,p.lambda.1=1,alpha.1)
    mxzmatKL.new[i,1]=DVKL(alpha.0,  alpha.1, Lambda.0,Lambda.1, Z.i,Z.a, Y.i,Y.ia)

    dis.new=sum(mxzmatKL.new[i,])/sum(mxzmatKL.new[i-1,])

    i=i+1
  }
  list(lsmhEZ.i=Z.i,lsmhEZ.a=Z.a,lsmhVZ.0=Lambda.0,lsmhVZ.1=Lambda.1,lsmhAlpha.0=alpha.0,lsmhAlpha.1=alpha.1,
       lsmhKL=mxzmatKL.new[i-1,1])
}

