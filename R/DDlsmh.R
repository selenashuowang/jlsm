DDlsmh<-function(Niter,Y.i,Y.ia,M,N,D,Zs){

  mxzmatKL.new=matrix(NA,nrow=Niter,ncol=1)
  mxzmatKL.new[1,]=10000000


  Lambda.0=diag(D)
  Lambda.1=diag(D)
  Z.i=Zs$Z.i
  Z.a=Zs$Z.a

  alpha.0=as.numeric(glm(c(Y.i)~c(as.matrix(dist(Z.i)^2)))$coeff[1])
  alpha.1=as.numeric(glm(c(Y.ia)~c(as.matrix(Zs$Z.i%*% t(Zs$Z.a))))$coeff[1])
  dis.new=.5
  i=2
  while( i<Niter ){

    B=solve(diag(D) + 2*Lambda.0 + 2 * Lambda.1)
    C=solve(diag(D)+4*Lambda.0)

    cont.i <- dist(Z.i %*% chol(C))^2


    alpha.0=DDupalpha.0(alpha.0, C, Y.i, cont.i)
    alpha.1=DDupalpha.1(alpha.1, B, Y.ia, Z.i, Z.a)
    Lambda.0=DDupLambda.0(alpha.0, alpha.1, C,B,  Y.i, Y.ia, cont.i, Z.i, Z.a, Lambda.0, Lambda.1 ,p.lambda.0=1)

   
    B=solve(diag(D) + 2*Lambda.0 + 2 * Lambda.1)

    Lambda.1=DDupLambda.1.ia( alpha.1, B, Y.ia,  Z.i, Z.a, Lambda.0, Lambda.1 ,p.lambda.1=1)

 

    Z.i=DDgetZ.i.t(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i, Z.a,alpha.0,p.lambda.0=1,alpha.1)

    times=1

    while (max(abs(range(Z.i)))>10 & times<20){
      Z.i=DDgetZ.i.t(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i, Z.a,alpha.0,p.lambda.0=1,alpha.1)
      times=times+1

    }

    Z.a=DDgetZ.a.t.i.ia(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,p.lambda.1=1,alpha.1)

    times=1

    while (max(abs(range(Z.a)))>10 & times<20){
      Z.a=DDgetZ.a.t.i.ia(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,p.lambda.1=1,alpha.1)
      times=times+1

    }


    mxzmatKL.new[i,1]=DDKL(alpha.0,  alpha.1, Lambda.0,Lambda.1, Z.i,Z.a, Y.i,Y.ia)
    dis.new=sum(mxzmatKL.new[i,1])/sum(mxzmatKL.new[i-1,1])
    i=i+1
  }
  list(lsmhEZ.i=Z.i,lsmhEZ.a=Z.a,lsmhVZ.0=Lambda.0,lsmhVZ.1=Lambda.1,lsmhAlpha.0=alpha.0,lsmhAlpha.1=alpha.1,
       lsmhKL=mxzmatKL.new[i-1,1])
}
