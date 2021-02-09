DVgetZ.a.t.i.ia<-function(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,p.lambda.1,alpha.1){

  B <- solve(diag(nrow(Lambda.1))- 2* sqrtm(Lambda.0) %*% sqrtm(Lambda.1))
  
  
  
  for (n in 1:nrow(Z.a)){
    
    
    Z.a.n=matrix(Z.a[n,],byrow=TRUE,nrow = nrow(Z.i),ncol = nrow(Lambda.1))
    
    dnj.ia <- t(sqrtm(Lambda.1) %*% t(Z.i) + sqrtm(Lambda.0) %*% t(Z.a.n))
    
    
    
    A.ia <- 1/sqrt(det(B)) * exp(- alpha.1) * exp( - apply(as.matrix(Z.i), 1, function(x) x %*% Z.a[n,])
                                                   - .5 * apply(dnj.ia, 1, function(x) x %*% B %*% x))   
    f1znTo.ia <-  colSums(Z.i / (1+A.ia)) + 
      sqrtm(Lambda.0) %*% (B + t(B)) %*% colSums(dnj.ia / (1+A.ia))
    
    f2znTo.ia <- .5 *  sqrtm(Lambda.0) %*% (B + t(B)) %*% sqrtm(Lambda.0) * sum(1 / (1 + A.ia))+
      .25 *  sqrtm(Lambda.0)%*%  (B + t(B)) %*% (t(dnj.ia / (2 + 1 / A.ia + A.ia)) %*% dnj.ia)  %*%  (B + t(B)) %*% sqrtm(Lambda.0) +
      .5 * (t(Z.i / (2 + 1 / A.ia + A.ia)) %*% dnj.ia) %*%  (B + t(B)) %*% sqrtm(Lambda.0) +
      .5 * sqrtm(Lambda.0)  %*% t((B + t(B)) %*% (t(dnj.ia / (2 + 1 / A.ia + A.ia)) %*% Z.i)) +
      (t(Z.i / (2 + 1 / A.ia + A.ia)) %*% Z.i)

    numZnT <-   t(Z.a[n, ]) %*% (.5*f2znTo.ia) +.5*colSums(as.matrix(Z.i)* (Y.ia[,n]))-.5 * t(f1znTo.ia)
    denZnT <- ( 1 / (2 * p.lambda.1^2)) * diag(nrow(Lambda.1)) +.5*f2znTo.ia
    
    Z.a[n, ]<- numZnT %*% solve(denZnT)
  }
  
  Z.a
}