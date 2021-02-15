
VVgetZ.i.t.vector<-function(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,alpha.0, p.lambda.0,alpha.1){
  
  
  B <- solve(diag(nrow(Lambda.1))- 2* sqrtm(Lambda.0) %*% sqrtm(Lambda.1))
  
  B.i <- solve(diag(nrow(Lambda.1))- 2* Lambda.0 )
  
  for (n in 1:nrow(Z.i)){
    
    
    Z.i.n.i=matrix(Z.i[n,] , byrow=TRUE, nrow=(nrow(Z.i)-1), ncol=nrow(Lambda.0))
    
    dnj.i <- t( t(Z.i.n.i) + t(Z.i[-n,]))
    
    
    
    A.i <- 1/sqrt(det(B.i)) * exp(- alpha.0 ) * exp( - apply(as.matrix(Z.i[-n,]), 1, function(x) x %*% Z.i[n,])
                                                     - .5 * apply(dnj.i, 1, function(x) x %*% sqrtm(Lambda.0) %*% B.i  %*% sqrtm(Lambda.0) %*% x))
    
    
    f1znTo.i <-  colSums(Z.i[-n,] / (1+A.i)) + 
      .5 *  sqrtm(Lambda.0) %*% (B.i + t(B.i)) %*% sqrtm(Lambda.0) %*% colSums(dnj.i / (1+A.i))
    
    
    f2znTo.i <- .5 *  sqrtm(Lambda.0) %*% (B.i + t(B.i)) %*% sqrtm(Lambda.0)  * sum(1 / (1 + A.i))+
      .25 * sqrtm(Lambda.0) %*% (B.i + t(B.i)) %*% sqrtm(Lambda.0) %*% (t(dnj.i / (2 + 1 / A.i + A.i)) %*% dnj.i)  %*% sqrtm(Lambda.0) %*% (B.i + t(B.i)) %*% sqrtm(Lambda.0) +
      .5 * sqrtm(Lambda.0) %*% (B.i + t(B.i)) %*% sqrtm(Lambda.0) %*% (t(dnj.i / (2 + 1 / A.i + A.i)) %*% Z.i[-n,]) +
      .5 *  (t( Z.i[-n,]/ (2 + 1 / A.i + A.i)) %*% dnj.i) %*% sqrtm(Lambda.0) %*% (B.i + t(B.i)) %*% sqrtm(Lambda.0) +
      (t(Z.i[-n,] / (2 + 1 / A.i + A.i)) %*% Z.i[-n,])
    
    
    
    
    
    Z.i.n=matrix(Z.i[n,] , byrow=TRUE, nrow=nrow(Z.a), ncol=nrow(Lambda.0))
    
    
    
    dnj.ia <- t(sqrtm(Lambda.1) %*% t(Z.i.n) + sqrtm(Lambda.0) %*% t(Z.a))
    
    
    A.ia <- 1/sqrt(det(B)) * exp(- alpha.1 ) * exp( - apply(as.matrix(Z.a), 1, function(x) x %*% Z.i[n,])
                                                    - .5 * apply(dnj.ia, 1, function(x) x %*% B %*% x))
    
    f1znTo.ia <-  colSums(Z.a / (1+A.ia)) + 
      .5 *  sqrtm(Lambda.1) %*% (B + t(B)) %*% colSums(dnj.ia / (1+A.ia))
    
    
    f2znTo.ia <- .5 *  sqrtm(Lambda.1) %*% (B + t(B)) %*% sqrtm(Lambda.1) * sum(1 / (1 + A.ia))+
      .25 *  sqrtm(Lambda.1)%*%  (B + t(B)) %*% (t(dnj.ia / (2 + 1 / A.ia + A.ia)) %*% dnj.ia)  %*%  (B + t(B)) %*% sqrtm(Lambda.1) +
      .5 * (t(Z.a / (2 + 1 / A.ia + A.ia)) %*% dnj.ia) %*%  (B + t(B)) %*% sqrtm(Lambda.1) +
      .5 * sqrtm(Lambda.1)  %*% t((B + t(B)) %*% (t(dnj.ia / (2 + 1 / A.ia + A.ia)) %*% Z.a)) +
      (t(Z.a / (2 + 1 / A.ia + A.ia)) %*% Z.a)
    
    
    
    
    numZnT <- .5 *colSums(as.matrix(Z.i[-n,]) * (Y.i[n,-n]+Y.i[-n,n])) - t(f1znTo.i)- .5 * t(f1znTo.ia) + Z.i[n, ] %*% (f2znTo.i+.5*f2znTo.ia)+ .5 * colSums(as.matrix(Z.a)* (Y.ia[n,]))
    denZnT <- ( 1 / (2 * p.lambda.0^2)) * diag(nrow(Lambda.0)) + f2znTo.i + f2znTo.ia *.5
    
    Z.i[n, ]<-  numZnT %*% solve(denZnT)
  }
  
  Z.i
  
  
}
