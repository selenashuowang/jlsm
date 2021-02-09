DVgetZ.i.t<-function(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,alpha.0, p.lambda.0,alpha.1){
  
  
  
  C <- solve(diag(nrow(Lambda.0)) + 4 * Lambda.0)
  
  B <- solve(diag(nrow(Lambda.1))- 2* sqrtm(Lambda.0) %*% sqrtm(Lambda.1))
  
  for (n in 1:nrow(Z.i)){
    
    dnj <- t(Z.i[n,] - t(as.matrix(Z.i[-n,])))
    
    A <-  1/sqrt(det(C)) * exp(- alpha.0) * exp(apply(dnj, 1, function(x) x %*% C %*% x))
    
    f1znTo <- - 2 * C %*% colSums(dnj / (1+A))
    f2znTo <- - 2 * sum(1 / (1 + A)) * C + 
      4 * C %*% (t(dnj / (2 + 1 / A + A)) %*% dnj) %*% C
    

    
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
    
    
    
    numZnT <- colSums(as.matrix(Z.i[-n,]) * (Y.i[n,-n]+Y.i[-n,n])) - t(f1znTo)- .5 * t(f1znTo.ia) + Z.i[n, ] %*% (f2znTo+.5*f2znTo.ia)+ .5 * colSums(as.matrix(Z.a)* (Y.ia[n,]))
    denZnT <- (sum(Y.i[n,-n] + Y.i[-n,n]) + 1 / (2 * p.lambda.0^2)) * diag(nrow(Lambda.0)) + f2znTo + f2znTo.ia *.5
    
    Z.i[n, ]<- numZnT %*% solve(denZnT)
  }
  
  Z.i
  
  
}
