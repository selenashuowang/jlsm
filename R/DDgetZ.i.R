DDgetZ.i.t<-function(Y.i,Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,alpha.0, p.lambda.0,alpha.1){
 
  C <- solve(diag(nrow(Lambda.0)) + 4 * Lambda.0)
  B <- solve(diag(nrow(Lambda.1)) + 2* Lambda.0 + 2 * Lambda.1)
  
  for (n in 1:nrow(Z.i)){
    
    dnj <- t(Z.i[n,] - t(as.matrix(Z.i[-n,])))
    
    A <-  1/sqrt(det(C)) * exp(- alpha.0) * exp(apply(dnj, 1, function(x) x %*% C %*% x))
    
    f1znTo <- - 2 * C %*% colSums(dnj / (1+A))
    f2znTo <- - 2 * sum(1 / (1 + A)) * C + 
      4 * C %*% (t(dnj / (2 + 1 / A + A)) %*% dnj) %*% C
    
    dnj.ia <- t(Z.i[n,] - t(as.matrix(Z.a)))
    
    A.ia <-  1/sqrt(det(B)) * exp(- alpha.1) * exp(apply(dnj.ia, 1, function(x) x %*% B %*% x))
    
    f1znTo.ia <- - 2 * B %*% colSums(dnj.ia / (1+A.ia))
    f2znTo.ia <- - 2 * sum(1 / (1 + A.ia)) * B + 
      4 * B %*% (t(dnj.ia / (2 + 1 / A.ia + A.ia)) %*% dnj.ia) %*% B
    
    
    numZnT <- colSums(as.matrix(Z.i[-n,]) * (Y.i[n,-n]+Y.i[-n,n])) - t(f1znTo) + Z.i[n, ] %*% (f2znTo+.5*f2znTo.ia)+  colSums(as.matrix(Z.a)* (Y.ia[n,]))- .5 * t(f1znTo.ia)
    denZnT <- (sum(Y.i[n,-n] + Y.i[-n,n]) + 1 / (2 * p.lambda.0^2) +   sum(Y.ia[n,])) * diag(nrow(Lambda.0)) + f2znTo + f2znTo.ia *.5
    
    Z.i[n, ]<- numZnT %*% solve(denZnT)
  }
  
  Z.i
  
  
}
