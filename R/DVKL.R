DVKL<-function(alpha.0,  alpha.1, Lambda.0,Lambda.1, Z.i,Z.a, Y.i,Y.ia){

  B <- solve(diag(nrow(Lambda.1))- 2* sqrtm(Lambda.0) %*% sqrtm(Lambda.1))
  
  
  C <- solve(diag(nrow(Lambda.0)) + 4 * Lambda.0)
  
  A <- log(1 + sqrt(det(C)) * exp(alpha.0) * exp(- dist(Z.i %*% chol(C))^2))
  
  
  rbZ.a=Z.a
  
  for(each in 2:nrow(Z.i)){rbZ.a=rbind(Z.a,rbZ.a)}
  
  rbZ.i=Z.i[rep(1:nrow(Z.i), times = rep(nrow(Z.a),nrow(Z.i))), ]
  dnj.ia <- t(sqrtm(Lambda.1)%*%t(rbZ.i)+sqrtm(Lambda.0)%*%t(rbZ.a))
  
  exp.ia = sapply(1:nrow(dnj.ia), function(i) 0.5 * dnj.ia[i,] %*% B %*% dnj.ia[i,] + rbZ.i[i,] %*% rbZ.a[i,])
  
  A.ia <- log(1 + sqrt(det(B)) * exp(alpha.1) * exp(exp.ia))
  
  
  
  KL=-(sum((alpha.0 - 2 * sum(diag(Lambda.0)) - as.matrix(dist(Z.i)^2)) * Y.i) - 2 * sum(A) +
         sum( Y.ia * (as.matrix(Z.i%*%t(Z.a)) + alpha.1) )  - sum(A.ia) )
  
  ifelse(KL==Inf,100000,KL)
  
}