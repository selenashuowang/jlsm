
DVupLambda.1.ia<-function(alpha.1, B, D, M, Z.i,Z.a, Lambda.0, Lambda.1 ,p.lambda.1){
  
  
  rbZ.a=Z.a
  
  for(each in 2:nrow(Z.i)){rbZ.a=rbind(Z.a,rbZ.a)}
  
  rbZ.i=Z.i[rep(1:nrow(Z.i), times = rep(nrow(Z.a),nrow(Z.i))), ]
  
  dnj.ia <- t(sqrtm(Lambda.1)%*%t(rbZ.i)+sqrtm(Lambda.0)%*%t(rbZ.a))
  
  
  
  temp.exp = sapply(1:nrow(rbZ.i), function(i) - rbZ.i[i,] %*% rbZ.a[i,] - .5 * dnj.ia[i,] %*% B %*% dnj.ia[i,])
  
  
  
  A.ia <- 1+  1/sqrt(det(B)) * exp(-alpha.1) * exp(temp.exp) 

  f1s2To.ia <- .25 * solve(sqrtm(Lambda.0)) %*% (t(dnj.ia / A.ia) %*% rbZ.a) %*% B +
    .25 * t(solve(sqrtm(Lambda.0)) %*% (t(dnj.ia / A.ia) %*% rbZ.a) %*% B) +
    .5* solve(sqrtm(Lambda.0)) %*% B %*% (t(dnj.ia / A.ia) %*% dnj.ia) %*% B %*% sqrtm(Lambda.1) +
    .5* sum(1 / A.ia) *sqrtm(Lambda.1) %*% B %*% solve(sqrtm(Lambda.0))
  
  
  
  
  solve((1 / p.lambda.1^2 * M / 2) * diag(D)  + f1s2To.ia) * M / 2
  
  
}
