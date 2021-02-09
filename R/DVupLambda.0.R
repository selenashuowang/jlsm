
DVupLambda.0<-function(alpha.0, alpha.1, C,B,  Y.i, cont.i, Z.i,Z.a, Lambda.0, Lambda.1 ,p.lambda.0){
  
  D<-nrow(C) 
  N<-nrow(Y.i) 
  
  
  A <- 1 + 1 / sqrt(det(C)) * exp(- alpha.0) * exp(cont.i) 
  
  cbZ.i <- combn(nrow(Z.i), 2) 
  dij <- as.matrix(Z.i[cbZ.i[1, ], ] - Z.i[cbZ.i[2, ], ]) 
  
  f1s2To <- 8 * C %*% (t(dij / A) %*% dij) %*% C - 4 * sum(1 / A) * C 
  
  
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
  
  
  
  
  solve((1 / p.lambda.0^2 * N / 2 +  sum(Y.i) * 2 ) * diag(D) +  f1s2To  + f1s2To.ia) * N / 2
  
}
