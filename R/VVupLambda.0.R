
VVupLambda.0.vector<-function(alpha.0, alpha.1, D,N, Y.i, Z.i,Z.a, Lambda.0, Lambda.1 ,p.lambda.0){
  

  B.i <- solve(diag(nrow(Lambda.1))- 2* Lambda.0 )
  B <- solve(diag(nrow(Lambda.1))-2*sqrtm(Lambda.0) %*% sqrtm(Lambda.1))
  
  cbZ.i <- combn(nrow(Z.i), 2) 
  dij <- as.matrix(Z.i[cbZ.i[1, ], ] + Z.i[cbZ.i[2, ], ]) 
  

  dnj.i <- t(sqrtm(Lambda.0) %*% t(Z.i[cbZ.i[1, ], ])+sqrtm(Lambda.0) %*% t(Z.i[cbZ.i[2, ], ]))
  
  
  temp.exp = sapply(1:nrow(dnj.i), function(i) - Z.i[cbZ.i[1, ], ][i,] %*% Z.i[cbZ.i[2, ],][i,] - .5 * dnj.i[i,] %*% B.i %*% dnj.i[i,])
  
  A.i <- 1+  1/sqrt(det(B.i)) * exp(-alpha.0) * exp(temp.exp) 
  
  
  
  f1s2To.i <- .25 * solve(sqrtm(Lambda.0)) %*% B.i %*% (t(dij / A.i) %*% dij) %*% sqrtm(Lambda.0) +
    .25 * t(solve(sqrtm(Lambda.0)) %*% B.i %*% (t(dij / A.i) %*% dij) %*% sqrtm(Lambda.0)) +
    sqrtm(Lambda.0) %*% B.i %*% (t(dij / A.i) %*% dij) %*% B.i %*% sqrtm(Lambda.0) +
    .5* sum(1 / A.i) * B.i 
  
  
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

  solve((1 / p.lambda.0^2 * N / 2  ) * diag(D) +  f1s2To.i*2  + f1s2To.ia) * N / 2
  
}

