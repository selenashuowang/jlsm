
DDupLambda.0<-function(alpha.0, alpha.1, C,B,  Y.i, Y.ia, cont.i, Z.i,Z.a, Lambda.0, Lambda.1 ,p.lambda.0){

  D<-nrow(C)
  N<-nrow(Y.i)


  A <- 1 + 1 / sqrt(det(C)) * exp(- alpha.0) * exp(cont.i)

  cbZ.i <- combn(nrow(Z.i), 2)
  dij <- as.matrix(Z.i[cbZ.i[1, ], ] - Z.i[cbZ.i[2, ], ])

  f1s2To <- 8 * C %*% (t(dij / A) %*% dij) %*% C - 4 * sum(1 / A) * C


  rbZ.a=Z.a

  for(each in 2:nrow(Z.i)){rbZ.a=rbind(Z.a,rbZ.a)}

  rbZ.i=Z.i[rep(1:nrow(Z.i), times = rep(nrow(Z.a),nrow(Z.i))), ]



  dnj.ia <- rbZ.a - rbZ.i

  A.ia <-  1 + 1 / sqrt(det(B)) * exp(- alpha.1) * exp(apply(dnj.ia, 1, function(x) x %*% B %*% x))


  f1s2To.ia <- 2 * B %*% (t( dnj.ia  / A.ia) %*% dnj.ia) %*% B -   sum(1 / A.ia) * B


  solve((1 / p.lambda.0^2 * N / 2 +  sum(Y.i) * 2 + sum(Y.ia)) * diag(D) +  f1s2To  + f1s2To.ia) * N / 2



}
