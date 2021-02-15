VVKL<-function(alpha.0,  alpha.1, Lambda.0,Lambda.1, Z.i,Z.a, Y.i,Y.ia){

  B <- solve(diag(nrow(Lambda.1))-2*sqrtm(Lambda.0)%*%sqrtm(Lambda.1))
  B.i <- solve(diag(nrow(Lambda.1))- 2* Lambda.0 )

  rbZ.second=Z.i[-1,]
  
  for(n in 2:nrow(Z.i)){rbZ.second=rbind(rbZ.second,Z.i[-n,])}
  
  rbZ.first=Z.i[rep(1:nrow(Z.i), times = rep((nrow(Z.i)-1),nrow(Z.i))), ]

  dnj.i <- t(sqrtm(Lambda.0) %*% t(rbZ.first) + sqrtm(Lambda.0) %*% t(rbZ.second))
  
  
  exp.i = sapply(1:nrow(dnj.i), function(i) 0.5 * dnj.i[i,] %*% B.i %*% dnj.i[i,] + rbZ.first[i,] %*% rbZ.second[i,])
  
  A.i <- log(1 + sqrt(det(B.i)) * exp(alpha.0) * exp(exp.i))

  rbZ.a=Z.a
  
  for(each in 2:nrow(Z.i)){rbZ.a=rbind(Z.a,rbZ.a)}
  
  rbZ.i=Z.i[rep(1:nrow(Z.i), times = rep(nrow(Z.a),nrow(Z.i))), ]

  dnj.ia <- t(sqrtm(Lambda.1)%*%t(rbZ.i)+sqrtm(Lambda.0)%*%t(rbZ.a))
  
  exp.ia = sapply(1:nrow(dnj.ia), function(i) 0.5 * dnj.ia[i,] %*% B %*% dnj.ia[i,] + rbZ.i[i,] %*% rbZ.a[i,])
  
  A.ia <- log(1 + sqrt(det(B)) * exp(alpha.1) * exp(exp.ia))

  KL=-(sum(Y.ia * (as.matrix(Z.i%*%t(Z.a)) + alpha.1)) - sum(A.ia) +
         sum(Y.i * (as.matrix(Z.i%*%t(Z.i)) + alpha.0)) - sum( A.i))
  
  ifelse(KL==Inf,100000000000000,KL)
  

}
