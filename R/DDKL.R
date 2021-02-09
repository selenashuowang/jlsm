DDKL<-function(alpha.0,  alpha.1, Lambda.0,Lambda.1, Z.i,Z.a, Y.i,Y.ia){

  M=ncol(Y.ia)
  N=nrow(Y.ia)


  B <- solve(diag(nrow(Lambda.1))+ 2*Lambda.0 + 2 * Lambda.1)


  C <- solve(diag(nrow(Lambda.0)) + 4 * Lambda.0)

  A <- log(1 + sqrt(det(C)) * exp(alpha.0) * exp(- dist(Z.i %*% chol(C))^2))

  rbZ.a=Z.a

  for(each in 2:nrow(Z.i)){rbZ.a=rbind(Z.a,rbZ.a)}

  rbZ.i=Z.i[rep(1:nrow(Z.i), times = rep(nrow(Z.a),nrow(Z.i))), ]



  dnj.ia <- rbZ.a - rbZ.i

  A.ia <-  log(1 +  sqrt(det(B)) * exp( alpha.1) * exp(- apply(dnj.ia, 1, function(x) x %*% B %*% x)))




  KL=-(sum((alpha.0 - 2 * sum(diag(Lambda.0)) - as.matrix(dist(Z.i)^2)) * Y.i) - 2 * sum(A) +
         sum(Y.ia * (alpha.1 -  sum(diag(Lambda.0)) - sum(diag(Lambda.1))
                     - matrix(as.numeric(apply(dnj.ia, 1, function(x) x %*% x)),byrow = TRUE, nrow=N, ncol=M)))
       - sum(A.ia) )

  ifelse(KL==Inf,100000000000000,KL)


}
