

DKL<-function( alpha.1, Lambda.0,Lambda.1, Z.i,Z.a, Y.ia){

 N=nrow(Y.ia)
 M=ncol(Y.ia)


  B <- solve(diag(nrow(Lambda.1))+ 2*Lambda.0 + 2 * Lambda.1)


  rbZ.a=Z.a

  for(each in 2:nrow(Z.i)){rbZ.a=rbind(Z.a,rbZ.a)}

  rbZ.i=Z.i[rep(1:nrow(Z.i), times = rep(nrow(Z.a),nrow(Z.i))), ]



  dnj.ia <- rbZ.a - rbZ.i

  A.ia <-  log(1 +  sqrt(det(B)) * exp( alpha.1) * exp(- apply(dnj.ia, 1, function(x) x %*% B %*% x)))


  KL=-( sum(Y.ia * (alpha.1 -  sum(diag(Lambda.0)) - sum(diag(Lambda.1))
                    - matrix(as.numeric(apply(dnj.ia, 1, function(x) x %*% x)),byrow = TRUE, nrow=N, ncol=M)))
        - sum(A.ia) )

  ifelse(KL==Inf,100000000000000,KL)

}
