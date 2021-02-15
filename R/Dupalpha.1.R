
Dupalpha.1<-function(alpha.1, B, Y.ia, Z.i,Z.a){

  rbZ.a=Z.a

  for(each in 2:nrow(Z.i)){rbZ.a=rbind(Z.a,rbZ.a)}

  rbZ.i=Z.i[rep(1:nrow(Z.i), times = rep(nrow(Z.a),nrow(Z.i))), ]


  dnj.ia <- rbZ.a - rbZ.i

  A.ia <-  sqrt(det(B)) * exp(-apply(dnj.ia, 1, function(x) x %*% B %*% x))

  ( sum(Y.ia) -  sum(1 / (1 + exp(-alpha.1) / A.ia)) + alpha.1 *  sum(1 / (2 + exp(-alpha.1) / A.ia+ exp(alpha.1) * A.ia))) /
    (   sum(1 / (2 + exp(-alpha.1) / A.ia + exp(alpha.1) * A.ia)))


}


