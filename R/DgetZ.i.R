DgetZ.i<-function(Y.ia, Lambda.0,Lambda.1,Z.i,Z.a, p.lambda.0,alpha.1){


  B <- solve(diag(nrow(Lambda.1)) + 2* Lambda.0 + 2 * Lambda.1)

  for (n in 1:nrow(Z.i)){


    dnj.ia <- t(Z.i[n,] - t(as.matrix(Z.a)))

    A.ia <-  1/sqrt(det(B)) * exp(- alpha.1) * exp(apply(dnj.ia, 1, function(x) x %*% B %*% x))

    f1znTo.ia <- - 2 * B %*% colSums(dnj.ia / (1+A.ia))
    f2znTo.ia <- - 2 * sum(1 / (1 + A.ia)) * B +
      4 * B %*% (t(dnj.ia / (2 + 1 / A.ia + A.ia)) %*% dnj.ia) %*% B


    numZnT <-   .5* Z.i[n, ] %*% f2znTo.ia +  colSums(as.matrix(Z.a)* (Y.ia[n,]))- .5 * t(f1znTo.ia)
    denZnT <- ( 1 / (2 * p.lambda.0^2) +   sum(Y.ia[n,])) * diag(nrow(Lambda.0)) +  f2znTo.ia *.5

    Z.i[n, ]<- numZnT %*% solve(denZnT)
  }

  Z.i


}
