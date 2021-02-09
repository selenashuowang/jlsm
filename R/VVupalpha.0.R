VVupalpha.0.vector<-function(alpha.0,  Y.i, Z.i, D,Lambda.0){ 
  
  B.i <- solve(diag(D)- 2* Lambda.0 )
  
  rbZ.i=Z.i[-1,]
  
  for(n in 2:nrow(Z.i)){rbZ.i=rbind(rbZ.i,Z.i[-n,])}
  
  rbZ.i.1=Z.i[rep(1:nrow(Z.i), times = rep((nrow(Z.i)-1),nrow(Z.i))), ]
  

  dnj.i <- t(sqrtm(Lambda.0) %*% t(rbZ.i.1)+sqrtm(Lambda.0) %*% t(rbZ.i))
  

  t.i=matrix(sapply(1:nrow(rbZ.i.1), function(i) t(rbZ.i.1[i,]) %*% rbZ.i[i,]), 
             nrow = nrow(Z.i), byrow = TRUE)

  temp.exp.i= - t.i - .5 * matrix(
    apply(dnj.i, 1, function(x) t(x %*% t(B.i) %*% x)) 
    , nrow = nrow(Z.i), byrow=TRUE)

  A.i <-  sqrt( det(B.i)) * exp(-temp.exp.i )
  
  
  
  (sum(Y.i) -  sum((1 / (1 + exp(-alpha.0) / A.i))) + alpha.0 *  sum(1 / (2 + exp(-alpha.0) / A.i+ exp(alpha.0) * A.i))) /
    (  sum((1 / (2 + exp(-alpha.0) / A.i + exp(alpha.0) * A.i))))
  
} 

