
VVupalpha.1<-function(alpha.1, B, Y.ia, dnj.ia,Z.i,Z.a){ 
  
  
  A.ia <-  sqrt(det(B)) *exp( Z.i%*%t(Z.a)+ .5*matrix(
    apply(dnj.ia, 1, function(x) t(x %*% t(B) %*% x))
    ,nrow = nrow(Z.i),byrow=TRUE))
  
  
  (sum(Y.ia) -  sum((1 / (1 + exp(-alpha.1) / A.ia))) + alpha.1 * sum(1 / (2 + exp(-alpha.1) / A.ia+ exp(alpha.1) * A.ia))) /
    (  sum((1 / (2 + exp(-alpha.1) / A.ia + exp(alpha.1) * A.ia))))
  
  
}