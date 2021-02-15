DDupalpha.0<-function(alpha.0, C, Y.i, cont.i){ 
  
  A <- sqrt(det(C)) * exp(-cont.i) 
  
  ( sum(Y.i) - 2 * sum(1 / (1 + exp(-alpha.0) / A)) + alpha.0 * 2 * sum(1 / (2 + exp(-alpha.0) / A+ exp(alpha.0) * A))) / (  2 * sum(1 / (2 + exp(-alpha.0) / A + exp(alpha.0) * A))) 
  
} 
