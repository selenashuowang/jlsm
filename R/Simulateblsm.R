#' @title Simulate from the BLSM model
#'
#' @description function to simulate networks from the BLSM
#'
#' @param model object of class BLSM
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{Y.ia} (\code{N} x \code{M}) matrix containing the simulated Y.ia
#'  }
#'
#' @export
#'
#' @examples
#' attach(french)
#' a=blsm(Niter=5,Y.ia,D=2)
#' Simulateblsm(a)


Simulateblsm<-function(model){

  est.alpha.1 = model$lsmbAlpha.1
  Z.i = model$lsmbEZ.i
  Z.a = model$lsmbEZ.a

  D=nrow(model$lsmbVZ.1)
  N=nrow(Z.i)
  M=nrow(Z.a)


  est.P.ia=Predictblsm(model)



  Y.ia<- matrix(rbinom(N * M, 1, est.P.ia), N, M)

  return(Y.ia)



}

