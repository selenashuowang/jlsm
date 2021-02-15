#' @title Simulate from the APLSM
#'
#' @description function to simulate networks from the APLSM
#'
#' @param type character indicating the types of model. It could be "DD", distance by distance model, "DV", distance by vector model,
#'   "VV", vector by vector model
#' @param model object of class APlsm
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{Y.i} (\code{N} x \code{N}) matrix containing the simulated Y.i
#'  \item \code{Y.ia} (\code{N} x \code{M}) matrix containing the simulated Y.ia
#'  }
#' @export
#'
#' @examples
#' attach(french)
#' b=aplsm(Niter=3,Y.i, Y.ia,D=2, type="DD")
#' Simulateaplsm(b,"DD")

Simulateaplsm<-function(model,type){

  est.alpha.0 = model$lsmhAlpha.0
  est.alpha.1 = model$lsmhAlpha.1
  Z.i = model$lsmhEZ.i
  Z.a = model$lsmhEZ.a

  D=nrow(model$lsmhVZ.1)
  N=nrow(Z.i)
  M=nrow(Z.a)

  Ps=Predictaplsm(model,type)

  Y.i <- matrix(rbinom(N * N, 1, Ps[[1]]), N, N)
  Y.ia<- matrix(rbinom(N * M, 1, Ps[[2]]), N, M)

  diag(Y.i)=0

  return(list(Y.i,Y.ia))



}

