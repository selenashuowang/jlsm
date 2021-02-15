#' @title Predict from the APLSM
#'
#' @description This function allows you to obtain the posterior edge values based on the APLSM
#'
#' @param type character indicating the types of model. It could be "DD", distance by distance model, "DV", distance by vector model,
#'  "VV", vector by vector model
#' @param model object of class the APLSM
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{est.P.i} (\code{N} x \code{N}) matrix containing the predicted probabilities of an edge
#'  \item \code{est.P.ia} (\code{N} x \code{M}) matrix containing the predicted probabilities of an edge
#'  }
#'
#' @export
#'
#' @examples
#' attach(french)
#' b=aplsm(Niter=3,Y.i, Y.ia,D=2, type="DD")
#' Predictaplsm(b,"DD")



Predictaplsm<-function(model,type){

  est.alpha.0 = model$lsmhAlpha.0
  est.alpha.1 = model$lsmhAlpha.1
  Z.i = model$lsmhEZ.i
  Z.a = model$lsmhEZ.a


  D=nrow(model$lsmhVZ.1)
  N=nrow(Z.i)
  M=nrow(Z.a)


  if(type == "DD"){

    Ps= DDgetEstP(est.alpha.0,est.alpha.1,Z.i,Z.a,D,M,N)

  }

  if(type == "DV"){
    Ps= DVgetEstP(est.alpha.0,est.alpha.1,Z.i,Z.a,D,M,N)


  }
  if(type == "VV"){
    Ps= VVgetEstP(est.alpha.0,est.alpha.1,Z.i,Z.a,D,M,N)

  }

  return(Ps)



}

