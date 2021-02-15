#' @title Predict from BLSM model
#'
#' @description This function allows you to obtain the posterior mean of the edges from the BLSM model
#'
#'
#' @param model object of class BLSM
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{est.P.ia} (\code{N} x \code{M}) matrix containing the predicted probabilities of an edge
#'  }
#'
#' @export
#'
#' @examples
#' attach(french)
#' a=blsm(Niter=5,Y.ia,D=2)
#' Predictblsm(a)


Predictblsm<-function(model){

  est.alpha.1 = model$lsmbAlpha.1
  Z.i = model$lsmbEZ.i
  Z.a = model$lsmbEZ.a


  D=nrow(model$lsmbVZ.1)
  N=nrow(Z.i)
  M=nrow(Z.a)

  est.P.ia=matrix(NA,N,M)




  for(i in 1:N){
    for(a in 1:M){
      est.P.ia[i,a]=inv.logit(est.alpha.1-sum((Z.i[i,] - Z.a[a,]) ^ 2))
    }
  }

  return(est.P.ia)



}
