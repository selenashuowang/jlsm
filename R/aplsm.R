#' @title The Attribute Person Latent Space model
#'
#' @description Jointly model social network with multivariate attributes
#'
#' @param Niter number of iterations
#' @param Y.i N by N matrix containing the binary social network
#' @param Y.ia N by M matrix containing the binary multivariate attributes
#' @param D number of dimensions in the data
#' @param type character indicating the types of model. It could be "DD", distance by distance model, "DV", distance by vector model,
#' "VV", vector by vector model
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{lsmhEZ.i} (\code{N} x \code{D}) matrix containing the posterior means of the latent person positions
#'  \item \code{lsmhEZ.a} (\code{M} x \code{D}) matrix containing the posterior means of the latent item positions
#'  \item \code{lsmhVZ.0} (\code{D} x \code{D}) matrix containing the posterior variance of the latent person positions
#'  \item \code{lsmhVZ.1} (\code{D} x \code{D}) matrix containing the posterior variance of the latent item positions
#'  \item \code{lsmhAlpha.0} scaler of mean of the posterior distributions of \eqn{\alpha.0}
#'  \item \code{lsmhAlpha.1} scaler of mean of the posterior distributions of \eqn{\alpha.1}
#'  \item \code{lsmhKL} expected log-likelihood
#'  }
#'
#' @export
#' @examples
#' attach(french)
#' a=aplsm(Niter=10,Y.i, Y.ia, D=2, type="DD")


aplsm=function (Niter, Y.i, Y.ia, D, type) 
{
  M = ncol(Y.ia)
  N = ncol(Y.i)

  Z.a = cmdscale(dist(t(Y.ia)), eig = TRUE, k = D)$points
  

  Z.i = cmdscale(dist(Y.ia), eig = TRUE, k = D)$points
  
  
  ZsMat = list(Z.i = Z.i, Z.a = Z.a)
  if (type == "DD") {
    lsmhMat <- DDlsmh(Niter, Y.i, Y.ia, M, N, D, ZsMat)
  }
  if (type == "DV") {
    lsmhMat <- DVlsmh(Niter, Y.i, Y.ia, M, N, D, ZsMat)
  }
  if (type == "VV") {
    lsmhMat <- VVlsmh(Niter, Y.i, Y.ia, M, N, D, ZsMat)
  }
  return(lsmhMat)
}

