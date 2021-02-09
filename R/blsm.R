#' @title The Bipartite Latent Space Model
#'
#' @description Function to fit the bipartite latent space model (BLSM) outlined in Wang et al. (2021)
#'
#' @param Niter number of iterations
#' @param Y.ia N by M matrix containing the binary multivariate attributes
#' @param D number of dimensions in the data
#'
#' @return list containing:
#'  \itemize{
#'  \item \code{lsmhEZ.i} (\code{N} x \code{D}) matrix containing the posterior means of the latent person positions
#'  \item \code{lsmhEZ.a} (\code{M} x \code{D}) matrix containing the posterior means of the latent item positions
#'  \item \code{lsmhVZ.0} (\code{D} x \code{D}) matrix containing the posterior variance of the latent person positions
#'  \item \code{lsmhVZ.1} (\code{D} x \code{D}) matrix containing the posterior variance of the latent item positions
#'  \item \code{lsmhAlpha.1} scaler of mean of the posterior distributions of \eqn{\alpha.1}
#'  \item \code{lsmhKL} expected log-likelihood
#'  }
#'
#' @export
#'
#' @examples
#' attach(french)
#' a=blsm(Niter=10,Y.ia,D=2)

blsm<-function(Niter,Y.ia,D){
  M=ncol(Y.ia)
  N=nrow(Y.ia)


  mxzmatKL.new=matrix(NA,nrow=Niter,ncol=1)
  mxzmatKL.new[1,]=10000000


  Lambda.0=diag(D)*.01
  Lambda.1=diag(D)*.01

  d <- dist(t(Y.ia))
  fit <- cmdscale(d,eig=TRUE, k=D)
  Z.a= fit$points

  d <- dist(Y.ia)
  fit <- cmdscale(d,eig=TRUE, k=D)
  Z.i= fit$points


  alpha.1=as.numeric(glm(c(Y.ia)~c(as.matrix(Z.i%*% t(Z.a))))$coeff[1])

  i=2
  while( i<Niter ){


    B=solve(diag(D) + 2*Lambda.0 + 2 * Lambda.1)

    alpha.1=Dupalpha.1(alpha.1, B, Y.ia, Z.i, Z.a)
    Lambda.0=DupLambda.0( alpha.1, B,  Y.ia, Z.i, Z.a, Lambda.0, Lambda.1 ,p.lambda.0=1)

    Lambda.0[upper.tri(Lambda.0)] <- t(Lambda.0)[upper.tri(t(Lambda.0))]

    if(!is.positive.definite(Lambda.0)){
      Lambda.0= nearPD(Lambda.0)$mat
    }


    Lambda.1=DupLambda.1( alpha.1, B, Y.ia, Z.i, Z.a, Lambda.0, Lambda.1 ,p.lambda.1=1)

    Lambda.1[upper.tri(Lambda.1)] <- t(Lambda.1)[upper.tri(t(Lambda.1))]

    if(!is.positive.definite(Lambda.1)){
      Lambda.1= nearPD(Lambda.1)$mat
    }


    Z.i=DgetZ.i(Y.ia, Lambda.0,Lambda.1,Z.i,Z.a, p.lambda.0=1,alpha.1)

    if(max(abs(range(Z.i)))>5){
      Z.i=DgetZ.i(Y.ia, Lambda.0,Lambda.1,Z.i,Z.a, p.lambda.0=1,alpha.1)

    }

    Z.a=DgetZ.a(Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,p.lambda.1=1,alpha.1)

    if(max(abs(range(Z.a)))>5){
      Z.a=DgetZ.a(Y.ia, Lambda.0,Lambda.1,Z.i,Z.a,p.lambda.1=1,alpha.1)

    }


    mxzmatKL.new[i,1]=DKL(  alpha.1, Lambda.0,Lambda.1, Z.i,Z.a, Y.ia)

    dis.new=sum(mxzmatKL.new[i,1])/sum(mxzmatKL.new[i-1,1])
    i=i+1
  }
  list(lsmbEZ.i=Z.i,lsmbEZ.a=Z.a,lsmbVZ.0=Lambda.0,lsmbVZ.1=Lambda.1,lsmbAlpha.1=alpha.1,
       lsmbKL=mxzmatKL.new[i-1,1])
}
