#' @title Two dimensional plot of the Bipartite Latent Space Model
#'
#' @description plot the latent space with two types of nodes and one type of relations
#'
#' @param model model output from BLSM
#' @param Y.ia N by M matrix containing the binary multivariate attributes
#' @param labels vector of characters containing the item names
#' @param plotedges TRUE or FALSE, whether the bipartite edges should be plotted
#' @param edgecolor color of the edge. Default \code{edgecolor = "black"}
#' @param xlab name of the x axis
#' @param ylab name of the y axis
#' @param colEll.i  \code{col} for the ellipses of persons. Default \code{rgb(.6, .6 ,.6 , alpha=.1)}
#' @param colEll.ia  \code{col} for the ellipses of attributes Default \code{rgb(1, .6 ,.6 , alpha=.1)}
#' @param LEVEL levels of confidence bounds shown when plotting the ellipses. Default \code{LEVEL = .95}
#' @param pchplot Default \code{pchplot = 20}
#' @param pchEll \code{pch} for the ellipses. Default \code{pchEll = 19}
#' @param pchPl \code{pch} for the points representing the nodes. Default \code{pchPl = 19}
#' @param cexPl \code{cex} for the points representing the nodes. Default \code{cexPl = 1.1}
#' @param arrowhead logical, if the arrowed are to be plotted. Default \code{arrowhead = FALSE}
#' @param curve curvature of edges. Default \code{curve = 0}
#' @param lwdLine lwd of edges. Default \code{lwdLine = .3}
#' @param xlim range for x
#' @param ylim range for y
#' @param ... Arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).
#' @return plot
#' @export
#'
#' @examples
#' attach(french)
#' a=blsm(Niter=3,Y.ia,D=2)
#' Plotblsm(Y.ia,  a)

Plotblsm<-function(Y.ia,  model,  labels = NULL, xlab = "", ylab = "", plotedges = TRUE, edgecolor = "black",
                     colEll.i = rgb(.6, .6 ,.6 , alpha=.1), colEll.ia = rgb(1, .6 ,.6 , alpha=.1),
                     LEVEL = .80,
                     pchplot = 20, pchEll = 19, pchPl = 19, cexPl = 1.1,
                     arrowhead = FALSE, curve = 0, xlim = c(-2,2), ylim = c(-2,2), lwdLine = .001, ...){

  Z.i=model$lsmbEZ.i
  Z.a=model$lsmbEZ.a

  D=nrow(model$lsmbVZ.0)
  M=nrow(model$lsmbEZ.a)
  N=nrow(model$lsmbEZ.i)

  VZ1=model$lsmbVZ.1
  VZ0=model$lsmbVZ.0
  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1

  par(mar=c(2,2,2,2))

  plot(NA, xlim = xlim, ylim = ylim,
       ylab=ylab, xlab = xlab,...)

  for(n in 1:N)
  {
    coEl<-ellipse::ellipse(VZ0,centre = Z.i[n,],level=LEVEL,pch=pchEll)
    polygon(coEl[,1], coEl[,2], col = colEll.i, border=colEll.i)
  }

  if (plotedges == TRUE) {
    for(i in 1:N){
      for(a in 1:M){
        if(Y.ia[i,a] == 1)
          network.arrow(Z.i[i,1], Z.i[i,2], Z.a[a,1], Z.a[a,2], col = edgecolor,
                        border = rgb(0, 0, 0, alpha =.02), lwd = .005,
                        arrowhead = arrowhead, curve = curve)
      }
    }
  }



  for(n in 1:M)
  {
    coEl<-ellipse::ellipse(VZ1,centre = Z.a[n,],level=LEVEL,pch=pchEll)
    polygon(coEl[,1], coEl[,2], col = colEll.ia, border=colEll.ia)
  }

  if(!is.null(labels)){
    text(Z.a, labels = labels, col = "black", font = 2, cex=1)
  }



}


