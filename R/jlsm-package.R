#' @title Create Joint Latent Space Model for Social networks and Multivariate Attributes
#'
#' @description \code{jlsm} provides a set of latent space models for jointly modeling
#' unipartite social networks with bipartite attribute networks. The latent space models are implemented using the
#' variational inference approach.
#'
#' @details Latent space models for bipartite networks: the function \code{\link{blsm}} implements the bipartite latent space model (BLSM) outlined in  Wang et al. (2021) using variational inference and squared Euclidian distance; the function
#' \code{\link{aplsm}} implements person and attribute latent space model (APLSM) introduced by
#' Wang et.al (2021).
#' These models assume that the person and attribute information can be summarized by latent person and attribute variables.
#' Both the Euclidean distances and the vector distances are used to describe relationships among persons and between persons and attributes.
#'
#' @references Wang, S. S., Paul, S., Logan, J., & De Boeck, P. (2019). Joint analysis of social and item response networks with latent space models. arXiv preprint arXiv:1910.12128.
#' @name jlsm-package
#' @aliases jlsm
#' @import MASS
#' @importFrom lvm4net lsm
#' @importFrom stats as.dist cmdscale dist glm rnorm rbinom
#' @importFrom expm sqrtm
#' @importFrom utils combn
#' @importFrom graphics abline boxplot legend lines matlines matplot matpoints mtext par plot points polygon text
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvnorm
#' @importFrom matrixcalc is.positive.definite
#' @importFrom boot inv.logit
#' @importFrom pROC roc
#' @importFrom grDevices rgb
#' @importFrom network network.arrow
#' @importFrom Matrix nearPD
#' @importFrom igraph layout.fruchterman.reingold graph.adjacency graph_from_incidence_matrix
#' @docType package
NULL
