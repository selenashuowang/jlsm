% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gofaplsm.R
\name{GOFaplsm}
\alias{GOFaplsm}
\title{Assess the fit of the APLSM}
\usage{
GOFaplsm(model, type, Y.i, Y.ia)
}
\arguments{
\item{model}{object of class the APLSM}

\item{type}{character indicating the types of model. It could be "DD", distance by distance model, "DV", distance by vector model,
"VV", vector by vector model}

\item{Y.i}{N by N matrix containing the binary social network}

\item{Y.ia}{N by M matrix containing the binary multivariate attributes}
}
\value{
list containing:
 \itemize{
 \item \code{Yi.auc} scaler of the area under the curve for the social network
 \item \code{Ya.auc} scaler of the area under the curve for the multivariate covariates
 }
}
\description{
assess the fit of the model using ROC curves and auc values
}
\examples{
attach(french)
b=aplsm(Niter=3,Y.i, Y.ia,D=2, type="DD")
GOFaplsm(b, "DD",Y.i, Y.ia)
}
