\name{tide.ci.zero.rho}
\alias{tide.ci.zero.rho}
\title{
Bootstrap samples of TIDE with zero correlation
}
\description{
Generate bootstrap samples of the total indirect effect (TIDE) when the correlation coefficient is zero
}
\usage{
tide.ci.zero.rho(y, M, tr, x = NULL, w = NULL, n.boot=2000)
}
\arguments{
  \item{y}{Vector of continuous outcomes}
  \item{M}{Matrix of compositional data}
  \item{tr}{Vector of continuous or binary treatments}
  \item{x}{Matrix of covariates}
  \item{w}{Vector of weights on samples}
  \item{n.boot}{Number of bootstrap samples}
}
\value{bootstrap samples of TIDE}
\references{
Sohn, M.B. and Li, H. (2017). Compositional Mediation Analysis for Microbiome Studies (AOAS: In revision)
}
\author{
Michael B. Sohn

Maintainer: Michael B. Sohn <msohn@mail.med.upenn.edu>
}
\examples{
# Load test data
data(ccmm_test_data);
outcome <- ccmm_test_data[,1];
treatment <- ccmm_test_data[,2];
mediators <- as.matrix(ccmm_test_data[,3:22]);
covariates <- as.matrix(ccmm_test_data[,23:24]);

cisa <- tide.ci.zero.rho(outcome, mediators, treatment, covariates, n.boot=200)
}