\name{ccmm.sensitivity}
\alias{ccmm.sensitivity}
\title{
Sensitivity analysis
}
\description{
Estimate the total indirect effect (TIDE) given a correlation coefficient
}
\usage{
ccmm.sensitivity(rh, y, M, tr, x = NULL, w = NULL)
}
\arguments{
  \item{rh}{Correlation coefficient}
  \item{y}{Vector of continuous outcomes}
  \item{M}{Matrix of compositional data}
  \item{tr}{Vector of continuous or binary treatments}
  \item{x}{Matrix of covariates}
  \item{w}{Vector of weights on samples}
}
\value{Estimated TIDE given a correlation coefficient}
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

ccmm.sensitivity(rh=0, outcome, mediators, treatment, covariates);
}