\name{ccmm.sa}
\alias{ccmm.sa}
\title{
Sensitivity analysis
}
\description{
Estimated total indirect effects (TIDE) given correlation coefficients (rho)
}
\usage{
ccmm.sa(y, M, tr, x = NULL, w = NULL, stp = 0.01)
}
\arguments{
  \item{y}{Vector of continuous outcomes}
  \item{M}{Matrix of compositional data}
  \item{tr}{Vector of continuous or binary treatments}
  \item{x}{Matrix of covariates}
  \item{w}{Vector of weights on samples}
  \item{stp}{Increment of the correlation coefficient}
}
\value{Matrix of rho and TIDE}
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

rslt.sa <- ccmm.sa(outcome, mediators, treatment, covariates);
}