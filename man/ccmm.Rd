\name{ccmm}
\alias{ccmm}
\title{
Causal Compositional Mediation Model
}
\description{
Estimate the direct and indirect (mediation) effects of treatment on the outcome when intermediate variables (mediators) are compositional and high-dimensional.
}
\usage{
ccmm(y, M, tr, x = NULL, w = NULL, method.est.cov = "bootstrap", n.boot = 2000,
	sig.level = 0.05, tol = 1e-06, max.iter = 5000)
}
\arguments{
  \item{y}{Vector of continuous outcomes}
  \item{M}{Matrix of compositional data}
  \item{tr}{Vector of continuous or binary treatments}
  \item{x}{Matrix of covariates}
  \item{w}{Vector of weights on samples}
  \item{method.est.cov}{One of two options ("bootstrap", "normal") to estimate the variance of indirect effects}
  \item{n.boot}{Number of bootstrap samples}
  \item{sig.level}{Significance level to estimate bootstrap confidence intervals for direct and indirect effects of treatment}
  \item{tol}{Error tolerance}
  \item{max.iter}{Maximum number of iteration in a debias procedure}
}
\value{
If method.est.cov is "bootstrap",
  \item{DE}{Direct effect of treatment on an outcome}
  \item{DE.CI}{Bootstrap confidence interval for the direct effect}
  \item{TIDE}{Total indirect effect of treatment on an outcome}
  \item{TIDE.CI}{Bootstrap confidence interval for the indirect effect}
  \item{IDEs}{Component-wise indirect effects of treatment on an outcome}
  \item{IDE.CIs}{Bootstrap confidence intervals for the component-wise indirect effects}
If method.est.cov is "normal",
  \item{DE}{Direct effect of treatment on an outcome}
  \item{Var.DE}{Variance of the direct effect}
  \item{TIDE}{Total indirect effect of treatment on an outcome}
  \item{Var.TIDE}{Variance of the indirect effect}
  \item{IDEs}{Component-wise indirect effects of treatment on an outcome}
  \item{Var.IDEs}{Variances of the component-wise indirect effects}
}
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

# Run CCMM
rslt.ccmm <- ccmm(outcome, mediators, treatment, covariates);
}