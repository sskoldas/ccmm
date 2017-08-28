\name{ccmm-package}
\alias{ccmm-package}
\docType{package}
\title{
Causal Compositional Mediation Model
}
\description{
Estimate the direct and indirect (mediation) effects of treatment on the outcome when intermediate variables (mediators) are compositional and high-dimensional.
}
\author{
Michael B. Sohn

Maintainer: Michael B. Sohn <msohn@mail.med.upenn.edu>}
\references{
Sohn, M.B. and Li, H. (2017). Compositional Mediation Analysis for Microbiome Studies. (AOAS: In revision)
}
\examples{
\dontrun{
# Load test data
data(ccmm_test_data);
head(ccmm_test_data);
outcome <- ccmm_test_data[,1];
treatment <- ccmm_test_data[,2];
mediators <- as.matrix(ccmm_test_data[,3:22]);
covariates <- as.matrix(ccmm_test_data[,23:24]);

# Run CCMM
rslt.ccmm <- ccmm(outcome, mediators, treatment, covariates);

# Sensitivity analysis
rslt.sa <- ccmm.sa(outcome, mediators, treatment, covariates);
plot(rslt.sa, type="l", xlab=expression(rho), ylab="TIDE")
abline(h=rslt.ccmm$TIDE, lty=2)
abline(h=0, lty=3)
cisa <- tide.ci.zero.rho(outcome, mediators, treatment, covariates)
csqs <- quantile(cisa, c(0.025, 0.975))
segments(0, csqs[1], 0, csqs[2])
}
}