\encoding{utf8} 
\name{linERR-package}
\alias{linERR-package}
\alias{linERR}
\docType{package}
\title{Fits the linear excess relative risk model
}
\description{
Usual approaches to the analysis of cohort and case control data often follow from risk-set sampling designs, where at each failure time a new risk set is defined, including the index case and all the controls that were at risk at that time. That kind of sampling designs are usually related to
the Cox proportional hazards model, available in most standard statistical packages but limited to log-linear models (except \emph{Epicure}, (Preston et al., 1993)) of the form \eqn{log(\phi(z, \beta)) = \beta_1 \cdot z_1 + \ldots \beta_k \cdot z_k}, where \eqn{z} is a vector of explanatory variables and \eqn{\phi} is the rate ratio. This implies exponential dose-response trends and multiplicative interactions, which may not be the best exposure-response representation in some cases, such as radiation exposures. One model of particular interest, especially in radiation environmental and occupational epidemiology is the ERR model,
 \eqn{\phi(z, \beta) = 1 + \alpha \cdot f(dose)}.
The ERR model represents the excess relative rate per unit of exposure and \eqn{z_1, \ldots, z_k} are covariates. Estimation of a dose-response trend under a linear relative rate model implies that for every 1-unit increase in the exposure metric, the rate of disease increases (or decreases) in an
additive fashion. The modification of the effect of exposure in linear relative rate models by a study covariate \eqn{m} can be assessed by including a log-linear subterm for the linear exposure effect (Preston et al., 2003; Ron et al., 1995), implying a model of the form \eqn{\phi(z, \beta) = e^{\beta_0 + \beta_1 \cdot z_1 + \ldots + \beta_k \cdot z_k} (1 + \alpha \cdot f(dose))}.
}

\details{
\tabular{ll}{
Package: \tab linERR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-02-23\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}

\author{
David Moriña, ISGlobal, Centre for Research in Environmental Epidemiology (CREAL)

Maintainer: David Moriña <david.morina@uab.cat>
}
\references{
B. Langholz and D. B. Richardson. Fitting general relative risk models for survival time and matched
case-control analysis. American journal of epidemiology, 171(3):377-383, 2010.
D. L. Preston, J. H. Lubin, D. A. Pierce, and M. E. McConney. Epicure: User's Guide. HiroSoft
International Corporation, Seattle, WA, 1993.
E. Ron, J. H. Lubin, R. E. Shore, K. Mabuchi, B. Modan, L. M. Pottern, A. B. Schneider, M. A.
Tucker, and J. D. Boice Jr. Thyroid Cancer after Exposure to External Radiation: A Pooled Analysis
of Seven Studies. Radiation Research, 141(3):259-277, 1995.
}
\keyword{ linERR }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
\seealso{
\code{\link{fit.linERR}}, \code{\link{ERRci}}
}