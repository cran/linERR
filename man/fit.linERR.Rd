\encoding{utf8} 
\name{fit.linERR}
\alias{fit.linERR}
\title{Fits linear ERR model
}
\description{
Usual approaches to the analysis of cohort and case control data often follow from risk-set sampling designs, where at each failure time a new risk set is defined, including the index case and all the controls that were at risk at that time. That kind of sampling designs are usually related to
the Cox proportional hazards model, available in most standard statistical packages but limited to log-linear models (except \emph{Epicure}, (Preston et al., 1993)) of the form \eqn{log(\phi(z, \beta)) = \beta_1 \cdot z_1 + \ldots \beta_k \cdot z_k}, where \eqn{z} is a vector of explanatory variables and \eqn{\phi} is the rate ratio. This implies exponential dose-response trends and multiplicative interactions, which may not be the best exposure-response representation in some cases, such as radiation exposures. One model of particular interest, especially in radiation environmental and occupational epidemiology is the ERR model,
 \eqn{\phi(z, \beta) = 1 + \alpha \cdot f(dose)}.
The ERR model represents the excess relative rate per unit of exposure and \eqn{z_1, \ldots, z_k} are covariates. Estimation of a dose-response trend under a linear relative rate model implies that for every 1-unit increase in the exposure metric, the rate of disease increases (or decreases) in an
additive fashion. The modification of the effect of exposure in linear relative rate models by a study covariate \eqn{m} can be assessed by including a log-linear subterm for the linear exposure effect (Preston et al., 2003; Ron et al., 1995), implying a model of the form \eqn{\phi(z, \beta) = e^{\beta_0 + \beta_1 \cdot z_1 + \ldots + \beta_k \cdot z_k} (1 + \alpha \cdot f(dose))}.
}
\usage{
fit.linERR(formula, beta = NULL, data, ages, lag = 0)
}
\arguments{
  \item{formula}{
An object of class \code{formula} (or one that can be coerced to that class), i.e. a symbolic description of the model to be fitted. The response must be a survival object as returned by the \code{Surv()} function, and the log-linear and linear terms are separated by the character ``|''. Stratum are defined using the \code{strata()} function.
}
  \item{beta}{
Starting values for parameter estimates. Its default value is \code{NULL}.
}
  \item{data}{
Data frame that contains the cohort.
}
  \item{ages}{
Age at each exposure.
}
  \item{lag}{
Lag to be applied. Its default value is zero.
}

}

\value{
An object of class \code{fit.linERR}, essentially a named list. The elements of this list are detailed below
\item{lowb }{Low boundary of the parameter in the linear part.}
\item{beta }{Initial values for the estimates.}
\item{max.exp }{Maximum number of exposures.}
\item{covariates1 }{Covariates in the loglinear part.}
\item{data_2 }{Original data reestructured as a list.}
\item{rsets_2 }{Risk sets reestructured as a list.}
\item{doses_2 }{Doses at each exposure reestructured as a list.}
\item{ages_2 }{Ages at each exposure reestructured as a list.}
\item{vcov }{Variance-covariance matrix.}
\item{aic }{Akaike's Information Criteria.}
\item{Call}{Call to the function.}
\item{llike}{Maximum log-likelihood.}
\item{deviance}{Deviance of the model.}
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
\author{
David Moriña, ISGlobal, Centre for Research in Environmental Epidemiology (CREAL)
}

\seealso{
\code{\link{ERRci}}, \code{\link{linERR-package}}
}
\examples{
  data(cohort1) 
  fit.1 <- fit.linERR(Surv(entryage, exitage, leu)~sex|dose1+dose2+dose3+dose4+dose5+dose6+
                      dose7+dose8+dose9+dose10+dose11+dose12+dose13+dose14+dose15+dose16+
                      dose17+dose18+dose19+dose20+dose21+dose22+dose23+dose24+dose25+dose26+
                      dose27+dose28+dose29+dose30+dose31+dose32, data=cohort1, beta=NULL, 
                      ages=cohort1[, 7:38], lag=2)
}

\keyword{ linERR }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
