\encoding{utf8} 
\name{ERRci}
\alias{ERRci}

\title{Profile likelihood based confidence intervals
}
\description{
The standard procedure for computing a confidence interval for a parameter \eqn{\beta} (Wald-type CI), based on \eqn{\hat{\beta} \pm z_{1-\frac{\alpha}{2}} SE(\hat{\beta})} may work poorly if the distribution of the parameter estimator is markedly skewed or if the standard error is a poor
estimate of the standard deviation of the estimator. Profile likelihood confidence intervals doesn't assume normality of
the estimator and perform better for small sample sizes or skewed estimates than Wald-type confidence intervals.
}
\usage{
ERRci(object, prob=0.95)
}

\arguments{
  \item{object}{
An object of class \code{fit.linERR}.
}
  \item{prob}{
Level of confidence, defaults to 0.95.
}
}
\value{
A numeric vector containing the \eqn{prob} profile likelihood based confidence interval.
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
  ERRci(fit.1, prob=0.9)
}


\keyword{ linERR }
\keyword{ models }
\keyword{ regression }
\keyword{ survival }
