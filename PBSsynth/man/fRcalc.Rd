\name{fRcalc}
\alias{fRcalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Fraction Recruitment
}
\description{
  Calculate the fraction (proportion) of \eqn{R_0}{R0} that is allocated to
  each subarea in a multi-area SS3 model.
}
\usage{
fRcalc(d.mcmc, narea)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d.mcmc}{\code{numeric} -- data frame of MCMC parameter estimates}
  \item{narea}{\code{narea} -- number of subareas used in multi-area model}
}
\details{
  If recruitment distribution (\code{log Rdist}) is estimated,
  one of the subareas is held constant at zero while the other
  \code{(narea-1)} areas are estimated.
  The fractions are calculated using:
  \deqn{f_{\alpha} = e^{p_{\alpha}} / ( e^{p_{\alpha=1}} + e^{p_{\alpha=2}} + e^{p_{\alpha=3}} )}{%
  f_a = exp(p_a) / ( exp(p_1) + exp(p_2) + exp(p_3) )}
  where \eqn{\alpha \in{1,2,3}}{a in 1,2,3}.
  Annual time-varying recruitment proportions are calculated by adjusting the estimated proportions in log space:
  \deqn{p_{t,\alpha} = p_{\alpha} + \nu_{t,\alpha} \zeta_{\alpha}}{%
  p_(t,a) = p_a + nu_(t,a) * zeta_a}
  where \eqn{\nu_{t,\alpha}}{nu_(t,a)} = annual deviance, and
  \eqn{\zeta_{\alpha}}{zeta_a} = standard error of the deviance.
}
\value{
  Data frame of MCMC samples by row with fields (columns): 
  \itemize{
    \item \code{Rdist_a} -- (or \eqn{p_{\alpha}}{p_a}) the log-space estimates of recruitment distribution for subareas \code{a}
    \item \code{RdevSE_a} -- (or \eqn{\zeta_{\alpha}}{zeta_a}) the standard errors of deviance estimated/specified for subareas \code{a}
    \item \code{Rdev_(t,a)} -- (or \eqn{\nu_{t,\alpha}}{nu_(t,a)}) the deviance estimates for years \code{t} and subareas \code{a}
    \item \code{fR_(t,a)} -- (or \eqn{f_{t,\alpha}}{f_(t,a)}) the calculated fraction of coastwide \eqn{R_0}{R0} allocated to each year \code{t} and subarea \code{a}
  }
}
\references{
Methot, R.D., Wetzel, C.R., Taylor, I.G., Doering, K.L., Perl, E.F. and Johnson, K.F. (2024)
Stock Synthesis: User Manual, version 3.30.23. Tagged release, NOAA Fisheries, Seattle WA, USA,
October 31, 2024
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-07-31}
}
\note{
  Autocorrelation of recruitment distribution deviations can also be 
  estimated but this function does not incorporate it.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.pmcmc}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{sumting}
% \concept{ ~cpt2 }
% Use only one concept per line.
