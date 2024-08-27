\name{small.funs}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{allEqual}
\alias{cquantile.vec}
\alias{getNpan}
\alias{getYrIdx}
\alias{is.numStr}
\alias{med5.95}
\alias{ptab}
\alias{qtab}
\alias{stab}
\title{
  Small Helper Functions
}
\description{
  Small helper routines to facilitate transactions.
}
\usage{
allEqual(x)
cquantile.vec(z, prob)
getNpan()
getYrIdx(yrNames, mod=5)
is.numStr(x)
med5.95(xx.MCMC, dig=0, quants3=tcall(quants3))
ptab(xx)
qtab(xx.MCMC, dig=0, quants3=tcall(quants3))
stab(xx.MCMC, dig=3, quants3=tcall(quants3), print=TRUE)
}
\arguments{
  \item{x}{\code{numeric} -- vector of values (can sometimes be \code{character} also)}
  \item{z}{\code{numeric} -- cumulative quantile of vector}
  \item{prob}{\code{numeric} -- scalar probability}
  \item{yrNames}{\code{character} -- years as characters (e.g., \code{'1935'})}
  \item{mod}{\code{numeric} -- used in modulo operation: \code{yrVals \%\% mod==0}}
  \item{xx}{\code{numeric} -- vector of values}
  \item{xx.MCMC}{\code{numeric} -- vector of values}
  \item{dig}{\code{numeric} -- number of decimal places or significant digits, depending on function}
  \item{quant3}{\code{numeric} -- vector of three quantiles, usually (0.05, 0.5, 0.95)}
  \item{print}{\code{logical} -- if \code{TRUE}, print results to current device}
}
%- maybe also 'usage' for other objects documented here.
\details{
  \tabular{ll}{
    \code{allEqual}      \tab determine if all elements in vector equal the first element\cr
    \code{cquantile.vec} \tab calculate cumulative quantile as a vector\cr
    \code{getNpan}       \tab get panel number when inside a multi-panel plot\cr
    \code{getYrIdx}      \tab select years for plotting\cr
    \code{is.numStr}     \tab check if strings can be converted to numerics\cr
    \code{med5.95}       \tab print median and credibility envelope, e.g., (0.05, 0.95), to text\cr
    \code{ptab}          \tab prior tabulation (LaTeX) for lines in table of priors\cr
    \code{qtab}          \tab quantile tabulation summary using decimal places\cr
    \code{stab}          \tab quantile tabulation summary using significant digits
  }
}
\author{
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew Edwards}, Quantitative Assessment Methods, PBS, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-21}
}
\keyword{manip}
\concept{smallFuns.r}

