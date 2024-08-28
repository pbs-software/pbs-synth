\name{findTarget}
\alias{findTarget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Find Target Probability
}
\description{
  For MCMC samples, find probability of some metric (e.g., \code{Bcurr}) occurring
  above target refence points (e.g., \code{Bmsy}, \code{B0}).
  For exploitation target, find probability of occurring below the target.
}
\usage{
findTarget(Vmat, yrU=as.numeric(dimnames(Vmat)[[2]]), yrG=90, ratio=0.5,
   target=B0.MCMC, conf=0.95, plotit=FALSE, retVal="N", op=">")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Vmat}{\code{matrix} -- projected biomass values \eqn{B_{Nt}}{B[N,t]}, \cr
    where \eqn{N} = number of MCMCs and \eqn{t} = projection year.}
  \item{yrU}{\code{numeric} -- user-specified projection years.}
  \item{yrG}{\code{numeric} -- number of years \eqn{G} for a moving target window (\emph{e.g.}, 
     3 YMR generations = 90y); might not work for all possibilities.}
  \item{ratio}{\code{numeric} -- recovery target ratio \eqn{R}.}
  \item{target}{\code{numeric} -- recovery target values \eqn{T_N}{T[N]} \cr
     = \code{B0.MCMC} for ratios of \eqn{B_0}{B0}; \cr
     = \code{Bmsy.MCMC} for ratios of \eqn{B_\text{MSY}}{Bmsy}; \cr
     = \code{Bt.MCMC} for moving window of \eqn{B_{N,t-G}}{B[N,t-G]}. }
  \item{conf}{\code{numeric} -- confidence level \eqn{C} required.}
  \item{plotit}{\code{logical} -- if \code{TRUE}, plot the probability \eqn{p_t}{p[t]} of exceeding target reference point.}
  \item{retVal}{\code{character} -- name of object to return: \cr
    \code{retVal="N"} : creates global object \code{"Ttab"} (see below); \cr
    \code{retVal="p.hi"} : creates global object \code{"Ptab"} (see below). }
  \item{op}{\code{character} -- binary operator for comparing values in atomic vectors.}
}
\details{
  As this function uses Bayesian output, there are \eqn{N} (\emph{e.g.}, 1,000)
  values of some target \eqn{T_N}{T[N}, which can remain fixed (\eqn{B_0}{B0}, \eqn{B_\text{MSY}}{Bmsy})
  or move forward in time \eqn{G} years before the projection year \eqn{t} 
  (that is \eqn{T_{N,t-G}}{T[N,t-G]}). For simplification, we'll just call all targets \eqn{T_N}{T[N]}.

  The probability of exceeding a target ratio \eqn{R} is:
  
  \deqn{p_t = \frac{1}{N} \sum^N  \left[ \frac{B_{Nt}}{T_N} > R \right], }{%
        p[t] = (1/N) \Sigma_N [ B[N,t]/T[N] > R ],} 
  
  where \eqn{R} = target ratio of the reference point 
  (\emph{e.g.}, 0.4\eqn{B_\text{MSY}}{Bmsy} (\eqn{R}=0.4), 0.2\eqn{B_0}{B0} (\eqn{R}=0.2), 0.5\eqn{B_{t-G}}{B[t-G]} (\eqn{R}=0.5) ).
  
  At a glance, we can see for any given projection year \eqn{t} whether the probability 
  of achieving a target ratio is greater than the confidence required:
  
  \deqn{p_t \ge C,}{p[t] \ge C,}
  
  where \eqn{C} is the confidence level acceptable.
}
\value{
  If \code{retVal="N"} then the function returns a data frame object 
  called \code{"Ttab"} in the user's global environment. 
  This table reports the number of years to achieve the target reference
  point at various catch levels with a specified confidence.

  If \code{retVal="p.hi"} then the function returns a list object
  called \code{"Ptab"} in the user's global environment. 
  This list contains data frames (tables) that report the probability
  of achieving various reference points at specified catch levels.

  Any other \code{retVal} will return a list of the specified object,
  if it exists in the function.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-17}
}
\note{
  This function originally appeared in package \pkg{PBSawatea}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{agileDT}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{logic}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
