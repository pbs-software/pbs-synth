\name{compTS}
\alias{compTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compare Time Series
}
\description{
  Compare MCMC medians of sensitivity run time series (trajectories).
}
\usage{
compTS(x, runs=c(24,39), val="Rtdev", type="bars")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{list} -- a \code{'senso'} list object created by the function \code{'gatherMCMC'}}
  \item{runs}{\code{numeric} -- sensitivity run numbers to compare}
  \item{val}{\code{character} -- derived quantity (parameter) for comparison}
  \item{type}{\code{character} -- type of plot (currently only \code{'bars'} is supported)}
}
\details{
  A very basic barplot comparison of two sensitivities.
  Function retained for future development, should the need arise.
}
\value{
  Nothing in particular.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-22}
}
\note{
  Only used for 2022 Canary Rockfish
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{gatherMCMC}},
  \code{\link[PBSsynth]{plotSS.comparisons}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSgetMCMC}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{plotFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
