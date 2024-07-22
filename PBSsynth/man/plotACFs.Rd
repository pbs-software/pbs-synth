\name{plotACFs}
\alias{plotACFs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Autocorrelations for MCMC Parameters
}
\description{
  Plot ACFs for a set of estimated parameters from MCMC samples.
}
\usage{
plotACFs(mcmc, lag.max=60, lang="e")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc}{\code{data.frame} -- MCMC samples (by row) of parameters (by column).}
  \item{lag.max}{\code{numeric} -- maximum lag at which to calculate the autocorelation function (see \code{\link[stats]{acf}}).}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
}
\details{
  Uses the estimated parameters from the data frame \code{currentMCMC$P} and 
  displays the first \code{lag.max} lags using the \code{acf} function from the \code{stats} package.
}
\value{
  Returns nothing in particular.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-16}
}
\note{
  \emph{Caution} -- (from \code{plot.acf}) The confidence interval plotted in \code{plot.acf}
  is based on an \emph{uncorrelated} series and should be treated with appropriate caution. 
  Using \code{ci.type = "ma"} may be less potentially misleading. 
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{panelFuns}},
  \code{\link[PBSsynth]{plotTraj}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
\concept{plotFuns.r}

