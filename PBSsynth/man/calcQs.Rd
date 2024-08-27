\name{calcQs}
\alias{calcQs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Quantiles
}
\description{
  Calculate quantiles from MCMC posteriors.
}
\usage{
calcQs(dat, ivec, ovec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{\code{matrix} -- MCMC posterior samples (by row) of parameter estimates or annual derived quantities (by column)}
  \item{ivec}{\code{character} -- vector of \code{'run.rwt.ver'} with elements for every record in \code{'dat'}}
  \item{ovec}{\code{numeric} -- vector of model runs (base and sensitivities), one value for each model run}
}
\details{
	Calculates quantiles (\code{0.05, 0.25, 0.5, 0.75, 0.95}) for each model run
	in the matrix \code{'dat'}.
	The data matrix contains \eqn{x} samples (e.g., \eqn{x}=2000) for each model run 1:\eqn{n},
	and so \code{'dat'} will have \eqn{nx} rows.
}
\value{
  Named list object of matrices; the list is named by a 
  model run descriptor (\code{'run.rwt.ver'}).
  Each element of the list object comprises a matrix of five quantiles
  (by row) for each parameter estimate or annual derived quantity (by column).
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-20}
}
\note{
  Originally a subfunction in appendix functions; 
  promoted to a new function for the 2024 Yellowtail Rockfish assessment.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{mergPA}},
  \code{\link[PBSsynth]{plotSS.senso}},
  \code{\link[PBSsynth]{tabSS.compo}},
  \code{\link[PBSsynth]{tabSS.senso}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
