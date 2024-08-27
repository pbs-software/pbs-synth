\name{mergePA}
\alias{mergePA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Agglomerate Parameters
}
\description{
  Agglomerate parameters from model sensitivity runs into common and infrequent groups.
}
\usage{
mergePA(PA, good, bad)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PA}{\code{numeric} -- matrix array of parameter MPDs or MCMC 
    posteriors from multiple sensitivity model runs}
  \item{good}{\code{character} -- text patterns to identify \emph{good} parameters:
    those that are common to all (most) sensitivity runs}
  \item{bad}{\code{character} -- text patterns to identify \emph{bad} parameters:
    those that occur in only one (a few) sensitivity runs}
}
\details{
  Separate parameter columns into two groups (common and infrequent) for processing
  (e.g., calculating quantiles and creating two tables).
  Most sensitivity runs have the same parameters in common; however,
  some sensitivity runs use parameters other than the common set.
  Additionally, fleets dropped or added create fleet numbering offsets
  and, consequently, parameter mismatches
  (e.g., \eqn{\mu_{\text{QCS(2)}}}{mu_QCS(2)} vs. \eqn{\mu_{\text{QCS(3)}}}{mu_QCS(3)}).
}
\value{
  List object:
  \tabular{ll}{
    \code{PA}   \tab matrix of parameters that occur in all or most sensitivity runs\cr
    \code{exPA} \tab matrix of parameters that occur in one or a few sensitivity runs
  }
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-20}
}
\note{
  New function introduced for 2024 Yellowtail Rockfish assessment.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{calcQs}},
  \code{\link[PBSsynth]{plotSS.senso}}
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
