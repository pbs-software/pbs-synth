\name{getSS.rdevs}
\alias{getSS.rdevs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get Recruitment Deviations
}
\description{
  Get recruitment deviations from an SS3 results file.
}
\usage{
getSS.rdevs(replist, forecast=FALSE, minyr=-Inf, maxyr=Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- report file create by a call to SS3 and 
    read into R by \code{r4ss::SS_output()}}
  \item{forecast}{\code{logical} -- if \code{TRUE}, also grab the forecast
    recruitment deviations; if \code{FALSE}, only one forecast Rdev is 
    returned (for current model year)}
  \item{minyr}{\code{numeric} -- minimum year of Rdev extraction}
  \item{maxyr}{\code{numeric} -- maximum year of Rdev extraction}
}
\details{
  Essentially an extraction from \code{'replist'} that mirrors code
  found in one of the \code{r4ss} functions.
}
\value{
  A data frame of Early, Main, Late, and Forecast recruitment deviations (identified by row name).
  Three columns (fields) are retrieved: \code{'Yr'}, \code{'Value'}, and \code{'Parm_StDev'}.
}
\references{
  Taylor, I.G., Doering, K.L., Johnson, K.F., Wetzel, C.R., Stewart, I.J., 2021.
  \href{https://doi.org/10.1016/j.fishres.2021.105924}{Beyond visualizing catch-at-age models: Lessons learned from the r4ss package about software to support stock assessments}.
  \emph{Fisheries Research}, \bold{239}:105924 
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-20}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.selex}},
  \code{\link[PBSsynth]{plt.cohortResids}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{data}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{r4ssFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
