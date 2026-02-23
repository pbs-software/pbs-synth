\name{compAE}
\alias{compAE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compare Ageing Error
}
\description{
  Compare ageing error vectors used in SS3 model.
}
\usage{
compAE(fnlen, fnage, fnmod=NULL, CASAL=0.1, maxage=45, strSpp="417",
   tabs=TRUE, png=FALSE, pngres=400, PIN=c(10,7), lang=c("f","e"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnlen}{\code{character} -- file name containing CVs-at-length}
  \item{fnage}{\code{character} -- file name containing CVs-at-age}
  \item{fnmod}{\code{character} -- file name containing CVs-at-age from Punt's AgeingError R package}
  \item{CASAL}{\code{numeric} -- constant CV adopted by the CASAL model}
  \item{maxage}{\code{numeric} -- maximum age used in SS3 (plus group)}
  \item{strSpp}{\code{character} -- string species Hart code}
  \item{tabs}{\code{logical} -- if \code{TRUE}, save compiled data to a CSV file}
  \item{png}{\code{logical} -- if \code{TRUE}, send figure to PNG device (file)}
  \item{pngres}{\code{numeric} -- resoltion (pixels/inch) of PNG file}
  \item{PIN}{\code{numeric} -- width and height of PNG file}
  \item{lang}{\code{character} -- language to render in figures: \code{'e'}=english, \code{'f'}=french}
}
\details{
  Plot standard deviations at age used for AE vectors in SS3.
}
\value{
  Nada but tables and figures.
}
\references{
Punt, A.E., Smith, D.C., KrusicGolub, K. and Robertson, S. (2008)
\href{https://doi.org/10.1139/F08-111}{%
Quantifying age-reading error for use in fisheries stock assessments, with application to species in Australia's southern and eastern scalefish and shark fishery}.
Can. J. Fish. Aquat. Sci. 65. 1991--2005.

Richards, L.J., Schnute, J.T., Kronlund, A.R. and Beamish, R.J. (1992)
\href{https://doi.org/10.1139/f92-200}{%
Statistical models for the analysis of ageing error}.
Can. J. Fish. Aquat. Sci. 49(9). 1801--1815.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-02-23}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{compWT}}\cr
In package \pkg{PBStools}:\cr
  \code{\link[PBStools]{compAF}},
  \code{\link[PBStools]{compLen}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{plotFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
