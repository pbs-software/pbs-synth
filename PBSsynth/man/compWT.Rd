\name{compWT}
\alias{compWT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Compare Weights
}
\description{
  Compare weights (or lengths) for females (or males) by survey across all years.
}
\usage{
compWT(dat, ssid, sex=2, wfld="len", ylim=c(0, 60), Nmin=10, 
   good=c(1, 3, 4, 16, 21, 79), bad=c(22, 36), outnam, strSpp="417", 
   lang="e", png=FALSE, pngres=400, PIN=c(9, 8))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{\code{composite} -- data frame of biological records (from \code{'gfb_bio.sql'})}
  \item{ssid}{\code{numeric} -- surveys series IDs}
  \item{sex}{\code{numeric} -- sex code: 1=males, 2=females}
  \item{wfld}{\code{character} -- weight field (can be other morphometrics like \code{'len'})}
  \item{ylim}{\code{numeric} -- limits of the y-axis}
  \item{Nmin}{\code{numeric} -- minimum number of \code{wfld} values for survey to be included}
  \item{good}{\code{numeric} -- survey series IDs that will have green boxplots}
  \item{bad}{\code{numeric} -- survey series IDs that will have red boxplots}
  \item{outnam}{\code{character} -- name of output file}
  \item{strSpp}{\code{character} -- string species Hart code}
  \item{lang}{\code{character} -- language to render in figures: \code{'e'}=english, \code{'f'}=french}
  \item{png}{\code{logical} -- if \code{TRUE}, send figure to PNG device (file)}
  \item{pngres}{\code{numeric} -- resoltion (pixels/inch) of PNG file}
  \item{PIN}{\code{numeric} -- width and height of PNG file}
}
\details{
  Compares distributions of a morphometric feature like weight or length
  for a selected sex among available surveys in the biological data file.
}
\value{
  Nada. Plots a set of boxplots comparing the selected morphometric.
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
  \code{\link[PBSsynth]{compTS}}\cr
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
