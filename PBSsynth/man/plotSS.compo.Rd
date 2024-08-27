\name{plotSS.compo}
\alias{plotSS.compo}
\alias{plotSS.senso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create Figures for Base and Sensitivity Runs
}
\description{
  Create figures for a base run (or a composite base case) and for sensitivity runs.
}
\usage{
plotSS.compo(compo, spp.code="YTR", istock="YTR", 
   subset=NULL, redo.figs=FALSE, redo.panels=FALSE, 
   ptypes, lang, pngres=400, PIN=c(9,9))

plotSS.senso(senso, spp.code="YTR", istock="YTR",
   subset=NULL, redo.figs=FALSE, redo.panels=FALSE,
   ptypes, lang, pngres=400, PIN=c(9,9), useRlow=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compo}{\code{list} -- data object containing results from MCMC base run or composite base case runs
    (see function \code{'gatherMCMC'})}
  \item{spp.code}{\code{character} -- three- or four-letter code for species}
  \item{istock}{\code{character} -- three- or four-letter code for species}
  \item{subset}{\code{NA} -- not currently used}
  \item{redo.figs}{\code{logical} -- if \code{TRUE}, re-create the figure PNG files}
  \item{redo.panels}{\code{logical} -- if \code{TRUE}, re-create the panel figures}
  \item{ptypes}{\code{character} -- output types for plots (e.g., \code{'png'})}
  \item{lang}{\code{character} -- language of text in figures: \code{'e'}=english, \code{'f'}=french}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for PNG files}
  \item{PIN}{\code{numeric} -- dimensions (width, height) for PNG files (inches)}
  \item{senso}{\code{list} -- data object containing results from MCMC sensitivity runs
    (see function \code{'gatherMCMC'})}
  \item{useRlow}{\code{logical} -- if \code{TRUE}, assume input data has been reduced to
    low-recruitment samples only (somewhat deprecated, needs attention)}
}
\details{
  Creates PNG files for use later by Sweave routine.
}
\value{
  Nothing in particular.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-16}
}
\note{
  Various input files of MCMC results are characterised by their name:\cr
  \code{'compo'} -- base run (or composite base case)\cr
  \code{'senso'} -- sensitivity runs\cr
  \code{'penso'} -- single-area model MCMCs if main model (\code{'compo'}) is multi-area\cr
  \code{'renso'} -- low-recruitment MCMC samples only
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{gatherMCMC}},
  \code{\link[PBSsynth]{tabSS.compo}},
  \code{\link[PBSsynth]{tabSS.senso}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{hplot}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{appFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
