\name{plotSS.pairs}
\alias{plotSS.pairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Parameter Pairs
}
\description{
  Plot pairs of parameters to view relationships between them.
}
\usage{
plotSS.pairs(P.mpd, P.mcmc, type="image", ptypes, 
   lang=c("e","f"), pngres=400, PIN=c(10,10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P.mpd}{\code{numeric} -- vector of MPD values for model parameters}
  \item{P.mcmc}{\code{numeric} -- matrix of MCMC values for model parameters}
  \item{type}{\code{character} -- type of plot, either \code{'image'} or \code{'contour'}}
  \item{ptypes}{\code{character} -- vector of plot types (\code{'win'}, \code{'png'}, \code{'eps'})}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
  \item{pngres}{\code{numeric} -- resolution for figure files (\code{'png'}, \code{'eps'})}
  \item{PIN}{\code{numeric} -- plot size for output figures in inches (width, height)}
}
\details{
  Figures appear in top right corner, parameter names appear along the diagonal,
  and correlation values appear in lower left corner.
}
\value{
  Nada, but a symphony of tie-dye patterns explodes one's mind.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-18}
}
\note{
  Used primarily for diagnostics in MCMC model run summaries.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.pars}},
  \code{\link[PBSsynth]{plotSS.rmcmc}}
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
