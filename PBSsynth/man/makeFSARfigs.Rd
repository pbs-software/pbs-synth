\name{makeFSARfigs}
\alias{makeFSARfigs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make FSAR Figures
}
\description{
  Make figures for Fishery Science Advisory Report (FSAR) documents.
}
\usage{
makeFSARfigs(xTS, xRP, xPJ, years=1935:2024, TAC, 
   png=F, pngres=400, PIN=c(10, 7), lang="e", 
   fig.p4=TRUE, fig.snail=FALSE, fig.catch=FALSE, fig.hbar=FALSE, fig.proj=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xTS}{\code{data.frame} -- time series object created by function \code{'gatherMCMC'};
    object can be \code{'xavgTS'} for multi-area models or \code{'avgTS'} from single-area models}
  \item{xRP}{\code{data.frame} -- reference parameter object created by function \code{'gatherMCMC'};
    object can be \code{'xavgRP'} for multi-area models or \code{'avgRP'} from single-area models}
  \item{xPJ}{\code{data.frame} -- projection object created by function \code{'gatherMCMC'};
    object can be \code{'xavgPJ'} for multi-area models or \code{'avgPJ'} from single-area models}
  \item{years}{\code{numeric} -- vector of years to appear in trajectory plots}
  \item{TAC}{\code{numeric} -- data frame of annual Total Allowable Catch (TAC) values by region/fishery}
  \item{png}{\code{logical} -- if \code{TRUE}, send plot(s) to PNG file(s)}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for PNG files}
  \item{PIN}{\code{numeric} -- dimensions (width, height) for PNG files (inches)}
  \item{lang}{\code{character} -- language of text in figures: \code{'e'}=english, \code{'f'}=french}
  \item{fig.p4}{\code{logical} -- if \code{TRUE}, plot the standard 4-panel plot
    (catch, spawning biomass, exploitation, recruitment)}
  \item{fig.snail}{\code{logical} -- if \code{TRUE}, construct phase plot of exploitation vs. spawning biomass}
  \item{fig.catch}{\code{logical} -- if \code{TRUE}, plot regional components of the catch history}
  \item{fig.hbar}{\code{logical} -- if \code{TRUE}, plot the current stock status}
  \item{fig.proj}{\code{logical} -- if \code{TRUE}, plot regional spawnining biomass with various projections}
}
\details{
  Originally designed for 2023 POP, which was a multi-area model with three subareas.
  Retrofitting to single-area models was partially done for 2024 Yellowtail.
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
  Will be an evolving function.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{gatherMCMC}},
  \code{\link[PBSsynth]{load_extra_mcmc}},
  \code{\link[PBSsynth]{plotSS.compo}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{appFuns.r}
\concept{FSAR}
% Use only one concept per line.
