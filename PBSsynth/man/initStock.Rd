\name{initStock}
\alias{initStock}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Initialize Stock Objects
}
\description{
  Initialize objects like fleets based on the stock (species and year)
  and the run number used in SS3.
}
\usage{
initStock(strSpp, assyr, run, rwt=NULL, ver=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strSpp}{\code{character} -- species Hart code (for groundfish) or 
    alpha-numeric code for invertebrates (not used in GF assessments)}
  \item{assyr}{\code{numeric} -- year of stock assessment}
  \item{run}{\code{numeric} -- SS3 model run number}
  \item{rwt}{\code{numeric} -- SS3 model reweight number (only for specific circumstances);
     default is \code{NULL}}
  \item{ver}{\code{numeric|character} -- SS3 model run version for specific combinations;
    MPD versions are numeric (e.g., \code{1}) but can be strings (e.g., \code{"1"}); 
    MCMC versions are strings (e.g., \code{"1a"});
    default is \code{NULL}}
}
\details{
  Based on stock model run number (and sometimes reweight and/or version),
  the code returns a list of defined objects.
}
\value{
  List object:
    \item{\code{strSpp}}{input value of \code{strSpp}}
    \item{\code{assyr}}{input value of \code{assyr}}
    \item{\code{run}}{input value of \code{run}}
    \item{\code{rwt}}{input value of \code{rwt}}
    \item{\code{ver}}{input value of \code{ver}}
    \item{\code{species.code}}{species 3-letter code}
    \item{\code{species.name}}{species proper name}
    \item{\code{d.base}}{directory where SS3 runs occur}
    \item{\code{fleets.all}}{all fleets for stock}
    \item{\code{fleets.use}}{fleets used in run}
    \item{\code{fleets.lab}}{fleet labels}
    \item{\code{nfleet}}{number of fleets in run}
    \item{\code{fleets.idx}}{fleet index in run}
    \item{\code{fleets.sel}}{fleets to display selectivity}
    \item{\code{fleets.af}}{fleets with AF data in run}
    \item{\code{fleets.gear}}{fleets with commercial gear}
    \item{\code{gear.names}}{gear names per fleet}
    \item{\code{ngear}}{number of gear types used}
    \item{\code{area.names}}{subarea names}
    \item{\code{narea}}{number of subareas in run}
    \item{\code{col.idx}}{line colors by subarea for abundance indices}
    \item{\code{bg.idx}}{backgound colors by subarea for index points}
    \item{\code{maxage.sel}}{maximum age in selectivity plots}
    \item{\code{assYrs}}{final model years from previous stock assessments}
For some of the earlier stocks, not all objects will appear in the returned list object.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-07-14}
}
\note{
  This function replaces command line code previously housed in \code{'initialise.r'}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{prepMPD}},
  \code{\link[PBSsynth]{runSweave}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
