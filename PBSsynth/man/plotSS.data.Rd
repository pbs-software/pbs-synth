\name{plotSS.data}
\alias{plotSS.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Timeline of Data -- Presence/Absence by Type, Year, and Fleet
}
\description{
  Plot shows graphical display of what data is being used in the model. 
  Some data types may not yet be included. 
  Note, this is based on output from the model, not the input data file.
}
\usage{
plotSS.data(replist, plot=TRUE, print=FALSE, plotdir=".", subplots=1:2, 
   fleetcol="default", datatypes="all", fleets="all", fleetnames="default", 
   ghost=FALSE, pwidth=10, pheight=7.5, punits="in", res=400, ptsize=10, 
   cex.main=1, margins=c(3,1,1,10), cex=2, lwd=12, maxsize=1, alphasize=1,
   mainTitle=FALSE, verbose=TRUE, subplot=lifecycle::deprecated(), lang="e")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- list object created by \code{SS_output()}}
  \item{plot}{\code{logical} -- if \code{TRUE}, send figure to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, send figure to \code{PNG} file}
  \item{plotdir}{\code{character} -- directory where \code{PNG} files will be written}
  \item{subplots}{\code{numeric} -- vector controlling which subplots to create;
    currently there are only two subplots:
    \itemize{
      \item 1 = equal size points showing presence/absence of data type by year/fleet
      \item 2 = points scaled to indicate quantity or precision of data
    }}
  \item{fleetcol}{\code{character|numeric} -- either the string `default' or a vector
    of colors to use for each fleet; If tagging data or environmental data are included,
    an additional color needs to be added for the tag releases which are not assigned to a fleet}
  \item{datatypes}{\code{character} -- either the string `all', or a vector including
    some subset of the following: `catch', `cpue', `lendbase', `sizedbase', `agedbase',
    `condbase', `ghostagedbase', `ghostcondbase', `ghostlendbase', `ladbase', `wadbase',
    `mnwgt', `discard', `tagrelease', `tagdbase1', and `morphcompdbase'}
  \item{fleets}{\code{character|numeric} -- either the string `all', or a vector of 
    numerical values, like c(1,3), listing fleets or surveys to be included in the plot}
  \item{fleetnames}{\code{character} -- optional replacement for fleetnames used in data file}
  \item{ghost}{\code{logical} -- if \code{TRUE}, show presence of composition data from 
    ghost fleets (data for which fit is shown, but not included in the likelihood calculations)}
  \item{pwidth}{\code{numeric} -- width of plots printed to files in units of \code{'punits'}}
  \item{pheight}{\code{numeric} -- height of plots printed to \code{png} files in units
    of \code{'punits'}; default is designed to allow two plots per page, with \code{'pheight_tall'}
    used for plots that work best with a taller format and a single plot per page}
  \item{punits}{\code{character} -- units for pwidth and pheight; can be `px' (pixels),
    `in' (inches), `cm' (centimeters), or `mm' (millimeters); default is \code{punits="in"}}
  \item{res}{\code{numeric} -- resolution of plots printed to files; default is \code{res=400} pixels/inch}
  \item{ptsize}{\code{numeric} -- point size for plotted text in plots printed to files; 
    type \code{help("png")} on the R command line for details}
  \item{cex.main}{\code{numeric} -- character expansion for plot titles; default is \code{cex.main=1}}
  \item{margins}{\code{numeric} -- margins of plot (passed to \code{par()} function), 
    which may need to be increased if fleet names run off the right-hand margin}
  \item{cex}{\code{numeric} -- character expansion for points showing isolated years of data}
  \item{lwd}{\code{numeric} -- line width for plot elements}
  \item{maxsize}{\code{numeric} -- size (cex) of the largest bubble in the datasize plot}
  \item{alphasize}{\code{numeric} -- transparency of bubbles in datasize plot; 
    defaults to 1 (no transparency); useful for models with lots of overlapping points}
  \item{mainTitle}{\code{logical} -- if \code{TRUE}, include a title at the top
    (not yet implemented for all plots)}
  \item{verbose}{\code{logical} -- if \code{TRUE}, print output to the screen}
  \item{subplot}{\code{deprecated} -- use \code{subplots} instead}
  \item{lang}{\code{character} -- language (\code{'e'}=english, \code{'f'}=french) to use in plots}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the non-fungible function \code{'SSplotData'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific purpose.
}
\author{
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA\cr
	Chantel Wetzel and Cole Monnahan (both from NOAA)

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-06-10}
}
\note{
  Function was retro-fitted to yield french-language plots.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotData}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{r4ssFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
