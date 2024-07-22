\name{plotSS.selex}
\alias{plotSS.selex}
\alias{plt.selectivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Selectivity
}
\description{
  Plot selectivity in the style to which we have become accustomed.
}
\usage{
plotSS.selex(replist, infotable=NULL, fleets="all", fleetnames="default",
   sizefactors=c("Lsel"), agefactors=c("Asel", "Asel2"), 
   years="endyr", minyr=-Inf, maxyr=Inf, maxage=50, season=1, 
   sexes="all", selexlines=1:6, subplot=101, skipAgeSelex10=TRUE,
   xlim, plot=TRUE, print=FALSE, add=FALSE, 
   labels=c("Length (cm)","Age (yr)","Year","Selectivity","Retention","Discard mortality"),
   col1="red", col2="blue", lwd=2, spacepoints=5, staggerpoints=1,
   legendloc="bottomright", pwidth=7, pheight=7, punits="in", res=400,
   ptsize=12, cex.main=1, showmain=TRUE, plotdir="default", verbose=TRUE,
   debug=FALSE, sobj=NULL, lang=c("e", "f"))

plt.selectivity <- function( obj, sobj=NULL, mainTitle="Rockfish", maxage,
   ptypes="win", pngres=400, PIN=c(7,8), lang=c("e","f"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- object created by \code{SS_output()}}
  \item{infotable}{\code{data.frame} -- optional table of information 
    controlling appearance of plot and legend; produced as output and 
    can be modified and entered as input}
  \item{fleets}{\code{numeric} -- optional vector of fleets}
  \item{fleetnames}{\code{character} -- optional replacement for fleetnames}
  \item{sizefactors}{\code{numeric|character} -- vector of factors within column 
    \code{'SIZE_SELEX'} to include in plot of selectivity across multiple fleets}
  \item{agefactors}{\code{numeric|character} -- vector of factors within column 
    \code{'AGE_SELEX'} to include in plot of selectivity across multiple fleets}
  \item{years}{\code{numeric} -- years for selectivity, shown in multi-line plot
    (default = last year of model)}
  \item{minyr}{\code{numeric} -- minimum year to show in plots (optional)}
  \item{maxyr}{\code{numeric} -- maximum year to show in plots (optional)}
  \item{maxage}{\code{numeric} -- maximum age to show in plots}
  \item{season}{\code{numeric} -- if seasonal model, season to show in plots (default=1)}
  \item{sexes}{\code{numeric} -- vector of sexes (1=female, 2=male)}
  \item{selexlines}{\code{numeric} -- selection of lines to plot (1=selectivity,
    2=retention, 3=discard mortality, and 4=keep)}
  \item{subplot}{\code{numeric} -- vector of subplot types (see \code{\link[r4ss]{SSplotSelex}});
    additionally:\cr
  \code{101} -- Awatea-type selectivity plot with maturity curve underlay;\cr
  \code{102} -- Like \code{101} but add Awatea's estimated selectivity for comparison.}
  \item{skipAgeSelex10}{\code{numeric} -- exclude plots for age selectivity type 10}
  \item{xlim}{\code{numeric} -- limits of the x-axis}
  \item{plot}{\code{logical} -- if \code{TRUE}, send plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, send plot a PNG file}
  \item{add}{\code{logical} -- if \code{TRUE}, add to existing plot (not yet implemented)}
  \item{labels}{\code{character} -- vector of labels for plots (titles and axis labels)}
  \item{col1}{\code{numeric|character} -- colour for female growth curve}
  \item{col2}{\code{numeric|character} -- colour for male growth curve}
  \item{lwd}{\code{numeric} -- line widths for plots}
  \item{spacepoints}{\code{numeric} -- number of years between points shown on top of lines
    (for long timeseries, points every year get mashed together)}
  \item{staggerpoints}{\code{numeric} -- number of years to stagger the first point
    (if \code{spacepoints > 1}) for each line (so that adjacent lines have points in
    different years)}
  \item{legendloc}{\code{numeric|character} -- location of legend}
  \item{pwidth}{\code{numeric} -- default width of plots printed to files in units of
    \code{punits}; default is \code{pwidth=6.5}}
  \item{pheight}{\code{numeric} -- height of plots printed to \code{'.png'} files
    in units of \code{punits}; default is designed to allow two plots per page, 
    with \code{pheight_tall} used for plots that work best with a taller format and
    a single plot per page}
  \item{punits}{\code{character} -- units for \code{pwidth} and \code{pheight}.
    Can be \code{'px'} (pixels), \code{'in'} (inches), \code{'cm'} (centimeters),
    or \code{'mm'} (millimeters); the default is \code{punits="in"}}
  \item{res}{\code{numeric} -- resolution of plots printed to files;
    the default is \code{res = 300}}
  \item{ptsize}{\code{numeric} -- point size for plotted text in plots printed to files
    (see \code{help("png")} in R for details)}
  \item{cex.main}{\code{numeric} -- character expansion for plot titles}
  \item{showmain}{\code{logical} -- if \code{TRUE}, show the main title;
    note: more recent versions of \code{SSplotSelex} have deprecated this argument in 
    favour of a new argument called \code{'mainTitle'}}
  \item{plotdir}{\code{character} -- directory where PNG files will be written; 
    by default it will be the directory where the model was run}
  \item{verbose}{\code{logical} -- if \code{TRUE}, print output to the screen}
  \item{debug}{\code{logical} -- if \code{TRUE}, call subfunction \code{'sumtingwong'}
    to display any error messages that might occur}
  \item{sobj}{\code{numeric} -- second object, usually an Awatea selectivity ogive}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french.}
%% arguments in plt.selectivity not already specified above
  \item{mainTitle}{\code{character} -- main title for plot, usually the species name}
  \item{maxage}{\code{numeric} -- maximum age for the y-axis}
  \item{ptypes}{\code{character} -- vector of plot types (\code{'win'}, \code{'png'}, \code{'eps'})}
  \item{pngres}{\code{numeric} -- resolution for figure files (\code{'png'}, \code{'eps'})}
  \item{PIN}{\code{numeric} -- plot size for output figures in inches (width, height)}
}
\details{
  Function has been borrowed from the package \pkg{r4ss}.
  Specifically, the function \code{'SSplotSelex'} has been modified to
  give additional subplots:\cr
  \code{101} -- Awatea-type selectivity plot with maturity curve underlay;\cr
  \code{102} -- Like \code{101} but add Awatea's estimated selectivity for comparison.
}
\value{
  Invisibly returns a list of \code{'infotable2'} and \code{'plotinfo'}.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-16}
}
\note{
  Thank you Ian Stewart and Ian Taylor.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.stdres}}, and other \code{plotSS} functions.
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotSelex}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{borrow}
\concept{plotFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
