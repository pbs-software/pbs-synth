\name{plotSS.index}
\alias{plotSS.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Indices of Abundance
}
\description{
  Plot indices of abundance with or without model fit as well as other diagnostic
  plots such as observed vs. expected index and plots related to time-varying
  catchability (if present).
}
\usage{
plotSS.index(replist, subplots=c(1:10,12), plot=TRUE, print=FALSE, 
  fleets="all", fleetnames="default", smooth=TRUE, add=FALSE, datplot=TRUE, 
  labels=list("Year", "Index", "Observed index", "Expected index", 
  "Log index", "Log observed index", "Log expected index", "Standardized index", 
  "Catchability (Q)", "Time-varying catchability", "Vulnerable biomass", 
  "Catchability vs. vulnerable biomass", "Residual", "Deviation"), 
  col1="default", col2="default", col3="blue", col4="red", pch1=21, pch2=16,
  cex=1, bg="white", legend=TRUE, legendloc="topright", seasnames=NULL, 
  pwidth=9, pheight=7, punits="in", res=400, ptsize=10, PIN=c(9, 9), 
  cex.main=1, mainTitle=FALSE, plotdir="default", minyr=NULL, maxyr=NULL, 
  maximum_ymax_ratio=Inf, show_input_uncertainty=TRUE, verbose=TRUE, 
  onepage=FALSE, outnam, lang="e", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- list object created by \code{'SS_output()'}}
  \item{subplots}{\code{numeric} -- vector controlling which subplots to create;
    numbering of subplots is as follows, where subplot 9 (comparison of 
    all indices) is provided first:
  \itemize{
    \item 1  index data by fleet
    \item 2  index data with fit by fleet
    \item 3  observed vs expected index values with smoother
    \item 4  index data by fleet on a log scale (lognormal error only)
    \item 5  index data with fit by fleet on a log scale (lognormal error only)
    \item 6  log(observed) vs log(expected) with smoother (lognormal error only)
    \item 7  time series of time-varying catchability (only if actually time-varying)
    \item 8  catchability vs. vulnerable biomass (if catchability is not constant)
    \item 9  comparison of all indices
    \item 10  index residuals based on total uncertainty
    \item 11  index residuals based on input uncertainty (not currently provided)
    \item 12  index deviations (independent of index uncertainty)
  }}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, print to PNG files}
  \item{fleets}{\code{numeric} -- optional vector to subset fleets for which plots will be made}
  \item{fleetnames}{\code{character} -- optional replacement for \code{'fleenames'}
    used in data file}
  \item{smooth}{\code{logical} -- if \code{TRUE}, add smoothed line to plots of
    observed vs. expected sample sizes}
  \item{add}{\code{logical} -- if \code{TRUE}, add to existing plot
    (\emph{not yet implemented})}
  \item{datplot}{\code{logical} -- if \code{TRUE}, make plot of data only}
  \item{labels}{\code{character} -- vector of labels for plots (titles and axis labels)}
  \item{col1}{\code{character} -- vector of colors for points in each season for 
    time series plot; default is red for single season models and a rainbow using
    the \code{'rich.colors.short'} function for multiple seasons}
  \item{col2}{\code{character} -- vector of colors for points in each season for
    obs. vs. exp. plot; default is blue for single season models and a rainbow using
    the \code{'rich.colors.short'} function for multiple seasons}
  \item{col3}{\code{character} -- color of line showing expected index in time
    series plot; default is blue}
  \item{col4}{\code{character} -- color of smoother shown in obs. vs. exp. plots;
    default is red}
  \item{pch1}{\code{numeric} -- single value or vector of plotting characters
    (\code{pch} parameter) for time-series plots of index fit; default=21}
  \item{pch2}{\code{numeric} -- single value or vector of plotting characters
    (\code{pch} parameter) for sample size plots of index fit; default=16}
  \item{cex}{\code{numeric} -- character expansion factor for points showing
    observed values; default=1}
  \item{bg}{\code{character} -- background color for points with \code{pch}=21}
  \item{legend}{\code{logical} -- if \code{TRUE}, add a legend to seasonal colors
    (only for seasonal models)}
  \item{legendloc}{\code{character} -- legend location for seasonal colors
    (default is \code{'topright'})}
  \item{seasnames}{\code{character} -- optional vector of names for each season
    to replace defaults if a legend is used}
  \item{pwidth}{\code{numeric} -- width of plot}
  \item{pheight}{\code{numeric} -- height of plot}
  \item{punits}{\code{character} -- units for PNG file}
  \item{res}{\code{numeric} -- resolution of plots printed to files;
    default is \code{res=400}}
  \item{ptsize}{\code{numeric} -- point size for PNG file}
  \item{PIN}{\code{numeric} -- dimensions (width, height in inches) for PNG file}
  \item{cex.main}{\code{numeric} -- character expansion for plot titles}
  \item{mainTitle}{\code{logical} -- if \code{TRUE}, add main title to plot}
  \item{plotdir}{\code{character} -- directory where PNG files will be written;
    by default it will be the directory where the model was run}
  \item{minyr}{\code{numeric} -- first year to show in plot (for zooming in
    on a subset of values)}
  \item{maxyr}{\code{numeric} -- last year to show in plot (for zooming in
    on a subset of values)}
  \item{maximum_ymax_ratio}{\code{numeric} -- maximum allowed value for \code{'ymax'}
    (specified as ratio of \code{y}), which overrides any value of \code{'ymax'}
    that is greater (default=\code{Inf})}
  \item{show_input_uncertainty}{\code{logical} -- if \code{TRUE}, switch 
    controlling whether to add thicker uncertainty interval lines indicating 
    the input uncertainty relative to the total uncertainty which may result
    from estimating a parameter for extra standard deviations; this is only 
    added for the plots with index fit included (the data-only plots only 
    show the input uncertainty)}
  \item{verbose}{\code{logical} -- if \code{TRUE}, report progress to R GUI}
  \item{onepage}{\code{logical} -- if \code{TRUE}, plot all figures on one page}
  \item{outnam}{\code{character} -- output name for plot (and other data, if applicable)}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
  \item{\dots}{\code{dots} -- extra arguments to pass to calls to \code{'plot'}}
}
\details{
  Index plots emulate \code{PBSawatea} index plots, which emulate PJS index plots:
  red squares with black error bars.
}
\value{
  Nothing in particular.
}
\references{
  Taylor, I.G., Doering, K.L., Johnson, K.F., Wetzel, C.R., Stewart, I.J., 2021.
  \href{https://doi.org/10.1016/j.fishres.2021.105924}{Beyond visualizing catch-at-age models: Lessons learned from the r4ss package about software to support stock assessments}.
  \emph{Fisheries Research}, \bold{239}:105924 
}
\author{
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-23}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotComps'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.comps}},
  \code{\link[PBSsynth]{plotSS.francis}},
  \code{\link[PBSsynth]{plotSS.pars}},
  \code{\link[PBSsynth]{plotSS.profile}},
  \code{\link[PBSsynth]{plotSS.rdevs}},
  \code{\link[PBSsynth]{plotSS.selex}},
  \code{\link[PBSsynth]{plotSS.stdres}},
  \code{\link[PBSsynth]{plotSS.stock.recruit}},
  \code{\link[PBSsynth]{plotSS.ts}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotIndices}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{r4ssFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
