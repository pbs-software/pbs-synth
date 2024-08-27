% Funcion is borrowed from 'r4ss::SSplotComparisons'
\name{plotSS.comparisons}
\alias{plotSS.comparisons}
\title{
  Plot Retrospective Model Comparisons
}
\description{
  Plot retrospective model results in same plotting space to compare trajectories.
}
\usage{
plotSS.comparisons(
  summaryoutput,  subplots=1:20, plot=TRUE, print=FALSE, png=print,
  pdf=FALSE, models="all", endyrvec=NULL, indexfleets=NULL, 
  indexUncertainty=TRUE, indexQlabel=TRUE, indexQdigits=4, 
  indexSEvec=NULL, indexPlotEach=FALSE, labels=c("Year", 
  "Spawning biomass (t)", "Fraction of unfished", 
  "Age-0 recruits (1,000s)", "Recruitment deviations", "Index", 
  "Log index", "SPR-related quantity", "Density", "Management target",
  "Minimum stock size threshold", "Spawning output", "Harvest rate"),
  col=NULL, shadecol=NULL, pch=NULL, lty=1, lwd=2, spacepoints=10,
  staggerpoints=1, initpoint=0, tickEndYr=TRUE, shadeForecast=TRUE,
  xlim=NULL, ylimAdj=1.05, xaxs="i", yaxs="i", type="o", 
  uncertainty=TRUE, shadealpha=0.1, legend=TRUE, legendlabels=NULL, 
  legendloc="topright", legendorder=NULL, legendncol=1, sprtarg=NULL, 
  btarg=NULL, minbthresh=NULL, pwidth=8.5, pheight=6, punits="in", 
  res=400, ptsize=12, plotdir=NULL, filenameprefix="", 
  densitynames=c("SSB_Virgin", "R0"), densityxlabs=NULL, rescale=TRUE,
  densityscalex=1, densityscaley=1, densityadjust=1, 
  densitysymbols=TRUE, densitytails=TRUE, densitymiddle=FALSE, 
  densitylwd=1, fix0=TRUE, new=TRUE, add=FALSE, 
  par=list(mar=c(5,4,1,1)+0.1), verbose=TRUE, mcmcVec=FALSE, 
  show_equilibrium=TRUE )
}
\arguments{
  \item{summaryoutput}{\code{list} -- created by function \code{'SSsummarize'}}
  \item{subplots}{\code{numeric} -- vector of subplots to be created;\cr
    Numbering of subplots is as follows:
    \itemize{
      \item 1  spawning biomass
      \item 2  spawning biomass with uncertainty intervals
      \item 3  biomass ratio (hopefully equal to fraction of unfished)
      \item 4  biomass ratio with uncertainty
      \item 5  SPR ratio
      \item 6  SPR ratio with uncertainty
      \item 7  F value
      \item 8  F value with uncertainty
      \item 9  recruits
      \item 10  recruits with uncertainty
      \item 11  recruit devs
      \item 12  recruit devs with uncertainty
      \item 13  index fits
      \item 14  index fits on a log scale
      \item 15  phase plot
      \item 16  densities
      \item 17  cumulative densities
  }}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, send plots to PNG files in 
    directory specified by \code{'plotdir'}}
  \item{png}{\code{logical} -- if \code{TRUE}, send plots to PNG files in 
    directory specified by \code{'plotdir'} (has same result as \code{'print'},
    included for consistency with \code{SS_plots})}
  \item{pdf}{\code{logical} -- if \code{TRUE}, write output to PDF file; 
    cannot be used in conjunction with \code{'png'} or \code{'print'}}
  \item{models}{\code{character} -- optional subset of the models described in \code{summaryoutput};
    either \code{'all'} or a vector of numbers indicating columns in summary tables}
  \item{endyrvec}{\code{numeric} -- optional single year or vector of years representing the final year of
    values to show for each model; by default it is set to the ending year specified in each model}
  \item{indexfleets}{\code{numeric} -- fleet numbers for each model to compare indices of abundance;
    can take different forms:
    \itemize{
      \item NULL: (default) create a separate plot for each index as long as the fleet
      numbering is the same across all models
      \item integer: create a single comparison plot for the chosen index
      \item vector of length equal to number of models: a single fleet number
      for each model to be compared in a single plot
      \item list: list of fleet numbers associated with indices within each
      model to be compared, where the list elements are each a vector of the
     same length but the names of the list elements don't matter and can be absent
    }}
  \item{indexUncertainty}{\code{logical} -- if \code{TRUE}, show uncertainty intervals on index data;
    default=\code{FALSE} because if models have any extra standard deviations added,
    these intervals may differ across models}
  \item{indexQlabel}{\code{logical} -- if \code{TRUE}, add catchability to legend in plot of index fits}
  \item{indexQdigits}{\code{numeric} -- number of significant digits for catchability in legend
    (if \code{indexQlabel=TRUE})}
  \item{indexSEvec}{\code{numeric} -- optional replacement for the SE values in 
    \code{summaryoutput[["indices"]]} to deal with the issue of differing uncertainty by models described above}
  \item{indexPlotEach}{\code{logical} -- if \code{TRUE}, plot the observed index for each model with
    colors; if \code{FALSE}, plot observed once in black dots}
  \item{labels}{\code{character} -- vector of labels for plots (titles and axis labels)}
  \item{col}{\code{character|numeric} -- optional vector of colors to be used for lines;
    input \code{NULL} makes use of \code{rich.colors.short} function.}
  \item{shadecol}{\code{character|numeric} -- optional vector of colors to be used for shading
    uncertainty intervals; the default (\code{NULL}) defaults to using the same colors provided by
    \code{col} (either the default or a user-chosen input) and make them more transparent
    by applying the \code{shadealpha} input as an alpha transparency value (using the \code{adjustcolor()} function)}
  \item{pch}{\code{character|numeric} -- optional vector of plot character values}
  \item{lty}{\code{character|numeric} -- optional vector of line types}
  \item{lwd}{\code{numeric} -- optional vector of line widths}
  \item{spacepoints}{\code{numeric} -- number of years between points shown on top of lines
    (for long timeseries, points every year get mashed together)}
  \item{staggerpoints}{\code{numeric} -- number of years to stagger the first point
    (if \code{spacepoints > 1}) for each line (so that adjacent lines have points in different years)}
  \item{initpoint}{\code{numeric} -- year value for first point to be added to lines;
    points added to plots are those that satisfy 
    (Yr-initpoint)\%\%spacepoints == (staggerpoints*iline)\%\%spacepoints}
  \item{tickEndYr}{\code{logical} -- if \code{TRUE}, switch to turn on extra axis mark at final
    year in timeseries plots}
  \item{shadeForecast}{\code{logical} -- if \code{TRUE}, switch to turn on shading of years beyond
    the maximum ending year of the models}
  \item{xlim}{\code{numeric} -- optional x-axis limits}
  \item{ylimAdj}{\code{numeric} -- multiplier for ylim parameter; 
    allows additional white space to fit legend if necessary; default=1.05}
  \item{xaxs}{\code{character} -- choice of xaxs parameter (type \code{'?par'} for more info)}
  \item{yaxs}{\code{character} -- choice of yaxs parameter (type \code{'?par'} for more info)}
  \item{type}{\code{character} -- type parameter passed to points (default \code{'o'}
    overplots points on top of lines)}
  \item{uncertainty}{\code{logical} -- if \code{TRUE}, show plots with uncertainty intervals;
    either a scalar or vector for each model; or a set of integers corresponding to the choice of models}
  \item{shadealpha}{\code{numeric} -- transparency adjustment used to make default 
    \code{shadecol} values (implemented as \code{adjustcolor(col=col, alpha.f=shadealpha)})}
  \item{legend}{\code{logical} -- if \code{TRUE}, add a legend}
  \item{legendlabels}{\code{character} -- optional vector of labels to include in legend;
    default is \code{'model1'}, \code{'model2'},etc.}
  \item{legendloc}{\code{character|numeric} -- location of legend; 
    either a string like \code{"topleft"} or a vector of two numeric values representing 
    the fraction of the maximum in the x and y dimensions, respectively;
    Type \code{'?legend'} for more info on the string options}
  \item{legendorder}{\code{numeric} -- optional vector of model numbers that can be used
    to have the legend display the model names in an order that is different than that
    which is represented in the summary input object}
  \item{legendncol}{\code{numeric} -- number of columns for the legend}
  \item{sprtarg}{\code{numeric} -- target value for SPR-ratio where line is drawn in the
    SPR plots and phase plot}
  \item{btarg}{\code{numeric} -- target biomass value at which to show a line (set to 0 to remove)}
  \item{minbthresh}{\code{numeric} -- minimum biomass threshold at which to show a line (set to 0 to remove)}
  \item{pwidth}{\code{numeric} -- width of plot}
  \item{pheight}{\code{numeric} -- height of plot}
  \item{punits}{\code{character} -- units for PNG file}
  \item{res}{\code{numeric} -- resolution of plots printed to files;
    the default is \code{res=300}}
  \item{ptsize}{\code{numeric} -- point size for PNG file}
  \item{plotdir}{\code{character} -- directory where PNG or PDF files will be written;
    by default it will be the directory where the model was run}
  \item{filenameprefix}{\code{character} -- additional text to append to PNG or PDF file names;
    it will be separated from default name by an underscore}
  \item{densitynames}{\code{character} -- vector of names (or subset of names) of parameters or
    derived quantities contained in \code{summaryoutput[["pars"]][["Label"]]} or\cr
    \code{summaryoutput[["quants"]][["Label"]]} for which to make density plots}
  \item{densityxlabs}{\code{character} -- optional vector of x-axis labels to use in the density
    plots (must be equal in length to the printed vector of quantities that
    match the \code{densitynames} input)}
  \item{rescale}{\code{logical} -- if \code{TRUE}, implement automatic rescaling of units into
    thousands, millions, or billions}
  \item{densityscalex}{\code{numeric} -- scalar for upper x-limit in density plots (values below 1 
    will cut off the right tail to provide better contrast among narrower distributions}
  \item{densityscaley}{\code{numeric} -- scalar for upper y-limit in density plots (values below 1 
    will cut off top of highest peaks to provide better contrast among broader distributions}
  \item{densityadjust}{\code{numeric} -- multiplier on bandwidth of kernel in density function
    used for smoothing MCMC posteriors; type \code{'?density'} and see 'adjust' for details}
  \item{densitysymbols}{\code{logical} -- if \code{TRUE}, add symbols along lines
    in density plots. Quantiles are \code{c(0.025,0.1,0.25,0.5,0.75,0.9,0.975)}}
  \item{densitytails}{\code{logical} -- if \code{TRUE}, shade tails outside of 95\%
    interval darker in density plots}
  \item{densitymiddle}{\code{logical} -- if \code{TRUE}, shade middle inside of 95\%
    interval darker in density plots}
  \item{densitylwd}{\code{numeric} -- line width for density plots}
  \item{fix0}{\code{logical} -- if \code{TRUE}, always include 0 in the density plots}
  \item{new}{\code{logical} -- if \code{TRUE}, create new empty plot window}
  \item{add}{\code{logical} -- if \code{TRUE}, allow single plot to be added to existing figure;
    this needs to be combined with specific 'subplots' input to make sure only one thing gets added}
  \item{par}{\code{list} -- graphics parameter values passed to the \code{par} function}
  \item{verbose}{\code{logical} -- if \code{TRUE}, report progress to R GUI}
  \item{mcmcVec}{\code{logical} -- if \code{TRUE} (scalar or vector), indicate whether input values are
    from MCMC or to use normal distribution around MLE}
  \item{show_equilibrium}{\code{logical} -- if \code{TRUE}, show the equilibrium values for SSB;
    for some model comparisons, these might not be comparable and thus useful to turn off}
}
\references{
  Taylor, I.G., Doering, K.L., Johnson, K.F., Wetzel, C.R., Stewart, I.J., 2021.
  \href{https://doi.org/10.1016/j.fishres.2021.105924}{Beyond visualizing catch-at-age models: Lessons learned from the r4ss package about software to support stock assessments}.
  \emph{Fisheries Research}, \bold{239}:105924 
}
\author{
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA\cr
  \href{mailto:john.wallace@noaa.gov}{John Wallace}, Northwest Fisheries Science Center, NOAA, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-21}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotComparisons'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{doRetros}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{retro}},
  \code{\link[r4ss]{SSplotComparisons}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{overFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
