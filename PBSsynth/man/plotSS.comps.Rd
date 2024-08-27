\name{plotSS.comps}
\alias{plotSS.comps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Composition Data
}
\description{
  Plot age frequency/proportion (AF) composition data and fits.
}
\usage{
plotSS.comps(replist, subplots=c(1:21, 24), kind="LEN", sizemethod=1, 
  aalyear=-1, aalbin=-1, plot=TRUE, print=FALSE, fleets="all", 
  fleetnames="default", sexes="all", yupper=0.4, datonly=FALSE, 
  samplesizeplots=TRUE, compresidplots=TRUE, bub=FALSE, showyears=TRUE, 
  showsampsize=TRUE, showeffN=TRUE, aggregates_by_mkt=FALSE, 
  sampsizeline=FALSE, effNline=FALSE, minnbubble=3, pntscalar=NULL, 
  scalebubbles=FALSE, cexZ1=1.5, bublegend=TRUE, 
  colsex=c("orange","limegreen","slategray3"), colvec=lucent(colsex, 0.75), 
  linescol=rep(lucent("black",0.75),3), xlas=0, ylas=NULL, axis1=NULL,
  axis2=NULL, axis1labs=NULL, sizebinlabs=NULL, blue=lucent("limegreen",0.75),
  red=lucent("purple",0.75), pwidth=6.5, pheight=5, punits="in", 
  ptsize=10, res=400, plotdir="default", cex.main=1, linepos=1, 
  fitbar=FALSE, do.sqrt=TRUE, smooth=TRUE, cohortlines=c(), 
  labels=c("Length (cm)", "Age (yr)", "Year", "Observed sample size", 
  "Effective sample size", "Proportion", "cm", "Frequency", "Weight", 
  "Length", "(mt)", "(numbers x1000)", "Stdev (Age)", 
  "Conditional AAL plot, ", "Size bin"), printmkt=TRUE, printsex=TRUE, 
  maxrows=4, maxcols=4, maxrows2=2, maxcols2=4, rows=1, cols=1, 
  andre_oma=c(3,0,3,0), andrerows=3, fixdims=TRUE, fixdims2=FALSE, 
  maxneff=5000, verbose=TRUE, scalebins=FALSE, addMeans=TRUE, 
  mainTitle=FALSE, outnam, lang="e", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- list object created by \code{'SS_output()'}}
  \item{subplots}{\code{numeric} -- vector controlling which subplots to create;
    numbering of subplots is as follows, where subplots 21 to 24
    (aggregated across years) are provided first, and subplots
    1 to 10 are all repeated for each fleet:
  \itemize{
    \item 1 multipanel (by year) fits to age proportions for multiple fleets
    \item 2 single panel bubble plot for numbers at length or age
    \item 3 multi-panel bubble plots for conditional age-at-length
    \item 4 multi-panel plot of fit to conditional age-at-length for specific years
    \item 5 Pearson residuals for A-L key
    \item 6 multi-panel plot of point and line fit to conditional
            age-at-length for specific length bins
    \item 7 sample size plot
    \item 8 TA1.8 Francis weighting plot
    \item 9 TA1.8 Francis weighting plot for conditional data
    \item 10 Andre's mean age and std. dev. in conditional AAL
    \item 21 composition by fleet aggregating across years
    \item 22 composition by fleet aggregating across years within each season
    \item 23 composition by fleet aggregating across seasons within a year
    \item 24 bubble plot comparison of length or age residuals
  }}
  \item{kind}{\code{character} -- indicator of type of plot can be 
    \code{'LEN'}, \code{'SIZE'}", \code{'AGE'}, \code{'cond'}, 
    \code{'GSTAGE'}, \code{'GSTLEN'}, \code{'L@A'}, or \code{'W@A'}}
  \item{sizemethod}{\code{numeric} -- if \code{kind='SIZE'}, this switch
    chooses which of the generalized size bin methods will be plotted}
  \item{aalyear}{\code{numeric} -- years to plot multi-panel conditional 
    age-at-length fits for all length bins; must be in a \code{'c(YYYY,YYYY)'}
    format; useful for checking the fit of a dominant year class, 
    critical time period, etc.; default=-1}
  \item{aalbin}{\code{numeric} -- length bin for which multi-panel plots
    of the fit to conditional age-at-length data will be produced for all years;
    useful to see if growth curves are `ok', or to see the information on
    year classes move through the conditional data; default=-1}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, print to PNG file(s)}
  \item{fleets}{\code{character} -- optional vector to subset fleets for
    which plots will be made}
  \item{fleetnames}{\code{character} -- optional vector of fleet names to put in the labels}
  \item{sexes}{\code{character} -- sexes to show plots for; default="all" which will include
    males, females, and unsexed; this option is not fully implemented for all plots}
  \item{yupper}{\code{numeric} -- upper limit on ymax for polygon/histogram composition plots}
  \item{datonly}{\code{logical} -- if \code{TRUE}, make plots of data without fits}
  \item{samplesizeplots}{\code{logical} -- if \code{TRUE}, make sample size plots}
  \item{compresidplots}{\code{logical} -- if \code{TRUE}, make plots of residuals
    for fit to composition data}
  \item{bub}{\code{logical} -- if \code{TRUE}, make bubble plot for numbers at age or size}
  \item{showyears}{\code{logical} -- if \code{TRUE}, add labels for years to sample size plots}
  \item{showsampsize}{\code{logical} -- if \code{TRUE}, add sample sizes to plot}
  \item{showeffN}{\code{logical} -- if \code{TRUE}, add effective sample sizes to plot}
  \item{aggregates_by_mkt}{\code{logical} -- if \code{TRUE}, separate plots of aggregates
    across years into different plots for each market category (retained, discarded)}
  \item{sampsizeline}{\code{logical} -- if \code{TRUE}, show line for input sample
    sizes on top of conditional age-at-length plots (\emph{still in development})}
  \item{effNline}{\code{logical} -- if \code{TRUE}, show line for effective sample
    sizes on top of conditional age-at-length plots (\emph{still in development})}
  \item{minnbubble}{\code{numeric} -- number of unique x values before adding 
    buffer; type \code{'?bubble3'} for more info}
  \item{pntscalar}{\code{numeric} -- scalar defines the maximum bubble size for bubble
    plots; this option is still available but a better choice is to use \code{'cexZ1'}
    which allow the same scaling throughout all plots}
  \item{scalebubbles}{\code{logical} -- if \code{TRUE}, scale data-only bubbles
    by sample size, not just proportion within sample}
  \item{cexZ1}{\code{numeric} -- character expansion (\code{cex}) for point 
    associated with value of 1}
  \item{bublegend}{\code{logical} -- if \code{TRUE}, add legend with example 
    bubble sizes to bubble plots}
  \item{colsex}{\code{character} -- colours for females, males, and both/unknown/unsexed}
  \item{colvec}{\code{character} -- vector of length 3 with colors for females, males, unsexed fish}
  \item{linescol}{\code{character} -- color for lines on top of polygons}
  \item{xlas}{\code{numeric} -- label style (\code{las}) input for x-axis; 
    default 0 has horizontal labels, input 2 provides vertical labels}
  \item{ylas}{\code{numeric} -- label style (\code{las}) input for y-axis;
    default \code{NULL} has horizontal labels when all labels have fewer 
    than 6 characters and vertical otherwise; input 0 would force vertical
    labels, and 1 would force horizontal ones}
  \item{axis1}{\code{numeric} -- optional position of bottom axis values}
  \item{axis2}{\code{numeric} -- optional position of left size axis values}
  \item{axis1labs}{\code{character} -- optional vector of labels for \code{'axis1'}
    (either \code{NULL} or needs to match length of \code{axis1})}
  \item{sizebinlabs}{\code{character} -- vector of size bin labels corresponding
    to the generalized size frequency method}
  \item{blue}{\code{character} -- color to use for males in bubble plots
    (default was slightly transparent blue, \emph{but RH modified it})}
  \item{red}{\code{character} -- color to use for females in bubble plots
    (default was slightly transparent red, \emph{but RH modified it})}
  \item{pwidth}{\code{numeric} -- width of plots printed to PNG files in 
    units of \code{'punits'}; default is \code{pwidth=6.5}}
  \item{pheight}{\code{numeric} -- height of plots printed to PNG files in 
    units of \code{'punits'}; default is designed to allow two plots per page,
    with \code{'pheight_tall'} used for plots that work best with a taller 
    format and a single plot per page}
  \item{punits}{\code{character} -- units for \code{'pwidth'} and \code{'pheight'};
    can be \code{'px'} (pixels), \code{'in'} (inches), \code{'cm'} (centimeters),
    or \code{'mm'} (millimeters); default is \code{punits='in'}}
  \item{ptsize}{\code{numeric} -- point size for plotted text in plots printed
    to files (see \code{help('png')} in R for details)}
  \item{res}{\code{numeric} -- resolution of plots printed to files;
    default is \code{res=400}}
  \item{plotdir}{\code{character} -- directory where PNG files will be written;
    by default it will be the directory where the model was run}
  \item{cex.main}{\code{numeric} -- character expansion parameter for plot titles}
  \item{linepos}{\code{numeric} -- lines added before points (\code{linepos=1})
    or after (\code{linepos=2})}
  \item{fitbar}{\code{logical} -- if \code{TRUE}, show fit to bars instead of points}
  \item{do.sqrt}{\code{logical} -- if \code{TRUE}, scale bubbles based on \code{sqrt}
    of size vector; type \code{'?bubble3'} for more info}
  \item{smooth}{\code{logical} -- if \code{TRUE}, add loess smoother to observed vs.
    expected index plots and input vs. effective sample size}
  \item{cohortlines}{\code{numeric} -- optional vector of birth years for cohorts
    to add growth curves to numbers at length bubble plots}
  \item{labels}{\code{character} -- vector of labels for plots (titles and axis labels)}
  \item{printmkt}{\code{logical} -- if \code{TRUE}, show market categories in plot titles}
  \item{printsex}{\code{logical} -- if \code{TRUE}, show sex in plot titles}
  \item{maxrows}{\code{numeric} -- maximum (or fixed) number or rows of panels in the plot}
  \item{maxcols}{\code{numeric} -- maximum (or fixed) number or columns of panels in the plot}
  \item{maxrows2}{\code{numeric} -- maximum number of rows for conditional age at length plots}
  \item{maxcols2}{\code{numeric} -- maximum number of columns for conditional age at length plots}
  \item{rows}{\code{numeric} -- number or rows to return to as default for next 
    plots to come or for single plots}
  \item{cols}{\code{numeric} -- number or cols to return to as default for next 
    plots to come or for single plots}
  \item{andre_oma}{\code{numeric} -- outer margins passed to Andre's multi-panel
    conditional age-at-length plots}
  \item{andrerows}{\code{numeric} -- number of rows of Andre's conditional 
    age-at-length plots within each page; default=3}
  \item{fixdims}{\code{logical} -- if \code{TRUE}, fix the dimensions at 
    \code{'maxrows'} by \code{'maxcols'} or resize based on number of years of data}
  \item{fixdims2}{\code{logical} -- if \code{TRUE}, fix the dimensions at 
    \code{'maxrows'} by \code{'maxcols'} in aggregate plots or resize based on number of fleets}
  \item{maxneff}{\code{numeric} -- maximum value to include on plots of input
    and effective sample size; occasionally a calculation of effective N blows
    up to very large numbers, rendering it impossible to observe the relationship
    for other data; default=5000}
  \item{verbose}{\code{logical} -- if \code{TRUE}, return updates of function progress to the R GUI}
  \item{scalebins}{\code{logical} -- if \code{TRUE}, rescale expected and observed proportions
    by dividing by bin width for models where bins have different widths;
    \emph{caution: may not work correctly in all cases}}
  \item{addMeans}{\code{logical} -- if \code{TRUE}, add parameter means in addition to medians
    for MCMC posterior distributions in which the median and mean differ}
  \item{mainTitle}{\code{logical} -- if \code{TRUE}, produce a title for the plot}
  \item{outnam}{\code{character} -- output name for plot (and other data, if applicable)}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
  \item{\dots}{\code{dots} -- additional arguments that will be passed to
    the \code{par} command in the \code{'make.multifig()'} function}
}
\details{
  Plot composition data and fits from Stock Synthesis output.
  Multi-figure plots depend on \code{'make.multifig'} (modified version of \code{'make_multifig'}).
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
  \code{\link[PBSsynth]{make.multifig}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{make_multifig}},
  \code{\link[r4ss]{SSplotComps}}
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
