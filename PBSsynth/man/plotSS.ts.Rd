%% Borrowed from 'r4ss::SSplotTimeseries'
\name{plotSS.ts}
\alias{plotSS.ts}
\title{
  Plot Time Series Data
}
\description{
  Plot time series data contained in \code{'TIME_SERIES'} output from a 
  Stock Synthesis report file. Some values have optional uncertainty intervals.
}
\usage{
plotSS.ts(replist, subplot, add=FALSE, areas="all", areacols="default",
  areanames="default", forecastplot=TRUE, uncertainty=TRUE,
  bioscale=1, minyr=-Inf, maxyr=Inf, plot=TRUE, print=FALSE,
  plotdir="default", verbose=TRUE, btarg="default", minbthresh="default",
  xlab="Year", labels=NULL, punits="in", ptsize=10, cex.main=1, sobj=NULL,
  res=400, outnam, PIN=c(9,9), lang="e")
}
\arguments{
  \item{replist}{\code{list} -- list object created by \code{'SS_output()'}}
  \item{subplot}{\code{numeric} -- number controlling which subplot to create;
    numbering of subplots is as follows, where the spawning biomass plots
    (7 to 10) are provided first when this function is called by \code{'SS_plots()'}:
    \itemize{
      \item 1 Total biomass (tonnes) with forecast
      \item 2 Total biomass by area (spatial models only)
      \item 3 Total biomass (t) at beginning of spawning season with forecast
      \item 4 Summary biomass (t) with forecast
      \item 5 Summary biomass (t) by area (spatial models only)
      \item 6 Summary biomass (t) at beginning of season 1 with forecast
      \item 7 Spawning output with forecast with ~95\% asymptotic intervals
      \item 8 Spawning output by area (spatial models only)
      \item 9 Relative spawning output with forecast with ~95\% asymptotic intervals
      \item 10 Relative spawning output by area (spatial models only)
      \item 11 Age-0 recruits (1,000s) with forecast with ~95\% asymptotic intervals
      \item 12 Age-0 recruits by area (spatial models only)
      \item 13 Fraction of recruits by area (spatial models only)
      \item 14 Age-0 recruits (1,000s) by birth season with forecast
      \item 15 Fraction of total age-0 recruits by birth season with forecast
      \item 101 Spawning biomass, vulnerable male/female biomass, total biomass
      \item 102 Vulnerable biomass
      \item 103 Harvest rate
    }}
  \item{add}{\code{logical} -- if \code{TRUE}, add to existing plot; (not yet implemented)}
  \item{areas}{\code{character} -- optional subset of areas to plot for spatial models}
  \item{areacols}{\code{character} -- vector of colours by area; default uses \code{'rich.colors'} by Arni Magnusson}
  \item{areanames}{\code{character} -- names for areas; default is to use \code{'Area1'}, \code{'Area2'}, ...}
  \item{forecastplot}{\code{logical} -- if \code{TRUE}, add points from forecast years}
  \item{uncertainty}{\code{logical} -- if \code{TRUE}, add intervals around quantities for which uncertainty is available}
  \item{bioscale}{\code{numeric} -- scaling for spawning biomass; default = 1;
    previously this was set to 0.5 for single-sex models, and 1.0 for all others, 
    but now single-sex models are assumed to use the -1 option for \code{'Nsexes'}
    in the data file so the scaling is done automatically by SS3}
  \item{minyr}{\code{numeric} -- optional input for minimum year to show in plots}
  \item{maxyr}{\code{numeric} -- optional input for maximum year to show in plots}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, print to PNG files}
  \item{plotdir}{\code{character} -- directory where PNG or PDF files will be written;
    by default it will be the directory where the model was run}
  \item{verbose}{\code{logical} -- if \code{TRUE}, report progress to R GUI}
  \item{btarg}{\code{numeric} -- target depletion to be used in plots showing depletion;
    may be omitted by setting to 0; value \code{'default'} chooses value based on modeloutput}
  \item{minbthresh}{\code{numeric} -- threshold depletion to be used in plots showing depletion;
    may be omitted by setting to 0; \code{'default'} assumes 0.25 unless \code{'btarg'} in model
    output is 0.25 in which case \code{minbthresh=0.125} (US west coast flatfish)}
  \item{xlab}{\code{character} -- x-axis label for all plots}
  \item{labels}{\code{character} -- vector of labels for plots (titles and axis labels)}
  \item{punits}{\code{character} -- units for PNG file}
  \item{ptsize}{\code{numeric} -- point size for PNG file}
  \item{cex.main}{\code{numeric} -- character expansion for plot titles}
  \item{sobj}{\code{character} -- originally a selectivity object from Awatea to compare with SS3 output;
    now can be a variety of objects from Awatea}
  \item{res}{\code{numeric} -- resolution (pixels/inch) of plots printed to files}
  \item{outnam}{\code{character} -- output name for figures and/or tables}
  \item{PIN}{\code{numeric} -- dimensions (width, height) for PNG file (inches)}
  \item{lang}{\code{character} -- language of text in figures: \code{'e'}=english, \code{'f'}=french}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotTimeseries'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\references{
  Taylor, I.G., Doering, K.L., Johnson, K.F., Wetzel, C.R., Stewart, I.J., 2021.
  \href{https://doi.org/10.1016/j.fishres.2021.105924}{Beyond visualizing catch-at-age models: Lessons learned from the r4ss package about software to support stock assessments}.
  \emph{Fisheries Research}, \bold{239}:105924 
}
\author{
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA\cr
  \href{mailto:ian@iphc.int}{Ian Stewart}, International Pacific Halibut Commission, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-21}
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.selex}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotTimeseries}}
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
