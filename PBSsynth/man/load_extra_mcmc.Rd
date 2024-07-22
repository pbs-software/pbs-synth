\name{load_extra_mcmc}
\alias{load_extra_mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Load Extra MCMC Data
}
\description{
  Load extra MCMC information from Report files generated for every sample.
  Only used for area-based models at this point.
}
\usage{
load_extra_mcmc(dir.mcmc=".", dir.extra="./sso", quants5=c(0.05,0.25,0.5,0.75,0.95),
   RC=c(TRUE,FALSE), loadCP=FALSE, startyr=1935, run="17v17a", 
   areas=c("5ABC","3CD","5DE"), Fmethod=3, plot=TRUE, png=FALSE, pngres=400,
   PIN=c(9,9), lang="e", show.combo=TRUE, vertical=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.mcmc}{\code{character} -- directory where MCMC results are located;
    can be relative to \code{getwd()}}
  \item{dir.extra}{\code{character} -- directory where extra MCMC results are located;
    can be relative to \code{getwd()}}
  \item{quants5}{\code{numeric} -- vector of five quantiles (usually 0.05, 0.25, 0.5, 0.75, and 0.95)}
  \item{RC}{\code{logical} -- (two-element vector) if \code{c(TRUE,TRUE)}, read in \code{'Results*.sso'} files and \code{'CompResults*.sso'} files}
  \item{loadCP}{\code{logical} -- if \code{TRUE}, load catch policies instead of regular MCMC runs}
  \item{startyr}{\code{numeric} -- start year of the model}
  \item{run}{\code{character} -- MCMC run number (versions amended with lower-case letters}
  \item{areas}{\code{character} -- vector of subareas used in the model}
  \item{Fmethod}{\code{numeric} -- method used to estimate F, where 1=Pope's approximation, 2=Baranov equation, 3=hybrid}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot the various time series}
  \item{png}{\code{logical} -- if \code{TRUE}, send the figures to \code{'.png'} files}
  \item{pngres}{\code{numeric} -- figure resolution in pixels per inch}
  \item{PIN}{\code{numeric} -- plot size for output figures in inches (width, height)}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
  \item{show.combo}{\code{logical} -- if \code{TRUE}, plot all area trajectories in an additional figure panel}
  \item{vertical}{\code{logical} -- if \code{TRUE}, plot all area trajectories in a vertical strip (single column)}
}
\details{
  Extracts data from thousands of Report files (one for each MCMC sample) and,
  if needed, thousands of CompReport files. The sheer magnitude of reading in
  so many files can be achieved through the wizardry of Kevin Ushey's R package
  called \pkg{Kmisc}. Thank you Kevin!
  
  The function then extracts time series, subsets them for desired 
  derived parameters, and plots the results.
}
\value{
  Nothing is returned explicitly; however, a number of binary \code{.rda} files
  are saved (assume 2,000 Report files for 2,000 MCMC samples):
\tabular{rl}{
  \code{reps.rda}        \tab the Report set, but only on the Linux Server;\cr
                         \tab not transferred to PC because it's humungous (ditto \code{comps.rda}).\cr
  \code{mcmc.ts.rda}     \tab \code{list} of time series data frames, one for each MCMC sample;\cr
                         \tab each data frame row comprises subareas and model years with eras \cr
                         \tab \code{'VIRG'}, \code{'INIT'}, \code{'TIME'}, and \code{'FORE'};\cr
                         \tab columns comprise whatever derived parameters SS3 spits out.\cr
  \code{mcmc.ts.sub.rda} \tab \code{list} of data frames that have been subset from \code{mcmc.ts};\cr
                         \tab each list element contains an extracted derived parameter (e.g, \code{'SpawnBio'})\cr
                         \tab as a data frame with the same row elements as those for data frames in \code{mcmc.ts};\cr
                         \tab however, the columns comprise \code{'Area'}, \code{'Yr'}, \code{'Era'} and \code{'Seas'},\cr
                         \tab along with 2,000 columns, each comprising an MCMC sample (e.g. \code{'s0001'}).\cr
  \code{mcmc.posts.rda}  \tab a collection of data frames with suffix \code{'.mcmc'};\cr
                         \tab prefixes: \code{'B'}, \code{'B0'}, \code{'BtBmsy'}, \code{'BtLRP'}, \code{'BtUSR'}, \code{'D'},\cr
                         \tab \code{'fR'}, \code{'MSY'}, \code{'pVB'}, \code{'R'}, \code{'u'}, \code{'utumsy'}, \code{'V'}, \code{'V0'}\cr
  \code{mcmc.quants.rda} \tab sumting\cr
}
}
\references{
Ushey K (2013)
  \href{https://github.com/kevinushey/Kmisc}{Kmisc: Kevin Miscellaneous}.
  R package version 0.5.1,\cr commit 0a35b51944528ed4adb24d3c6687c8c7cc9fe549
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-22}
}
\note{
  Best run on Linux Server to extract the information from thousands of
  Report files using \code{Kmisc::read} (Ushey 2013).
  Once objects have been extracted and saved as binary \code{.rda} files,
  the function can be run on individuals' work stations (in Windows).

  \bold{Important}: MSY values are not broken down by subarea so the code uses
  the observed ratios of \eqn{B_0}{B0} by subarea to partition MSY values.
  Not ideal, but then nothing about using extra MCMC files is ideal.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.pmcmc}},
  \code{\link[PBSsynth]{plotSS.ts}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSgetMCMC}},
  \code{\link[r4ss]{SS_output}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{hplot}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{appFuns.r}
\concept{Linux}
% Use only one concept per line.
