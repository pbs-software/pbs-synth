\name{plotSS.stdres}
\alias{plotSS.stdres}
\alias{plt.ageResids}
\alias{plt.yearResids}
\alias{plt.cohortResids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Studentised Residuals
}
\description{
  Plot Studentised (Pearson) residuals reported in SS3 output.
}
\usage{
plotSS.stdres(replist, kind="AGE", fleets="all", fleetnames="default", 
   sexes="both", datonly=FALSE, aggregates_by_mkt=FALSE, 
   labels=c("Length (cm)", "Age (yr)", "Year", "Observed sample size", 
   "Effective sample size", "Proportion", "cm", "Frequency", "Weight", 
   "Length", "(mt)", "(numbers x1000)", "Stdev (Age)", 
   "Conditional AAL plot, ", "Size bin"), plot=TRUE, print=FALSE, 
   type="Multinomial", useOSA=FALSE, ptypes="png", pngres=400, 
   PIN=c(7, 9), outnam, lang="e", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- report file create by a call to SS3 and 
    read into R by \code{r4ss::SS_output()} }
  \item{kind}{\code{charater} -- indicator of type of plot:
    \code{'LEN'}, \code{'SIZE'}, \code{'AGE'}, \code{'cond'}, \code{'GSTAGE'}, 
    \code{'GSTLEN'}, \code{'L@A'}, or \code{'W@A'}}
  \item{fleets}{\code{numeric} -- fleet numbers that have AF data}
  \item{fleetnames}{\code{character} -- fleet names}
  \item{sexes}{\code{character} -- label that describes sexes 
    (e.g., \code{'both'}, \code{'female'}, \code{'male'})}
  \item{datonly}{\code{logical} -- if \code{TRUE}, make plots of data without fits}
  \item{aggregates_by_mkt}{\code{logical} -- if \code{TRUE}, make separate plots of aggregates
    across years into different plots for each market category (retained, discarded)}
  \item{labels}{\code{character} -- vector of labels for plots (titles and axis labels)}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, print to PNG files}
  \item{type}{\code{character} -- type of distribution to fit AF data; 
    if \code{'Fournier'} or \code{'Coleraine'}, use function \code{'calcStdRes'} to 
    calculate Pearson residuals, else use SS3 calculation of Pearson residuals}
  \item{useOSA}{\code{logical} -- if \code{TRUE}, calculate one-step-ahead residuals 
    using function\cr \code{compResidual::resMulti}}
  \item{ptypes}{\code{character} -- plot types (\code{'png'}, \code{'eps'})}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for \code{.png} output}
  \item{PIN}{\code{numeric} -- plot device size (inches) for width and height}
  \item{outnam}{\code{character} -- file output name}
  \item{lang}{\code{character} -- figure language (\code{'e'} for english, \code{'f'} for french)}
  \item{\dots}{\code{dots} -- additional arguments to pass to functions:
    \code{'plt.ageResids'}, \code{'plt.yearResids'}, and \code{'plt.cohortResids'}}
}
\details{
  The function was adapted from \pkg{PBSawatea} to plot a triptych of residuals by:
  (a) age class using \code{'plt.ageResids'}, (b) year using \code{'plt.yearResids'}, and
  (c) year of birth using \code{'plt.cohortResids'}. The latter three functions were 
  originally created by Andrew M. Edwards for the POP stock assessment in 2010.
  Rowan Haigh revised these in years following.
}
\value{
  Function returns nothing; creates a plot.
}
\author{
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew Edwards}, Quantitative Assessment Methods\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-08 (eclipse)}
}
\note{
  Added functionality to use one-step-ahead residuals by Vanessa Trijoulet and Anders Nielsen.\cr
  Trijoulet, V. and Nielsen, A. (2022)\cr
  \href{https://github.com/fishfollower/compResidual}{compResidual: Residual calculation for compositional observations}\cr
  R package version 0.0.1, commit d4c74845089074d8016454f235044c5d13ded3a5
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{calcStdRes}}\cr
In package \pkg{compResidual}:\cr
  \code{\link[compResidual]{resMulti}}
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
