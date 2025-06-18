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
plotSS.stdres(replist, kind="AGE", fleets="all",
   fleetnames="default", sexes="both", datonly=FALSE, aggregates_by_mkt = FALSE,
   labels=c("Length (cm)", "Age (yr)", "Year", "Observed sample size",
   "Effective sample size", "Proportion", "cm", "Frequency", "Weight", "Length",
   "(mt)", "(numbers x1000)", "Stdev (Age)", "Conditional AAL plot, ", "Size bin"),
   plot=TRUE, print=FALSE, type="Multinomial", useOSA=TRUE, usePearson=FALSE,
   ptypes="png", pngres=400, PIN=c(7,9), outnam, lang="e", ...)

plt.ageResids (obj, ages=NULL, main=NULL, lang="e", yrs, resfld="Pearson", ...)

plt.yearResids (obj, ages=NULL, main=NULL, fill.in=TRUE, lang="e", resfld="Pearson", ...)

plt.cohortResids (obj, ages=NULL, main=NULL, lang="e", use.rdevs=T, resfld="Pearson", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- report file create by a call to SS3 and 
    read into R by \code{r4ss::SS_output()} (or \code{getSS.output()}, if available) }
  \item{kind}{\code{charater} -- indicator of type of plot:
    \code{'LEN'}, \code{'SIZE'}, \code{'AGE'}, \code{'cond'}, \code{'GSTAGE'}, 
    \code{'GSTLEN'}, \code{'L@A'}, or \code{'W@A'}}
  \item{fleets}{\code{numeric} -- fleet numbers that have AF data}
  \item{fleetnames}{\code{character} -- fleet names}
  \item{sexes}{\code{character} -- label that describes sexes 
    (e.g., \code{'female'}, \code{'male'}, \code{'both'})}
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
  \item{usePearson}{\code{logical} -- if \code{TRUE}, display Pearson residuals used
    by stock assessments before OSA residuals became a thing}
  \item{ptypes}{\code{character} -- plot types (\code{'png'}, \code{'eps'})}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for \code{.png} output}
  \item{PIN}{\code{numeric} -- plot device size (inches) for width and height}
  \item{outnam}{\code{character} -- file output name}
  \item{lang}{\code{character} -- figure language (\code{'e'} for english, \code{'f'} for french)}
  \item{\dots}{\code{dots} -- additional arguments to pass to functions:
    \code{'plt.ageResids'}, \code{'plt.yearResids'}, and \code{'plt.cohortResids'}}
  %% Auxiliary functions
  \bold{Additional arguments:}
  \item{obj}{\code{data.frame} -- age database extracted from \code{replist$agedbase}}
  \item{ages}{\code{numeric} -- subset of available ages}
  \item{main}{\code{character} -- optional title for the plot}
  \item{yrs}{\code{NA} -- not used for anything apparently}
  \item{resfld}{\code{character} -- field name of residuals desired; 
    \code{'OSA'} triggers a cascade of Trijoulets}
  \item{use.rdevs}{\code{logical} -- if \code{TRUE}, call function \code{'getSS.rdevs()'}
    to get recruitment deviations (for colouring cohort boxes)}
  \item{fill.in}{\code{logical} -- if \code{TRUE}, fill the boxes using \code{'pars=tcall(boxpars)'};
    seems to be done regardless for both \code{'plt.ageResids'} and \code{'plt.yearResids'}}
}
\details{
  The function was adapted from the \code{'PBSawatea'} package to plot a triptych of residuals by:
  (a) age class using \code{'plt.ageResids'}, (b) year using \code{'plt.yearResids'}, and
  (c) year of birth using \code{'plt.cohortResids'}. The latter three functions were 
  originally created by Andrew M. Edwards for the POP stock assessment in 2010.

  When SS3 was adopted, the function \code{'r4ss::SSplotComps'} was commandeered to provide
  a wrapper to the AME \code{'plt'} functions.
}
\value{
  Nothing in particular.
}
\references{
  Taylor, I.G., Doering, K.L., Johnson, K.F., Wetzel, C.R., Stewart, I.J., 2021.
  \href{https://doi.org/10.1016/j.fishres.2021.105924}{Beyond visualizing catch-at-age models: Lessons learned from the r4ss package about software to support stock assessments}.
  \emph{Fisheries Research}, \bold{239}:105924 

  Trijoulet, V., Albertsen, C.M., Kristensen, K., Legault, C.M., Miller, T.J. and Nielsen, A. 2023.
  \href{https://doi.org/10.1016/j.fishres.2022.106487}{Model validation for compositional data in stock assessment models: calculating residuals with correct properties}.
  \emph{Fish. Res.} \bold{257}(106487). 1--10.
}
\author{
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA\cr
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew Edwards}, Quantitative Assessment Methods, PBS, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-01-02}
}
\note{
  In 2024, functionality was updated to use one-step-ahead (OSA) residuals by Vanessa Trijoulet and Anders Nielsen.\cr
  Trijoulet, V. and Nielsen, A. (2022)
  \href{https://github.com/fishfollower/compResidual}{compResidual: Residual calculation for compositional observations}
  R package version 0.0.1, commit d4c74845089074d8016454f235044c5d13ded3a5
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotComps'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{calcStdRes}},
  \code{\link[PBSsynth]{getSS.rdevs}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotComps}}\cr
In package \pkg{compResidual}:\cr
  \code{\link[compResidual]{resMulti}}
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
