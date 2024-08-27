\name{plotSS.francis}
\alias{plotSS.francis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Mean Ages \ifelse{latex}{\out{\`{a}}}{\ifelse{html}{\out{&#224;}}{a}} la Francis
}
\description{
  Plot mean age fits using Francis (2011) methodology.
}
\usage{
plotSS.francis(fit, type, fleet, part=0:2, sexes=0:3, seas=NULL, 
   method=NULL, plotit=TRUE, printit=TRUE, datonly=FALSE, 
   plotadj=!datonly, maxpanel=1000, fleetnames=NULL, label.part=TRUE, 
   label.sex=TRUE, set.pars=TRUE, col.obs=c("green3","green"), 
   col.fit=lucent(c("blue2","cyan"), 0.5), png=FALSE, pngres=400, 
   PIN=c(8,9), outnam, lang=c("e","f") )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{\code{list} -- SS3 report list from a call to function \code{'r4ss::SS_output'}}
  \item{type}{\code{character} -- either \code{'len'} (for length composition data),
    \code{'size'} (for generalized size composition data),
    \code{'age'} (for age composition data), or
    \code{'con'} (for conditional age-at-length data)}
  \item{fleet}{\code{numeric} -- vector of one or more fleet numbers whose data are to
    be analysed simultaneously (the output N multiplier applies to all fleets combined)}
  \item{part}{\code{numeric} -- vector of one or more partition values; analysis is restricted
    to composition data with one of these partition values;
    default is to include all partition values (0,1,2)}
  \item{sexes}{\code{numeric} -- vector of one or more values for \code{Sexes};
    analysis is restricted to composition data with one of these
    \code{Sexes} values; ignored if \code{type=='con'}}
  \item{seas}{\code{character} -- how to treat data from multiple seasons:
    \code{'comb'} -- combine seasonal data for each year and plot against \code{Yr};
    \code{'sep'} -- treat seasons separately, plotting against \code{Yr}; 
    if \code{is.null(seas)}, it is assumed that there is only one season in
    the selected data (a warning is output if this is not true) and
    option \code{'comb'} is used}
  \item{method}{\code{character} -- vector of one or more size-frequency method numbers
    (ignored unless \code{type='size'}); if \code{!is.null(method)}, analysis is restricted
    to size-frequency methods in this vector; NB comps are separated by method}
  \item{plotit}{\code{logical} -- if \code{TRUE}, make an illustrative plot like one or more
    panels of Fig.4 in Francis (2011)}
  \item{printit}{\code{logical} -- if \code{TRUE}, print results to R console}
  \item{datonly}{\code{logical} -- if \code{TRUE}, don't show the model expectations}
  \item{plotadj}{\code{logical} -- if \code{TRUE}, plot the confidence intervals associated with
    the adjusted sample sizes (\code{TRUE} by default unless \code{datonly=TRUE})}
  \item{maxpanel}{\code{numeric} -- maximum number of panels within a plot}
  \item{fleetnames}{\code{character} -- vector of alternative fleet names to draw
    from for plot titles and captions; it should have length equal to the number
    of fleets in the model, not the number of fleets considered in this function}
  \item{label.part}{\code{logical} -- if \code{TRUE}, include labels indicating which partitions are included}
  \item{label.sex}{\code{logical} -- if \code{TRUE}, include labels indicating which sexes are included}
  \item{set.pars}{\code{logical} -- if \code{TRUE}, set the graphical parameters
    such as \code{'mar'} and \code{'mfrow'}; can be set to \code{FALSE} in order to 
    add plots form multiple calls to this function as separate panels in one larger figure}
  \item{col.obs}{\code{character} -- vector of colours for observed mean ages}
  \item{col.fit}{\code{character} -- vector of colours for fitted mean ages}
  \item{png}{\code{logical} -- if \code{TRUE}, send figure to a PNG file}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for PNG file}
  \item{PIN}{\code{numeric} -- dimensions (width, height) for PNG file (inches)}
  \item{outnam}{\code{character} -- output name for figures and/or tables}
  \item{lang}{\code{character} -- language of text in figures: \code{'e'}=english, \code{'f'}=french}
}
\details{
	Calculate Francis (2011) mean-age weights for the next reweighting of AF data.
	If \code{'plotit=TRUE'}, create a figure showing curve of mean-age fit by year
	overtop the distribution of input (observed) ages.
}
\value{
List of objects:
\preformatted{list(agedat=pldat, w=Nmult, lo=confint[1,], hi=confint[2,], w.francis=wj)}
}
\references{
  Francis, R.I.C.C. 2011. 
  \href{https://doi.org/10.1139/f2011-025}{Data weighting in statistical fisheries stock assessment models}.
  \emph{Can. J. Fish. Aquat. Sci.} \bold{68(6)}. 1124--1138.

  Taylor, I.G., Doering, K.L., Johnson, K.F., Wetzel, C.R., and Stewart, I.J. 2021.
  \href{https://doi.org/10.1016/j.fishres.2021.105924}{Beyond visualizing catch-at-age models: Lessons learned from the r4ss package about software to support stock assessments}.
  \emph{Fisheries Research}, \bold{239}:105924
}
\author{
  \href{mailto:francisc@niwa.co.nz}{Chris Francis}, National Institute of Water & Atmospheric Research, Wellington NZ\cr
  \href{mailto:aepunt@uw.edu}{Andre Punt}, School of Aquatic & Fishery Sciences, Univ. of Washington, Seattle WA\cr
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-21}
}
\note{
  Years ago, \code{'r4ss::SSMethod.TA1.8'} appeared to have a bug when processing multiple fleets.
  This was fixed in \code{'plotSS.francis'}.
  The issue might have been fixed in the most recent version of \code{r4ss}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSMethod.TA1.8'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{calcMA}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSMethod.TA1.8}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{hplot}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{overFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
