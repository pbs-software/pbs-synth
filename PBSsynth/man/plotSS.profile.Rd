\name{plotSS.profile}
\alias{plotSS.profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Likelihood Profiles
}
\description{
  Makes a plot of change in negative-log-likelihood for each likelihood
  component that contributes more than some minimum fraction of change in total.
}
\usage{
plotSS.profile(summaryoutput, plot=TRUE, print=FALSE, models="all", 
  profile.string="steep", profile.label="Spawner-recruit steepness (h)", 
  exact=FALSE, ylab="Change in -log-likelihood", components=c("TOTAL", 
  "Catch", "Equil_catch", "Survey", "Discard", "Mean_body_wt", 
  "Length_comp", "Age_comp", "Size_at_age", "SizeFreq", "Morphcomp", 
  "Tag_comp", "Tag_negbin", "Recruitment", "InitEQ_Regime", 
  "Forecast_Recruitment", "Parm_priors", "Parm_softbounds", "Parm_devs", 
  "F_Ballpark", "Crash_Pen"), component.labels=c("Total", "Catch", 
  "Equilibrium catch", "Index data", "Discard", "Mean body weight", 
  "Length data", "Age data", "Size-at-age data", "Generalized size data", 
  "Morph composition data", "Tag recapture distribution", 
  "Tag recapture total", "Recruitment", "Initital equilibrium recruitment", 
  "Forecast recruitment", "Priors", "Soft bounds", "Parameter deviations", 
  "F Ballpark", "Crash penalty"), minfraction=0.01, sort.by.max.change=TRUE, 
  col="default", pch="default", lty=1, lty.total=1, lwd=2, lwd.total=3, 
  cex=1, cex.total=1.5, xlim="default", ymax="default", xaxs="r", yaxs="r", 
  type="o", legend=TRUE, legendloc="topright", pwidth=8, pheight=6, 
  punits="in", res=400, ptsize=10, cex.main=1, plotdir=NULL, add_cutoff=FALSE, 
  cutoff_prob=0.95, verbose=TRUE, lang=c("f", "e"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summaryoutput}{\code{list} -- list created by the function \code{'SSsummarize()'}}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, print to PNG files}
  \item{models}{\code{numeric} -- optional subset of the models described in
    \code{'summaryoutput'}; either \code{'all'} or a vector of numbers indicating
    columns in summary tables}
  \item{profile.string}{\code{character} -- string used to find parameter
    over which the  profile was conducted; if \code{'exact=FALSE'}, 
    this can be a substring of one of the SS parameter labels found in 
    the `Report.sso' file; for instance, the default input \code{'steep'}
    matches the parameter \code{'SR_BH_steep'}; if \code{'exact=TRUE'},
    then profile.string needs to be an exact match to the parameter label}
  \item{profile.label}{\code{character} -- label for x-axis describing 
    the parameter over which the profile was conducted}
  \item{exact}{\code{logical} -- if \code{TRUE}, match the 
    \code{'profile.string'} parameter label exactly}
  \item{ylab}{\code{character} -- label for y-axis; default is 
    \code{'Change in -log-likelihood'}}
  \item{components}{\code{character} -- vector of likelihood components
    that may be included in plot; list is further refined by any 
    components that are not present in model or have little change over
    range of profile (based on limit \code{'minfraction'})}
  \item{component.labels}{\code{character} -- vector of labels for use
    in the legend that matches the vector in \code{'components'}}
  \item{minfraction}{\code{numeric} -- minimum change in likelihood
    (over range considered) as a fraction of change in total likelihood
    for a component to be included in the figure}
  \item{sort.by.max.change}{\code{logical} -- if \code{TRUE}, switch 
    giving option to sort components in legend in order of maximum 
    amount of change in likelihood (over range considered)}
  \item{col}{\code{character} -- optional vector of colors for each line}
  \item{pch}{\code{character|numeric} -- optional vector of plot 
    characters for the points}
  \item{lty}{\code{character|numeric} -- line types for the likelihood components}
  \item{lty.total}{\code{character|numeric} -- line type for the total likelihood}
  \item{lwd}{\code{numeric} -- line width for the likelihood components}
  \item{lwd.total}{\code{numeric} -- line width for the total likelihood}
  \item{cex}{\code{numeric} -- character expansion for the points representing
    the likelihood components}
  \item{cex.total}{\code{numeric} -- character expansion for the points
    representing the total likelihood}
  \item{xlim}{\code{numeric} -- range for x-axis; change in likelihood
    is calculated relative to values within this range}
  \item{ymax}{\code{numeric} -- maximum y-value; default is 10\%
    greater than largest value plotted}
  \item{xaxs}{\code{character} -- style of axis interval calculation
    to be used for the x-axis (type \code{'?par'} for more info)}
  \item{yaxs}{\code{character} -- style of axis interval calculation
    to be used for the y-axis (type \code{'?par'} for more info)}
  \item{type}{\code{character} -- line type (type \code{'?plot'} for more info)}
  \item{legend}{\code{logical} -- if \code{TRUE}, include legend}
  \item{legendloc}{\code{character} -- location of legend (see ?legend for more info)}
  \item{pwidth}{\code{numeric} -- width of plot}
  \item{pheight}{\code{character} -- height of plot}
  \item{punits}{\code{character} -- units for PNG file}
  \item{res}{\code{numeric} -- resolution of plots printed to files;
    default is \code{res=400}}
  \item{ptsize}{\code{numeric} -- point size for PNG file}
  \item{cex.main}{\code{numeric} -- character expansion for plot titles}
  \item{plotdir}{\code{character} -- directory where PNG files will be
    written; by default it will be the directory where the model was run}
  \item{add_cutoff}{\code{logical} -- if \code{TRUE}, add dashed line
    at ~1.92 to indicate 95\% confidence interval based on common 
    cutoff of half of chi-squared of p=0.95 with 1 degree of freedom:
    \code{0.5*qchisq(p=cutoff_prob, df=1)}; the probability value
    can be adjusted using the \code{'cutoff_prob'} below}
  \item{cutoff_prob}{\code{numeric} -- probability associated with
    \code{'add_cutoff'} above}
  \item{verbose}{\code{logical} -- if \code{TRUE}, return updates of
    function progress to the R GUI (\emph{doesn't do anything yet})}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
  \item{\dots}{\code{dots} -- additional arguments passed to \code{'plot()'}}
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
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA\cr
  \href{mailto:ian@iphc.int}{Ian Stewart}, International Pacific Halibut Commission, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-23}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotProfile'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.comps}},
  \code{\link[PBSsynth]{plotSS.francis}},
  \code{\link[PBSsynth]{plotSS.index}},
  \code{\link[PBSsynth]{plotSS.pars}},
  \code{\link[PBSsynth]{plotSS.rdevs}},
  \code{\link[PBSsynth]{plotSS.selex}},
  \code{\link[PBSsynth]{plotSS.stdres}},
  \code{\link[PBSsynth]{plotSS.stock.recruit}},
  \code{\link[PBSsynth]{plotSS.ts}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotProfile}}
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
