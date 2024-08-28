\name{plotSS.stock.recruit}
\alias{plotSS.stock.recruit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Spawner-Recruit Curve
}
\description{
  Plot spawner-recruit curve based on output from Stock Synthesis model.
}
\usage{
plotSS.stock.recruit(replist, subplot=1:3, add=FALSE, plot=TRUE, 
  print=FALSE, xlim=NULL, ylim=NULL, labels=c("Spawning biomass (mt)", 
  "Recruitment (1,000s)", "Spawning output", 
  expression(paste("Spawning output (relative to ", italic(B)[0], ")")), 
  expression(paste("Recruitment (relative to  ", italic(R)[0], ")")), 
  "Log recruitment deviation"), bioscale="default", plotdir="default", 
  pwidth=6.5, pheight=6.5, punits="in", res=300, ptsize=10, verbose=TRUE, 
  colvec=c("blue", "black", "black", gray(0,0.7)), ltyvec=c(1,2,1,NA), 
  ptcol="default", legend=TRUE, legendloc=NULL, minyr="default", 
  textmindev=0.5, relative=FALSE, expected=TRUE, estimated=TRUE, 
  bias_adjusted=TRUE, show_env=TRUE, virg=TRUE, init=TRUE, forecast=FALSE, 
  lang=c("f", "e"), ptypes="win", pngres=400, PIN=c(8, 7), outnam="stockRecruit")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- list object created by \code{'SS_output()'}}
  \item{subplot}{\code{numeric} -- vector of subplots to show;
    \code{1}=plot without labels, \code{2}=plot with year labels}
  \item{add}{\code{logical} -- if \code{TRUE}, add to existing plot}
  \item{plot}{\code{numeric} -- plot to active plot device}
  \item{print}{\code{numeric} -- print to PNG files}
  \item{xlim}{\code{numeric} -- optional control of x range}
  \item{ylim}{\code{numeric} -- optional control of y range}
  \item{labels}{\code{character} -- vector containing x-axis label for 
    models with spawning biomass in metric tons, y-axis label, and 
    alternative x-axis for models with a fecundity relationship 
    making spawning output not equal to spawning biomass}
  \item{bioscale}{\code{numeric} -- multiplier on spawning biomass;
    set to 0.5 for single-sex models}
  \item{plotdir}{\code{character} -- directory where PNG files will 
    be written; by default it will be the directory where the model was run}
  \item{pwidth}{\code{numeric} -- width of plots printed to files in
    units of \code{'punits'}; default is \code{pwidth=6.5}}
  \item{pheight}{\code{numeric} -- height of plots printed to png files
    in units of \code{'punits'}; default is designed to allow two plots
    per page, with \code{'pheight_tall'} used for plots that work best
    with a taller format and a single plot per page}
  \item{punits}{\code{character} -- units for \code{'pwidth'} and 
    \code{'pheight'}; can be \code{'px'} (pixels), \code{'in'} (inches),
    \code{'cm'} (centimeters), or \code{'mm'} (millimeters); 
    default is \code{punits='in'}}
  \item{res}{\code{numeric} -- resolution of plots printed to files;
    default is \code{res=400}}
  \item{ptsize}{\code{numeric} -- point size for plotted text in 
  plots printed to files (see \code{help('png')} in R for details)}
  \item{verbose}{\code{logical} -- if \code{TRUE}, report progress to R GUI}
  \item{colvec}{\code{character} -- vector of length 4 with colors for
    3 lines and 1 set of points (where the 4th value for the points is the
    color of the circle around the background color provided by \code{'ptcol'}}
  \item{ltyvec}{\code{numeric} -- vector of length 4 with line types for the
    3 lines and 1 set of points, where the points are disconnected (lty=NA) by default}
  \item{ptcol}{\code{character} -- vector or single value for the color
    of the points; \code{'default'} will by replaced by a vector of 
    colors of length equal to \code{nrow(replist[['recruit']])}}
  \item{legend}{\code{logical} -- if \code{TRUE}, add a legend to the figure}
  \item{legendloc}{\code{character} -- location of legend; by default
    it is chosen as the first value in the set of \code{'topleft'}, 
    \code{'topright'}, \code{'bottomright'} that results in no
    overlap with the points in the plot, but the user can override this
    with their choice of location; see ?legend for more info on the options}
  \item{minyr}{\code{numeric} -- minimum year of recruitment deviation to show in plot}
  \item{textmindev}{\code{numeric} -- minimum recruitment deviation for
    label to be added so only extreme devs are labeled (labels are added
    to first and last years as well); default=0.7}
  \item{relative}{\code{logical} -- if \code{TRUE}, scale both axes so
    that \code{B0} and \code{R0} are at 1 to show spawning output and
    recruitment relative to the equilibrium}
  \item{expected}{\code{logical} -- if \code{TRUE}, show line for 
    expected recruitment (stock-recruit curve)}
  \item{estimated}{s\code{logical} -- if \code{TRUE}, how points for
    estimated recruitment values (including deviations)}
  \item{bias_adjusted}{\code{logical} -- if \code{TRUE}, show lines
    for bias adjusted expected recruitment}
  \item{show_env}{\code{logical} -- if \code{TRUE}, add line for 
    expected recruitment with environmental variability}
  \item{virg}{\code{logical} -- if \code{TRUE}, add point for 
    equilibrium conditions (\code{x=B0}, \code{y=R0})}
  \item{init}{\code{logical} -- if \code{TRUE}, add point for initial
    conditions (\code{x=B1}, \code{y=R1}); only appears if this point
    differs from virgin values}
  \item{forecast}{\code{logical} -- if \code{TRUE}, include forecast years in the curve}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
  \item{ptypes}{\code{character} -- output types for plots (e.g., \code{'png'})}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) of the PNG file}
  \item{PIN}{\code{numeric} -- plot dimensions (width, height) for PNG file(s)}
  \item{outnam}{\code{character} -- output name for plot (and other data, if applicable)}
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
  \href{mailto:ian@iphc.int}{Ian Stewart}, International Pacific Halibut Commission, Seattle WA\cr
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-26}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotSpawnrecruit'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.comps}},
  \code{\link[PBSsynth]{plotSS.francis}},
  \code{\link[PBSsynth]{plotSS.index}},
  \code{\link[PBSsynth]{plotSS.pars}},
  \code{\link[PBSsynth]{plotSS.profile}},
  \code{\link[PBSsynth]{plotSS.rdevs}},
  \code{\link[PBSsynth]{plotSS.selex}},
  \code{\link[PBSsynth]{plotSS.stdres}},
  \code{\link[PBSsynth]{plotSS.ts}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotSpawnrecruit}}
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
