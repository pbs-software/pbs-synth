\name{plotSS.pars}
\alias{plotSS.pars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Parameter Distributions
}
\description{
  Make multi-figure plots of parameter distributions: priors, MPD/MLE estimates,
  and MCMC posteriors.
  Function can work without MCMC analysis.
}
\usage{
plotSS.pars(replist, plotdir=NULL, xlab="Parameter value", 
  ylab="Density", showmle=TRUE, showpost=TRUE, showprior=TRUE, 
  showinit=TRUE, showdev=FALSE, showlegend=TRUE, fitrange=FALSE, 
  fitnudge=0, xaxs="i", xlim=NULL, ylim=NULL, verbose=TRUE, 
  debug=FALSE, nrows=3, ncols=3, ltyvec=c(1,1,3,4), 
  colvec=c("blue","red","black","green", rgb(0,0,0,0.5)), 
  add=FALSE, plot=TRUE, print=FALSE, punits="in", ptsize=10, 
  strings=NULL, exact=FALSE, newheaders=NULL, 
  outnam, res=400, PIN=c(8, 8), lang="e")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- list object created by \code{'SS_output()'}}
  \item{plotdir}{\code{character} -- path to the folder where the plots 
    will be saved; default is \code{NULL}, which leads to the plots being 
    created in the folder that contains the results}
  \item{xlab}{\code{character} -- label on horizontal axis }
  \item{ylab}{\code{character} -- label on vertical axis}
  \item{showmle}{\code{logical} -- if \code{TRUE}, show MLE estimate and 
    asymptotic variance estimate with blue lines}
  \item{showpost}{\code{logical} -- if \code{TRUE}, show posterior distribution
    as bar graph if MCMC results are available in \code{'replist'}}
  \item{showprior}{\code{logical} -- if \code{TRUE}, show prior distribution as black line}
  \item{showinit}{\code{logical} -- if \code{TRUE}, show initial value as red triangle}
  \item{showdev}{\code{logical} -- if \code{TRUE}, include devs in the plot}
  \item{showlegend}{\code{logical} -- if \code{TRUE}, show the legend}
  \item{fitrange}{\code{logical} -- if \code{TRUE}, fit range tightly around MLE 
    and posterior distributions instead of using full parameter range}
  \item{fitnudge}{\code{numeric} -- nudge factor to fit around central tendency
    (\emph{RH foofaraw})}
  \item{xaxs}{\code{character} -- parameter input for x-axis; type \code{'?par'} for more info}
  \item{xlim}{\code{numeric} -- optional x-axis limits to be applied to all plots;
    otherwise, limits are based on the model results}
  \item{ylim}{\code{numeric} -- optional y-axis limits to be applied to all plots;
    otherwise, limits are based on the model results}
  \item{verbose}{\code{logical} -- if \code{TRUE}, send output text (if any) to command line}
  \item{debug}{\code{logical} -- if \code{TRUE}, provide additional messages
    to help with debugging when the function fails}
  \item{nrows}{\code{numeric} -- number of rows in multi-figure plot}
  \item{ncols}{\code{numeric} -- number of columns in multi-figure plot}
  \item{ltyvec}{\code{numeric} -- vector of line types used for lines showing MLE
    and prior distributions and the median of the posterior distribution}
  \item{colvec}{\code{character} -- vector of colors used for lines and polygons
    showing MLE, initial value, prior, posterior, and median of the posterior}
  \item{add}{\code{logical} -- if \code{TRUE}, add to existing plot}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, print to PNG files}
  \item{punits}{\code{character} -- units for \code{'PIN'}; can be \code{'px'}
    (pixels), \code{'in'} (inches), \code{'cm'} or \code{'mm'}; default=\code{'in'}}
  \item{ptsize}{\code{numeric} -- point size for plotted text in plots printed to files
    (see \code{help('png')} in R for details); default=12}
  \item{strings}{\code{character} -- subset parameters included in the plot using
    substring from parameter names (i.e. \code{'SR'} will get \code{'SR_LN(R0)'}
    and \code{'SR_steep'} if they are both estimated quantities in this model)}
  \item{exact}{\code{logical} -- if \code{TRUE}, \code{'strings'} match input 
    parameter names exactly; otherwise substrings are allowed}
  \item{newheaders}{\code{character} -- optional vector of headers for each panel
    to replace the parameter names}
  \item{outnam}{\code{character} -- output name for plot (and other data, if applicable)}
  \item{res}{\code{numeric} -- resolution of plots printed to files;
    default is \code{res=400}}
  \item{PIN}{\code{numeric} -- plot dimensions (width, height) for PNG file(s)}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
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
  \href{mailto:cole.monnahan@noaa.gov}{Cole Monnahan}, Alaska Fisheries Science Center, NOAA, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-23}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotPars'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.comps}},
  \code{\link[PBSsynth]{plotSS.francis}},
  \code{\link[PBSsynth]{plotSS.index}},
  \code{\link[PBSsynth]{plotSS.profile}},
  \code{\link[PBSsynth]{plotSS.rdevs}},
  \code{\link[PBSsynth]{plotSS.selex}},
  \code{\link[PBSsynth]{plotSS.stdres}},
  \code{\link[PBSsynth]{plotSS.stock.recruit}},
  \code{\link[PBSsynth]{plotSS.ts}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotPars}}
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
