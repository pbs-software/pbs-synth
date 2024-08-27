\name{plotSS.rdevs}
\alias{plotSS.rdevs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Recruitment Deviations
}
\description{
  Plot recruitment deviations and associated quantities including
  derived measures related to bias adjustment.
}
\usage{
plotSS.rdevs(replist, subplots=1:3, plot=TRUE, print=FALSE, add=FALSE, 
  uncertainty=TRUE, minyr=-Inf, maxyr=Inf, forecastplot=FALSE, 
  col1="black", col2="blue", col3="green3", col4="red", 
  legendloc="topleft", labels=c("Year", 
  "Asymptotic standard error estimate", "Log recruitment deviation", 
  "Bias adjustment fraction, 1 - stddev^2 / sigmaR^2"), 
  pwidth=8, pheight=6, punits="in", res=400, ptsize=10, cex.main=1, 
  plotdir="default", verbose=TRUE, outnam, lang="e")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- list object created by \code{'SS_output()'}}
  \item{subplots}{\code{numeric} -- vector controlling which subplots to create}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active plot device}
  \item{print}{\code{logical} -- if \code{TRUE}, print to PNG files}
  \item{add}{\code{logical} -- if \code{TRUE}, add to existing plot
    (\emph{not yet implemented})}
  \item{uncertainty}{\code{logical} -- if \code{TRUE}, include plots showing uncertainty}
  \item{minyr}{\code{numeric} -- optional input for minimum year to show in plots}
  \item{maxyr}{\code{numeric} -- optional input for maximum year to show in plots}
  \item{forecastplot}{\code{logical} -- if \code{TRUE}, include points from forecast years}
  \item{col1}{\code{character} -- first color used}
  \item{col2}{\code{character} -- second color used}
  \item{col3}{\code{character} -- third color used}
  \item{col4}{\code{character} -- fourth color used}
  \item{legendloc}{\code{character} -- location of legend;
    see ?legend for more info}
  \item{labels}{\code{character} -- vector of labels for plots
    (titles and axis labels)}
  \item{pwidth}{\code{numeric} -- width of plot}
  \item{pheight}{\code{numeric} -- height of plot}
  \item{punits}{\code{character} -- units for PNG file}
  \item{res}{\code{numeric} -- resolution of plots printed to files;
    default is \code{res=400}}
  \item{ptsize}{\code{numeric} -- point size for PNG file}
  \item{cex.main}{\code{numeric} -- character expansion for plot titles}
  \item{plotdir}{\code{character} -- directory where PNG files will be
    written; by default it will be the directory where the model was run}
  \item{verbose}{\code{logical} -- if \code{TRUE}, report progress to R GUI}
  \item{outnam}{\code{character} -- output name for plot (and other data, if applicable)}
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
  \href{mailto:ian@iphc.int}{Ian Stewart}, International Pacific Halibut Commission, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-23}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotRecdevs'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.comps}},
  \code{\link[PBSsynth]{plotSS.francis}},
  \code{\link[PBSsynth]{plotSS.index}},
  \code{\link[PBSsynth]{plotSS.pars}},
  \code{\link[PBSsynth]{plotSS.profile}},
  \code{\link[PBSsynth]{plotSS.selex}},
  \code{\link[PBSsynth]{plotSS.stdres}},
  \code{\link[PBSsynth]{plotSS.stock.recruit}},
  \code{\link[PBSsynth]{plotSS.ts}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotRecdevs}}
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
