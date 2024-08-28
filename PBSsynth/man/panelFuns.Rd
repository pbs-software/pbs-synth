\name{panelFuns}
\alias{panelBoxes}
\alias{panelChains}
\alias{panelTraces}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Something in Panels
}
\description{
  Plot the following items from MCMC samples for base or sensitivity runs:\cr
  \code{Boxes} -- quantile boxes;\cr
  \code{Chains} -- cumulative frequency of \code{n} chains by partitioning one set of samples;\cr
  \code{Traces} -- traces (set of samples).
}
\usage{
panelBoxes(mcmc, nchains=9, pdisc=0, axes=FALSE, same.limits=FALSE,
   between=list(x=axes,y=axes), div=1, log=FALSE, base=10, main=NULL,
   xlab=NULL, ylab=NULL, xlim=c(0.25,9.75), 
   boxfill=paste0(rep(c("cyan","green","coral"),each=3),rep(1:3,3)),
   cex.main=1.2, cex.lab=1.2, cex.strip=0.9, cex.axis=0.9, las=0,
   tck=0.4, tick.number=5, xfac=paste0("B",1:nchains), outline=TRUE, 
   lang="e", ...)

panelChains(mcmc, nchains=3, pdisc=0.1, axes=FALSE, same.limits=FALSE,
   between=list(x=axes,y=axes), div=1, span=1/4, log=FALSE,
   base=10, main=NULL, xlab=NULL, ylab=NULL, cex.main=1.2,
   cex.lab=1, cex.strip=0.8, cex.axis=0.8,
   las=0, tck=0.4, tick.number=5, lty.trace=1, lwd.trace=1,
   col.trace="grey", lty.median=1, lwd.median=1,
   col.median="black", lty.quant=2, lwd.quant=1,
   col.quant="black", plot=TRUE, probs=tcall(quants3), lang="e", ...)

panelTraces(mcmc, mpd=mcmc[1,], nchains=1, pdisc=0, axes=FALSE, 
   same.limits=FALSE, between=list(x=axes,y=axes), div=1, log=FALSE,
   base=10, main=NULL, xlab=NULL, ylab=NULL, cex.main=1.2, cex.lab=1.2,
   cex.strip=0.9, cex.axis=0.9, las=0, tck=0.4, tick.number=5,
   xfac=NULL, lang="e", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc}{\code{data.frame} -- MCMC samples (by row) of parameters or derived quantities (by column)}
  \item{nchains}{\code{numeric} -- number of chains to create from one trace of MCMC samples}
  \item{pdisc}{\code{numeric} -- proportion of the initial trace to discard before creating chains}
  \item{axes}{\code{logical} -- if \code{TRUE}, plot axis values}
  \item{same.limits}{\code{logical} -- if \code{TRUE}, panels should have same x-axis limits}
  \item{between}{\code{deprecate} -- list with \code{x} and \code{y} indicating panel spacing (\emph{not used})}
  \item{div}{\code{numeric} -- denominator to shorten values on the y axis}
  \item{log}{\code{logical} -- if \code{TRUE}, log-transform values}
  \item{base}{\code{numeric} -- logarithm base}
  \item{main}{\code{character} -- main title}
  \item{xlab}{\code{character} -- x-axis title}
  \item{ylab}{\code{character} -- y-axis title}
  \item{xlim}{\code{numeric} -- limits of the x-axis (mainly to add outer spacing to the outer boxes)}
  \item{boxfill}{\code{character} -- vector of colours to fill the quantile boxes}
  \item{cex.main}{\code{numeric} -- size of main title}
  \item{cex.lab}{\code{numeric} -- size of axis labels}
  \item{cex.strip}{\code{numeric} -- size of strip labels}
  \item{cex.axis}{\code{numeric} -- size of tick labels}
  \item{las}{\code{numeric} -- orientation of tick labels: 0=parallel, 1=horizontal,
    2=perpendicular, 3=vertical}
  \item{tck}{\code{numeric} -- tick mark length}
  \item{tick.number}{\code{numeric} -- number of tick marks}
  \item{xfac}{\code{character} -- vector of labels for x-axis ticks (under quantile boxes)}
  \item{outline}{\code{logical} -- if \code{TRUE}, add outliers to the quantile plots}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french)}
  \item{\dots}{\code{dots} -- passed to panel's subfunction}
%% Specific to plotChains
  \item{span}{\code{numeric} -- parameter which controls the degree of smoothing}
  \item{lty.trace}{\code{numeric|character} -- line type of trace}
  \item{lwd.trace}{\code{numeric} -- line width of trace}
  \item{col.trace}{\code{numeric|character} -- colour of trace}
  \item{lty.median}{\code{numeric|character} -- line type of median}
  \item{lwd.median}{\code{numeric} -- line width of median}
  \item{col.median}{\code{numeric|character} -- colour of median}
  \item{lty.quant}{\code{numeric|character} -- line type of quantile trace}
  \item{lwd.quant}{\code{numeric} -- line width of quantile trace}
  \item{col.quant}{\code{numeric|character} -- colour of quantile trace}
  \item{plot}{\code{logical} -- if \code{TRUE}, draw the plot}
  \item{probs}{\code{numeric} -- quantile values for quantile trace}
%% Specific to plotTraces
  \item{mpd}{\code{numeric} -- MPD estimate, defaults to first sample in MCMC set (i.e., search strats at the MPD)}
}
\details{
  Primarily used as diagnostics to visualise how parameters or 
  derived quantities compare across base runs or sensitivity runs
  when building the \sQuote{Model Results} appendix.
}
\value{
  Invisibly returns the data used by the function \code{'mochaLatte'} to plot the panels.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-16}
}
\note{
  The interaction between panel functions and the function \code{'mochaLatte'} can be confusing.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{mochaLatte}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\concept{plotFuns.r}


