\name{plot_cres}
\alias{plot_cres}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot OSA Residuals
}
\description{
  Create a 5-panel plot of the one-step-ahead (OSA) residuals from an MLE | MPD.
}
\usage{
plot_cres(x, pick_one=NULL, maxLag=NULL, main="", 
   png=FALSE, pngres=400, PIN=c(10,9), lang=c("f","e"), 
   outnam=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{numeric} -- residual object as returned from one of the residual functions}
  \item{pick_one}{\code{numeric} -- number of the plot to provide if one of the 4 plot types printed by default needs to be extracted}
  \item{maxLag}{\code{numeric} -- maximum lag to use for ACF plot}
  \item{main}{\code{character} -- main title for the 5-panel figure}
  \item{png}{\code{logical} -- if \code{TRUE}, send the figure to a \code{PNG} device}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) of the \code{PNG} file}
  \item{PIN}{\code{numeric} -- width and height (inches) of the \code{PNG} file}
  \item{lang}{\code{character} -- language to use when constructing the figure 
    (either \code{'e'}=english or \code{'f'}=french)}
  \item{outnam}{\code{character} -- name of output file (without the suffix `\code{.png}')}
  \item{\dots}{\code{dots} -- extra arguments}
}
\details{
Adapted from \code{compResidual:::plot.cres}

Produces a 5-panel overview of the residuals by default or 
one specific plot if the \code{pick_one} argument is provided.

\itemize{
  \item Panel 1 (upper left) : bubble plot of the OSA residuals
    (red=negative, blue=positive; scale at top of plot)
  \item Panel 2 (lower left) : autocorrelation and cross-correlation plot for 
    all dimensions of the residual matrix (rows, columns and diagonal) with 95\% confidence limits
  \item Panel 3 (upper right) : Q-Q plot of the residuals, each color is a row (age)
    of the residual matrix (i.e., compositional group) 
  \item Panel 4 (middle right) : compositional plot of the residuals by column (year)
    with \ifelse{latex}{\out{$\pm$}}{\ifelse{html}{\out{&plusmn;}}{+/-}}2 standard
    deviations indicated by horizontal dashed lines and medians as dots
  \item Panel 5 (lower right) : compositional plot of the residuals by row (age)
    with \ifelse{latex}{\out{$\pm$}}{\ifelse{html}{\out{&plusmn;}}{+/-}}2 standard
    deviations indicated by horizontal dashed lines and medians as dots
}
}
\value{
	Nothing in particular -- no values returned.
}
\references{
Trijoulet, V. and Nielsen, A. (2022). 
compResidual: Residual calculation for compositional observations.
R package version 0.0.1, commit d4c74845089074d8016454f235044c5d13ded3a5,\cr
  <https://github.com/fishfollower/compResidual>.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-02-20}
}
\note{
Thank you Vanessa and Anders.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.stdres}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{plotFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
