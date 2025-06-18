\name{calcRhat}
\alias{calcRhat}
\title{
  Calculate Split R-hat Statistic
}
\description{
  Calculate scale-reduction split-\eqn{\widehat{R}}{Rhat} statistic
  and Effective Sample Size (ESS) for MCMC model parameters.
}
\usage{
calcRhat(dir=".", nchains=8, parpos, rhat.only=FALSE, only.bad=FALSE, 
   badhat=1.01, offset=0.0025, recdevs=FALSE, zero.bar=TRUE, 
   barcols=.colBlind, xlim=NULL, png=FALSE, pngres=400, PIN=c(10,7.5), 
   lang="e", outnam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{\code{character} -- file directory or folder}
  \item{nchains}{\code{numeric} -- number of chains used in the MCMC analysis}
  \item{parpos}{\code{numeric} -- vector of parameter positions (columns) in the data object}
  \item{rhat.only}{\code{logical} -- if \code{TRUE}, only plot R-hat figure(s)}
  \item{only.bad}{\code{logical} -- if \code{TRUE}, only show parameters that don't meet the criteria}
  \item{badhat}{\code{numeric} -- upper threshold of R-hat for parameter acceptance}
  \item{offset}{\code{numeric} -- offset plot buffer for R-hat values falling below 1}
  \item{recdevs}{\code{logical} -- if \code{TRUE}, plot the R-hats for the recruitment deviation parameters}
  \item{zero.bar}{\code{logical} -- if \code{TRUE}, plot R-hat lines starting at zero}
  \item{barcols}{\code{character} -- colours for histogram bars stacked by chain}
  \item{xlim}{\code{numeric} -- limits of the x-axis}
  \item{png}{\code{logical} -- if \code{TRUE}, send plot to a PNG file}
  \item{pngres}{\code{numeric} -- resolution of PNG file (pixels/inch)}
  \item{PIN}{\code{numeric} -- dimensions (width, height) of PNG file (inches)}
  \item{lang}{\code{character} -- language of text in PNG file ('e'=english, 'f'=french)}
  \item{outnam}{\code{character} -- output name of PNG file}
}
\details{
  Creates various plots:
  (i) \eqn{\widehat{R}}{R-hat} + ESS for main parameters;
  (ii) \eqn{\widehat{R}}{R-hat} for recruitment deviations;
  (iii) stacked histogram of parameter posteriors, colour-coded by chain;
  (iv) Vehtari \eqn{\widehat{R}}{R-hat} vs. Langley \eqn{\widehat{R}}{R-hat} (simple).
}
\value{
  Nothing in particular.
}
\references{
Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and B\ifelse{latex}{\out{\"{u}}}{\ifelse{html}{\out{&uuml;}}{u}}rkner, P.C. (2021)
Rank-normalization, folding, and localization: an improved \eqn{\widehat{R}}{R-hat} for assessing convergence of MCMC (with discussion). 
\emph{Bayesian Analysis} \bold{16(2)}. 667--718.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-18}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{panelChains}}
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
