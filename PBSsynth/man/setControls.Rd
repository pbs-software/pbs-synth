\name{setControls}
\alias{setControls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set Controls for Model Results Appendix
}
\description{
  Set various default values for the stock assessment's model results appendix.
}
\usage{
setControls(fsa="SGR2025", stospp=c("SGR","SGR"), stolab=c("3area","coast"), 
   assyr=2025, appG.dir, clean=FALSE, lang="e")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fsa}{\code{character} -- fisheries stock assessment label}
  \item{stospp}{\code{character} -- three-letter species code for stock under consideration}
  \item{stolab}{\code{character} -- either a model label or a substock label}
  \item{assyr}{\code{numeric} -- major assessment year}
  \item{appG.dir}{\code{character} -- full system path to the directory containing model runs}
  \item{clean}{\code{logical} -- if \code{TRUE}, clean the user's work environment before setting controls}
  \item{lang}{\code{character} -- language identifier for plots (either \code{'e'} or \code{'f'})}
}
\details{
  The primary task of the function is to populate a list called
  \code{'stock'}, which contains information specific to each stock,
  e.g., model runs, model labels, etc.
}
\value{
List object containing:
  \item{saveDate}{\code{character} -- date the list object was created}
  \item{stospp}{\code{character} -- three-letter species code for each stock}
  \item{stolab}{\code{character} -- label used to identify each stock (or model variant)}
  \item{assyr}{\code{numeric} -- year when the most of the models were run}
  \item{appG.dir}{\code{character} -- full system directory where the model results reside}
  \item{Nstock}{\code{numeric} --number of stocks (or model variants)}
  \item{runA}{\code{numeric} -- array of model runs comprising axis of uncertainty for first stock}
  \item{runB}{\code{numeric} -- array of model runs comprising axis of uncertainty for second stock}
  \item{A.base}{\code{numeric} -- specifics (run, reweight, version) for first axis of uncertainty}
  \item{B.base}{\code{numeric} -- specifics (run, reweight, version) for second axis of uncertainty}
  \item{stock}{\code{list} -- main set of controls for building model results appendix}
  \item{strSpp}{\code{character} -- Hart code for assessed species}
  \item{spp.code}{\code{character} -- 3-letter code for assessed species}
  \item{spp.name}{\code{character} -- common (proper) name of assessed species}
  \item{latin.name}{\code{character} -- latin name of assessed species}
  \item{spp.year}{\code{character} --  assessment name (species and year)}
  \item{the.stocks}{\code{character} -- stock names}
  \item{ptype}{\code{character} -- plot type, usually \code{'png'}}
  \item{pngres}{\code{numeric} -- figure file output resolution (pixels/inch)}
  \item{sigdig}{\code{numeric} -- number of significant digits for output tables}
  \item{decdig}{\code{numeric} -- number of decimal places for output tables}
  \item{quants3}{\code{numeric} -- three quantile probabilities (median and 90\% CI)}
  \item{ quants5}{\code{numeric} -- five quantile probabilities (median, 50\% and 90\% CI)}
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-07-24}
}
\note{
  Function \code{setControls} was created from R code called \code{'set.controls'},
  which was modified and sourced for each fisheries stock assessment (FSA).
  In the Silvergray Rockfish (SGR) FSA 2026, the model results appendix for 
  Offshore Rockfish was changed from Appendix F to Appendix G
  (Appendix E: Equations \ifelse{latex}{\out{$\rightarrow$}}{\ifelse{html}{\out{&rArr;}}{->}} Ecosystem;
  Appendix F: Model Results \ifelse{latex}{\out{$\rightarrow$}}{\ifelse{html}{\out{&rArr;}}{->}} Model Platform).

  The first use of \code{'set.controls'} occurred in 2020 for the
  Rougheye/Blackspotted Rockfish (REBS) FSA.
  Prior to 2020, the \code{'stock'} object was built by Sweave files.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{initStock}}, etc.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
