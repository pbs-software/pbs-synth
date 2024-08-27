\name{calcMA}
\alias{calcMA}
\title{
  Calculate Mean Age
}
\description{
  Calculate various measures of age, including observed and expected mean age,
  by year for fleets with AF data.
}
\usage{
calcMA(runs=1, rwts=0:1, vers=1, overwrite=FALSE, fleets.lab, fleets.af, 
   cwd="C:/Users/haighr/Files/GFish/PSARC24/YTR/Data/SS3/YTR2024",
   use.cwd.as.is=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{runs}{\code{numeric} -- model run number(s)}
  \item{rwts}{\code{numeric} -- model reweight number(s)}
  \item{vers}{\code{numeric} -- model run version}
  \item{overwrite}{\code{logical} -- if {TRUE}, overwrite existing csv file(s)}
  \item{fleets.lab}{\code{character} -- fleet labels}
  \item{fleets.af}{\code{numeric} -- fleet numbers of fleets with AF data}
  \item{cwd}{\code{character} -- currrent working directory where model runs are located;
    usually with a base name comprising three-letter species code and assessment year}
  \item{use.cwd.as.is}{\code{logical} -- if {TRUE}, use the \code{'cwd'} value as is
    (user may specify an MPD run explicitly that does not conform to the specified run, rwt, and ver)}
}
\details{
  Taps into the Francis (2011) routine TA1.8 by calling the function
  \code{'plotSS.francis'}, which is borrowed and modified from the function
  \code{'r4ss::SSMethod.TA1.8'}.
}
\value{
  A data frame (also sent to a CSV file), with rows specifying fleet and calendar year, and columns:
  \tabular{ll}{
    \code{Obsmn}    \tab observed mean age\cr
    \code{Obslo}    \tab observed minimum age\cr
    \code{Obshi}    \tab observed maximum age\cr
    \code{semn}     \tab standard error of the mean age\cr
    \code{Expmn}    \tab expected (model-derived) mean age\cr
    \code{Vexp}     \tab variance of the expected mean age\cr
    \code{N}        \tab number of observations per year\cr
    \code{Std.res}  \tab standardised residual of mean age\cr
    \code{ObsloAdj} \tab adjusted observed minimum age\cr
    \code{ObshiAdj} \tab adjusted observed maximum age\cr
    \code{Fleet}    \tab SS3 fleet number\cr
    \code{Yr}       \tab model year (usually calendar)\cr
  }
}
\references{
  Francis, R.I.C.C. 2011. 
  \href{https://doi.org/10.1139/f2011-025}{Data weighting in statistical fisheries stock assessment models}.
  \emph{Can. J. Fish. Aquat. Sci.} \bold{68(6)}. 1124--1138.}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-19}
}
\note{
  Years ago, \code{'r4ss::SSMethod.TA1.8'} appeared to have a bug when processing multiple fleets.
  This was fixed in \code{'plotSS.francis'}.
  The issue might have been fixed in the most recent version of \code{r4ss}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.francis}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSMethod.TA1.8}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
