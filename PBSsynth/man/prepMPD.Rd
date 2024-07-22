\name{prepMPD}
\alias{prepMPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prepare MPD Run
}
\description{
  Prepare a new MPD run based on a previous run.
}
\usage{
prepMPD(run.rwt.ver,
   d.base="C:/Users/haighr/Files/GFish/PSARC24/YTR/Data/SS3/YTR2024",
   w=NULL, cvpro=NULL, modify=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{run.rwt.ver}{\code{numeric} -- six-element vector of:\cr
    new run, new reweight, new version, old run, old reweight, old version}
  \item{d.base}{\code{character} -- directory of model runs for current stock assessment}
  \item{w}{\code{numeric} -- vector of AF weights to apply to AF sample size, named by fleet number}
  \item{cvpro}{\code{numeric} -- vector of process CVs to add to observed CVs, named by fleet number}
  \item{modify}{\code{logical} -- if \code{TRUE}, prompts user to modify the copied \code{'.ss'} files;
    if \code{FALSE}, copy model input and output of previous run [NOT fully tested yet]}
}
\details{
  Essentially copies the four \code{'.ss'} input files from a previous run
  to a newly created run, and sets to the new working directory.
  User must modify one or more input files in the new run
  before calling the \code{'SS'} executable.
}
\value{
  Four copied \code{'.ss'} input files in a new model run directory.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-11}
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.francis}}, 
  \code{\link[PBSsynth]{runSweave}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{IO}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
