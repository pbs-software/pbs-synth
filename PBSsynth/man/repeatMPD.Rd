\name{repeatMPD}
\alias{repeatMPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Repeat MPD Runs
}
\description{
  Repeat MPD runs for likelihood profile analysis.
}
\usage{
repeatMPD(Pfix = list(FFF=seq(0.05, 0.06, 0.01), 
   MMM=seq(0.05, 0.06, 0.01), RRR=NULL, SSS=NULL), A=60, 
   dir.mpd=getwd(), dir.par, prefix="control.MMM.", 
   clean=FALSE, strSpp="CAR", argsMPD="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pfix}{\code{numeric} -- list of sequential fixed values for model parameters:
    \itemize{
      \item \code{'FFF'} female natural mortality (\eqn{M_1}{M female})
      \item \code{'MMM'} male natural mortality (\eqn{M_2}{M male})
      \item \code{'RRR'} log intial recruitment (\eqn{\log R_0}{log R0})
      \item \code{'SSS'} standard deviation of recruitment deviations (\eqn{\sigma_R}{sigmaR})
  }}
  \item{A}{\code{numeric} -- maximum age in model (plus group)}
  \item{dir.mpd}{\code{character} -- directory containing MPD results of model run}
  \item{dir.par}{\code{character} -- subdirectory where parameter runs 
    will be placed (e.g., \code{'LL.Mf'})}
  \item{prefix}{\code{character} -- prefix for \code{'control.ss'} file names}
  \item{clean}{\code{logical} -- if \code{TRUE}, clean model output from previous runs}
  \item{strSpp}{\code{character} -- three- or four-letter code for rockfish species}
  \item{argsMPD}{\code{character} -- additional arguments for the likelihood model runs}
}
\details{
  Creates a subdirectory called \code{'NLL.<run>.<rwt>.<ver>'} under 
  the model run directory (parallel with the MPD and MCMC model runs). 
  Under the \code{'NLL'} directory, subdirectories are created (depending
  on the value of \code{'dir.par'}) for whichever parameter is held constant.
  In the parameter subdirectory, MPD runs are executed for each 
  incremental fixed value of the parameter.

  Once all model runs for the parameter have finished, log likelihood
  profiles can be visualised using the function \code{'plotSS.profile'}.
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
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-26}
}
\note{
  Order of function exection:\cr
  \code{'repeatMPD'} \eqn{\rightarrow}{-->} \code{'SSgetoutput'}
  \eqn{\rightarrow}{-->} \code{'SSsummarize'}  \eqn{\rightarrow}{-->}
  \code{'plotSS.profile'}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.profile}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{profile}},
  \code{\link[r4ss]{SSgetoutput}},
  \code{\link[r4ss]{SSsummarize}},
  \code{\link[r4ss]{SSplotProfile}},
  \code{\link[r4ss]{PinerPlot}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{model}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
