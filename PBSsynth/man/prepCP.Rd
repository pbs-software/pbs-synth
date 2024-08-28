\name{prepCP}
\alias{prepCP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prepare Catch Policies
}
\description{
  Prepare catch policies by evaluating fitted models using a selection
  of catch policies.
}
\usage{
prepCP(run.rwt, ver="", cp=list(BC=4000), d.cp="CC", tag="", 
   fyrs=2025:2035, season=1, fleet=1, 
   d.base="C:/Users/haighr/Files/GFish/PSARC24/YTR/Data/SS3/YTR2024", 
   w=NULL, cvpro=NULL, linux=TRUE, onlinux=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{run.rwt}{\code{numeric} -- vector of length two giving run and reweight numbers}
  \item{ver}{\code{character} -- version of \code{'run.rwt'} adopted for the base run}
  \item{cp}{\code{numeric} -- list of catch policies (\code{'CC'}) for
    each region or subarea; if stock is coastwide, list will be length one
    using name \code{'BC'} (or \code{'CST'} in past assessments)}
  \item{d.cp}{\code{character} -- subdirectory name for the set of catch policies;
    use \code{'CC'} for constant catch and \code{'HR'} for harvest rate policies
    (\code{'HR'} \emph{policies not implemented yet})}
  \item{tag}{\code{charaacter} -- optional tag string in case MCMC model directory 
    name has been qualified by some arbitrary text (e.g., \code{'.nuts4K'})}
  \item{fyrs}{\code{numeric} -- vector of forecast years}
  \item{season}{\code{numeric} -- number of seasons within a year
    (\emph{offshore rockfish assessments never use more than one season})}
  \item{fleet}{\code{numeric} -- vector of commercial fleets that catch target species}
  \item{d.base}{\code{character} -- base directory name where models 
    runs are located for the target species}
  \item{w}{\code{numeric} -- vector of weights to apply to the composition data
    (\emph{currently not used})}
  \item{cvpro}{\code{numeric} -- vector of process errors to add to the abundance data
    (\emph{currently not used})}
  \item{linux}{\code{logical} -- if \code{TRUE}, MCMCs were created on 
    Linux and transferred to Windows}
  \item{onlinux}{\code{logical} -- if \code{TRUE}, the R code is to be
    run on the Linux server}
}
\details{
  The code creates directories under the specified model run (usually
  the base run) to accomomdate the various catch policies, copies the
  relevant model run files to each directory, and isses the command:\cr
  \code{'/bin/ss3/build/ss3 -mceval'} on Linux, or\cr
  \code{'ss -mceval'} on Windows (assumes user has \code{'ss.exe'} on the path).
}
\value{
  Nothing in particular
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
  Much quicker to run \code{'-mceval'} on Linux server.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{agileDT}},
  \code{\link[PBSsynth]{findTarget}},
  \code{\link[PBSsynth]{gatherMCMC}},
  \code{\link[PBSsynth]{load_extra_mcmc}},
  \code{\link[PBSsynth]{tabSS.decision}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{run}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{datagen}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
