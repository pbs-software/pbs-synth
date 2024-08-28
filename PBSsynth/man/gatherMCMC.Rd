\name{gatherMCMC}
\alias{gatherMCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Gather MCMCs from Model Runs
}
\description{
  Gather MCMC samples from a variety of model runs (e.g., composite base runs, sensitivity runs).
}
\usage{
gatherMCMC(mcdir=".", type="compo", strSpp="YTR",
   basedir="C:/Users/haighr/Files/GFish/PSARC24/YTR/Data/SS3/YTR2024",
   ryrs=1935:2025, pyrs=2026:2035, valTS=c("SSB", "F", "Recr", "RecrDev"),
   valRP=c("SSB_MSY", "annF_MSY", "Dead_Catch_MSY"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcdir}{\code{character} -- name of directory containing MCMC results}
  \item{type}{\code{character} -- type of runs (\code{'compo'}, \code{'senso'}, etc.)}
  \item{strSpp}{\code{character} -- three-letter code of species}
  \item{basedir}{\code{character} -- directory where SS3 runs are stored for species \code{strSpp}}
  \item{ryrs}{\code{numeric} -- reconstruction years for the population model}
  \item{pyrs}{\code{numeric} -- projection years for the population model}
  \item{valTS}{\code{character} -- Array dimension names for time series (TS) output}
  \item{valRP}{\code{character} -- Array dimension names for reference point (RP) output}
}
\details{
  Collects a variety of arrays for use later, 
  primarily in the processing of Appendix F (Model Results).
}
\value{
A list of arrays are produced for single-area models:
\tabular{rcl}{
  \code{ampdPA} \tab Dim \tab MPD parameter estimates for each run:\cr
                \tab   1 \tab run designator (run.rwt), where rwt=reweight number\cr
                \tab   2 \tab parameter names\cr
  \code{avgLL}  \tab  -- \tab log likelihood (LL) values for each run:\cr
                \tab   1 \tab log-likelihood components\cr
                \tab   2 \tab run designator (run.rwt)\cr
  \code{avgPA}  \tab  -- \tab MCMC parameter estimates by sample:\cr
                \tab   1 \tab sample designator (run.rwt.sample)\cr
                \tab   2 \tab parameter names\cr
  \code{avgRP}  \tab  -- \tab MCMC reference points (derived parameters) by sample:\cr
                \tab   1 \tab sample designator (run.rwt.sample)\cr
                \tab   2 \tab reference point names\cr
  \code{avgTS}  \tab  -- \tab MCMC time series (reconstruction) by sample:\cr
                \tab   1 \tab sample designator (run.rwt.sample)\cr
                \tab   2 \tab model reconstruction year\cr
                \tab   3 \tab time series names\cr
  \code{avgPJ}  \tab  -- \tab MCMC projection series by sample:\cr
                \tab   1 \tab sample designator (run.rwt.sample)\cr
                \tab   2 \tab model projection year\cr
                \tab   3 \tab projection series names\cr
                \tab   4 \tab catch policy\cr
  \code{avgCP}  \tab  -- \tab MCMC projected catch policies:\cr
                \tab   1 \tab model projection year\cr
                \tab   2 \tab run designator (run.rwt)\cr
                \tab   3 \tab catch policy\cr
}
An additional list of arrays are produced for multi-area models (prefixed by \code{'x'} for extra):
\tabular{rcl}{
  \code{xavgRP} \tab Dim \tab MCMC reference points (derived parameters) by sample and subarea:\cr
                \tab   1 \tab sample designator (run.rwt.sample)\cr
                \tab   2 \tab reference point names\cr
                \tab   3 \tab subarea names\cr
  \code{xavgTS} \tab  -- \tab MCMC time series (reconstruction) by sample and subarea:\cr
                \tab   1 \tab sample designator (run.rwt.sample)\cr
                \tab   2 \tab model reconstruction year\cr
                \tab   3 \tab time series names\cr
                \tab   4 \tab subarea names\cr
  \code{xavgPJ} \tab  -- \tab MCMC projection series by sample and subarea:\cr
                \tab   1 \tab sample designator (run.rwt.sample)\cr
                \tab   2 \tab model projection year\cr
                \tab   3 \tab projection series names\cr
                \tab   4 \tab subarea names\cr
                \tab   5 \tab catch policy\cr
  \code{xavgCP} \tab  -- \tab MCMC projected catch policies by subarea:\cr
                \tab   1 \tab model projection year\cr
                \tab   2 \tab subarea names\cr
                \tab   3 \tab catch policy\cr
}
  Various output files of MCMC results are characterised by their name:\cr
  \code{'compo'} -- base run (or composite base case)\cr
  \code{'senso'} -- sensitivity runs\cr
  \code{'penso'} -- single-area model MCMCs if main model (\code{'compo'}) is multi-area\cr
  \code{'renso'} -- low-recruitment MCMC samples only
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-16}
}
\note{
  Perhaps output arrays should have a consistent order to their dimensions,
  e.g., (1) run/sample, (2) metric name, (3) year, (4) catch policy, (5) subarea;
  however, changes now will likely break code.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{load_extra_mcmc}},
  \code{\link[PBSsynth]{plotSS.compo}},
  \code{\link[PBSsynth]{plotSS.senso}},
  \code{\link[PBSsynth]{tabSS.compo}},
  \code{\link[PBSsynth]{tabSS.decision}},
  \code{\link[PBSsynth]{tabSS.senso}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{IO}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{appFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
