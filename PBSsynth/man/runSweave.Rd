\name{runSweave}
\alias{runSweave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Run Sweave
}
\description{
  Run Sweave code to collect SS3 model run output and create a PDF summary.
}
\usage{
runSweave(d.model=getwd(), d.sweave, type="MPD", figs.only=FALSE, debug=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d.model}{\code{character} -- directory containing SS3 model run output}
  \item{d.sweave}{\code{character} -- directory where master Sweave code is stored}
  \item{type}{\code{character} -- type of model output: \code{'MPD'} or \code{'MCMC'}}
  \item{figs.only}{\code{logical} -- if \code{TRUE}, select Sweave code that
    only collates the figures into a PDF file
    (\emph{has not been supported for some time, may or may not work})}
  \item{debug}{\code{logical} -- if \code{TRUE}, print various values 
    to the R console for debugging purposes}
}
\details{
  Creates a copy of \code{'sweaveMPD.Rnw'} or \code{'sweaveMCMC.Rnw'},
  depending on \code{'type'}, in the \code{'d.model'} directory, 
  and runs the Sweave code using R's \code{'Sweave'} function.
  The code will collect and process model output, create figures, and 
  compile tables and figures into a summary PDF file used for model
  run evaluation.
}
\value{
  Nothing in particular; however, the code saves a number of binary 
  output files for later use by Appendix F (model results) Sweave code.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-02-20}
}
\note{
  Used exclusivley to evaluate model runs.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{agileDT}},
  \code{\link[PBSsynth]{importCor}},
  \code{\link[PBSsynth]{makeFSARfigs}},
  \code{\link[PBSsynth]{plotSS.pairs}},
  \code{\link[PBSsynth]{plotSS.pmcmc}},
  \code{\link[PBSsynth]{plotSS.ts}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SS_plots}},
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
