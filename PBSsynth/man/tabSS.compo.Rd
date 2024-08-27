\name{tabSS.compo}
\alias{tabSS.compo}
\alias{tabSS.decision}
\alias{tabSS.senso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create Tables for Base and Sensitivity Runs
}
\description{
  Create tables for a base run (or a composite base case), including
  decision tables, and for sensitivity runs.
}
\usage{
tabSS.compo(istock="YTR", prefix="ytr.", compo, useRlow=FALSE, qRlow=0.25, sigdig=4)

tabSS.decision(istock="YTR", prefix="ytr.", compo, useRlow=FALSE, qRlow=0.25, 
   decdig=2, cp.type="CC", cp.num=c(1,3,6,9:15), use.agile=FALSE)

tabSS.senso(istock="YTR", prefix="ytr.", senso, sigdig=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{istock}{\code{character} -- three- or four-letter code for species}
  \item{prefix}{\code{character} -- prefix for table labels in LaTeX}
  \item{compo}{\code{list} -- data object containing results from MCMC base run or composite base case runs
    (see function \code{'gatherMCMC'})}
  \item{qRlow}{\code{numeric} -- lower quantile of recruitments to subset MCMC data}
  \item{sigdig}{\code{numeric} -- significant digits to display for table entries}
  \item{decdig}{\code{numeric} -- decimal places to display for table entries}
  \item{cp.type}{\code{character} -- catch policy type, either \code{'CC'} for constant catch
    or \code{'HR'} for harvest rate}
  \item{cp.num}{\code{numeric} -- subset of catch policies from an available set}
  \item{use.agile}{\code{logical} -- if \code{TRUE}, use the function \code{'agileDT'} to gather decision tables}
  \item{senso}{\code{list} -- data object containing results from MCMC sensitivity runs
    (see function \code{'gatherMCMC'})}
}
\details{
  Creates tables in LaTeX for use later by Sweave routine.
  Tables include quantile summaries of parameter estimates and derived quantities,
  and decision tables for catch projections that summarise the probability of 
  spawning biomass being above references points like LRP or USR 
  (or harvest rates being below reference points like RR).
}
\value{
  Binary \code{'.rda'} file containing LaTeX tables.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-16}
}
\note{
  LRP = Limit Reference Poin\cr
  USR = Upper Stock Reference\cr
  RR  = Removal Reference
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{gatherMCMC}},
  \code{\link[PBSsynth]{plotSS.compo}},
  \code{\link[PBSsynth]{plotSS.senso}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{print}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{appFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
