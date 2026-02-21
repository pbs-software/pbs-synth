\name{agileDT}
\alias{agileDT}
\alias{quickDT}
\alias{tabDQs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Agile Decision Tables
}
\description{
  Create a PDF of decision tables independently of pre-existing Sweave code.
}
\usage{
agileDT(compo, currYear=2026, projYear=2036, Ngen=3, gen1=25,
   area=c("5ABC","3CD","5DE"), cp=paste0("CC.",pad0(1:3, 2)),
   useRlow=FALSE, RPbase="B0", pdf=FALSE, onepage=TRUE, 
   tables.per.page=3, outnam, sigdig=2)

quickDT(compo, vals=c("BtLRP","BtUSR","BtBmsy","utumsy"),
   areas=c("5ABC","3CD","5DE"), cp=paste0("CC.", pad0(1:3,2)) )

tabDQs(xavgRP, xavgTS, csv=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compo}{\code{composite} -- MCMC output from one or more (composite) base runs;
    output includes constant catch and/or harvest rate catch policies;
    \code{'compo'} obect is created by a call to function \code{'gatherMCMC'}}
  \item{currYear}{\code{numeric} -- current year in model}
  \item{projYear}{\code{numeric} -- last projected year in model}
  \item{Ngen}{\code{numeric} -- number of generations}
  \item{gen1}{\code{numeric} -- length (years) of one generation}
  \item{area}{\code{character} -- area name(s), useful for multi-area models}
  \item{cp}{\code{character} -- names of catch policy dimensions in output arrays}
  \item{useRlow}{\code{logical} -- if \code{TRUE}, the MCMC will be reduced to records
    with the lowest recruitment}
  \item{RPbase}{\code{character} -- reference point base, either \code{'B0'} or \code{'BMSY'}}
  \item{pdf}{\code{logical} -- if \code{TRUE}, create PDF files of decision tables}
  \item{onepage}{\code{logical} -- if \code{TRUE}, render each decision table
    on a separate page of the PDF output}
  \item{tables.per.page}{\code{numeric} -- number of decision tables to place on one page}
  \item{outnam}{\code{character} -- output name for the PDF file}
  \item{sigdig}{\code{numeric} -- significant digits to use for table probabilities}
  \bold{Additional arguments:}
  \item{vals}{\code{character} -- field names of derived quantities in projection
    array object within \code{'compo'}}
  \item{areas}{\code{character} -- subareas or regions along the BC coast}
  \item{xavgRP}{\code{numeric} -- three-dimensional array of reference
    points from a multi-area model:\cr
    \code{'names(dimnames(xavgRP)) = c('mcmc','val','area')'}}
  \item{xavgTS}{\code{numeric} -- four-dimensional array of time series
    values from a multi-area model:\cr
    \code{'names(dimnames(xavgRP)) = c('mcmc','year','val','area')'}}
  \item{csv}{\code{logical} -- if \code{TRUE}, send tables to a CSV file}
}
\details{
  \code{'agileDT'} uses the function \code{'findTarget'} to calculate 
  probabilities, and arranges the information into decision tables. 
  These are compiled using \code{'pdflatex.exe'} and \code{'latexmk.exe'}.
}
\value{
  Save output (\code{'agile'}, \code{'is.good'}, and \code{'Nmcmc'})
  to a binary \code{'.rda'} file.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-02-20}
}
\note{
  Function \code{'quickDT'} offers a quick and dirty 
  method for extracting a set of calculated reference points 
  (e.g., \code{'BtLRP'}) from a single base run when using a multi-area
  model, and printing probabilties to decision tables.

  Function \code{'agileDT'} arose from \code{'quickDT'} for use in
  constructing the final set of decision tables for one base run or a 
  base case (composite set of base runs).
  
  Function \code{'tabDQs'} provides a quick extraction of derived quantities
  and their ratios with reference points (probably for PJS).
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{findTarget}}\cr
In package \pkg{PBStools}:\cr
  \code{\link[PBStools]{formatCatch}},
  \code{\link[PBStools]{texArray}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{models}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
