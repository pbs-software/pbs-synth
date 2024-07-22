\name{agileDT}
\alias{agileDT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Agile Decision Tables
}
\description{
  Create a PDF of decision tables independently of pre-existing Sweave code.
}
\usage{
agileDT(compo, currYear=2024, projYear=2034, Ngen=3, gen1=25, 
   area=c("5ABC","3CD","5DE"), cp=paste0("CC.", pad0(1:3, 2)), 
   onepage=TRUE, tables.per.page=3, outnam, sigdig=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compo}{\code{list} -- MCMC output from one or more (composite) base runs;
    output includes constant catch and/or harvest rate catch policies;
    \code{compo} obect is created by a call to function \code{'gatherMCMC'}}
  \item{currYear}{\code{numeric} -- current year in model}
  \item{projYear}{\code{numeric} -- last projected year in model}
  \item{Ngen}{\code{numeric} -- number of generations}
  \item{gen1}{\code{numeric} -- length (years) of one generation}
  \item{area}{\code{character} -- area name(s), useful for multi-area models}
  \item{cp}{\code{character} -- names of catch policy dimensions in output arrays}
  \item{onepage}{\code{logical} -- if \code{TRUE}, render each decision table
    on a separate page of the PDF output}
  \item{tables.per.page}{\code{numeric} -- number of decision tables to place on one page}
  \item{outnam}{\code{character} -- output name for the PDF file}
  \item{sigdig}{\code{numeric} -- significant digits to use for table probabilities}
}
\details{
  Uses the function \code{'findTarget'} to calculate the probabilities and
  arranges the information into decision tables. These are compiled using 
  \code{'pdflatex.exe'} and \code{'latexmk.exe'}.
}
\value{
  Save output (\code{'agile'}, \code{'is.good'}, and \code{'Nmcmc'}) to a binary \code{'.rda'} file.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-17}
}
\note{
  Useful for a quick and dirty set of decision tables when running
  exploratory analyses, but is also used for the final set of decision
  tables for one base run or a base case (composite set of base runs).
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
\concept{utilFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
