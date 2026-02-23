\name{copySSfiles}
\alias{copySSfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Copy SS Files to Git Archive
}
\description{
  Copy SS input files to an archive folder that can be copied later to a git repository.
}
\usage{
copySSfiles(strSpp="405", assyr=2025)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strSpp}{\code{character} -- string species Hart code}
  \item{assyr}{\code{numeric} -- year that \code{strSpp} stock was assessed}
}
\details{
  Uses R's \code{'file.copy'} function to copy SS files from model runs to 
  an input archive directory that can later be copied to a local git repository.
}
\value{
  A beautiful transfer. A transfer like nobody has ever seeen before.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-02-23}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{base}:\cr
  \code{\link[base]{dir.create}},
  \code{\link[base]{dir.exists}},
  \code{\link[base]{file.copy}},
  \code{\link[base]{file.path}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{IO}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
