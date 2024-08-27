\name{extract.between}
\alias{extract.between}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract Between
}
\description{
  Extract character strings between two delimiters.
}
\usage{
extract.between(x, open="{", close="}", first.only=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{character} -- string containing delimiters}
  \item{open}{\code{character} -- opening delimiter character}
  \item{close}{\code{character} -- closing delimiter character}
  \item{first.only}{\code{logical} -- if \code{TRUE}, extract all instances
    of text between the delimiters; if \code{FALSE}, extract only the first}
}
\details{
  Extracts text that occurs between specified delimiters.
  Involves functions \code{'regmatches'} and \code{'gregexpr'}
  suggested by users on stackoverflow (see References)
}
\value{
  snippets -- all instances of text between the user-specified delimiters
}
\references{
  \href{https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r}{%
  stackoverflow: Extract info inside all parenthesis in R (comment 77)}
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-22}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{allEqual}},
  \code{\link[PBSsynth]{getYrIdx}},
  \code{\link[PBSsynth]{is.numStr}},
  \code{\link[PBSsynth]{med5.95}}\cr
In package \pkg{base}:\cr
  \code{\link[base]{gregexpr}},
  \code{\link[base]{regmatches}}
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
