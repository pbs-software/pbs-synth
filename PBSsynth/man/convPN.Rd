\name{convPN}
\alias{convPN}
\title{
  Convert Parameter Names
}
\description{
  Convert parameter names from SS3's obtuse naming system to either
  Awatea equivalents or more concise SS names, depending on the 
  parameter overlap between the two platforms.
}
\usage{
convPN(pnams)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pnams}{\code{character} -- SS3 parameter names (use the active ones only).}
}
\details{
  Execute a series of \code{sub} statements followed by some text manipulation.
}
\value{
  Parameter names transformed to something more mathy.
  Use Awatea parameter names for SS parameters that parallel/mimic Awatea parameters.
}
\references{
\href{https://vlab.noaa.gov/web/stock-synthesis/document-library}{SS3 Documentation by Methot et. al}
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Regional Headquartes (RHQ), Vancouver BC\cr
  Last modified \code{Rd: 2022-11-01}
}
\seealso{
In \pkg{PBSsynth}:\cr
\code{\link[PBSsynth]{prepMPD}}, etc.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
