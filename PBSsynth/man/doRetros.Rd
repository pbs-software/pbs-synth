\name{doRetros}
\alias{doRetros}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Perform a Retrospective Analysis
}
\description{
  Run a retrospective analysis using the function \code{'retro'} from the package \pkg{r4ss}.
}
\usage{
doRetros(strSpp="418", assyr=2024, stock="BC", newsubdir="retros", 
   exe="C:/Users/haighr/Files/Archive/Bat/ss.exe")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strSpp}{\code{character} -- alpha-numeric Hart+ code}
  \item{assyr}{\code{numeric} -- assessment year}
  \item{stock}{\code{character} -- stock name(s)}
  \item{newsubdir}{\code{character} -- subdirectory where results will be placed}
  \item{exe}{\code{character} -- SS3 executable, including directory path}
}
\details{
  Function works backward through the years, successively removing data for the last year.
}
\value{
  A binary file called \code{'retroModels.rda'} is created and contains various objects:\cr
  \code{"retroModels","strSpp","assyr","stock","years","run","rwt","retdir"}
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
  Last modified \code{Rd: 2024-08-16}
}
\note{
  To plot retrospective results, run \code{'plotSS.comparisons'},
  which was \emph{borrowed} from function \code{'SSplotComparisons'}
  in package \pkg{r4ss}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.comparisons}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{retro}},
  \code{\link[r4ss]{SSplotComparisons}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{utilFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
