\name{getSS.output}
\alias{getSS.output}
\title{
  Create a List Object from SS3 Output
}
\usage{
getSS.output(dir="C:/myfiles/mymodels/myrun/", dir.mcmc=NULL,
  repfile="Report.sso", compfile="CompReport.sso", covarfile="covar.sso",
  forefile="Forecast-report.sso", wtfile="wtatage.ss_new",
  warnfile="warning.sso", ncols=lifecycle::deprecated(),
  forecast=TRUE, warn=TRUE, covar=TRUE, readwt=TRUE, verbose=TRUE,
  printstats=TRUE, hidewarn=FALSE, NoCompOK=TRUE, aalmaxbinrange=4)
}
\arguments{
  \item{dir}{\code{character} -- a file path to the directory of interest;
    the default value is \code{dir=NULL}, which leads to using the current working directory}
  \item{dir.mcmc}{\code{character} -- pptional directory containing MCMC output;
    this can either be relative to \code{dir}, such that \code{file.path(dir, dir.mcmc)}
    will end up in the right place, or an absolute path}
  \item{repfile}{\code{character} -- name of the big report file (could be renamed by user)}
  \item{compfile}{\code{character} -- name of the composition report file}
  \item{covarfile}{\code{character} -- name of the covariance output file}
  \item{forefile}{\code{character} -- name of the forecast file}
  \item{wtfile}{\code{character} -- name of the file containing weight at age data}
  \item{warnfile}{\code{character} -- name of the file containing warnings}
  \item{ncols}{\code{deprecated} -- this value is now calculated automatically}
  \item{forecast}{\code{logical} -- if \code{TRUE}, read the forecast-report file}
  \item{warn}{\code{logical} -- if \code{TRUE}, read the \code{Warning.sso} file}
  \item{covar}{\code{logical} -- if \code{TRUE}, read \code{covar.sso}}
  \item{readwt}{\code{logical} -- if \code{TRUE}, read the weight-at-age file}
  \item{verbose}{\code{logical} -- if \code{TRUE}, print output to screen}
  \item{printstats}{\code{logical} -- if \code{TRUE}, print summary statistics about 
    the output to the R GUI}
  \item{hidewarn}{\code{logical} -- if \code{TRUE}, hide warnings output from the R GUI}
  \item{NoCompOK}{\code{logical} -- if \code{TRUE}, allow the function to work without
    a \code{CompReport} file}
  \item{aalmaxbinrange}{\code{numeric} -- largest length bin range allowed for composition
    data to be considered as conditional age-at-length data}
}
\value{
  Many values are returned. Complete list would be quite long, but
  should probably be created at some point in the future.
}
\description{
  Read the \code{Report.sso} and (optionally) the \code{covar.sso}, 
  \code{CompReport.sso} and other files produced by Stock Synthesis and
  formats the important  content of these files into a list in the R workspace.
  A few statistics unavailable elsewhere are taken from the \code{.par} file.
  Summary information and statistics can be returned to the R console or just
  contained within the list produced by this function.
}
\examples{
\dontrun{
## read model output
myreplist <- getSS.output(dir="c:/SS/Simple/")

## make a bunch of plots
SS_plots(myreplist)

# read model output and also read MCMC results (if run), which in
# this case would be stored in c:/SS/Simple/mcmc/
myreplist <- getSS.output(dir="c:/SS/Simple/", dir.mcmc="mcmc")
}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SS_output'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{r4ss}:\cr
\code{\link[r4ss]{SS_output}}
}
\author{
  \href{mailto:ian@iphc.int}{Ian Stewart}, International Pacific Halibut Commission, Seattle WA\cr
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{IO}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{r4ssFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
