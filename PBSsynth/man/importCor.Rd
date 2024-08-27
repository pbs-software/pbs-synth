\name{importCor}
\alias{importCor}
\alias{importEva}
\alias{importPar}
\alias{importStd}
\title{
  Import SS3 Output Files
}
\description{
  Import SS3 output files: correlation, eigenvalues, parameters, and standard deviations.
}
\usage{
importCor(cor.file)
importEva(eva.file)
importPar(par.file) 
importStd(std.file, vnam="name")
}
\arguments{
  \item{cor.file}{\code{character} -- name of correlation file (e.g., \code{'ss3.cor'})}
  \item{eva.file}{\code{character} -- name of eigenvalue file (e.g., \code{'ss3.eva'})}
  \item{par.file}{\code{character} -- name of parameter file (e.g., \code{'ss3.par'})}
  \item{std.file}{\code{character} -- name of standard deviation file (e.g., \code{'ss3.std'})}
  \item{vnam}{\code{character} -- column name(s) to exclude from standard deviation table}
}
\details{
  Extracts data from character vectors and makes various data objects.
}
\value{
  The \bold{correlation} output is a list object with the following components:
  \item{cfile}{character vector representing the file line-by-line}
  \item{cmat}{qualified \code{cfile}}
  \item{cor}{data.frame of the correlation file}
  \item{cor.mat}{matrix representing the correlation matrix only}
  \item{index}{character vector used in the row and column names of \code{cor} and \code{cor.mat}}
  \item{cor.name}{Awatea parameter names in the correlation file}
  \item{cor.value}{Awatea parameter values in the correlation file}
  \item{cor.std.dev}{Awatea parameter standard deviations in the correlation file}
  \item{hessian_log_determinant}{log of the determinant of the hessian reported in the header}

  The \bold{eigenvalue} output is a list object with the following components:
  \item{efile}{character vector representing the file line-by-line}
  \item{eva}{data.frame of the eigenvalue file (colnames: unsorted, sorted}

  The \bold{parameter} output is a list object with the following components:
  \item{pfile}{character vector representing the file line-by-line}
  \item{npars}{number of parameters}
  \item{fval}{function value}
  \item{maxgrad}{maximum gradient}
  \item{dummy_parm}{number of dummy parameters}
  \item{<par_name>}{parameter values from 1 to \code{npars}}

  The \bold{standard deviation} output is a list object with the following components:
  \item{sfile}{character vector representing the file line-by-line}
  \item{std}{data.frame of the standard deviation file (colnames: index, name, value, std.dev}
  \item{<various>}{assorted quantities: \code{'recdev_early'}, \code{'recdev2'},
    \code{'Fcast_recruitments'}, \code{'SSB_std'}, \code{'recr_std'}, \code{'SPR_std'}, 
    \code{'F_std'}, \code{'depletion'}, \code{'Mgmt_quant'}, \code{'Extra_Std'}}
  \item{<par_name>}{standard deviation of selected parameters}
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-22}
}
\seealso{
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{read.admbFit}},
  \code{\link[r4ss]{SS_readpar_3.30}},
  \code{\link[r4ss]{getADMBHessian}}
}
\keyword{file}
\keyword{interface}
\concept{utilFuns.r}

