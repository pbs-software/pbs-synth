\name{calcStdRes}
\alias{calcStdRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Studentised Residuals
}
\description{
  Calculate standardised residuals for age proportions fit by the Multinomial or Robustified Normal likelihoods.
}
\usage{
calcStdRes(obj, trunc=3, myLab="Age Residuals", prt=TRUE, type="Multinomial", 
   afld="Bin", yfld="Yr", ofld="Obs", ffld="Exp", nfld="Nsamp_adj")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{\code{data.frame} -- database table of ages from \code{'replist'}}
  \item{trunc}{\code{numeric} -- maximum magnitude of standardised residual to use}
  \item{myLab}{\code{character} -- label for printed output}
  \item{prt}{\code{logical} -- if \code{TRUE}, concatenate standard deviation of standardised residuals to the console}
  \item{type}{\code{character} -- type of distribution (\code{'Multinomial'}, \code{'Fournier'}, or \code{'Coleraine'})}
  \item{afld}{\code{character} -- name of field containing age bins}
  \item{yfld}{\code{character} -- name of field containing years}
  \item{ofld}{\code{character} -- name of field containing observed values}
  \item{ffld}{\code{character} -- name of field containing expected values}
  \item{nfld}{\code{character} -- name of field containing sample sizes}
}
\details{
  The calculation of residuals depends on distribution.\cr
  Pearson residual = \code{(O-F)/std.dev(O)}\cr
  See CASAL Manual, Section 6.8 Residuals, where \code{std.dev(O)} is calculated as:
\preformatted{
## Fournier uses Fitted, Coleraine use Observed, Multinomial uses Fitted
Nprime  <- min( obj[,nfld][idx],1000)                         ## N prime
Fprime  <- obj[,ffld][idx]*(1.0-obj[,ffld][idx]) + 0.1/nage   ## F prime
Oprime  <- obj[,ofld][idx]*(1.0-obj[,ofld][idx]) + 0.1/nage   ## O prime
Mprime  <- Z(obj[,ffld][idx]) * (1-Z(obj[,ffld][idx]))        ## M prime
SD      <- sqrt(
  switch(type,'Multinomial'=Mprime,'Fournier'=Fprime,'Coleraine'=Oprime)/Nprime)
## Pearson residuals = Normalised residuals for normal error distributions:
result$stdRes[idx] <- res[idx]/SD
}
}
\value{
  Data frame with field name \code{'stdRes'} containing
  the Pearson residual based on distribution \code{type}.
}
\references{
Bull, B., Francis, R.I.C.C., Dunn, A., McKenzie, A., Gilbert, D.J. and Smith, M.H. (2005)
\href{https://docs.niwa.co.nz/library/public/NIWAtr127.pdf}{CASAL (C++ algorithmic stock assessment laboratory): CASAL User Manual v2.07-2005/08/21}.
NIWA Tech. Rep. 127. 274 p.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-04-11}
}
\note{
  No support for the Dirichlet-Multinomial (yet); however, see SS3 code \sQuote{SS_write_report.tpl} where:\cr

  \preformatted{show_Pearson = value((ocomp-ecomp) / sqrt(ecomp*(1.0-ecomp) / 
     nsamp*(nsamp+dirichlet_Parm) / (1.+dirichlet_Parm)));
   // Pearson for Dirichlet-multinomial using negative-exponential parameterization}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.stdres}}
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
