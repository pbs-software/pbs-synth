\name{weightAF}
\alias{weightAF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Weight Age Frequency Data
}
\description{
  Weight age frequencies using harmonic mean ratio method.
}
\usage{
weightAF(replist, fleets, abase="agedbase", afld="Nsamp_adj", 
   hfld="effN", rbase="recruit", rfld="pred_recr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- SS3 results read in by \code{'r4ss::SS_output()'}}
  \item{fleets}{\code{numeric} -- fleet number(s)}
  \item{abase}{\code{character} -- model output name for the age database}
  \item{afld}{\code{character} -- field name containing ages}
  \item{hfld}{\code{character} -- field name containing effective sample sizes}
  \item{rbase}{\code{character} -- model output name for the recruits database}
  \item{rfld}{\code{character} -- field name containing predicted recruitments}
}
\details{
  The harmonic mean ratio method of McAllister and Ianelli (1997) involves
  calculating  effective sample size:
  \deqn{\widetilde{n}_{tg} = \frac{\sum_{a=1}^A \widehat{p}_{atg} (1 - \widehat{p}_{atg})}{\sum_{a=1}^A (p_{atg} - \widehat{p}_{atg}} }{%
    n_tilde[tg] = ( sum [a=1:A] ( p_hat[atg] * (1 - p_hat[atg]) ) ) / ( sum [a=1:A] ( p[atg] - p_hat[atg] ) ) }
  and then calculating the harmonic mean ratio:
  \deqn{w_g = \frac{ \sum_t N_g / (1/\widetilde{n}_{tg}) }{ \sum_t n_{tg} / N_g } }{%
    w[g] = ( sum [t] N_g / (1/n_tilde[tg]) ) / ( sum [t] n[tg] / N[g] ) }
  where\cr
  \eqn{p_{atg}}{p[atg]} = observed proportions at age \eqn{a} in year \eqn{t} for fleet \eqn{g}\cr
  \eqn{\widehat{p}_{atg}}{p_hat[atg]} = expected proportions at age \eqn{a} in year \eqn{t} for fleet \eqn{g}\cr
  \eqn{\widetilde{n}_{tg}}{n_tilde[tg]} = effective sample size in year \eqn{t} for fleet \eqn{g}\cr
  \eqn{n_{tg}}{n[tg]} = observed sample size in year \eqn{t} for fleet \eqn{g}\cr
  \eqn{N_g}{N[g]} = number years with AF data for fleet \eqn{g}
}
\value{
  List object containing:
  \item{\code{agedat}}{age database qualified by selected fleet(s)}
  \item{\code{hmean}}{harmonic mean sample size for each fleet}
  \item{\code{amean}}{arithmetic mean sample size for each fleet}
  \item{\code{wmean}}{ratio of harmonic mean to arithmetic mean;
    use this value to adjust AF sample size (e.g., 2021 YMR)}
  \item{\code{Rmax}}{maximum recruitment estimate}
  \item{\code{M}}{estimated (or fixed) female natural mortality}
  \item{\code{wsub}}{some sort of weird-ass discount factor\cr
    \code{1/((Rmax/10^floor(log10(Rmax)))*M^0.7)}}
  \item{\code{wadj}}{ratio \code{'wmean'} adjusted by \code{'wsub'}; not used}
}
\references{
  McAllister, M.K. and Ianelli, J.N. 1997. 
  \href{https://doi.org/10.1139/f96-285}{%
  Bayesian stock assessment using catch-age data and the sampling-importance resampling algorithm}.
  \emph{Can. J. Fish. Aquat. Sci.} \bold{54}(2). 284--300.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-22}
}
\note{
  Calculation of harmonic mean weights analagous to the Francis mean-age weights.
  
  \emph{Extract from the SS3 (v.3.30.22.1) manual:}\cr
  The \code{'Length_Comp_Fit_Summary'} and \code{'Age_Comp_Fit_Summary'} sections in the
  Report file include information on the harmonic mean of the effective sample size and
  arithmetic mean of the input sample size used in this tuning method.
  In the \code{r4ss} package, these tables are returned by the \code{'SS_output'} function as
  \code{'$Length_comp_Eff_N_tuning_check'} and\cr \code{'$Age_comp_Eff_N_tuning_check'}.
  A convenient way to process these values into the format required by the control file
  is to use the function:\cr
  \code{SS_tune_comps(replist, option='MI')} 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.francis}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSMethod.TA1.8}},
  \code{\link[r4ss]{SS_tune_comps}}
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
