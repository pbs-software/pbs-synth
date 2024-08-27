\name{predictRec}
\alias{predictRec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Predict Recruitment
}
\description{
  Fit rockfish recruitment deviations to environmental indices.
}
\usage{
predictRec(rec, indices, mos=1:12, ryrs=NULL, rfun=median, ifun=mean, 
  nmcmc=2000, refit=TRUE, polyno=2, outnam, 
  png=FALSE, pngres=400, PIN=c(10, 8), lang=c("f","e"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rec}{\code{numeric} -- SS3 MCMC results data frame of recruitment values (e.g., deviations)}
  \item{indices}{\code{character} -- name of index and binary file RDA that contains the index anomalies}
  \item{mos}{\code{numeric} -- vector of months as numbers}
  \item{ryrs}{\code{numeric} -- vector of years to use from the recruitment data;
    e.g., if using recruitment deviations, use years for which the SS3 model 
    calculated MAIN recruitment deviations}
  \item{rfun}{\code{function} -- function used to transform recruitment values
    (e.g., log10 for recruitments, non for deviations)}
  \item{ifun}{\code{function} -- function used to summarise monthly environmental indices}
  \item{nmcmc}{\code{numeric} -- number of MCMC samples extracted from the SS3 model run}
  \item{refit}{\code{logical} -- if{TRUE}, refit the polynomial for every MCMC sample;
    set to FALSE if not changing the number of MCMC samples and plot tweaks are desired}
  \item{polyno}{\code{numeric} -- order of the polynomial to fit}
  \item{outnam}{\code{character} -- explicit output name desired}
  \item{png}{\code{logical} -- if{TRUE}, send the figure to a PNG file}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) of the PNG file}
  \item{PIN}{\code{numeric} -- dimensions (width, height) of the PNG (inches)}
  \item{lang}{\code{character} -- figure language (\code{'e'} for english, \code{'f'} for french)}
}
\details{
  Fit a first-order autoregressive polynomial model of recruitment (dependent)
  to environmental indices (independent)
  \deqn{\epsilon_t = \beta_0 + \beta_1 z_t + \beta_2 z_t^2 + \phi_1}{%
        epsilon[t] = beta[0] + beta[1]*z[t] + beta[2]*z[t]*z[t] + phi[1]}
  where\cr
  \eqn{\epsilon_t}{epsilon[t]} = recruitment deviation in year \eqn{t}\cr
  \eqn{z_t}{z[t]} = environmental index anomaly in year \eqn{t}\cr
  \eqn{\beta_0}{beta[0]} = intercept for the model\cr
  \eqn{\beta_{1,2}}{beta[1,2]} = parametric terms for the model\cr
  \eqn{\phi_1}{phi[1]} = autoregressive process of order 1
}
\value{
  Nothing in particular, although the function creates a wondrous plot
  of the fit for each of the \code{'nmcmc'} MCMC samples of recruitment deviations to
  the environmental index series (qualified and summarised by months selected).
  Annual posteriors of recruitment deviation are displayed as medians (circles
  with year labels) and 90\% credibility intervals (vertical bars). 
  Quantiles of predicted fit are displayed as black lines (solid = median, dashed=90\%
  credibility envelope, dotted=50\% credibility envelope).
  Negative index space is shaded blue, positive shaded orange. 
  Median (and 90\% credibility limits in parentheses) for the fitted parameters appear in the legend.
}
\references{
  Pinheiro, J., Bates, D. and R Core Team. 2023. 
  \href{https://cran.r-project.org/package=nlme}{nlme: Linear and Nonlinear Mixed Effects Models}.
  R package version 3.1-164.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-21}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{calcStdRes}},
  \code{\link[PBSsynth]{getSS.rdevs}}\cr
In package \pkg{nlme}:\cr
  \code{\link[nlme]{gls}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
\keyword{hplot}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{appFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
