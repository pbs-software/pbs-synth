\name{plotSS.pmcmc}
\alias{plotSS.pmcmc}
\alias{plotSS.dmcmc}
\alias{plotSS.rmcmc}
\title{
  Plot MCMC Time Series
}
\description{
  Plot either boxplots of quantiles or credibility envelopes 
  for various MCMC derived parameters.
}
\usage{
plotSS.pmcmc(obj, pqs=tcall(quants5), xyType="quantBox", lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL, userPrompt=FALSE, save=TRUE,
   tcl.val=-0.2, yrs, y0=TRUE, pyrs=NULL, LRP=NULL, USR=NULL, catpol=NULL,
   yaxis.by, yLab="Recruitment", outnam, lang=c("e","f"), ptypes="win",
   pngres=400, PIN=c(8, 6), ...)

plotSS.dmcmc(mcmcObj, mpdObj, ptypes, lang, pngres=400, PIN=c(9,9))

plotSS.rmcmc(mcmcObj, mpdObj, ptypes, lang, pngres=400, PIN=c(9,9))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{\code{numeric} -- data frame or matrix of MCMC 
    time series estimates (rows=samples, columns=years)}
  \item{pqs}{\code{numeric} -- vector of quantiles}
  \item{xyType}{\code{character} -- type of plot; choose between \code{'quantBox'} and \code{'envelope'}}
  \item{lineType}{\code{numeric} -- line types to use in \code{quantBox}}
  \item{refLines}{\code{numeric} -- reference lines to add to plot}
  \item{xLim}{\code{numeric} -- limits of the x-axis}
  \item{yLim}{\code{numeric} -- limits of the y-axis}
  \item{userPrompt}{\code{logical} -- \bold{not used}}
  \item{save}{\code{logical} -- \bold{not used}}
  \item{tcl.val}{\code{numeric} -- tick length}
  \item{yrs}{\code{numeric} -- reconstruction years to use along x-axis}
  \item{y0}{\code{logical} -- if \code{TRUE}, start the y-axis at zero}
  \item{pyrs}{\code{numeric} -- projection years to use along x-axis}
  \item{LRP}{\code{numeric} -- value of the limit reference point on x-axis
    (e.g., \code{0.2} for \eqn{B_t/B_0}{Bt/B0}, \code{0.4} for \eqn{B_t/B_\text{MSY}}{Bt/Bmsy})}
  \item{USR}{\code{numeric} -- value of the upper stock reference point on x-axis
    (e.g., \code{0.4} for \eqn{B_t/B_0}{Bt/B0}, \code{0.8} for \eqn{B_t/B_\text{MSY}}{Bt/Bmsy})}
  \item{catpol}{\code{numeric} -- catch policy (usually constant catch in tonnes)}
  \item{yaxis.by}{\code{numeric} -- increments along the y-axis to place tick marks}
  \item{yLab}{\code{character} -- label for the y-axis}
  \item{outnam}{\code{character} -- name for output file(s), usually a figure}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
  \item{ptypes}{\code{character} -- vector of plot types (\code{'win'}, \code{'png'}, \code{'eps'})}
  \item{pngres}{\code{numeric} -- resolution for figure files (\code{'png'}, \code{'eps'})}
  \item{PIN}{\code{numeric} -- plot size for output figures in inches (width, height)}
  \item{\dots}{\code{dots} -- additional arguments passed to subfunction \code{'plt.qB'}}
  \bold{Additional arguments:}
  \item{mcmcObj}{\code{numeric} -- list of data frames of MCMC 
    time series estimates (rows=samples, columns=years) for 
    B=spawning biomass, R=recruitment, and P=parameters}
  \item{mcmcObj}{\code{numeric} -- list of vectors of MPD 
    time series estimates (rows=samples, columns=years) for 
    B=spawning biomass, R=recruitment, and P=parameters}
}
\details{
  Used primarily for generating time series figures for model run summaries.

  \tabular{rrl}{
    \code{plotSS.pmcmc} \tab trajectories \tab recruits, exploit, deplete, sbiomass, recdev, boverbmsy, uoverumsy\cr
    \code{plotSS.dmcmc} \tab diagnostics  \tab traceBiomass, traceRecruits, traceParams, splitChain, paramACFs, rhat\cr
    \code{plotSS.rmcmc} \tab routine      \tab pairsPars, pdfParameters, snail
  }
}
\value{
  Makes beauteous figures.
}
\author{
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew Edwards}, Quantitative Assessment Methods, PBS, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-21}
}
\note{
  Modified from \pkg{PBSawatea}'s \code{'plotBmcmcPOP'} function by Andrew Edwards (PBS).
  Input comes from the SS3 output file 'derived.parameters.sso' (use \code{r4ss::SSgetMCMC})
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{gatherMCMC}},
  \code{\link[PBSsynth]{plotSS.compo}},
  \code{\link[PBSsynth]{plotSS.senso}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SS_output}},
  \code{\link[r4ss]{SSgetMCMC}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{plotFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
