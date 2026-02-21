\name{plotSS.yield}
\alias{plotSS.yield}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Yield and Surplus Production
}
\description{
  Plot yield and surplus production from Stock Synthesis output. 
  Surplus production is based on Walters et al. (2008).
}
\usage{
plotSS.yield(replist, subplots=1:5, 
   refpoints=c("MSY", "Btgt", "SPR", "Current", "LRP", "USR"), 
   add=FALSE, plot=TRUE, print=FALSE, labels=c("Fraction unfished",
   "Equilibrium yield (t)", "Total biomass (t)", "Surplus production (t)",
   "Yield per recruit (kg)", "Spawning output"), col="black", col2="black",
   lty=1, lwd=2, cex.main=1, pwidth=6.5, pheight=5, punits="in", 
   res=400, ptsize=10, plotdir="./english", verbose=TRUE, lang=c("f","e")) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{composite} -- list object created by \code{SS_output()}}
  \item{subplots}{\code{numeric} -- vector controlling which subplots to create\cr
    Numbering of subplots as follows:
    \enumerate{
      \item yield curve
      \item yield curve with reference points
      \item surplus production vs. total biomass plots (Walters et al. 2008)
      \item surplus production vs. spawning biomass plots (Forrest et al. 2023)
      \item yield per recruit plot
    } }
  \item{refpoints}{\code{character} -- vector of reference points to display 
    in subplot 2; options: \code{'MSY'}, \code{'Btgt'}, \code{'SPR'}, \code{'Current'},
    \code{'LRP'}, and \code{'USR'} }
  \item{add}{\code{logical} -- if \code{TRUE}, add something to current plot (\emph{pending})}
  \item{plot}{\code{logical} -- if \code{TRUE}, plot to active device}
  \item{print}{\code{logical} -- if \code{TRUE}, print to PNG files}
  \item{labels}{\code{character} -- vector of labels for plots (titles and axis labels)}
  \item{col}{\code{numeric|character} -- line colour for equilibrium plot\cr
    \emph{note} -- internal object \code{'colvec'} is used to control vertical line colours in equilibrium plot}
  \item{col2}{\code{numeric|character} -- line color for dynamic surplus production plot}
  \item{lty}{\code{numeric} -- line type (only applied to equilibrium yield curve)\cr
    \emph{note} -- internal object \code{'ltyvec'} is used to control vertical lines types in equilibrium plot}
  \item{lwd}{\code{numeric} -- line width for plot elements}
  \item{cex.main}{\code{numeric} -- character expansion for plot titles (default=1)}
  \item{pwidth}{\code{numeric} -- default width of plots printed to files in units of punits}
  \item{pheight}{\code{numeric} -- height of plots printed to PNG files in units of punits;
    default is designed to allow two plots per page, with \code{pheight\_tall} used for plots
    that work best with a taller format and a single plot per page}
  \item{punits}{\code{character} -- units for \code{pwidth} and \code{pheight};
    can be \code{'px'} (pixels), \code{'in'} (inches), \code{'cm'} (centimeters), 
    or \code{'mm'} (millimeters); the default is \code{'in'}}
  \item{res}{\code{numeric} -- resolution of plots printed to files}
  \item{ptsize}{\code{numeric} -- point size for text in plots printed to files;
    type \code{help("png")} in R's command line for details}
  \item{plotdir}{\code{character} -- directory where PNG files will be written}
  \item{verbose}{\code{logical} -- if \code{TRUE}, print output to the console}
  \item{lang}{\code{character} -- language for text in plots, 
    where \code{'e'} is english and \code{'f'} is french}
}
\details{
  Handy way to show \eqn{B_0}{B0} reference points on the equilibrium yield curve plotted in SSB space.
}
\value{
	Nada to console; gorgeous plots to device
}
\references{
Walters, C.J., Hilborn, R., and Christensen V. (2008)
\href{https://doi.org/10.1139/F08-170}{%
Surplus production dynamics in declining and recovering fish populations}.
\emph{Can. J. Fish. Aquat. Sci.} \bold{65}: 2536--2551. https://doi.org/10.1139/F08-170.

Forrest, R.E., Kronlund, A.R., Cleary, J.S., and Grinnell, M.H. (2023)
\href{https://doi.org/10.1139/cjfas-2022-0168.}{%
An evidence-based approach for selecting a limit reference point for Pacific herring
(\emph{Clupea pallasii}) stocks in British Columbia, Canada}.
\emph{Can. J. Fish. Aquat. Sci.} \bold{80}: 1071--1083.
}
\author{
  \href{mailto:ian@iphc.int}{Ian Stewart}, International Pacific Halibut Commission, Seattle WA\cr
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2026-02-18}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'SSplotYield'} from
  the R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.comps}},
  \code{\link[PBSsynth]{plotSS.francis}},
  \code{\link[PBSsynth]{plotSS.pars}},
  \code{\link[PBSsynth]{plotSS.profile}},
  \code{\link[PBSsynth]{plotSS.rdevs}},
  \code{\link[PBSsynth]{plotSS.selex}},
  \code{\link[PBSsynth]{plotSS.stdres}},
  \code{\link[PBSsynth]{plotSS.stock.recruit}},
  \code{\link[PBSsynth]{plotSS.ts}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSplotYield}}
}
\examples{
\dontrun{
## WWR R29v2 -- MPD rerun with benchmark reference point option = 1
##             (calculate F_spr, F_btgt, F_msy)
wwr.dir = "C:/Users/haighr/Files/GFish/PSARC25/SGR/Data/SS3/SGR2025/Run29/MPD.29.01.v2.Btgt"
replist = SS_output(www.dir)
plotSS.yield(replist, subplots=2, refpoints=c("MSY","Btgt","Current","LRP","USR"), print=F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{r4ssFuns}
% \concept{ ~cpt2 }
% Use only one concept per line.
