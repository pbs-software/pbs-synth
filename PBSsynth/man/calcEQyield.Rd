\name{calcEQyield}
\alias{calcEQyield}
\alias{calcMSY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Equilibrium Yield
}
\description{
  Calculate euilibrium yield from SS3 output.
}
\usage{
calcEQyield(replist, areas=c("5ABC","3CD","5DE"))

calcMSY(replist, strategy, method=1, proj_gears=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{replist}{\code{list} -- SS3 results read in by \code{'r4ss::SS_output()'}}
  \item{areas}{\code{character} -- subareas or regions (PMFC areas) along the BC coast}
  \bold{Additional arguments:}
  \item{strategy}{\code{list} -- catch/harvest policies:
    \itemize{
      \item \code{StrategyType} -- strategy type (1=constant catch; 2=harvest rate)
      \item \code{MSYmaxIter} -- maximum number of iterations
      \item \code{MSYtol} -- tolerance for convergence
      \item \code{StartStrategy} -- start strategy (lower bound of catch for projections)
      \item \code{EndStrategy} -- end strategy (upper bound of catch for projections)
      \item \code{StepStrategy} -- step strategy (interval to use for catch projections)
  }}
  \item{method}{\code{numeric} -- method of fishing 
    (e.g., trawl, other, hook and line, midwater trawl, etc.);
    depends on order of method in the model}
  \item{proj_gears}{\code{logical} -- if \code{TRUE}, use end-year harvest rates by method;
    if {FALSE}, use user-supplied harvest rate proportions}
}
\details{
  \code{calcEQyield}:\cr
  Ian Taylor (see section `Taylor, pers. comm. Apr 4, 2023') suggested using 
  \code{'replist$equil_yield'} to obtain MSY by area.

  Works for MPD but not MCMC samples because \code{'equil_yield'}
  table missing from individual sample replists; function rendered
  obsolete but interesting.

  \code{calcMYS}:\cr
  Attempt to replicate Awatea's method of brute-force calculations to estimate MSY.
  Attempt abandoned for SS3 because this model platform is far more complicated.
  Function retained in case a future fool wishes to resume the quest.
}
\value{
  \code{calcEQyield}: list object containing\cr
  \code{EQY_table, MSY_area, MSY_prop, Rdexpyr, Rpropyr}

  \code{calcMSY}: list object containing\cr
  \code{MSYtab, MSY, msy.idx}
}
\references{
  Taylor, I.G., Doering, K.L., Johnson, K.F., Wetzel, C.R., Stewart, I.J., 2021.
  \href{https://doi.org/10.1016/j.fishres.2021.105924}{%
  Beyond visualizing catch-at-age models: Lessons learned from the r4ss package about software to support stock assessments}.
  \emph{Fisheries Research}, \bold{239}:105924 
}
\author{
  \href{mailto:ian.taylor@noaa.gov}{Ian Taylor}, Northwest Fisheries Science Center, NOAA, Seattle WA

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-22}
}
\note{
  In \code{'Report.sso'}, look for \code{'SPR/YPR_Profile report:54'}.
  
  Conversation with PJS 230413: use proportions from area-specific
  \eqn{B_0}{B0} to allocate MSY, etc.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Taylor, pers. comm. Apr 4, 2023}{
  The \code{'$equil_yield'} element of the list returned by \code{'r4ss::SS_output()'}
  is from the table in `Report.sso' with header \code{'SPR/YPR_Profile'} which 
  has more detail than what is brought into R. 
  That table is created by having SS3 calculate various equilibrium quantities
  over a range of \code{'Fmult'} values scaling the catch up and down in order
  to calculate various quantities.
  All of the rows in the table are equilibrium values associated with different
  levels of \code{F}. 
  The iterations are described at the bottom of the table in `Report.sso' as follows:
  \preformatted{%
Finish SPR/YPR profile
#Profile 0 is descending additively from max possible F:  4.00001
#Profile 1 is descending multiplicatively half of max possible F
#Profile 2 is additive back to Fcrash: 3.58975
#value 3 uses endyr F, which has different fleet allocation than benchmark
#value 4 is Fspr: 0.339128
#value 5 is Fbtgt: 0.367703
#value 6 is Fmsy: 0.417663
#Profile 7 increases from Fmsy to Fcrash
#NOTE: meanage_of_catch_is_for_total_catch_of_fleet_type==1_or_bycatch_fleets_with_scaled_Hrate
}
  I think that if you take that table, filter for \code{'SPRloop != 3'} (since 
  profile 3 has different fleet allocation as noted in the pasted output above),
  and then find the row which has the maximum value of the \code{'Tot_Catch'}
  column, you have information associated with MSY.
  Adding up the values for the fleet-specific catches in each area which are 
  reported in the columns with header like \code{'[fleetname]Dead'}, should get
  you the area-specific fractions of the global MSY associated with the 
  specified benchmark fleet allocation and selectivities. 
  This is not the same as trying to maximize yield within an area, which I think
  is a more complex problem and hopefully not what your reference point is based on.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{gatherMCMC}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SS_output}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
