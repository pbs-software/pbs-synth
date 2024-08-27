%% Borrowed from 'r4ss::make_multifig'
\name{make.multifig}
\alias{make.multifig}
\title{Create Multi-figure Plots}
\description{
  Function created by Ian Taylor as an alternative to lattice package for
  multi-figure plots of composition data and fits from Stock Synthesis output.
}
\usage{
make.multifig(ptsx, ptsy, yr, linesx=0, linesy=0, ptsSD=0, 
   sampsize=0, effN=0, showsampsize=TRUE, showeffN=TRUE, 
   sampsize_label="N=", effN_label="effN=", sampsizeround=1, 
   maxrows=6, maxcols=6, rows=1, cols=1, fixdims=TRUE, 
   main="", cex.main=1, xlab="", ylab="", size=1, 
   cexZ1=1.5, bublegend=TRUE, maxsize=NULL, do.sqrt=TRUE, 
   minnbubble=8, allopen=TRUE, xbuffer=c(0.1,0.1), 
   ybuffer=c(0,0.15), yupper=NULL, ymin0=TRUE, xlas=0, ylas=NULL, 
   axis1=NULL, axis2=NULL, axis1labs=NULL, linepos=1, 
   type="o", polygons=TRUE, bars=FALSE, barwidth="default", 
   ptscex=1, ptscol=1, ptscol2=1, 
   colsex = c("orange","limegreen","slategray3"),
   colvec=c(rgb(1,0,0,0.7), rgb(0,0,1,0.7), rgb(0.1,0.1,0.1,0.7)),
   linescol=c(rgb(0,0.8,0,0.7), rgb(1,0,0,0.7), rgb(0,0,1,0.7)),
   lty=1, lwd=2, pch=1, nlegends=3, 
   legtext=list("yr", "sampsize", "effN"), legx="default", 
   legy="default", legadjx="default", legadjy="default", 
   legsize=c(1.2, 1), legfont=c(2, 1), venusmars=TRUE, 
   sampsizeline=FALSE, effNline=FALSE, sampsizemean=NULL, 
   effNmean=NULL, ipage=0, scalebins=FALSE, sexvec=NULL, 
   multifig_colpolygon=c("grey60","grey80","grey70"), 
   multifig_oma=NULL, lang="e", obsN=0, ...) 
}
\arguments{
  \item{ptsx}{\code{numeric} -- vector of x values for points or bars}
  \item{ptsy}{\code{numeric} -- vector of y values for points or bars of same length as ptsx}
  \item{yr}{\code{numeric} -- vector of category values (years) of same length as ptsx}
  \item{linesx}{\code{numeric} -- optional vector of x values for lines}
  \item{linesy}{\code{numeric} -- optional vector of y values for lines}
  \item{ptsSD}{\code{numeric} -- optional vector of standard deviations used to plot error bars
    on top of each point under the assumption of normally distributed error}
  \item{sampsize}{\code{numeric} -- optional sample size vector of same length as ptsx}
  \item{effN}{\code{numeric} -- optional effective sample size vector of same length as ptsx}
  \item{showsampsize}{\code{logical} -- if \code{TRUE}, show sample size values on plot}
  \item{showeffN}{\code{logical} -- if \code{TRUE}, show effective sample size values on plot}
  \item{sampsize_label}{\code{character} -- label on sampsize}
  \item{effN_label}{\code{character} -- label on effN}
  \item{sampsizeround}{\code{numeric} -- rounding level for sample size values}
  \item{maxrows}{\code{numeric} -- maximum (or fixed) number or rows of panels in the plot}
  \item{maxcols}{\code{numeric} -- maximum (or fixed) number or columns of panels in the plot}
  \item{rows}{\code{numeric} -- number or rows to return to as default for next plots
    to come or for single plots}
  \item{cols}{\code{numeric} -- number or cols to return to as default for next plots 
    to come or for single plots}
  \item{fixdims}{\code{logical} -- if \code{TRUE}, fix the dimensions at \code{'maxrows'}
    by \code{'maxcols'} or resize based on number of elements in \code{'yr'} input}
  \item{main}{\code{character} -- title of plot}
  \item{cex.main}{\code{numeric} -- character expansion for title}
  \item{xlab}{\code{character} -- x-axis label}
  \item{ylab}{\code{character} -- y-axis label}
  \item{size}{\code{numeric} -- vector of bubbles sizes if making a bubble plot}
  \item{cexZ1}{\code{numeric} -- character expansion (cex) for point associated with value of 1}
  \item{bublegend}{\code{logical} -- if \code{TRUE}, add legend with example bubble sizes to bubble plots}
  \item{maxsize}{\code{numeric} -- maximum size of bubbles}
  \item{do.sqrt}{\code{logical} -- if \code{TRUE}, scale bubbles based on sqrt of size vector;
    type \code{'?bubble3'} for more info}
  \item{minnbubble}{\code{numeric} -- number of unique x values before adding buffer;
    type \code{'?bubble3'} for more info}
  \item{allopen}{\code{logical} -- if \code{TRUE}, open all bubbles;
    type \code{'?bubble3'} for more info}
  \item{xbuffer}{\code{numeric} -- extra space around points on the left and right
    as fraction of total width of plot}
  \item{ybuffer}{\code{numeric} -- extra space around points on the bottom and top
    as fraction of total height of plot}
  \item{yupper}{\code{numeric} -- upper limit on ymax (applied before addition of ybuffer)}
  \item{ymin0}{\code{logical} -- if \code{TRUE}, fix minimum y-value at 0}
  \item{xlas}{\code{numeric} -- label style (las) input for x-axis; default 0 has horizontal
    labels, value 2 would provide vertical lables}
  \item{ylas}{\code{numeric} -- label style (las) input for y-axis; default NULL has horizontal
    labels when all labels have fewer than 6 characters and vertical otherwise;
    value 0 would force vertical labels, and value 1 would force horizontal}
  \item{axis1}{\code{numeric} -- optional position of bottom axis values}
  \item{axis2}{\code{numeric} -- optional position of left size axis values}
  \item{axis1labs}{\code{character} -- optional vector of labels for \code{'axis1'}
    (either NULL or needs to match length of \code{axis1})}
  \item{linepos}{\code{logical} -- if \code{TRUE}, add lines on top of points (linepos=1)
    or behind (linepos=2); a value of linepos=0 will result in no line}
  \item{type}{\code{character} -- type of line/points used for observed values
    (see \code{'type'} in \code{'?plot'} for details) on top of a grey polygon;
    default is \code{'o'} for overplotting points on lines}
  \item{polygons}{\code{logical} -- if \code{TRUE}, add polygons to the plot;
    (turning off is required for sex-ratio plot)}
  \item{bars}{\code{logical} -- if \code{TRUE}, render the ptsx/ptsy values as bars instead of points;
    \emph{not currently functional}}
  \item{barwidth}{\code{numeric} -- width of bars in barplot; default method
    chooses based on quick and dirty formula also, current method of 
    \code{plot(...type='h')} could be replaced with better approach}
  \item{ptscex}{\code{numeric} -- character expansion factor for points (default=1)}
  \item{ptscol}{\code{character} -- color for points/bars}
  \item{ptscol2}{\code{character} -- color for negative value points in bubble plots}
  \item{colsex}{\code{character} -- colour for females, males, and other (both or unknown)}%% (RH 240822) appears to be redundant
  \item{colvec}{\code{character} -- vector of length 3 with colors for females, males, unsexed fish}
  \item{linescol}{\code{character} -- color for lines}
  \item{lty}{\code{numeric} -- line type}
  \item{lwd}{\code{numeric} -- line width}
  \item{pch}{\code{numeric} -- point character type}
  \item{nlegends}{\code{numeric} -- number of lines of text to add as legends in each plot}
  \item{legtext}{\code{character} -- text in legend, a list of length=nlegends;
    values may be any of \code{'yr'}, \code{'sampsize'}, \code{'effN'},
    or a vector of length = \code{'ptsx'}.}
  \item{legx}{\code{character} -- vector of length=nlegends of x-values 
  of legends (default is first one on left, all after on right)}
  \item{legy}{\code{character} -- vector of length=nlegends of y-values
    of legends (default is top for all plots)}
  \item{legadjx}{\code{numeric} -- left/right adjustment of legends around legx}
  \item{legadjy}{\code{numeric} -- left/right adjustment of legends around legy}
  \item{legsize}{\code{numeric} -- font size for legends; default=c(1.2,1.0)
    (larger for year and normal for others)}
  \item{legfont}{\code{character} -- font type for legends, same as `font' under ?par}
  \item{venusmars}{\code{logical} -- if \code{TRUE}, label females and males with venus and mars symbols}
  \item{sampsizeline}{\code{logical} -- if \code{TRUE}, show line for input sample sizes
    on top of conditional age-at-length plots; (still in development)}
  \item{effNline}{\code{logical} -- if \code{TRUE}, show line for effective sample sizes
    on top of conditional age-at-length plots; (still in development)}
  \item{sampsizemean}{\code{numeric} -- mean input sample size value
    (used when \code{'sampsizeline=TRUE'})}
  \item{effNmean}{\code{numeric} -- mean effective sample size value
    (used when \code{'effNline=TRUE'})}
  \item{ipage}{\code{numeric} -- which page of plots when covering more
    than will fit within \code{'maxrows'} by \code{'maxcols'}}
  \item{scalebins}{\code{logical} -- if \code{TRUE}, rescale expected and 
    observed proportions by dividing by bin width for models where bins have
    different widths; \emph{caution!}: may not work correctly in all cases}
  \item{sexvec}{\code{numeric} -- vector of sex codes if more than one present (otherwise NULL)}
  \item{multifig_colpolygon}{\code{character} -- vector of polygon fill colors
    of length 3 (for females, males, and unsexed fish); can be input to 
    \code{'SS_plots'} and will be passed to this function via the \code{'...'} argument}
  \item{multifig_oma}{\code{numeric} -- vector of outer margins;
    can be input to \code{'SS_plots'} and will be passed to this 
    function via the \code{'...'} argument}
  \item{lang}{\code{character} -- language(s) to use in figures, 
    where \code{'e'} = english and \code{'f'} = french}
  \item{obsN}{\code{numeric} -- observed number (affects the labelling of sample numbers)}
  \item{\dots}{\code{dots} -- additional arguments passed to \code{par}}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Disclosure}{
  A cheeky person \emph{borrowed} the fungible function \code{'make_multifig'} from
  R-package \code{r4ss}, and modified it for a bespoke scientific use.
}
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{plotSS.comps}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{make_multifig}},
  \code{\link[r4ss]{SSplotComps}}
}
\references{
  Taylor, I.G., Doering, K.L., Johnson, K.F., Wetzel, C.R., Stewart, I.J., 2021.
  \href{https://doi.org/10.1016/j.fishres.2021.105924}{Beyond visualizing catch-at-age models: Lessons learned from the r4ss package about software to support stock assessments}.
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
  This is an unwieldy alternative to lattice plots!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{r4ssFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.

