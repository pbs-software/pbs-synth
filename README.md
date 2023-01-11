## PBSsynth: Tools for running Stock Synthesis and visualizing the results ##
&copy; Fisheries and Oceans Canada (2020-2022)

<b>PBSsynth</b> provides an R interface for running ADMB Stock Synthesis (SS) software, including its support packages (<b>r4ss</b>, <b>adnuts</b>). Currently, the repo is a collection of code to facilitate the Offshore Program's use of SS and should be considered extremely preliminary. The added functionality includes automation of:

1. reweighting abundance data (by weighting survey and commercial index CVs) and composition data (by weighting sample sizes <i>N</i> of proportion-at-ages) before choosing an MPD (mode of the posterior distribution) run;
2. launching an MCMC (Monte Carlo Markoff Chain) simulation using the <b>adnuts</b> package;
3. calculating <i>B</i><sub>MSY</sub> (biomass at maximum sustainable yield) and <i>u</i><sub>MSY</sub> (exploitation rate at MSY), and
   4. customising <tt>Sweave</tt> files for individual runs and reweightings from various master `Sweave` files.

The automation offers substantial time-saving when trying numerous model runs.

***** Modified to here.

<b>PBSsynth</b> requires the R packages <b>PBSmodelling</b>, <b>r4ss</b>, <b>adnuts</b>, <b>xtable</b>, <b>lattice</b>, <b>coda</b>, <b>gplots</b>, and <b>Hmisc</b> (all posted on <a href="https://cran.r-project.org/web/packages/">CRAN</a>).

<b>PBSsynth</b> borrows functionality from the <b>r4ss</b> package by adopting the code from a few functions and creating variants. We try to acknowledge the original source wherever possible.

<b>WARNING (check if TRUE):</b> The reliance on so many packages is not ideal because each relies on other packages (e.g., Frank Harrell's <b>Hmisc</b> taps into a cascade of dependencies that rivals the machinations of Wizard Wickham). This means that a user likely has to install various additional R packages, which is a nuisance. Below you can see this maze of cross-dependencies in the <b>scape</b> package. Given enough time, we could extract ourselves from this morass by designing our own functions; however, extrication is slow.

<b>PBSsynth</b> depends on and imports <b>PBSmodelling</b>. Additionally <b>PBSsynth</b> imports the following packages:

<b>methods</b>, <b>lattice</b>, <b>r4ss</b>, <b>adnuts</b>, and <b>xtable</b>.

<i>Needs updating when/if package is built:</i>

<b>PBSsynth</b> imports specific functions from the following:

<!--- importFrom("coda", "mcmc") --->

<!--- importFrom("gplots", "plotCI") --->

importFrom("graphics", "abline", "arrows", "axis", "barplot", "box", "boxplot", "legend", "lines", "mtext", "pairs", "par", "plot", "points", "polygon", "rect", "segments", "text")

importFrom("grDevices", "boxplot.stats", "colorRampPalette", "dev.cur", "dev.list", "dev.off", "extendrange", "pdf", "png", "postscript", "savePlot", "win.metafile", "windows")

importFrom("Hmisc", "Cbind", "panel.xYplot")

importFrom("stats", "acf", "cor", "lowess", "median", "na.pass", "qqnorm", "quantile", "rnorm", "runif", "sd", "var")

importFrom("utils", "Sweave", "data", "flush.console", "help", "installed.packages", "menu", "read.table", "tail", "write.table")

In the end, the target audience is very small, comprising only those people who actually use Stock Synthesis (in the Offshore Program, PBS).

<b>PBSsynth</b> is not available on <a href="https://cran.r-project.org/">CRAN</a> (Comprehensive R Archive Network); however, the source code is available on <a href="https://github.com/pbs-software/pbs-awatea">GitHub</a>. Additionally, both source tarball and Windows binary &#8212; built using R-devel and checked via `R CMD check --as-cran` routine on a <b>Windows 10</b> 64-bit system &#8212; are available on <a href="https://drive.google.com/drive/folders/0B2Bkic2Qu5LGOGx1WkRySVYxNFU?usp=sharing">Google Drive</a>. The best way to use this code is to source it directly when performing stock assessments, so think of the repo as a depot.

As with any freely available product, there is no warranty or promise that **PBSsynth** will perform adequately for all circumstances. Additionally, coding errors are possible, and users should contact the package maintainer if bugs are detected.

Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 


### PBSsynth dependencies ###
 Packages in <b>bold</b> depend on other packages, packages in <i>italics</i> do not. (needs something here, maybe)

