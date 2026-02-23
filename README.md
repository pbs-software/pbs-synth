## PBSsynth: Tools for Manipulating SS3 ##
&copy; Fisheries and Oceans Canada (2020-2026)

<b>PBSsynth</b> provides R-code support for running NOAA's Stock Synthesis 3 (SS3) software, including its complementary packages (<b>r4ss</b>, <b>adnuts</b>). The <b>PBSsynth</b> repo offers a collection of code to facilitate the Offshore Rockfish program's use of SS3, and should be considered a work in flux. The added functionality includes:

1. reweighting abundance data (by weighting survey and commercial index CVs directly in the 'data.ss' file) and composition data (by providing weights in the 'control.ss' file to adjust sample sizes of proportion-at-ages);
2. launching MCMC (Monte Carlo Markoff Chain) simulations using the <b>adnuts</b> package (switched to Chris Grandin's  code on a Linux server);
3. calculating <i>B</i><sub>MSY</sub> (biomass at maximum sustainable yield) and <i>u</i><sub>MSY</sub> (exploitation rate at MSY), and
4. customising <tt>Sweave</tt> files for individual runs and reweightings from various master `Sweave` files.

The semi-automation offers substantial time-saving when trying numerous model runs.

<b>PBSsynth</b> requires the R packages <b>PBSmodelling</b>, <b>PBSmapping</b>, <b>PBStools</b>, and <b>PBSdata</b>. 
It also needs the packages <b>r4ss</b>, <b>adnuts</b>, and <b>xtable</b>, although the dependencies have not yet been formalised.

<b>PBSsynth</b> borrows heavily the functionality from the <b>r4ss</b> package by adopting the code from a few functions and creating variants. We try to acknowledge the original source wherever possible.
<font color="red"><h3>Installation</h3></font>

Although **PBSsynth** is not available on <a href="https://cran.r-project.org/">CRAN</a> (Comprehensive R Archive Network), the source code appears on <a href="https://github.com/pbs-software/pbs-synth">GitHub</a> and can be built in R using:

`devtools::install_github("pbs-software/pbs-synth/PBSsynth")`

The source code has been checked via CRAN's `R CMD check --as-cran` routine using a recent R-devel installation on a **Windows 11** 64-bit system. Some flagrant violations are ignored (like the use of ':::' and 'attach'), but the package is not meant to pass through the pearly gates of CRAN. Consider it a miracle that it works at all.

<font color="red"><h3>Disclaimer</h3></font>

"Fisheries and Oceans Canada (DFO) GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. DFO relinquishes control of the information and assumes no responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against DFO stemming from the use of its GitHub project will be governed by all applicable Canadian Federal laws. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favouring by DFO. The Fisheries and Oceans Canada seal and logo, or the seal and logo of a DFO bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DFO or the Canadian Government.”

As with any freely available product, there is no warranty or promise that **PBSsynth** will perform adequately for all circumstances. Additionally, coding errors are inevitable, and users should contact the package maintainer if bugs are detected. Merci au revoir.


Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 

