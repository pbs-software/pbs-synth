## PBSsynth: Tools for Manipulating SS3 ##
&copy; Fisheries and Oceans Canada (2020-2026)

<b>PBSsynth</b> provides R-code support for running NOAA's Stock Synthesis 3 (SS3) software, including its complementary packages (<b>r4ss</b>, <b>adnuts</b>). Currently, the <b>PBSsynth</b> repo is a collection of code (*not an R package*) to facilitate the Offshore Rockfish program's use of SS3 and should be considered a work in flux. The added functionality includes:

1. reweighting abundance data (by weighting survey and commercial index CVs directly in the 'data.ss' file) and composition data (by providing weights in the 'control.ss' file to adjust sample sizes of proportion-at-ages);
2. launching MCMC (Monte Carlo Markoff Chain) simulations using the <b>adnuts</b> package (switched to Chris Grandin's  code on a Linux server);
3. calculating <i>B</i><sub>MSY</sub> (biomass at maximum sustainable yield) and <i>u</i><sub>MSY</sub> (exploitation rate at MSY), and
4. customising <tt>Sweave</tt> files for individual runs and reweightings from various master `Sweave` files.

The semi-automation offers substantial time-saving when trying numerous model runs.

<b>PBSsynth</b> requires the R packages <b>PBSmodelling</b>, <b>PBSmapping</b>, <b>PBStools</b>, and <b>PBSdata</b>. 
It also needs the packages <b>r4ss</b>, <b>adnuts</b>, and <b>xtable</b>, although the dependencies have not yet been formalised.

<b>PBSsynth</b> borrows heavily the functionality from the <b>r4ss</b> package by adopting the code from a few functions and creating variants. We try to acknowledge the original source wherever possible.


Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 

