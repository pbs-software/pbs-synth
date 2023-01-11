## PBSsynth: Tools for Manipulating Stock Synthesis 3 ##
&copy; Fisheries and Oceans Canada (2020-2023)

<b>PBSsynth</b> provides an R interface for running NOAA's Stock Synthesis 3 (SS3) software, including its support packages (<b>r4ss</b>, <b>adnuts</b>). Currently, the repo is a collection of code to facilitate the Offshore Program's use of SS3 and should be considered a work in flux. The added functionality includes automation of:

1. reweighting abundance data (by weighting survey and commercial index CVs) and composition data (by weighting sample sizes <i>N</i> of proportion-at-ages) before choosing an MPD (mode of the posterior distribution) run;
2. launching an MCMC (Monte Carlo Markoff Chain) simulation using the <b>adnuts</b> package;
3. calculating <i>B</i><sub>MSY</sub> (biomass at maximum sustainable yield) and <i>u</i><sub>MSY</sub> (exploitation rate at MSY), and
4. customising <tt>Sweave</tt> files for individual runs and reweightings from various master `Sweave` files.

The automation offers substantial time-saving when trying numerous model runs.

<b>PBSsynth</b> requires the R packages <b>PBSmodelling</b>, <b>PBSmapping</b>, <b>PBStools</b>, and <b>PBSdata</b>. 
It also needs the packages <b>r4ss</b>, <b>adnuts</b>, and <b>xtable</b>, although the dependencies have not yet been formalised.

<b>PBSsynth</b> borrows heavily the functionality from the <b>r4ss</b> package by adopting the code from a few functions and creating variants. We try to acknowledge the original source wherever possible.


Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 

