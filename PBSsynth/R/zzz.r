# Taking cue from Roger Bivand's maptools:
.PBSsynthEnv <- new.env(FALSE, parent=globalenv())  # be sure to exportPattern("^\\.PBS") in NAMESPACE

.onAttach <- function(lib, pkg)
{
	pkg_info = utils::sessionInfo( package="PBSsynth" )$otherPkgs$PBSsynth
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()

	userguide_path <- system.file( "doc/PBSsynthIntro.pdf", package = "PBSsynth" )
	year <- substring(date(),nchar(date())-3,nchar(date()))

	packageStartupMessage("
-----------------------------------------------------------
PBS Synthesis ", pkg_info$Version, " -- Copyright (C) 2020-",year," Fisheries and Oceans Canada

Functions to manipulate Stock Synthesis 3 input/output.

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
https://github.com/pbs-software
-----------------------------------------------------------

")
}
.onUnload <- function(libpath) {
	rm(.PBSsynthEnv)
}

# No Visible Bindings
# ===================
#if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
#	".",
#	"A",
#	"B",
#	"C",
#	"D",
#	"E",
#	"F",
#	"G",
#	"H",
#	"I",
#	"J",
#	"K",
#	"L",
#	"M",
#	"N",
#	"O",
#	"P",
#	"Q",
#	"R",
#	"S",
#	"T",
#	"U",
#	"V",
#	"W"
#	"X",
#	"Y",
#	"Z"
#	), package="PBSsynth")

## Force this globally (without have to make a documentation file)
do.call("assign", args=list(x="species.name", value="Canary Rockfish", envir=.PBSmodEnv))
do.call("assign", args=list(x="quants3", value=c(0.05,0.50,0.95), envir=.PBSmodEnv))
do.call("assign", args=list(x="quants5", value=c(0.05,0.25,0.50,0.75,0.95), envir=.PBSmodEnv))
do.call("assign", args=list(x="ptypes", value=c("win","png"), envir=.PBSmodEnv))
do.call("assign", args=list(x="pngres", value=400, envir=.PBSmodEnv))
do.call("assign", args=list(x="lang", value=c("e"), envir=.PBSmodEnv))
do.call("assign", args=list(x="boxpars", value=list(
	boxwex=0.8, boxfill="aliceblue", staplewex=0.5, 
	medlwd=2, medcol=lucent("slategrey",0.75), 
	outwex=0.5, whisklty=1, outpch=3, outcex=0.5, outcol=lucent("grey20",0.5)
	), envir=.PBSmodEnv))

## Flush the cat down the console (change to '.flash.cat' to avoid conflict with function in PBStools)
## Note: `.flush.cat' already in PBStools but this package may not be loaded.
.flash.cat = function(...) { cat(...); flush.console() }

