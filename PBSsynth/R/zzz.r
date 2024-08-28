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
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
	"%>%", ".findSquare",
	"0.2B0", "0.3Gen", "0.4B0", "0.4Bmsy", "0.5B0", "0.5Gen", "0.7B0", "0.8Bmsy",
	"B.index", "B.proj", "B0.MCMC", "Bmsy", "Bsens",
	"EndStrategy",
	"L1",
	"M", "MSY.mcmc", "MSYmaxIter", "MSYtol",
	"N2", "Narea", "Ngen", "Ninit", "Nmat", "NrefM", "Nsens",
	"P.mcmc", "P.mpd", "P.rc",
	"RecFraction", "Rtdev.mcmc",
	"SB0", "SSB", "StartStrategy", "StepStrategy", "StrategyType", 
#	"add_legend", "bubble3", "SSMethod.Cond.TA1.8", "SSMethod.TA1.8", "SS_output", "SS_readctl", "SS_readdat", "SSgetMCMC", "SSgetoutput",  ## r4ss
	"VB", "VB0",
	"acf2", "agile", "alpha", "ampdPA", "area.name", "area.names", "areaPJ", "areaTS", "assYrs", "avgCP",
		"avgLL", "avgPA", "avgPJ", "avgRP", "avgTS",
	"base.dir", "base.lab", "bf", "boxlim", "boxpars", "brm", "bubble3",
	"cattab", "ccex.strip", "central.mpd.dir", "central.run", "condbase", "control.file", "control_file_name",
		"corARMA", "currYear", "currYr",
	"data.compo.figs", "data_file_name", "derposts_file_name", "detectCores",
	"endyr", "exp.run.num", "exp.run.rwt",
	"fft", "fleets.lab", "forecast_file_name", "fval",
	"gear.names", "gen1", "get_model_executable", "ghostagedbase", "ghostlendbase", "gls", "gray", "green", "gseries",
##	"H",
	"intervals", "isnull", "istock",
##	"J",
	"kcol",
	"ladbase", "latex.bold", "lendbase",
	"m", "maxgrad", "mcmc.ts.sub", "modYrs", "mpd", "mutate_at",
	"name", "narea", "ngear", "nmcmc", "npars",
	"osares", "out",
	"pen.lab", "pen.long", "pen.run.num", "pen.run.rwt", "pen.rwt.num", "pgenYear", "posts_file_name", "prefix",
		"prevYear", "proYrs", "projYear",
	"qRlow", "quants3", "quants5",
	"read_admb_cov", "recentCatch", "red", "refCC.sentence", "refHR.sentence", "replist", "resMulti",
		"retro", "rich.colors.short", "rootd_models", "run.num",
	"sample_admb", "sample_inits", "sarima", "scol", "sen.lab", "sen.mcsubs", "senRP", "senTS", 
		"set_prior", "sizebinlabs", "sizedbase", "sizemethod", "slwd", "smpdPA", "so", "splitGear",
		"spp.code", "ss_executable", "startYear", "starter_file_name", "startyr", "stock", "strSpp",
		"strip.columns", "subplot", "survVec", "symbol", "system_",
	"tic", "time", "toc",
	"ucurr", "umsy", "use.num",
	"var", "ver",
	"wadbase", "weight_at_age_file_name",
	"xavgCP", "xavgRP"
##	"Y",
##	"Z"
	), package="PBSsynth")

## Force this globally (without have to make a documentation file)
#do.call("assign", args=list(x="control.file", value="rebsn_control.ss", envir=.PBSmodEnv))
#do.call("assign", args=list(x="species.name", value="Rougheye/Blackspotted Rockfish", envir=.PBSmodEnv))
#do.call("assign", args=list(x="species.name", value="Yellowmouth Rockfish", envir=.PBSmodEnv))
#do.call("assign", args=list(x="species.name", value="Canary Rockfish", envir=.PBSmodEnv))
#do.call("assign", args=list(x="species.name", value="Pacific Ocean Perch", envir=.PBSmodEnv))
do.call("assign", args=list(x="species.name", value="Yellowtail Rockfish", envir=.PBSmodEnv))
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
#.flash.cat = function(...) { cat(...); flush.console() }

