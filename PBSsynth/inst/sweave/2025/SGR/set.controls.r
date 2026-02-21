## --------------------------------------------SGR
## Controls for one or more SGR stocks
##  last modified: (RH 250924)
## -----------------------------------------------

## ===============================================
## Create common functions
## ===============================================
allEqual <- function(x) {
	result <- all( x==x[1] )
	result
}
divTwo <- function(x, dig=decdig) {  ## dig is number of dec places
	show0(round(x[1]/x[2],dig),dig)
}
findPos <- function(target, choices) {
	## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
	target = as.numeric(target)
	choices = as.numeric(choices)
	which(abs(choices-target)==min(abs(choices-target)))
}
gather.vectors <- function(run, rwt, ver, use) { ## gather vectors of base run components
	## Initialise objects
	run.num = rwt.num = numeric() ## RH 200508
	ver.num = character()         ## RH 230829
	use.num = logical()           ## RH 200623
	run.sub = 0
	for(k in 1:dim(run)[3]){
		for(i in 1:dim(run)[1]){
			for(j in 1:dim(run)[2]){
				run.sub = run.sub + 1
				if(is.na(run[i,j,k])) next
				else {
					run.tmp =run[i,j,k]; rwt.tmp=rwt[i,j,k]; ver.tmp=ver[i,j,k]; use.tmp=use[i,j,k]
					names(run.tmp) = names(rwt.tmp) = names(ver.tmp) = names(use.tmp) = run.sub
					run.num = c(run.num, run.tmp)
					rwt.num = c(rwt.num, rwt.tmp)
					ver.num = c(ver.num, ver.tmp)
					use.num = c(use.num, use.tmp)
#browser();return()
				}
			}
		}
	}   ## RH 200508|200623 (for subsets of Base runs)
	return(list(run.num=run.num, rwt.num=rwt.num, ver.num=ver.num, use.num=use.num))
}
getActDim <- function(A) { ## get number of active dimensions from the axis of uncertainty matrix
	if (prod(dim(A))==1)
		return(0)
	Nact = 0
	for (i in 1:length(dim(A))) {
#browser();return()
		Apop = apply(A,i,function(j){!is.na(j)})
		if (is.null(dim(Apop)) && sum(Apop)>1)
			return(1)
		if (sum(apply(Apop,2,any))>1) Nact = Nact + 1
	}
	return(Nact)
}
getAssYrs <- function(spp.code, area.name) { ## get modelled current years for previous stock assessments
	if (length(area.name)==1 && grepl("\\.|-|_|\\s+",area.name))
		area.name = strsplit(area.name,split="\\.|-|_|\\s+")[[1]]
	assYrs = NULL

	if (is.element(spp.code,c("SGR"))){ ## Silvergray
		if (any(is.element(area.name,c("CST","BC")))){
			assYrs = c(2014)  ## Modelled current year
		}
	} else if (is.element(spp.code,c("YTR"))){ ## Yellowtail
		if (any(is.element(area.name,c("CST","BC")))){
			assYrs = c(1996, 1997, 2015, 2025)  ## Modelled current year
		}
	} else if (is.element(spp.code,c("POP"))){
		if (any(is.element(area.name,c("CST","BC")))){
			assYrs = c(2001, 2010, 2012, 2017, 2024)
		} else if (any(is.element(area.name,c("5ABC")))){
			assYrs = c(2001, 2010, 2017, 2023)
		} else if (any(is.element(area.name,c("3CD","5DE")))){
			assYrs = c(2012, 2023)
		}
	} else if (is.element(spp.code,c("CAR"))){ ## Canary
		assYrs = c(1999, 2005, 2007, 2009, 2023)
	} else if (is.element(spp.code,c("YMR"))){ ## Yellowmouth
		assYrs = c(2011,2022)
	} else if (is.element(spp.code,c("BOR"))){ ## Bocaccio
		assYrs = c(2008, 2012, 2020, 2022, 2024)
	} else if (is.element(spp.code,c("WWR"))){ ## Widow
		assYrs = c(2019)
	} else if (is.element(spp.code,c("RSR"))){ ## Redstripe
		assYrs = c(2010, 2018)
	} else if (is.element(spp.code,c("WAP"))){ ## Walleye Pollock
		assYrs = c(2017)
	} else if (is.element(spp.code,c("SBF"))){  ## Sablefish
		assYrs = c(2016)
	} else if (is.element(spp.code,c("SGR","SST","YYR"))){
		assYrs = c(2015)
	} else if (is.element(spp.code,c("ARF","RBR"))){
		assYrs = c(2014)
	} else if (is.element(spp.code,c("ROL"))){
		assYrs = c(2013)
	} else {
		assYrs = NULL
	}
#browser();return()
	return(assYrs)
}
getMPD <- function(obj) {
	unpackList(obj$extra$general)
	unpackList(obj$extra$parameters)
	mpd = list()
	mpd[["R_0"]] = R0
	mpd[["R_avg"]] = avgR0
	mpd[["h"]]   = h
	for (i in 1:Nsexes) {
		mpd[[paste0("M_",i)]] = M1[i]
	}
	for (i in 1:Nsexes) {
		mpd[[paste0("M2_",i)]] = M2[i]
	}
	for (i in 1:Nsurveyindex) {
		mpd[[paste0("q_",i)]] = log_qsurvey[i]
	}
	for (j in 1:NCPUEindex) {
		jj = j + Nsurveyindex
		mpd[[paste0("q_",jj)]] = log_qCPUE[j]
	}
	for (i in 1:Nsurveyindex){
		mpd[[paste0("mu_",i)]]        = surveySfull[i]
		mpd[[paste0("Delta_",i)]]     = survey_SfullDelta[i]
		mpd[[paste0("log v_",i,"L")]] = log_surveyvarL[i]
		mpd[[paste0("log v_",i,"R")]] = log_surveyvarR[i]
	}
	for (j in 1:Nmethods) {
		jj = j + Nsurveyindex
		mpd[[paste0("mu_",jj)]]        = Sfullest[j]
		mpd[[paste0("Delta_",jj)]]     = SfullDelta[j]
		mpd[[paste0("log v_",jj,"L")]] = log_varLest[j]
		mpd[[paste0("log v_",jj,"R")]] = log_varRest[j]
	}
	mpd[["sigmaR"]] = obj$extra$residuals$p_log_RecDev[6]
	names(log_RecDev) = StartYear:EndYear
	mpd[["log R_dev"]] = log_RecDev
	return(mpd)
}
medCI <- function(x, dig=decdig, CI=quants3[c(1,3)]) { ## dig is number of dec places, CI = credible interval
	print(paste0(c(
	prettyNum(round(quantile(x, 0.50, na.rm=T), digits=dig), big.mark=","), "~(",
	prettyNum(round(quantile(x, CI[1], na.rm=T), digits=dig), big.mark=","), ",~",
	prettyNum(round(quantile(x, CI[2], na.rm=T), digits=dig), big.mark=","), ")"), collapse=""))
}
relabelTex <- function(texinput, prefix, caption) { ## (RH 190513)
	type = class(texinput)[1]
	if (type=="character")   tfile = readLines(paste0(texinput,".tex"))
	else if (type=="xtable"){tfile = attributes(texinput)$label; cfile=attributes(texinput)$caption}
	else                     tfile = texinput
	tfile = gsub("\\{tab:",paste0("{tab:", prefix), tfile)
	tfile = gsub("\\{fig:",paste0("{fig:", prefix), tfile)

	## Identify figure lines and slap an extra two arguments on the end
	isfig = grep("[one|two|three]fig", tfile)
	tfile[isfig] = paste0(tfile[isfig],"{",caption,"}{",prefix,"}")

	## Change the original call to a figure function to the new figure function (now use generic functions)
	#tfile = gsub("fig\\{",paste0("fig",substring(prefix,1,1),"{"), tfile)
	#tfile = gsub("figH\\{",paste0("figH",substring(prefix,1,1),"{"), tfile)
	#tfile = gsub("figWH\\{",paste0("figWH",substring(prefix,1,1),"{"), tfile)

	if (!missing(caption)) {
		if (type=="xtable")
			cfile = paste0(caption,cfile)
		else
			tfile = gsub("caption\\{",paste0("caption{",caption),tfile)
	}
	if (type=="character")
		writeLines(tfile,paste0(texinput,".relab.tex"))
	else if (type=="xtable") {
		tfile = gsub("^tab:",paste0("tab:", prefix), tfile)
		attr(texinput,"label")   = tfile
		attr(texinput,"caption") = cfile
		return (texinput)
	} else
		return(invisible(tfile))
}
## ===========================end common functions

## Start the main computations
## ===========================
#rm(list=setdiff(ls(all=T),c(".First","so","Rcode","qu","Scode")))  ## Start with clean slate
require(PBStools)
require(xtable)
require(r4ss)
options(scipen=10)

#d.awatea = "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/"
#R.awatea =c("plotFuns")
#for (i in R.awatea) source(paste0(d.awatea,"current/",i,".r"))

d.tools = "C:/Users/haighr/Files/Projects/R/Develop/PBStools/Authors/Rcode/"
r.tools = c("linguaFranca", "createFdir", "clearFiles", "extractAges", "calcStockArea", "compBmsy", "formatCatch", "texArray", "plotSnail")
scenario=0 ## to source compBmsy without running anything
for (i in r.tools) source(paste0(d.tools,"develop/",i,".r"))

d.synth = "C:/Users/haighr/Files/Projects/R/Develop/PBSsynth/Authors/Rcode/"
#R.synth = c("plotFuns","utilFuns","smallFuns","r4ssFuns","zzz.local")
R.synth = c("smallFuns","zzz.local") ## need these small functions and settings
r.synth = c("plotSS.pars", "plotSS.ts", "plotSS.index", "plotSS.francis", "plt.selectivity", "plotSS.selex", "make.multifig", "plotSS.comps", "plt.cohortResids", "plotSS.stdres", "plotSS.rdevs", "calcStdRes", "ptab", "convPN", "weightAF", "plotSS.pmcmc", "tabDQs", "mochaLatte", "panelBoxes", "plotTraj", "panelTraces", "panelChains", "plotACFs")
for (i in R.synth) source(paste0(d.synth,"current/",i,".r"))  ## this overwrites more recent changes (e.g., mergePA) so be careful
for (i in r.synth) source(paste0(d.synth,"develop/",i,".r"))

## For the Sweave, choose only one language even though some plots will be made in both.
lang = "e"; tput(lang)  ## 'e'=english, 'f'=francais
## Create a subdirectory called `french' for French-language figures if 'f' in lang
createFdir(lang)

appG.dir = "C:/Users/haighr/Files/GFish/PSARC25/SGR/Docs/WP/AppG_Results" #getwd()
##----------------------------
## Stock control list object
##----------------------------
stock     = list()
stolab    = c("3area", "coast")
stospp    = c("SGR", "SGR")  ## need this up front
assyr     = 2025
Nstock    = length(stolab)
saveDate  = sub("/0","/",sub("^0","",format.Date(Sys.time(),"%m/%d/%Y")))

## Axes of uncertainty is a ragged hot mess
## ----------------------------------------
if (any(stolab %in% c("BSR","RER")) && assyr==2020) {
	axis.dim = c(3,3,2)
	axisA.dimnames = list(M=c(0.035,0.045,0.055), cp=c(0.1,0.2759,0.4), AE=c(3,5))
	axisB.dimnames = list(M=c(0.035,0.045,0.055), cp=c(0.1,0.2529,0.4), AE=c(3,5))
	runA = rwtA = array(NA, dim=axis.dim, dimnames=axisA.dimnames)  ## empty array
	runB = rwtB = array(NA, dim=axis.dim, dimnames=axisB.dimnames)  ## empty array
	useA = array(FALSE, dim=axis.dim, dimnames=axisA.dimnames) ## logical array set to FALSE
	useB = array(FALSE, dim=axis.dim, dimnames=axisB.dimnames) ## logical array set to FALSE
	
	## Populate A arrays for first stock
	runA["0.035",,"3"] = 49:51; runA["0.045",,"3"] = c(47,46,48); runA["0.055",,"3"] = 52:54
	rwtA[,,"3"] = 1
	useA[,,"3"] = TRUE
	
	## Populate B arrays for second stock
	runB["0.035",,"3"] = c(18,12,15); runB["0.045",,"3"] = c(17,11,14); runB["0.055",,"3"] = c(19,13,16)
	rwtB[,,"3"] = 2; rwtB["0.055","0.1","3"] = 1
	useB["0.035",,"3"] = rep(TRUE,3); useB["0.045",,"3"] = c(F,T,T); useB["0.055",,"3"] = c(F,F,T) ## RPR subset
	
	#runB["0.045",,"3"] = c(NA,11,14) ## RPR subset
	#runB["0.055",,"3"] = c(NA,NA,16) ## RPR subset
	#runB["0.055",,"3"] = c(19,NA,NA) #,13,16)  ## remove B8 and B9

#	rwt.numA = rep(1,9)[as.numeric(names(run.numA))]           ## RH 200508 (for subsets of Base runs)
#	rwt.numB = c(rep(2,6),1,2,2)[as.numeric(names(run.numB))]  ## RH 200508 (for subsets of Base runs)
#	rwt.numB = rep(2,6)[as.numeric(names(run.numB))]  ## RH 200508 (for subsets of Base runs)

	A.base = gather.vectors(runA, rwtA, useA)
	B.base = gather.vectors(runB, rwtB, useB)
}
if (any(stolab %in% c("YMR")) && assyr==2021) {
	axis.dim = c(5,1,1)
	axisA.dimnames = list(M=c(0.04,0.045,0.05,0.055,0.06), cp=c(0.3296), RM=c(2))  ## RM = Reweight method: 1=Francis mean age, 2=Harmonic mean ratio
	runA = rwtA = array(NA, dim=axis.dim, dimnames=axisA.dimnames)  ## empty array
	useA = array(FALSE, dim=axis.dim, dimnames=axisA.dimnames) ## logical array set to FALSE
	
	## Populate A arrays for first stock
	runA[,1,1] = c(77,71,75,72,76)
	rwtA[,1,1] = 1
	useA[,1,1] = TRUE
	A.base = gather.vectors(runA, rwtA, useA)
}
if (any(stolab %in% c("CAR")) && assyr==2022) {
	axis.dim = c(1,1,1)
	axisA.dimnames = list(M="est", cp=c(0.178), RM=c(0))  ## RM = Reweight method: 0=None/D-M, 1=Francis mean age, 2=Harmonic mean ratio
	runA = rwtA = array(NA, dim=axis.dim, dimnames=axisA.dimnames)  ## empty array
	useA = array(FALSE, dim=axis.dim, dimnames=axisA.dimnames) ## logical array set to FALSE
	
	## Populate A arrays for first stock
	runA[1,1,1] = c(24)
	rwtA[1,1,1] = 1
	useA[1,1,1] = TRUE
	A.base = gather.vectors(runA, rwtA, useA)
}
if (any(stolab %in% c("POP")) && assyr==2023) {
	axis.dim = c(1,1,1)
	#axisA.dimnames = list(M="est", cp=c(0), RM=c(0))  ## RM = Reweight method: 0=None/D-M, 1=Francis mean age, 2=Harmonic mean ratio
	axisA.dimnames = list(M="est", cp=c(0), RM=c("D-M"))  ## RM = Reweight method: 'NO'=None, 'D-M'=Dirichlet-Multinomial, 'FMA'=Francis mean age, 'HMR'=Harmonic mean ratio ## RH 230829
	runA = rwtA = verA = array(NA, dim=axis.dim, dimnames=axisA.dimnames)  ## empty array
	useA = array(FALSE, dim=axis.dim, dimnames=axisA.dimnames)      ## logical array set to FALSE
	
	## Populate A arrays for first stock
	runA[1,1,1] = c(21)
	rwtA[1,1,1] = 1
	verA[1,1,1] = "3a"
	useA[1,1,1] = TRUE
	A.base = gather.vectors(runA, rwtA, verA, useA)
	B.base = A.base  ## there is no second stock
}
if (any(stolab %in% c("YTR")) && assyr==2024) {
	axis.dim = c(1,1,1)
	#axisA.dimnames = list(M="est", cp=c(0), RM=c(0))  ## RM = Reweight method: 0=None/D-M, 1=Francis mean age, 2=Harmonic mean ratio
	axisA.dimnames = list(M="est", cp=c(0), RM=c("FMA"))  ## RM = Reweight method: 'NO'=None, 'D-M'=Dirichlet-Multinomial, 'FMA'=Francis mean age, 'HMR'=Harmonic mean ratio ## RH 230829
	runA = rwtA = verA = array(NA, dim=axis.dim, dimnames=axisA.dimnames)  ## empty array
	useA = array(FALSE, dim=axis.dim, dimnames=axisA.dimnames)      ## logical array set to FALSE
	
	## Populate A arrays for first stock
	runA[1,1,1] = c(2)
	rwtA[1,1,1] = 1
	verA[1,1,1] = "1c"
	useA[1,1,1] = TRUE
	A.base = gather.vectors(runA, rwtA, verA, useA)
	B.base = A.base  ## there is no second stock
}
if (any(grepl("SGR",stospp)) && assyr==2025) {
	axis.dim = c(1,1,1)  ## no axes of uncertainty
	axisA.dimnames = list(M="3area", cp=c(0), RM=c("FMA"))  ## RM = Reweight method: 'NO'=None, 'D-M'=Dirichlet-Multinomial, 'FMA'=Francis mean age, 'HMR'=Harmonic mean ratio ## RH 230829
	axisB.dimnames = list(M="coast", cp=c(0), RM=c("FMA"))  ## RM = Reweight method: 'NO'=None, 'D-M'=Dirichlet-Multinomial, 'FMA'=Francis mean age, 'HMR'=Harmonic mean ratio ## RH 230829
	runA = rwtA = verA = array(NA, dim=axis.dim, dimnames=axisA.dimnames)  ## empty array
	runB = rwtB = verB = array(NA, dim=axis.dim, dimnames=axisB.dimnames)  ## empty array
	useA = array(FALSE, dim=axis.dim, dimnames=axisA.dimnames)      ## logical array set to FALSE
	useB = array(FALSE, dim=axis.dim, dimnames=axisB.dimnames)      ## logical array set to FALSE
	
	## Populate A arrays for first stock
	runA[1,1,1] = c(28)
	rwtA[1,1,1] = 1
	verA[1,1,1] = "2c"
	useA[1,1,1] = TRUE
	## Populate B arrays for second stock
	runB[1,1,1] = c(29)
	rwtB[1,1,1] = 1
	verB[1,1,1] = "2c"
	useB[1,1,1] = TRUE
	A.base = gather.vectors(runA, rwtA, verA, useA)
	B.base = gather.vectors(runB, rwtB, verB, useB)
}

for (i in 1:length(stolab)) {
	##----------------------------
	## Stock control list object (keep here temporarily until code for SS evolves)
	##----------------------------
	ii = stolab[i]
	jj = "Controls"
	stock[[ii]] = list()
	stock[[ii]][[jj]] = list()
	stock[[ii]][[jj]][["strSpp"]]        = switch(i, "405", "405")
	stock[[ii]][[jj]][["name"]]          = switch(i, paste0("SGR ",assyr," (3area)"), paste0("SGR ",assyr," (coast)")) ## Don't specify 3 stocks if using SS3's multi-area model
	stock[[ii]][[jj]][["spp.code"]]      = stospp[i]  ## defined above to use in conditions
	stock[[ii]][[jj]][["spp.name"]]      = "Silvergray Rockfish"
	stock[[ii]][[jj]][["latin.name"]]    = "Sebastes brevispinis"
	stock[[ii]][[jj]][["prefix"]]        = switch(i, "sgr.3area.", "sgr.coast.")
	stock[[ii]][[jj]][["model.dir"]]     = "C:/Users/haighr/Files/GFish/PSARC25/SGR/Data/SS3"
	stock[[ii]][[jj]][["base.dir"]]      = paste0(stock[[ii]][[jj]][["model.dir"]], switch(i, "/SGR2025","/SGR2025")) ## Top level directory for runs
	stock[[ii]][[jj]][["RPbase"]]        = "B0"   ## as oppposed to 'BMSY' or 'BHIS'
	stock[[ii]][[jj]][["virginYear"]]    = 1933
	stock[[ii]][[jj]][["startYear"]]     = 1935
	stock[[ii]][[jj]][["currYear"]]      = 2026
	stock[[ii]][[jj]][["prevYear"]]      = 2025
	stock[[ii]][[jj]][["projYear"]]      = 2036  ## 10-yr projection
	stock[[ii]][[jj]][["pgenYear"]]      = 3*25  ## 3G (G=25y)
	stock[[ii]][[jj]][["gen1"]]          = 25
	stock[[ii]][[jj]][["Ngen"]]          = 3
	stock[[ii]][[jj]][["Nsex"]]          = 2
	stock[[ii]][[jj]][["Nsurv"]]         = switch(i, 6, 6)
	stock[[ii]][[jj]][["iseries"]]       = c("QCS Synoptic", "WCHG Synoptic", "WCVI Synoptic", "HS Synoptic", "GIG Historical", "NMFS Triennial")  ## index series
	stock[[ii]][[jj]][["Ncpue"]]         = switch(i, 3, 3)
	stock[[ii]][[jj]][["useries"]]       = c("5ABC Fishery","5DE Fishery","3CD Fishery")  ## cpue series
	stock[[ii]][[jj]][["Ngear"]]         = switch(i, 3, 3)
	stock[[ii]][[jj]][["gseries"]]       = c("5ABC Trawl", "5DE Trawl", "3CD Trawl")  ## gear series
	stock[[ii]][[jj]][["Narea"]]         = switch(i, 1, 1)
	stock[[ii]][[jj]][["area.name"]]    = switch(i, "BC", "BC")
	stock[[ii]][[jj]][["Nsubarea"]]      = switch(i, 3, 1)
	stock[[ii]][[jj]][["area.names"]]    = switch(i, c("5ABC","5DE","3CD"), "BC")
	#stock[[ii]][[jj]][["areas"]]         = switch(i, c("5ABC","5DE","3CD"), "BC")
	stock[[ii]][[jj]][["Nfleet"]]        = switch(i, 9, 9)
	stock[[ii]][[jj]][["fleets"]]        = c("5ABC Trawl", "5DE Trawl", "3CD Trawl", "QCS Synoptic", "WCHG Synoptic", "WCVI Synoptic", "HS Synoptic", "GIG Historical", "NMFS Triennial", "HBLL North", "HBLL South")  ## fleet names
	stock[[ii]][[jj]][["mcsub"]]         = switch(i, rep(list(c(1,2000)),1), rep(list(c(1,2000)),1))
	stock[[ii]][[jj]][["run.num"]]       = switch(i, A.base$run.num, B.base$run.num)
	stock[[ii]][[jj]][["rwt.num"]]       = switch(i, A.base$rwt.num, B.base$rwt.num)
	stock[[ii]][[jj]][["ver.num"]]       = switch(i, A.base$ver.num, B.base$ver.num)
	stock[[ii]][[jj]][["use.num"]]       = switch(i, A.base$use.num, B.base$use.num)
	stock[[ii]][[jj]][["exp.run.num"]]   = switch(i, 28, 29) ## central run (originally called example run)
	stock[[ii]][[jj]][["exp.run.rwt"]]   = switch(i, "28.01.v2c", "29.01.v2c")
	stock[[ii]][[jj]][["recentCatch"]]   = switch(i, list(  ## RH 230829: change to a list object for multiple areas
		'5ABC'=c('2020'= 960, '2021'= 844, '2022'= 926, '2023'=1099, '2024'=1322), 
		'5DE' =c('2020'= 184, '2021'= 237, '2022'= 203, '2023'= 167, '2024'= 214),
		'3CD' =c('2020'= 334, '2021'= 249, '2022'= 265, '2023'= 345, '2024'= 312)
		), list(  ## not sure yet if coastwide uses this 
		'BC'  =c('2020'=1479, '2021'=1330, '2022'=1394, '2023'=1611, '2024'=1847) )
		) 
	## Axes of uncertainty (set to NULL if none)
	stock[[ii]][[jj]][["axisU"]]         = switch(i, runA, runB)
	stock[[ii]][[jj]][["base.lab"]]      = switch(i,
		paste0("B",names(A.base$run.num),": R",pad0(A.base$run.num,2), paste0(" (M=",rownames(runA),")")), 
		paste0("B",names(B.base$run.num),": R",pad0(B.base$run.num,2), paste0(" (M=",rownames(runB),")")) 
		)
	stock[[ii]][[jj]][["base.long"]]     = switch(i, "three-area model with three subareas", "coastwide model with three subareas" )
	stock[[ii]][[jj]][["Nsens"]]         = switch(i, 13, 11)
	stock[[ii]][[jj]][["sen.run.num"]]   = switch(i, c(32,21,30,seq(34,46,2),50,48:49), c(33,23,31,seq(35,47,2),51) )
	stock[[ii]][[jj]][["sen.rwt.num"]]   = switch(i, c(rep(1,13)), c(rep(1,11)) )  ## even D-M is reweighted once using CPUE cvpro
	stock[[ii]][[jj]][["sen.ver.num"]]   = switch(i, 
		paste0("v",c(1,3,2, rep(1,10)), c(letters[c(rep(2,5),rep(1,8))] ) ), 
		paste0("v",c(1,3,2, rep(1,8)), c(letters[c(rep(2,5),rep(1,6))] ) )
		)
	stock[[ii]][[jj]][["sen.run.rwt"]]   = paste(pad0(stock[[ii]][[jj]][["sen.run.num"]],2), pad0(stock[[ii]][[jj]][["sen.rwt.num"]],2), stock[[ii]][[jj]][["sen.ver.num"]], sep=".")
	stock[[ii]][[jj]][["sen.lab"]]       = c("drop_CPUE", "swept_area", "add_HBLL", "sigmaR=0.6", "sigmaR=1.2", "dirichlet", "no_age_error", "reduce_catch", "increase_catch", "drop_GIG_NMFS", "loosen_M_prior", if(i==1) c("fix_5ABC_Rdist", "fix_5DE_Rdist") )  ## almost the same for both stocks
	stock[[ii]][[jj]][["sen.long"]]      = c("drop CPUE index series", "use design-based series", "use HBLL surveys (North \\\\& South)", "reduce $\\\\sigma_R$ to 0.6", "increase $\\\\sigma_R$ to 1.2", "use Dirichlet-Mutinomial parameterisation", "apply no ageing error", "reduce commercial catch (1965-95) by 30\\\\pc{}", "increase commercial catch (1965-95) by 50\\\\pc{}", "drop GIG \\\\& NMFS surveys", "loosen M prior means: N(0.065,0.065)", if(i==1) c("fix 5ABC recruitment distribution parameter", "fix 5DE recruitment distribution parameter") )   ## almost the same for both stocks
	stock[[ii]][[jj]][["sseries"]]       = c( rep(list(c(4:9)),2), list(c(4:11)), rep(list(c(4:9)),6), list(c(4:7)), list(c(4:9)), if(i==1) rep(list(c(4:9)),2) )  ## survey series for comparing sens to base parameters (RH 240627: not sure if this is used anywhere, keep an eye open)
	## Add PJS sensitivities
	#stock[[ii]][[jj]][["pen.run.num"]]   = switch(i, c(24:26), c(0) )
	#stock[[ii]][[jj]][["pen.rwt.num"]]   = switch(i, rep(1,3), c(0) )
	#stock[[ii]][[jj]][["pen.run.rwt"]]   = switch(i, c(paste0(24:26,".01.v1a")),  "")
	#stock[[ii]][[jj]][["pen.lab"]]       = switch(i,
	#	c("single-area 5ABC", "single-area 3CD", "single-area 5DE"),
	#	c("lemony persnickets"))
	#stock[[ii]][[jj]][["pen.long"]]      = switch(i,
	#	c("single-area model for 5ABC", "single-area model for 3CD", "single-area model for 5DE"),
	#	c("lemony persnickets"))
	## stock[[ii]][[jj]][["dia.run.nums"]]= switch(i, c(23,24), c(24,25))
	## stock[[ii]][[jj]][["dia.rwt.nums"]]= switch(i, c(1,2),   c(3,3))
	## stock[[ii]][[jj]][["dia.lab"]]     = switch(i, c("Base + Fix_M","S4 + Fix_M"), c("Base + Fix_M","S4 + Fix_M"))
	stock[[ii]][[jj]][["sen.mcsubs"]]    = switch(i, rep(list(c(1,2000)),12), rep(list(c(1,2000)),10))
	#stock[[ii]][[jj]][["policies"]]      = switch(i, c(CC=800,HR=0.1), c(CC=300,HR=0.06))  ## set policies here but there could be a disconnect in policies determined by 'make.base.figs.r' depending on what policies are available.
}
##----------------------------
#browser();return()

##strSpp="418"; assyr=2024; species.code="YTR"; species.name="Yellowtail Rockfish"  ## also put this line in 'load.preview.r'

## Commonalities (use unique() rather than .su() because don't want alphabetic sorting)
## NOTE: probably only works for one species with multiple stocks
strSpp     = unique(sapply(stock,function(x){x$Controls$strSpp}))
spp.code   = unique(sapply(stock,function(x){x$Controls$spp.code})); species.code = spp.code
spp.name   = unique(sapply(stock,function(x){x$Controls$spp.name})); species.name = spp.name
latin.name = unique(sapply(stock,function(x){x$Controls$latin.name}))
spp.year   = paste(spp.name, assyr)
the.stocks = gsub(" ","~",unique(sapply(stock,function(x){x$Controls$name})))
ptype      = "png" ## make sure this is the same choice on line 19
pngres     = 400   ## pixels per inch
sigdig     = 3     ## Number of significant digits to output in tables
decdig     = 2     ## Number of decimal digits (e.g., may only want 2 for decision tables)
## Create quants and put into .PBSmodEnv just in case a PBSawatea plotting function might need either.
quants3    = c(0.05,0.50,0.95);           tput(quants3)
quants5    = c(0.05,0.25,0.50,0.75,0.95); tput(quants5)

## Unpack the years to get Sweave started in the wrapper Rnw
unpackList(stock[[1]][["Controls"]][c("startYear","currYear","prevYear","projYear","pgenYear")])


## Grab findSquare from PBSmodelling
#.findSquare <- PBSmodelling:::.findSquare



