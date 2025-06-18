## -----------------------------------------------
## Controls for one or more stocks (RH 2107263)
## Revised from Awatea for SS
## -----------------------------------------------
#rm(list=setdiff(ls(all=T),c(".First","so","Rcode","qu","Scode")))  ## Start with clean slate

require(PBStools)
require(xtable)
require(r4ss)
options(scipen=10)

d.awatea = "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop/"
r.awatea = c("plotFuns", "mochaLatte", "panelBoxes","plotTraj","panelTraces","panelChains","plotACFs")
for (i in r.awatea) source(paste0(d.awatea,i,".r"))

d.tools = "C:/Users/haighr/Files/Projects/R/Develop/PBStools/Authors/Rcode/develop/"
r.tools = c("linguaFranca", "clearFiles", "extractAges", "calcStockArea", "compBmsy", "formatCatch", "texArray", "plotSnail")
scenario=0 ## to source compBmsy without running anything
for (i in r.tools) source(paste0(d.tools,i,".r"))

d.synth = "C:/Users/haighr/Files/Projects/R/Develop/PBSsynth/Authors/Rcode/develop/YMR2021/"
r.synth = c("PBSsynth","plotFuns","utilFuns","plotSS.pars","plotSS.ts","plotSS.index","plotSS.francis", "plt.selectivity", "plotSS.selex", "make.multifig", "plotSS.comps", "plt.cohortResids", "plotSS.stdres", "plotSS.rdevs", "calcStdRes", "ptab", "convPN", "weightAF")
for (i in r.synth) source(paste0(d.synth,i,".r"))

## For the Sweave, choose only one language even though some plots will be made in both.
lang = "e"; tput(lang)  ## 'e'=english, 'f'=francais
## Create a subdirectory called `french' for French-language figures if 'f' in lang
createFdir(lang)

appF.dir = getwd()
##----------------------------
## Stock control list object
##----------------------------
stock     = list()
stolab    = c("YMR")
Nstock    = length(stolab)

## Gather vectors of base run components
gather.vectors = function(run, rwt, use) {
	## Initialise objects
	run.num = rwt.num = numeric() ## RH 200508
	use.num = logical() ## RH 200623
	run.sub = 0
	for(k in 1:dim(run)[3]){
		for(i in 1:dim(run)[1]){
			for(j in 1:dim(run)[2]){
				run.sub = run.sub + 1
				if(is.na(run[i,j,k])) next
				else {
					run.tmp =run[i,j,k]; rwt.tmp=rwt[i,j,k]; use.tmp=use[i,j,k]
					names(run.tmp) = names(rwt.tmp) = names(use.tmp) = run.sub
					run.num = c(run.num, run.tmp)
					rwt.num = c(rwt.num, rwt.tmp)
					use.num = c(use.num, use.tmp)
				}
			}
		}
	}   ## RH 200508|200623 (for subsets of Base runs)
	return(list(run.num=run.num, rwt.num=rwt.num, use.num=use.num))
}
## Axes of uncertainty is a ragged hot mess
if (any(stolab %in% c("BSR","RER"))) {
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
if (any(stolab %in% c("YMR"))) {
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

for (i in 1:length(stolab)) {
	##----------------------------
	## Stock control list object (keep here temporarily until code for SS evolves)
	##----------------------------
	ii = stolab[i]
	jj = "Controls"
	stock[[ii]] = list()
	stock[[ii]][[jj]] = list()
	stock[[ii]][[jj]][["name"]]            = switch(i, "YMR~BC", "")
	stock[[ii]][[jj]][["spp.code"]]        = "YMR"
	stock[[ii]][[jj]][["spp.name"]]        = "Yellowmouth Rockfish"
	stock[[ii]][[jj]][["latin.name"]]      = "Sebastes reedi"
	stock[[ii]][[jj]][["area.name"]]       = switch(i, "CST", "")
	stock[[ii]][[jj]][["prefix"]]          = switch(i, "ymr.","xxx.") ## use istock
	stock[[ii]][[jj]][["model.dir"]]       = "C:/Users/haighr/Files/GFish/PSARC21/YMR/Data/SS"
	stock[[ii]][[jj]][["base.dir"]]        = paste0(stock[[ii]][[jj]][["model.dir"]], switch(i, "/YMR2021","/nada")) ## Top level directory for runs
	stock[[ii]][[jj]][["startYear"]]       = 1935
	stock[[ii]][[jj]][["currYear"]]        = 2022
	stock[[ii]][[jj]][["prevYear"]]        = 2021
	stock[[ii]][[jj]][["projYear"]]        = 2032  ## 10-yr projection
	stock[[ii]][[jj]][["pgenYear"]]        = 2112  ## 90-yr projection  ## 3G (G=30y)
	stock[[ii]][[jj]][["gen1"]]            = 30
	stock[[ii]][[jj]][["Ngen"]]            = 3
	stock[[ii]][[jj]][["Nsex"]]            = 2
	stock[[ii]][[jj]][["Nsurv"]]           = switch(i, 4, 1)
	stock[[ii]][[jj]][["Ncpue"]]           = switch(i, 1, 1)
	stock[[ii]][[jj]][["Ngear"]]           = switch(i, 1, 1)
	stock[[ii]][[jj]][["mcsub"]]           = switch(i, rep(list(c(1,2000)),5), rep(list(c(1,2000)),1))
	stock[[ii]][[jj]][["run.num"]]         = switch(i, A.base$run.num, B.base$run.num)
	stock[[ii]][[jj]][["rwt.num"]]         = switch(i, A.base$rwt.num, B.base$rwt.num)
	stock[[ii]][[jj]][["use.num"]]         = switch(i, A.base$use.num, B.base$use.num)
	stock[[ii]][[jj]][["exp.run.rwt"]]     = switch(i, "75.01", "")
	stock[[ii]][[jj]][["recentCatch"]]     = switch(i, c('2016'=1162.4,'2017'=1404.2,'2018'=1216.3,'2019'=1520.7,'2020'=1056.9), 0) ## for now just provide the data
	## Axes of uncertainty (set to NULL if none)
	stock[[ii]][[jj]][["axisU"]]           = switch(i, runA, runB)
	stock[[ii]][[jj]][["iseries"]]         = switch(i, c("Bottom Trawl CPUE", "QCS Synoptic", "WCVI Synoptic", "WCHG Synoptic", "GIG Historical"), c( "nada") )  ## index series
	stock[[ii]][[jj]][["gseries"]]         = switch(i, c("Trawl+"), c("nada"))  ## gear series
	stock[[ii]][[jj]][["base.lab"]]        = switch(i, paste0("B",names(A.base$run.num),": R",A.base$run.num, paste0(" (M=",rownames(runA),")")), paste0("B",names(B.base$run.num)) )
	stock[[ii]][[jj]][["Nsens"]]           = switch(i, 14, 0)
	stock[[ii]][[jj]][["sseries"]]         = switch(i, c(rep(list(1:5),2),list(2:5),rep(list(1:5),11)), rep(list(1:5),14) )  ## survey series for comparing sens to base parameters
	stock[[ii]][[jj]][["sen.run.num"]]     = switch(i, c(78:88,91:93), c(0) )
	stock[[ii]][[jj]][["sen.rwt.num"]]     = switch(i, c(rep(1,14)), c(0) )
	stock[[ii]][[jj]][["sen.lab"]]         = switch(i,
		c("add_1997_WCHG_index", "estimate_M", "drop_CPUE", "Tweedie_CPUE", "sigmaR=0.6", "sigmaR=1.2", "reduce_catch_33%", "increase_catch_50%", "upweight_QCS_AF", "start_Rdevs_in_1970", "no_ageing_error", "steepness_h=0.5", "double_2021_catch", "AE_from_age_readers"),
		c("lemony persnickets"))
	## stock[[ii]][[jj]][["dia.run.nums"]] = switch(i, c(23,24), c(24,25))
	## stock[[ii]][[jj]][["dia.rwt.nums"]] = switch(i, c(1,2),   c(3,3))
	## stock[[ii]][[jj]][["dia.lab"]]      = switch(i, c("Base + Fix_M","S4 + Fix_M"), c("Base + Fix_M","S4 + Fix_M"))
	stock[[ii]][[jj]][["sen.mcsubs"]]      = switch(i, rep(list(c(1,2000)),14), rep(list(c(1,2000)),1))
	stock[[ii]][[jj]][["policies"]]        = switch(i, c(CC=1250,HR=0.1), c(CC=300,HR=0.06))  ## set policies here but there could be a disconnect in policies determined by 'make.base.figs.r' depending on what policies are available.
}
##----------------------------

## Commonalities
spp.name   = .su(sapply(stock,function(x){x$Controls$spp.name}))
latin.name = .su(sapply(stock,function(x){x$Controls$latin.name}))
spp.code   = .su(sapply(stock,function(x){x$Controls$spp.code}))
the.stocks = gsub(" ","~",.su(sapply(stock,function(x){x$Controls$name})))
ptype      = "png" ## make sure this is the same choice on line 19
pngres     = 400   ## pixels per inch
sigdig     = 3     ## Number of significant digits to output in tables
decdig     = 2     ## Number of decimal digits (e.g., may only want 2 for decision tables)
## Create quants and put into .PBSmodEnv just in case a PBSawatea plotting function might need either.
quants3    = c(0.05,0.50,0.95);           tput(quants3)
quants5    = c(0.05,0.25,0.50,0.75,0.95); tput(quants5)

## Unpack the years to get Sweave started in the wrapper Rnw
unpackList(stock[[1]][["Controls"]][c("startYear","currYear","prevYear","projYear","pgenYear")])

divTwo = function(x, dig=decdig)   ## dig is number of dec places
{
	show0(round(x[1]/x[2],dig),dig)
}
medCI = function(x, dig=decdig, CI=quants3[c(1,3)])  ## dig is number of dec places, CI = credible interval
{ print(paste0(c(
  prettyNum(round(quantile(x, 0.50), digits=dig), big.mark=","), "~(",
  prettyNum(round(quantile(x, CI[1]), digits=dig), big.mark=","), ",~",
  prettyNum(round(quantile(x, CI[2]), digits=dig), big.mark=","), ")"), collapse=""))
}
relabelTex = function(texinput, prefix, caption) ## (RH 190513)
{
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
findPos = function(target, choices)
{
	## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
	target = as.numeric(target)
	choices = as.numeric(choices)
	which(abs(choices-target)==min(abs(choices-target)))
}
allEqual <- function(x)
{
	result <- all( x==x[1] )
	result
}
## Get the number of active dimensions from the axis of uncertainty matrix
getActDim = function(A){
	Nact = 0
	for (i in 1:length(dim(A))) {
		Apop = apply(A,i,function(j){!is.na(j)})
		if (is.null(dim(Apop)) && sum(Apop)>1)
			return(1)
		if (sum(apply(Apop,2,any))>1) Nact = Nact + 1
	}
	return(Nact)
}

saveDate = sub("/0","/",sub("^0","",format.Date(Sys.time(),"%m/%d/%Y")))

getMPD = function(obj)
{
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

#browser();return()
	return(mpd)
}
#getMPD(currentRes.base)

