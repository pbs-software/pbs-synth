##==========================================================
## PBS Stock Synthesis small helper functions (some not used)
##----------------------------------------------------------
## allEqual.........Determine if all elements in vector equal the first element
## cquantile.vec....Calculate cumulative quantile as a vector
## divTwo...........Divide first element of vector by second element
## findPos..........Find nearest position of a number in a vector of numbers
## gather.vectors...Gather vectors of base run components
## getActDim........Get number of active dimensions from the axis of uncertainty matrix
## getAssYrs........Get modelled current years for previous stock assessments
## getMPD...........Get MPD estimated parameters
## getNpan..........Get panel number when inside a multi-panel plot
## getYrIdx.........Select years for plotting
## is.numStr........Check if strings can be converted to numerics
## med5.95..........Print median (0.05, 0.95) to text
## medCI............Print the median and the credible interval
## ptab.............Prior tabulation (LaTeX) for lines in table of priors
## qtab.............Quantile tabulation summary using decimal places
## relabelTex.......Relabel label and caption of table/figure from MPD/MCMC run
## stab.............Quantile tabulation summary using significant digits
##==========================================================

## -------------------------------------2010-10-20
##  Determine if all elements in vector equal the first element
## ---------------------------------------------RH
allEqual <- function(x)
{
  result <- all( x==x[1] )
  result
}

## cquantile.vec------------------------2010-10-20
##  Calculate cumulative quantile as a vector
##  AME doing this, just do one prob at a time 
##  (so it returns a vector not a matrix)
## --------------------------------------------AME
cquantile.vec <- function(z, prob)  # cumulative quantile of vector; prob is a single number
{
  cquant <- rep(NA, length(z))
  if(length(prob) != 1) stop("length prob should be 1")
  for (i in 1:length(z)) {
    cquant[i] <- quantile(z[1:i], probs = prob, names = FALSE)
  }
  return(cquant)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cquantile.vec

## divTwo ------------------------------2026-07-16
##  Divide first element of vector by second element
## ---------------------------------------------RH
divTwo <- function(x, dig=decdig)
{
	show0(round(x[1]/x[2],dig),dig)
}

## findPos -----------------------------2026-07-16
##  Find nearest position of a number in a vector of numbers
## ---------------------------------------------DM
findPos <- function(target, choices)
{
	## https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
	target = as.numeric(target)
	choices = as.numeric(choices)
	which(abs(choices-target)==min(abs(choices-target)))
}

## gather.vectors ----------------------2026-07-16
##  Gather vectors of base run components
## ---------------------------------------------RH
gather.vectors <- function(run, rwt, ver, use)
{
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
				}
			}
		}
	}   ## RH 200508|200623 (for subsets of Base runs)
	return(list(run.num=run.num, rwt.num=rwt.num, ver.num=ver.num, use.num=use.num))
}

## getActDim ---------------------------2026-07-16
##  Get number of active dimensions from the axis of uncertainty matrix
## ---------------------------------------------RH
getActDim <- function(A)
{
	if (prod(dim(A))==1)
		return(0)
	Nact = 0
	for (i in 1:length(dim(A))) {
		Apop = apply(A,i,function(j){!is.na(j)})
		if (is.null(dim(Apop)) && sum(Apop)>1)
			return(1)
		if (sum(apply(Apop,2,any))>1) Nact = Nact + 1
	}
	return(Nact)
}

## getAssYrs ---------------------------2026-07-16
##  Get modelled current years for previous stock assessments
## ---------------------------------------------RH
getAssYrs <- function(spp.code, area.name)
{
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

## getMPD ------------------------------2026-07-16
##  Get MPD estimated parameters
## ---------------------------------------------RH
getMPD <- function(obj)
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
	return(mpd)
}

## getNpan------------------------------2019-05-10
##  Get panel number when inside a multi-panel plot.
## ---------------------------------------------RH
getNpan <- function()
{
	mfg=par()$mfg
	mfg[2]+(mfg[1]-1)*mfg[4]
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getNpan

## getYrIdx-------------------------------2011-08-31
##  Purpose is to return selected years for plotting.
##  Default is to select 5 year increments.
##---------------------------------------------AME
getYrIdx <- function(yrNames, mod=5)
{
  ## Coerce to numeric and select the years modulo "mod".
  yrVals <- as.numeric( gsub("[^[:digit:]]","",yrNames) )
  idx <- yrVals %% mod==0

  ## Select years from character vector yrNames.
  result <- yrNames[ idx ]
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getYrIdx

## is.numStr----------------------------2020-09-20
##  Check if strings can be converted to numerics.
##----------------------------------------------RH
is.numStr <- function(x)
{
	out = sapply(x, function(xx) {
		xx = as.character(xx)
		all(grepl("[[:digit:]]|\\-|\\.|[eE]|[[:space:]]", strsplit(xx,"")[[1]]))
	})
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~is.numStr

## med5.95 -----------------------------2010-10-20
##  Print median (0.05, 0.95) to text
## -----------------------------------------AME|RH
med5.95 <- function(xx.MCMC, dig=0, quants3=tcall(quants3))
{
	## dig is number of dec places
	if (is.null(quants3)) quants3=c(0.05,0.5,0.95)
	big.mark = options()$big.mark; if (is.null(big.mark)) big.mark=","
	mess = paste0( c( prettyNum(round(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=big.mark),
		"~(", prettyNum(round(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=big.mark), 
		",\\,", prettyNum(round(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=big.mark), ")"), collapse="")
	print(mess)
}

## medCI -------------------------------2026-07-16
##  Print the median and the credible interval
## ---------------------------------------------RH
medCI <- function(x, dig=decdig, CI=quants3[c(1,3)])
{
	print(paste0(c(
	prettyNum(round(quantile(x, 0.50, na.rm=T), digits=dig), big.mark=","), "~(",
	prettyNum(round(quantile(x, CI[1], na.rm=T), digits=dig), big.mark=","), ",~",
	prettyNum(round(quantile(x, CI[2], na.rm=T), digits=dig), big.mark=","), ")"), collapse=""))
}

## ptab---------------------------------2025-04-08
##  Function to use for priors in table (adapted from PBSawatea).
## ---------------------------------------------RH
ptab <- function(xx)
{
	xx = sub("^\\s+", "", xx)  ## remove leading and trailing whitespace
	xlab = gsub("\\_+"," ",xx[1])
#browser();return()
	xnum =xx[-1]
	xnum[4] = switch(xnum[4], 'Normal'=6, 'No_prior'=0, 'Full_Beta'=2, 'Sym_Beta'=1)
	xnum = lapply(xnum,function(x){
		if(is.numStr(x))   as.numeric(x)
		else if (is.na(x)) "--"
		else x
	})
	xout = paste0(c(xlab, " & ", xnum[[1]], " & (", xnum[[2]], ", ", xnum[[3]], ") & ", xnum[[4]], " & (", xnum[[5]], ", ", xnum[[6]], ") & ", xnum[[7]], " & ", show0(round(xnum[[8]],3),3), " \\\\\\\\"), collapse="")
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ptab

## qtab---------------------------------2021-04-19
## Quantile tabulation summary using decimal places
## --------------------------------------------AME
qtab <- function(xx.MCMC, dig=0, quants3=tcall(quants3))
{  ## dig is number of dec places
	if (is.null(quants3)) quants3=c(0.05,0.5,0.95)
	big.mark = options()$big.mark; if (is.null(big.mark)) big.mark=","
	print(paste0( c( prettyNum(round(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=big.mark),
		" & ", prettyNum(round(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=big.mark),
		" & ", prettyNum(round(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=big.mark)), collapse=""))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~qtab

## relabelTex --------------------------2019-05-13
##  Relabel label and caption of table/figure from
##  MPD/MCMC run used to summarize model results
## ---------------------------------------------RH
relabelTex <- function(texinput, prefix, caption)
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

## stab---------------------------------2021-04-19
## Quantile tabulation summary using significant digits
## --------------------------------------------AME
stab <- function(xx.MCMC, dig=3, quants3=tcall(quants3), print=TRUE)
{  ## dig is number sig digits
	if (is.null(quants3)) quants3=c(0.05,0.5,0.95)
	big.mark = options()$big.mark; if (is.null(big.mark)) big.mark=","
	out = paste0( c( prettyNum(signif(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=big.mark), 
		" & ", prettyNum(signif(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=big.mark),
		" & ", prettyNum(signif(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=big.mark)), collapse="")
	if (print) print(out)
	invisible(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~stab
