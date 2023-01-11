##==========================================================
## PBS Stock Synthesis utility functions:
## --------------------------------------
## calcStdRes............Calculate standardised residuals for robustified normal likelihood
## convPN................Convert parameter names from SS to Awatea
## cquantile.vec.........Calculate cumulative quantile as a vector
## extract.between.......Extract character strings between two delimiters.
## findTarget............Derive decision tables for reference points, etc.
## gatherMCMC............Gather and combine MCMC (and MPD) values from compo and/or senso runs
## getNpan...............Get panel number when inside a multi-panel plot.
## getSS.control.........Get select values out of the control file.
## importCor.............Import SS parameter correlations (mod.PBSawatea)
## importEva.............Import SS Hessian eigenvlaues (mod.PBSawatea)
## importPar.............Import all SS parameters (mod.PBSawatea).
## is.numStr.............Check if strings can be converted to numerics.
## prepCP................Prepare Catch Policies -- 'CC'=constant catch
## prepMPD...............Prepare MPD runs for likelihood analysis
## ptab..................Function to use for priors in table (adapted from PBSawatea).
## qtab..................Quantile tabulation summary using decimal places
## repeatMPD ............Repeat MPDs for axes of uncertainty to visualise likelihood density
## runADnuts.............Based on Chris Grandis' R code that was used to manually run adnuts
## runSweave.............Run Sweave code to build pdfs for MPD and MCMC runs (not appendix)
## stab..................Quantile tabulation summary using significant digits
## weightAF..............Weight age frequencies using harmonic mean ratio method
##
## Hidden functions
## ----------------
## .flash.cat............Flush the cat down the console (also appears in PBStools)
##==========================================================

allEqual <- function(x)
{
  result <- all( x==x[1] )
  result
}

## calcStdRes---------------------------2022-11-01
##  This implements standardised residuals for the Awatea
##  implementation of the Fournier robustified normal likelihood
##  for proportions at length. Based on PJS summary of CASAL doc and ACH change to length.
##------------------------------------------AME|RH
calcStdRes <- function( obj, trunc=3, myLab="Age Residuals", prt=TRUE, type="Multinomial",
   afld="Bin", yfld="Yr", ofld="Obs", ffld="Exp", nfld="Nsamp_adj" )
{
	# Make a column for the standardised residuals.
	result <- cbind( obj,stdRes=rep(NA,nrow(obj)) )

	# Raw residuals.
	res <- obj[,ofld] - obj[,ffld]  ## (O-F)

	# Number of age bins.
	# QUESTION: Should this be from FIRST age to plus group?
	nage <- length( unique(obj[,afld]) )

	# Kludgy, but loop through years.
	# Could reformat into matrix and vectorize.

	## Robustifying function
	Z = function (x, r=0.001){
		is.soft = x < r
		if (any(is.soft))
			x[is.soft] = r / (2 - x[is.soft]/r)
		return(x)
	}

	yrList <- sort( unique( obj[,yfld] ) )
	for ( i in 1:length( yrList ) )
	{
		idx <- yrList[i]==obj[,yfld]
		## Pearson residual = (O-F)/std.dev(O)  : see CASAL Manual, Section 6.8 Residuals
		## where std.dev(O) is calculated as:
		Nprime  <- min( obj[,nfld][idx],1000)                         ## N prime
		Fprime  <- obj[,ffld][idx]*(1.0-obj[,ffld][idx]) + 0.1/nage   ## F prime (Fournier uses Fitted)
		Oprime  <- obj[,ofld][idx]*(1.0-obj[,ofld][idx]) + 0.1/nage   ## O prime (Coleraine uses Observed)
		Mprime  <- Z(obj[,ffld][idx]) * (1-Z(obj[,ffld][idx]))        ## M prime (Multinomial uses Fitted)
		SD      <- sqrt(
			switch(type, 'Multinomial'=Mprime, 'Fournier'=Fprime, 'Coleraine'=Oprime) / Nprime )  
		result$stdRes[idx] <- res[idx]/SD                             ## Pearson residuals = Normalised residuals for normal error distributions
#browser();return()
		## For Dirichlet-Multinomial, SS code "SS_write_report.tpl" shows:
		##   show_Pearson = value((ocomp - ecomp) / sqrt(ecomp * (1.0 - ecomp) / nsamp * (nsamp + dirichlet_Parm) / (1. + dirichlet_Parm))); // Pearson for Dirichlet-multinomial using negative-exponential parameterization
		## But this is too complicated to replicated here (for now)
	}
	if ( prt ) {
		sdRes <- sqrt( var( result$stdRes,na.rm=TRUE ) )
		sdTrunc <- ifelse( result$stdRes > trunc, trunc, result$stdRes )
		sdTrunc <- ifelse( result$stdRes < -trunc, -trunc, result$stdRes )
		sdResTrunc <- sqrt( var( sdTrunc,na.rm=TRUE ) )
		cat( "\n",myLab,"\n" )
		cat( "\n     Std. Dev. of standardised Residuals=",sdRes,"\n" )
		cat( "     Std. Dev. of Truncated standardised Residuals=",sdResTrunc,"\n" )
	}
	return(result)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcStdRes


## convPN-------------------------------2022-06-15
##  Convert parameter names from SS to Awatea
## ---------------------------------------------RH
convPN = function(pnams) {
	cnams =
		sub("_steep","_h",
		sub("NatM_p_1_Fem_GP_1|NatM_uniform_Fem_GP_1|NatM_break_(\\d+)?_Fem_GP_1","M\\1_Female",
		sub("NatM_p_1_Mal_GP_1|NatM_uniform_Mal_GP_1|NatM_break_(\\d+)?_Mal_GP_1","M\\1_Male",
		sub("Early_RecrDev|Main_RecrDev|Late_RecrDev|ForeRecr","RecrDev",
		sub("FISHERY","", 
		sub("SYNOPTIC","",
		sub("HISTORIC(AL)?","",
		sub("TRIENNIAL","",
		sub("HBLL_NORTH","HBLLN_",
		sub("HBLL_SOUTH","HBLLS_",
		sub("^AgeSel_(\\d+)?(Male|Fem)_Scale","delta5(\\1)",
		sub("^AgeSel_(\\d+)?(Male|Fem)_Final","delta4(\\1)",
		sub("^AgeSel_(\\d+)?(Male|Fem)_Descend","delta3(\\1)",
		sub("^AgeSel_(\\d+)?(Male|Fem)_Ascend","delta2(\\1)",
		sub("^AgeSel_(\\d+)?(Male|Fem)_Peak","delta1(\\1)",
		sub("^Age_DblN_end_logit","beta6",
		sub("^Age_DblN_top_logit","beta2",
		sub("^Age_DblN_descend_se","varR",
		sub("^Age_DblN_ascend_se","varL",
		sub("^Age_DblN_peak","mu",
		sub("^SR_","",
		pnams)))))))))))))))))))))
		onams   = sapply(strsplit(cnams,"_"), function(x){
		if(length(x)==3 && x[1] %in% c("mu","beta2","varL","varR","beta6")) {
#browser();return()
			paste0(x[1],x[3],"_",x[2])
		} else {
			paste0(x,collapse="_")
		}
	})
	onams = sub("^M1_", "M_", onams)
	return(onams)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~convPN


## cquantile.vec------------------------2010-10-20
##  Calculate cumulative quantile as a vector
##  AME doing this, just do one prob at a time 
##  (so it returns a vector not a matrix)
## --------------------------------------------AME
cquantile.vec <- function(z, prob)  # cumulative quantile of vector
{                                   #  prob is a single number
  cquant <- rep(NA, length(z))
  if(length(prob) != 1) stop("length prob should be 1")
  for (i in 1:length(z))
    {
    cquant[i] <- quantile(z[1:i], probs = prob, names = FALSE)
    }
  return(cquant)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cquantile.vec


## extract.between----------------------2020-09-17
##  Extract character strings between two delimiters; based on:
##  https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
## ---------------------------------------------RH
extract.between = function(x, open="{", close="}", first.only=TRUE)
{
	metas = c(".", "\\", "|", "(", ")", "[", "{", "^", "$", "*", "+", "?")
	if (open  %in% metas) open  = paste0("\\\\",open)
	if (close %in% metas) close = paste0("\\\\",close)
	x.names = names(x)
	snippets = lapply(1:length(x), function(i){
		mess  = paste0("gsub(\"[", open, close, "]\", \"\", regmatches(\"", x[i], "\", gregexpr(\"", open, ".*?", close, "\", \"", x[i], "\"))[[1]])")
		snip  = eval(parse(text=mess))
		return(snip)
	})
	if (first.only)
		snippets = sapply(snippets,function(s){return(s[1])})
	names(snippets) = x.names
	return(snippets)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~extract.between


## findTarget --------------------------2019-12-04
##  To derive decision tables for reference points, moving windows,
##  and to find # yrs to achieve recovery with given confidence.
##   Vmat   = matrix of projected B-values (MCMC projections x Year)
##   yrP    = user-specified projection years
##   yrG    = number of years for moving target window (e.g, 90y=3 YMR generations). Might not work for all possibilities.
##   ratio  = recovery target ratio
##   target = recovery target values (e.g., B0, Bmsy).
##          = B0.MCMC for ratios of B0
##          = Bmsy.MCMC for ratios of Bmsy
##          = Bt.MCMC for moving window
##   conf   = confidence level required
##   plotit = logical to plot the probability of Bt/target
##   retVal = character name of object to return
##          = "N", look for the global object "Ttab" (number of years to acheive target)
#           = "p.hi" gives global object "Ptab", a list of decision tables where row is the catch option and column is the year.
##  Values are probabilities of acheiving target.
##  (2012-02-20) 'xhi=x>=r' change to 'xhi=x>r'
## ---------------------------------------------RH
findTarget=function(Vmat, yrU=as.numeric(dimnames(Vmat)[[2]]), yrG=90, 
   ratio=0.5, target=B0.MCMC, conf=0.95, plotit=FALSE, retVal="N", op=">")
{
	## oldpar=par(no.readonly=TRUE);  on.exit(par(oldpar))
	yrA   =as.numeric(dimnames(Vmat)[[2]])   ## years available
	yrP   =sort(intersect(yrA,yrU))          ## years for proj
	yr0   =yrP[1]; yrN=rev(yrP)[1]
#browser();return()

	vmat=Vmat[,is.element(dimnames(Vmat)[[2]],as.character(yrP))]             ## include only yrP years
#browser();return()

	## Check if COSEWIC reference criterion
	if (is.data.frame(target) || is.matrix(target)) {
		yrM  = yrP - yrG                                                       ## moving target years
		yrM1 = intersect(as.numeric(dimnames(target)[[2]]),yrM)                ## available target years from MCMC
		if (length(yrM1)==0) {                                                 ## projection not long enough for any overlap with 3 generations
			if (retVal=="N") return(NA)
			else {p.hi=rep(NA,length(yrP)); names(p.hi)=yrP }; return(p.hi) }
		yrMr =range(yrM1)                                                      ## range of years to use from MCMC
		targM=target[,as.character(yrM1)]                                      ## target data from MCMC
		yrM2 =setdiff(yrM,yrM1)                                                ## missing target years (can occur before and after the MCMC years)
#browser();return()

		if (length(yrM2)>0) {
			nrow=dim(target)[1]
			if (any(yrM2<yrMr[1])) {
				yrMo =yrM2[yrM2<yrMr[1]]                                         ## years of data older than MCMCs
				ncol =length(yrMo)
				targ0=matrix(rep(target[,as.character(yrM1[1])],ncol),
					nrow=nrow, ncol=ncol, dimnames=list(1:nrow,yrMo))             ## repeat B0 (first column)
				targM=cbind(as.data.frame(targ0),targM)                          ## moving target
			}
			if (any(yrM2>yrMr[2])) {
				yrMn =yrM2[yrM2>yrMr[2]]                                         ## years of data newer than MCMCs
				ncol =length(yrMn)
				targN=vmat[,as.character(yrMn)]                                  ## start using projections
				targM=cbind(targM,targN)                                         ## moving target
			}
		}
		rats=vmat/targM                                                        ## matrix of ratios Bt/ moving target
	}
	else    ## if it's a vector, so no moving window
		rats=apply(vmat,2,function(x,targ){x/targ},targ=target)                ## matrix of ratios Bt/ target (B0 or Bmsy)

	#p.hi=apply(rats,2,function(x,r){xhi=x>r; sum(xhi)/length(xhi)},r=ratio)  ## vector of probabilities Bt/B0 > target ratio for each year.

	## vector of probabilities Bt/B0 op (>|<) target ratio for each year.
	p.hi=apply(rats,2,function(x,r){xhi= eval(call(op,x,r)); sum(xhi)/length(xhi)}, r=ratio)

	## p.hi can become each row of a decision table (AME checked)
	##  the numbers for 0.4 Bmsy match my existing
	##  independent calculations). Need to save this for moving window.

	z.hi = p.hi >= conf                               ## logical: is p.hi >= confidence limit specified
#browser();return()

	if (all(z.hi) ||                                  ## all p.hi exceed the confidence level, also check:
		(all(rats[,1]==1) && all(z.hi[-1]))) yrT=yr0   ## if values in first year of projection = values in target (e.g., Bcurr) before proceeding
	else if (!any(z.hi)) yrT=yrN                      ## no  p.hi exceed the confidence level
	else {
		pdif = round(diff(p.hi),5)                     ## one-year change in trend
		z1 = pdif >= 0                                 ## logical: trend increasing? -- sometimes it just remains flat
		## z2 = c(pdif[-1],rev(pdif)[1]) >= 0          ## logical: trend one period later increasing? (or flat) -- PJS does not like this requirement
		z3 = z.hi[-1]                                  ## does the probability of equalling or exceeding the target ratio exceed the confidence level?
		## z  = z1 & z2 & z3                           ## logical: potential years when target reached
		z  = z1 & z3                                   ## logical: potential years when target reached
		if (!any(z)) yrT=yrN                           ## target not reached within the projection period
		else {
			yrT=as.numeric(names(z)[z][1])              ## first year when target reached
			
		}
	}
	N=yrT - yr0                                       ## number of years to reach target
	if (plotit) {
		par(mar=c(4,5,0.5,0.5))
		#ylim=c(0, max(p.hi,ratio))
		ylim=c(min(p.hi,0.5),1)
		plot(yr0:yrN,p.hi,type="n",ylim=ylim,ylab="",mgp=c(2.0,0.75,0))
		lines(yr0:yrN,p.hi,col="grey")
		points(yr0:yrN,p.hi,pch=20,col="orange",cex=1.2)
		mtext(text=expression(p~~frac(B[t],B[Target]) ), side=2, line=1.5, cex=1.5)
		abline(h=conf,v=yrT,col="dodgerblue")
	}
#browser();return()
	eval(parse(text=paste("return(",retVal,")",sep=""))) 
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~findTarget


## gatherMCMC---------------------------2022-06-23
## All values occur at start of year so B2022 is B at end of 2021 or start of 2022;
##   u2022 must be u in 2021 after final catch was specified for 2021 in the data file.
## Ian Taylor (Aug 9, 2021):
##   Everything is start of year throughout.
##   Setting 2021 as the final year of the model allows 2021 catch.
##   The values for derived quantities like SSB_2021 and Bratio_2021 represent the beginning of the year values at the start of 2021 in the main time series.
##   The SSB_2022 is one year later but will only depend on catches up through 2021, so changes in 2022 catch won't impact the 2022 quantities.
## ---------------------------------------------RH
gatherMCMC = function( mcdir=".", type="compo", strSpp="CAR",
   basedir="C:/Users/haighr/Files/GFish/PSARC22/CAR/Data/SS/YMR2022",
   ryrs=1935:2023, pyrs=2024:2033, valTS=c("SSB","F","Recr","RecrDev"), valRP=c("SSB_MSY","annF_MSY","Dead_Catch_MSY"))
{
	ayrs  = c(ryrs, pyrs); nyrs = length(ayrs)
	currYr = rev(ryrs)[1]; byr = cyr = as.character(currYr)
	#prevYr = currYr-1;     byr = as.character(prevYr) ## see note above
	mclst = Nmcmc = rowN = valPAs = valLLs = list()

	for (m in 1:length(mcdir)){
		mdir = file.path( sub("/$","",basedir), sub("/$","",mcdir[m]) )
		if(!file.exists(mdir)) next
		mm    = substring(basename(mdir),6,10)
		if (penso) {
			## Special versions for PJS
			mm = substring(mdir,regexpr("Run",mdir)+3,regexpr("Run",mdir)+7)
		}
#browser();return()
		mclst[[mm]] = list()
		cp = "AC.00"

		## Gather MPD information
		d.mpd = sub("\\.mh.+","",sub("\\.hm.+","",sub("\\.nuts.+","",sub("MCMC","MPD",mdir))))
		mpd = SS_output(dir=d.mpd, verbose=F, printstats=F)

		## Grab likelihood components
		run = as.numeric(strsplit(mm,"\\.")[[1]][1]); names(run)="Run"
		LL.fleet  = mpd$likelihoods_by_fleet
		LL.fleet.ex = LL.fleet[is.element(LL.fleet$Label,"Surv_like"),grep("OTHER",mpd$FleetNames,invert=T,value=T)]
		names(LL.fleet.ex) = sapply(strsplit(names(LL.fleet.ex),"_"),function(x){return(paste0(x[1],"_",substring(x[2],1,3))) })
		#names(LL.fleet.ex) = sub("TRAWL","CPUE",names(LL.fleet.ex))
		names(LL.fleet.ex)[grep("TRAWL",names(LL.fleet.ex))]="CPUE_BT"
		LL.used   = mpd$likelihoods_used
		LL.used.ex  = LL.used[c("Survey","Age_comp","Recruitment","TOTAL"),"values"]
		names(LL.used.ex) = c("Index","AF","Recruit","Total")
		LL.compo = c(run, LL.fleet.ex, LL.used.ex)
		LL.compo = unlist(LL.compo)
		mclst[[mm]][["LL"]] = LL.compo
		valLLs[[mm]] = names(LL.compo)  ## collect all names
#browser();return()

		## Grab parameter estimates
		parameters = mpd$parameters
		pactive = parameters[!is.na(parameters$Active_Cnt) & parameters$Phase>0 & !is.element(parameters$Pr_type,"dev"),]
		pactive$Label = convPN(pactive$Label)
		P.mpd = pactive$Value; names(P.mpd)=pactive$Label
		valPA = names(P.mpd)
		valPAs[[mm]] = valPA  ## collect all
		if (m==1 && strSpp %in% c("YMR")) {
			valPA = setdiff(valPA,c("M_Female","M_Male"))  ## just in case R79 comes first
			if (any(grepl("Run79",mcdir))) ## Specific to YMR sensitivity run with natural mortality
				valPA = c(valPA, "M_Female", "M_Male")
		}
		## Check for catch policies
		CP = c("AC", c("CC","HR")[file.exists(file.path(mdir,c("CC","HR")))])
		for (n in 1:length(CP)){
			nn = CP[n]
			if (nn=="AC") next
			ndir = file.path(mdir,nn)
			ncps = dir (ndir)  ## name of catch policy subdirectories
			cp = c(cp, paste0(nn,".",sub("^[[:alpha:]]+","",ncps)))
			mdir = c(mdir, file.path(ndir,ncps))
		}
		for (n in 1:length(mdir)) {
			mmm   = mdir[n]
			ccc   = cp[n]
			.flush.cat(paste0("Getting MCMC: Run ", mm, "; catch policy: ", ccc), "\n")
			mmc   = SSgetMCMC(mmm, verbose=F)
#browser();return()
			colnames(mmc) = convPN(colnames(mmc))
			msid  = mmc$Iter
			ii    = as.character(msid)
			nmc   = length(msid)
			#run   = as.numeric(strsplit(mm,"\\.")[[1]][1])
			run   = mm
			if (m==1){
				Nmcmc[[n]] = nmc
				rowN[[n]]  = paste0(run,".", pad0(msid,3))
			} else {
#browser();return()
				Nmcmc[[n]] = Nmcmc[[n]] + nmc
				rowN[[n]]  = c(rowN[[n]], paste0(run,".", pad0(msid,3)))
			}
#if (m==2) {browser();return()}

			if (n==1) { ## only need to look at the base catch policy for PA and RP
			## Parameters matrix
				matPA = array(NA, dim=c(nmc, length(valPA)), dimnames=list(sid=msid, val=valPA))
				for (k in 1:length(valPA)) {
					kk = valPA[k]
					vecPA = mmc[[kk]]
					if (is.null(vecPA)) next
					matPA[ii,kk] = vecPA
				}
				mclst[[mm]][["PA"]] = matPA

				## Reference points matrix
				matRP = array(NA, dim=c(nmc, length(valRP)), dimnames=list(sid=msid, val=valRP))
				for (k in 1:length(valRP)) {
					kk = valRP[k]
					vecRP = mmc[[kk]]
					if (is.null(vecRP)) next
					matRP[ii,kk] = vecRP
				}
				mclst[[mm]][["RP"]] = matRP

#browser();return()
				## Time series matrix of population reconstruction
				matTS = array(NA, dim=c(nmc, length(ryrs), length(valTS)), dimnames=list(sid=msid, yr=ryrs, val=valTS))
				for ( k in 1:length(valTS) ) {
					kk = valTS[k]
					for (j in ryrs) {
						jj = as.character(j)
						vecTS = mmc[[paste0(kk,"_",jj)]]
						if (is.null(vecTS)) next
						matTS[ii,jj,kk] = vecTS
					}
				}
				mclst[[mm]][["TS"]] = matTS
			}
			## Projection series matrix (average catch and whatever catch policies are available)
			matPJ = array(NA, dim=c(nmc, length(pyrs), length(valTS)), dimnames=list(sid=msid, yr=pyrs, val=valTS))
			for ( k in 1:length(valTS) ) {
				kk = valTS[k]
				for (j in pyrs) {
					jj = as.character(j)
					vecPJ = mmc[[paste0(kk,"_",jj)]]
					if (is.null(vecPJ)) next
					matPJ[ii,jj,kk] = vecPJ
				}
			}
#if (n==2){browser();return()}
			mclst[[mm]][["PJ"]][[ccc]] = matPJ

			## Catch policy matrix
			matCP = array(NA, dim=c(nmc, length(pyrs)), dimnames=list(sid=msid, yr=pyrs))
			for (j in pyrs) {
				jj = as.character(j)
				vecCP = mmc[[paste0("ForeCatch_",jj)]]
				if (is.null(vecCP)) next
				matCP[ii,jj] = vecCP
			}
#browser();return()
			#if ( all(apply(matCP,2,function(x){length(.su(x))}) %in% c(0,1)) )
			#	mclst[[mm]][["CP"]][[ccc]] = apply(matCP,2,unique)
			#else
			#	mclst[[mm]][["CP"]][[ccc]] = matCP
#if (ccc=="CC.02") {browser();return()}
			## Note; 0-catch policy is sometimes set to a very small amount; sometimes other catch policies have small amounts added. WTF?
			## PJS thinks it may be a difference between v.2.30.16 and 3.30.17
			#mclst[[mm]][["CP"]][[ccc]] = apply(matCP,2,function(x){xx=mean(x,na.rm=T); xx[is.na(xx) | !is.finite(xx)]=NA; return(xx)})
			## CAUTION: Mode may not be sufficiently robust at higher CPs
			mclst[[mm]][["CP"]][[ccc]] = apply(round(matCP),2,function(x){xx=as.numeric(names(rev(sort(table(x))))[1]); xx[is.na(xx) | !is.finite(xx)]=NA; return(xx)})
		}
		mclst[[mm]][["MPD"]] = P.mpd
	}
	if (strSpp %in% "CAR") {
		valPA = unique(unlist(valPAs))
		valPA = c(grep("R0",valPA,value=T), grep("theta|R0",valPA,invert=T,value=T), grep("theta",valPA,value=T))
		valLL = unique(unlist(valLLs))
		chunk = "Index|AF|Recruit|Total"
		valLL = c(grep(chunk,valLL,invert=T,value=T), grep(chunk,valLL,value=T))
	}
	mpdPA = list()
	#avgPA = array(NA, dim=c(Nmcmc[[1]], length(P.mpd)), dimnames=list(mcmc=rowN[[1]], par=names(P.mpd)))
	avgPA = array(NA, dim=c(Nmcmc[[1]], length(valPA)), dimnames=list(mcmc=rowN[[1]], par=valPA))
	vRP   = c("Bcurr","B0","20B0","40B0","Fcurr","ucurr","MSY","Bmsy","LRP","USR","Fmsy","umsy")
	avgRP = array(NA, dim=c(Nmcmc[[1]], length(vRP)), dimnames=list(mcmc=rowN[[1]], val=vRP))
	vTS   = c("Bt","BtB0","BtBmsy","Ft","FtFmsy","ut","utumsy","Rt","Rtdev")
	avgTS = array(NA, dim=c(Nmcmc[[1]], length(ryrs), length(vTS)), dimnames=list(mcmc=rowN[[1]], yr=ryrs, val=vTS))
	avgPJ = array(NA, dim=c(Nmcmc[[1]], length(pyrs), length(vTS), length(cp)), dimnames=list(mcmc=rowN[[1]], yr=pyrs, val=vTS, proj=cp))
	avgCP = array(NA, dim=c(length(pyrs), length(mclst), length(cp)), dimnames=list(yr=pyrs, run=names(mclst), proj=cp))

	avgLL =array(NA, dim=c(length(valLL), length(mclst)), dimnames=list(ll=valLL, run=names(mclst)))
#browser();return()


	## Build the composite run ('model average')
	for (m in 1:length(mclst)){
		mm    = names(mclst)[m]
		matrix(mclst[[mm]][["MPD"]][valPA],nrow=1)
		mpdPA = rbind(mpdPA, matrix(mclst[[mm]][["MPD"]][valPA], nrow=1))
		mcPA  = mclst[[mm]][["PA"]]
		mcRP  = mclst[[mm]][["RP"]]
		mcTS  = mclst[[mm]][["TS"]]
		mcPJ  = mclst[[mm]][["PJ"]]
		mcCP  = mclst[[mm]][["CP"]]
		#run   = as.numeric(strsplit(mm,"\\.")[[1]][1])
		run   = mm
		iii   = paste0(run,".",pad0(as.numeric(rownames(mcTS)),3))
		#iii   = (lastmc+1):(lastmc+dim(mcTS)[1])

		## Populate likelihood array
		for (n in 1:length(mclst[[mm]][["LL"]])){
			lll = names(mclst[[mm]][["LL"]])[n]
			avgLL[lll, mm] = mclst[[mm]][["LL"]][lll]
		}
		## Populate parameter array
		for (n in 1:length(mclst[[mm]][["MPD"]])){
			ppp = names(mclst[[mm]][["MPD"]])[n]
			avgPA[iii,ppp] = mcPA[,ppp]
		}
#if (m==2) {browser();return()}
		## Populate reference point array
		avgRP[iii,"Bcurr"] = mcTS[,cyr,"SSB"]
		avgRP[iii,"B0"]    = mcTS[,1,"SSB"]
		avgRP[iii,"20B0"]  = 0.2 * mcTS[,1,"SSB"]
		avgRP[iii,"40B0"]  = 0.4 * mcTS[,1,"SSB"]
		avgRP[iii,"Fcurr"] = mcTS[,byr,"F"]
		avgRP[iii,"ucurr"] = 1 - exp(-mcTS[,byr,"F"])
		avgRP[iii,"MSY"]   = mcRP[,"Dead_Catch_MSY"]
		avgRP[iii,"Bmsy"]  = mcRP[,"SSB_MSY"]
		avgRP[iii,"LRP"]   = 0.4 * mcRP[,"SSB_MSY"]
		avgRP[iii,"USR"]   = 0.8 * mcRP[,"SSB_MSY"]
		avgRP[iii,"Fmsy"]  = mcRP[,"annF_MSY"]
		avgRP[iii,"umsy"]  = 1-exp(-mcRP[,"annF_MSY"])

		## Populate time series array
		avgTS[iii,dimnames(mcTS)$yr,"Bt"]     = mcTS[,,"SSB"]
		avgTS[iii,dimnames(mcTS)$yr,"BtB0"]   = mcTS[,,"SSB"] / avgRP[iii,"B0"]
		avgTS[iii,dimnames(mcTS)$yr,"BtBmsy"] = mcTS[,,"SSB"] / avgRP[iii,"Bmsy"]
		avgTS[iii,dimnames(mcTS)$yr,"Ft"]     = mcTS[,,"F"]
		avgTS[iii,dimnames(mcTS)$yr,"FtFmsy"] = mcTS[,,"F"] / avgRP[iii,"Fmsy"]
		avgTS[iii,dimnames(mcTS)$yr,"ut"]     = 1 - exp(-mcTS[,,"F"])
		avgTS[iii,dimnames(mcTS)$yr,"utumsy"] = (1 - exp(-mcTS[,,"F"])) / avgRP[iii,"umsy"]
		avgTS[iii,dimnames(mcTS)$yr,"Rt"]     = mcTS[,,"Recr"]
		avgTS[iii,dimnames(mcTS)$yr,"Rtdev"]  = mcTS[,,"RecrDev"]
#browser();return()

		for (n in 1:length(cp)){
			ccc = cp[n]
#if (n==2){browser();return()}
			avgCP[as.character(pyrs),mm,ccc] =  mcCP[[ccc]]

			## Populate projection array
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"Bt",ccc]     = mcPJ[[ccc]][,,"SSB"]
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"BtB0",ccc]   = mcPJ[[ccc]][,,"SSB"] / avgRP[iii,"B0"]
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"BtBmsy",ccc] = mcPJ[[ccc]][,,"SSB"] / avgRP[iii,"Bmsy"]
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"Ft",ccc]     = mcPJ[[ccc]][,,"F"]
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"FtFmsy",ccc] = mcPJ[[ccc]][,,"F"] / avgRP[iii,"Fmsy"]
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"ut",ccc]     = 1 - exp(-mcPJ[[ccc]][,,"F"])
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"utumsy",ccc] = (1 - exp(-mcPJ[[ccc]][,,"F"])) / avgRP[iii,"umsy"]
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"Rt",ccc]     = mcPJ[[ccc]][,,"Recr"]
			avgPJ[iii,dimnames(mcPJ[[ccc]])$yr,"Rtdev",ccc]  = mcPJ[[ccc]][,,"RecrDev"]
		}
	}
	storage.mode(mpdPA)="double"
	rownames(mpdPA) = names(mclst); colnames(mpdPA) = valPA  ## because of occasional extra parameters
	if (type=="senso")
		out = list(senPA=avgPA, senRP=avgRP, senTS=avgTS, senPJ=avgPJ, senCP=avgCP, smpdPA=mpdPA, senLL=avgLL)
	else
		out = list(avgPA=avgPA, avgRP=avgRP, avgTS=avgTS, avgPJ=avgPJ, avgCP=avgCP, ampdPA=mpdPA, avgLL=avgLL)
#browser();return()
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~gatherMCMC


## getNpan------------------------------2019-05-10
##  Get panel number when inside a multi-panel plot.
## ---------------------------------------------RH
getNpan = function()
{
	mfg=par()$mfg
	mfg[2]+(mfg[1]-1)*mfg[4]
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getNpan


## getSS.control------------------------2020-11-10
##  Get select values out of the control file.
## (May or may not continue using this fn)
## ---------------------------------------------RH
getSS.control = function(fnam="control.ss")
{
	octl = list()
	if (!file.exists(fnam))
		stop (paste0("Control file ", fnam, " does not exist!"))
	rctl = readLines(fnam)
	## Maturity
	i  = grep("maturity_option",rctl)
	ii = as.numeric(substring(rctl[i],1,1))
	octl[["maturity_option"]] = ii
	if (ii==3) {
		iii = strsplit(rctl[i+1],split=" +")[[1]]
		octl[["maturity_ogive"]] = as.numeric(grep("^#",iii,invert=TRUE,value=TRUE))
	}
	return(octl)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getSS.control
#junk = getSS.control("./BSR/models/zanadu/rebsn.control.46.01.ss")


## importCor----------------------------2020-09-28
## Import SS parameter correlations.
##  (Adapted from 'importCor' in PBSawatea.)
## ---------------------------------------------RH
importCor = function(cor.file) {
	cfile  = readLines(cor.file)
	header = cfile[1]
	cmat   = cfile[-1]
	hessian = as.numeric(strsplit(header,split=" = ")[[1]][2])
	cmat = sub("std dev","std.dev",cmat)
	cmat = gsub("]","",gsub("\\[",".",cmat))
	nI = length(cmat)-1 # number of indices
	for (i in 1:nI) {
		ii = i + 1
		cmat[ii] = paste(c(cmat[ii],rep("NA",nI-i)),collapse=" ")
	}
	writeLines(cmat,"cor.tmp")
	cor = read.table("cor.tmp",header=TRUE,sep="")
	iii = paste("i",pad0(1:nI,n=ceiling(log10(nI))),sep="")
	colnames(cor)[grep("X",colnames(cor))] = iii
	row.names(cor) = iii
	zi = upper.tri(cor[iii,iii])
	cor[iii,iii][zi] = t(cor[iii,iii])[zi] # populate upper right triangle with lower left values
	cor.mat = as.matrix(cor[iii,iii])
	cor.name = cor$name
	cor.value = cor$value
	cor.std.dev = cor$std.dev
#browser();return()
	out = list(cfile=cfile, cmat=cmat, cor=cor, cor.mat=cor.mat, index=iii, cor.name=cor.name, cor.value=cor.value, cor.std.dev=cor.std.dev, hessian_log_determinant=hessian)
	file.remove("cor.tmp")
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importCor
#test=importCor("ss.cor")


## importEva----------------------------2020-09-28
## Import SS Hessian eigenvlaues.
##  (Adapted from 'importEva' in PBSawatea.)
## ---------------------------------------------RH
importEva = function(eva.file)
{
	out   = list()
	efile = readLines(eva.file)
	out[["efile"]] = efile
	eline = sub(":\\t","",sub("^ +","",efile))
	elist = strsplit(eline,split="[[:space:]]")
	names(elist) = sapply(elist,head,1)
	edata = as.data.frame(sapply(elist,function(x){as.numeric(x[-1])}))
	colnames(edata) = names(elist)
	#efile = as.numeric(read.table(eva.file,header=FALSE,sep=""))
	out[["eva"]] = edata
#browser();return()
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importEva
#test=importEva("ss.eva")


## importLik----------------------------2012-08-08
## Import Awatea likelihoods.
## ---------------------------------------------RH
#importLik = function(lik.file)
#{
#	lfile = readLines(lik.file)
#	out   = list(lfile=lfile)
#	liks  = lfile[!is.element(lfile,c("","**Likelihoods**"))]
#	liks  = gsub("@","",gsub("  "," ",gsub("   "," ",liks)))
#	liksplit = strsplit(liks,split=" ")
#	liknams  = sapply(liksplit,function(x){x[1]})
#	likvals  = sapply(liksplit,function(x,nf){
#		mess  = paste("assign(\"",x[1],"\",c(",paste(x[-1],collapse=","),"),envir=sys.frame(which=nf))",sep="")
#		eval(parse(text=mess))}, nf=sys.nframe())
#	for (i in liknams) out[[i]] = get(i)
#	return(out)
#}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importLik
#test=importLik("likelihood.dat")


## importPar----------------------------2020-09-24
##  Import all SS parameters.
##  (Adapted from 'importPar' in PBSawatea.)
## ---------------------------------------------RH
importPar = function(par.file) 
{
	pfile  = readLines(par.file)
	out    = list(pfile=pfile)
	header = pfile[1]
	pfile  = pfile[-1]
	mess = gsub(" +","_",gsub(" += +","=",sub(" +Max",";Max",sub(" +Obj",";Obj",sub("#[_ +]","",header)))))
	mess = sub("Maximum_gradient_component","maxgrad",mess)
	mess = sub("Number_of_parameters","npars",mess)
	mess = sub("Objective_function_value","fval",mess)
	eval(parse(text=mess))
	out = c(out, list(npars=npars,fval=fval,maxgrad=maxgrad))
	vnampos = grep("#",pfile)
	vvalbeg = vnampos + 1
	vvalend = vnampos+c(diff(vnampos)-1,length(pfile)-rev(vnampos)[1])
	pfile = gsub("]","",gsub("\\[",".",pfile))
#browser();return()
	pfile = gsub(":","",gsub("# ","",pfile))
	pfile = gsub(" ",",",gsub("^ ","",pfile))
	vals  = sapply(1:length(vnampos),function(x,f,n,v1,v2){
		#mess = paste("out[\"",f[n[x]],"\"]= c(", paste(f[v1[x]:v2[x]],collapse=","),")",sep="")
		mess = paste(f[n[x]]," = c(", paste(f[v1[x]:v2[x]],collapse=","),")",sep="")
		eval(parse(text=mess))
	},f=pfile,n=vnampos,v1=vvalbeg,v2=vvalend)
	names(vals)=pfile[vnampos]
	out=c(out,vals)
#browser();return()
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importPar
#test=importPar("ss.par")


## importStd----------------------------2020-09-28
## Import SS table of standard deviations.
##  (Adapted from 'importStd' in PBSawatea.)
## ---------------------------------------------RH
importStd = function(std.file, vnam="name") 
{
	out   = list()
	sfile = readLines(std.file)
	out[["sfile"]] = sfile
	sfile = sub("std dev","std.dev",sfile)
	sfile = gsub("]","",gsub("\\[",".",sfile))
	writeLines(sfile,"std.tmp")
	std   = read.table("std.tmp",header=TRUE,sep="")
	out[["std"]] = std
	onam  = setdiff(names(std),vnam)
	vars  = unique(std[,vnam])
	for (v in vars) {
		mess = paste(v,"=std[is.element(std[,vnam],\"",v,"\"),onam,drop=FALSE]",sep="")
		mess = c(mess, paste("out[\"",v,"\"] = list(",v,"=",v,")",sep=""))
		eval(parse(text=paste(mess,collapse=";")))
	}
	file.remove("std.tmp")
#browser();return()
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~importStd
#test=importStd("ss.std")


## is.numStr----------------------------2020-09-20
##  Check if strings can be converted to numerics.
##----------------------------------------------RH
is.numStr = function(x)
{
	out = sapply(x, function(xx) {
		xx = as.character(xx)
		all(grepl("[[:digit:]]|\\-|\\.|[eE]|[[:space:]]", strsplit(xx,"")[[1]]))
	})
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~is.numStr

## to give median (5%-95%) to put in text.
med5.95 = function(xx.MCMC, dig=0, quants3=tcall(quants3)){  ## dig is number of dec places
	print(paste0( c( prettyNum(round(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=options()$big.mark),
		"~(", prettyNum(round(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=options()$big.mark), 
		"-", prettyNum(round(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=options()$big.mark), ")"), collapse=""))
}


## prepCP-------------------------------2022-06-30
##  Prepare Catch Policies -- 'CC'=constant catch, 'HR'=harvest rate (not implemented yet)
##  fyrs includes the current year (e.g., 2022=beginning of 2022, end of 2021), which is treated as a projection
##  Update for Canary Rockfish 2022
## ---------------------------------------------RH
prepCP = function(run.rwt, cp=list('1'=800), d.cp="CC", tag="", #tag=".nuts4K", 
	fyrs=2023:2033, season=1, fleet=1:2, ## will need to alter if more than one fleet
	d.base = "C:/Users/haighr/Files/GFish/PSARC22/CAR/Data/SS/CAR2022",
	w=NULL, cvpro=NULL)
{
	padded = pad0(run.rwt,2)
	t.run=padded[1]; t.rwt=padded[2]
	tar.run.rwt = paste(t.run,t.rwt,sep=".")
	d.target = file.path(d.base, paste0("Run", t.run), paste0("MCMC.",tar.run.rwt,tag))
	poop = function(){ setwd(d.target); gdump = gc(verbose=FALSE) }
	on.exit(poop())
	
	## Check for same number of catch policies by fleet
	n.pol = sort(unique(sapply(cp,length)))
	if (length(n.pol)>1) stop("Revise catch policy input list to have equal numbers by fleet")

	need.files = c("starter.ss","forecast.ss",paste0(c("control","data"),".",tar.run.rwt,".ss"), "ss.par","ss.cor","ss.psv", paste0("admodel.",c("hes","cov")))
	for (i in 1:n.pol) {
		ii = sapply(cp,function(x,npol){return(x[i])},npol=i)
		iii = pad0(i,2)
		#d.catch.policy = file.path(d.target, d.cp, pad0(i,4))
		d.policy.type = file.path(d.target, d.cp)
		if (!file.exists(d.policy.type))
			dir.create(d.policy.type)
		d.catch.policy = file.path(d.policy.type, paste0("CP",iii))
		if (!file.exists(d.catch.policy)){
			dir.create(d.catch.policy)
		} else {
			clearFiles( setdiff(list.files(d.catch.policy,full.names=T),list.dirs(d.catch.policy)) )  ## need full names and explicit exclusion of directories
		}
		pee = file.copy(from=file.path(d.target,need.files), to=d.catch.policy, copy.date=T, overwrite=T)
		
		## Modify forecast file to reflect new catch policies
		forecast = readLines(file.path(d.catch.policy,"forecast.ss"))
#browser();return()

		fline    = grep("#_fcast_nyrs",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,(regexpr(" #_fcast_nyrs",forecast[fline])-1)), length(fyrs), forecast[fline])
		fline    = grep("#_ctl_rule_ul",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,(regexpr(" #_ctl_rule_ul",forecast[fline])-1)), 0.002, forecast[fline])
		fline    = grep("#_ctl_rule_ll",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,(regexpr(" #_ctl_rule_ll",forecast[fline])-1)), 0.001, forecast[fline])
#_ctl_rule_ul
		fline    = grep("#_fcast_yr1",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,4), fyrs[2] + 100, forecast[fline])
		f1     = grep("#_Yr Seas Fleet Catch",forecast)
		f2     = grep("#_end_catpols",forecast)
		nyrs   = length(fyrs)
		nfleet = length(fleet)
		newpol = data.frame(Yr=rep(fyrs,nfleet), Seas=rep(season,nyrs*length(fleet)), Fleet=rep(fleet,each=nyrs), Catch=rep(ii,each=nyrs))
		foopol = apply(newpol,1,paste0,collapse=" ")
		forecast = c(forecast[1:f1], foopol, forecast[f2:length(forecast)])
		writeLines(forecast, con=file.path(d.catch.policy,"forecast.ss"))
#		};
#browser();return()
#	{
		.flush.cat(paste0("CP = ", i, " -- run mceval to generate catch policies.\n"))
		setwd(d.catch.policy)
		gdump = gc(verbose=FALSE)
		syssy = system("ss -mceval", intern=TRUE)
		writeLines(syssy,paste0("./sys", pad0(i,4), ".txt"))
	}
	.flush.cat("C'est toute, folks\n")
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~prepCP


## prepMPD------------------------------2022-02-17
##  Prepare MPD runs for likelihood analysis
## ---------------------------------------------RH
prepMPD = function(run.rwt, #n.run, n.rwt, l.run, l.rwt,
	d.base = "C:/Users/haighr/Files/GFish/PSARC22/CAR/Data/SS/CAR2022",
	w=NULL, cvpro=NULL)
{
	padded = pad0(run.rwt,2)
	n.run=padded[1]; n.rwt=padded[2]; l.run=padded[3]; l.rwt=padded[4]
	pre.run.rwt = paste(l.run,l.rwt,sep=".")
	new.run.rwt = paste(n.run,n.rwt,sep=".")
	## Destination directories
	#d.cwd = getwd(); on.exit(setwd(d.cwd))
	d.run = file.path(d.base,paste0("Run",n.run))
	d.mpd = file.path(d.run,paste("MPD",new.run.rwt,sep="."))

	## Previous run from which to poach ssfiles
	if (run.rwt[3]==0 && run.rwt[4]==0) {
		## Preliminary run so get files from d.base
		p.run = p.mpd = d.base
	} else {
		p.run = file.path(d.base,paste0("Run",l.run))
		p.mpd = file.path(p.run,paste("MPD",pre.run.rwt,sep="."))
	}
	if (!all(file.exists(c(p.run,p.mpd))))
		stop("Previous Run directory and/or MPD directory does not exist")
	if (!file.exists(d.run)) dir.create(d.run)
	if (!file.exists(d.mpd)) dir.create(d.mpd)
#browser();return()
	p.ss  = setdiff(list.files(p.mpd, pattern="\\.ss$"),"runnumber.ss")
	pee   = file.copy(from=file.path(p.mpd,p.ss), to=d.mpd, copy.date=T, overwrite=T)
	poo   = file.rename(from=file.path(d.mpd,paste(c("data","control"),pre.run.rwt,"ss",sep=".")), to=file.path(d.mpd,paste(c("data","control"),new.run.rwt,"ss",sep=".")))

	## Modify starter file to reflect new run
	starter = readLines(file.path(d.mpd,"starter.ss"))
	dline   = grep(paste("data",pre.run.rwt,sep="."),starter)
	cline   = grep(paste("control",pre.run.rwt,sep="."),starter)
	starter[dline] = sub(paste0("data\\.",l.run,"\\.",l.rwt),paste0("data.",new.run.rwt),starter[dline])
	starter[cline] = sub(paste0("control\\.",l.run,"\\.",l.rwt),paste0("control.",new.run.rwt),starter[cline])
#browser();return()
	writeLines(starter, con=file.path(d.mpd,"starter.ss"))

	## Modify data file to add cvpro using the Francis method
	data = readLines(file.path(d.mpd,paste0("data.",new.run.rwt,".ss")))
	if (!is.null(cvpro)){
		vline = intersect(grep("_index$",data), grep("^#",data,invert=TRUE))
		vbits = strsplit(data[vline], split=" +")
		if (length(.su(sapply(vbits,function(x){x[3]}))) != length(cvpro)){
			.flush.cat("User inputs for cvpro do not match indices in data file\n")#; browser();return()
		}
		vdump = sapply(1:length(vbits),function(i){vbits[[i]][5] <<- sqrt(as.numeric(vbits[[i]][5])^2 + cvpro[vbits[[i]][3]]^2) })
		vnew  = sapply(vbits,function(x){ paste0(x,collapse=" ")})
		data[vline] = vnew
		writeLines(data, con=file.path(d.mpd,paste0("data.",new.run.rwt,".ss")))
	}
	## Modify control file to add Francis reweights
	control = readLines(file.path(d.mpd,paste0("control.",new.run.rwt,".ss")))
	if (!is.null(w)){
		fline = intersect(grep("vadj_af",control), grep("^#",control,invert=TRUE))
#browser();return()
		if (length(fline)!=length(w)){
			.flush.cat("Francis reweights do not match control lines\n"); browser();return()
		}
		fbits = strsplit(control[fline], split=" +")
		wdump = sapply(1:length(fbits),function(i){fbits[[i]][3] <<- w[fbits[[i]][2]] }) #w[i]})
		fnew  = sapply(fbits,function(x){ paste0(x,collapse=" ")})
		control[fline] = fnew
		writeLines(control, con=file.path(d.mpd,paste0("control.",new.run.rwt,".ss")))
	}
#	if (!is.null(w) || !is.null(cvpro)){
#	}
	setwd(d.mpd)
	.flush.cat(paste0("***Edit the ss files in\n\t",d.mpd,"\nbefore proceeding with an MPD fit in SS."),"\n\n")
	return(invisible(list(starter=starter, control=control)))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~prepMPD


## ptab---------------------------------2022-02-18
##  Function to use for priors in table (adapted from PBSawatea).
## ---------------------------------------------RH
ptab = function(xx) {
	xx = sub("^\\s+", "", xx)  ## remove leading and trailing whitespace
	xlab = gsub("\\_+"," ",xx[1])
#browser();return()
	xnum =xx[-1]
	xnum[4] = switch(xnum[4], 'Normal'=6, 'No_prior'=0, 'Full_Beta'=2)
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
qtab = function(xx.MCMC, dig=0, quants3=tcall(quants3)) {  ## dig is number of dec places
	print(paste0( c( prettyNum(round(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=options()$big.mark),
		" & ", prettyNum(round(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=options()$big.mark),
		" & ", prettyNum(round(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=options()$big.mark)), collapse=""))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~qtab


## repeatMPD----------------------------2022-05-26
## Repeat MPDs for axes of uncertainty to visualise likelihood density.
## Need to run from a specific Likelihood directory withing RunNN.
## ---------------------------------------------RH
repeatMPD = function(Pfix=list(FFF=seq(0.05,0.06,0.01), MMM=seq(0.05,0.06,0.01), RRR=NULL), A=60, 
   dir=getwd(), dir.pro, prefix="control.MMM.", clean=FALSE, 
   strSpp="CAR", argsMPD="")
{
	## Start subfunctions
	## Determine number of decimal places with non-trailing zeroes
	## See user 'darocsig' @ https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r
	decimalplaces <- function(x) {
		dc = as.character(); dp = as.numeric()
		for (xx in x) {
#if(xx==0.10) {browser();return()}
			if (abs(xx - round(xx)) > .Machine$double.eps^0.5) {
				dc = c(dc, strsplit(sub('0+$', '', as.character(xx)), ".", fixed = TRUE)[[1]][2] )
				dp = c(dp, nchar(strsplit(sub('0+$', '', as.character(xx)), ".", fixed = TRUE)[[1]][2]) )
			} else {
				dc = c(dc, NA)
				dp = c(dp, 0)
			}
		}
		## Automatically determine minimum number of decimal places to yield unique values
		udp = max(dp) + 1
		for (i in max(dp):1){
			#print(length(unique(round(x,i)))==length(x))
			if (length(unique(round(x,i))) == length(x) && !any(round(x,i)==0))
				udp = udp - 1
		}
#browser() ;return()
		return(list(dc=dc, dp=dp, mindp=udp))
	}
	cleanup = function(){
		junkpat = c("^admodel","^ss\\.","\\.log$","\\.sso$","\\.ss_new$","\\.png$","\\.bak$","\\.old$","runnumber\\.ss")
		junkit  = sapply(junkpat,function(x){list.files(pattern=x)})
		junkit  = sapply(junkit,setdiff,"SS.exe")
		junk    = sapply(junkit,function(x){ if (length(x)>0) for (i in x) if (file.exists(i)) file.remove(i)})
	}
	## End subfunctions

	ncA   = max(nchar(A))
	Pfix = Pfix[ sapply(Pfix,function(x){!is.null(x) & all(!is.na(x)) }) ]
	if ( !all(diff(sapply(Pfix,length))==0) )
		stop("Pfix list elements must be vectors of equal length")

	Y = Ydec = Yval = Pfix
	Yvar  = names(Pfix)
	Ypref = paste0(substring(Yvar,1,2),"\\.")
	dcY   = lapply(Pfix, decimalplaces)
	dpY   = sapply(dcY,function(x){x$mindp})
	ipY   = sapply(Y,function(x){ max(0,ceiling(log10(x + 0.00001))) })  ## integer places

#browser() ;return()

	#if (!is.null(M)) {
	#	Y = Ydec = Yval = M; Ypref="MM\\."; Yvar="MMM"
	#} else if (!is.null(R0)) {  ## R0 expressed as log(R0)
	#	Y = R0; Ydec=Y/100; Yval=exp(Y); Ypref="RR\\."; Yvar="RRR"
	#}
	#
	#dcY = decimalplaces(Ydec)
	#dpY = dcY$mindp

	keep  = c("admodel.hes",paste0(c("Report","CompReport","covar","warning"),".sso"), paste("ss",c("cor","eva","par","rep","std","mpd"),sep="."), "data.ss_new") ## MPD files to save
	Nprof = length(Y[[1]]) * length(A)
	Iprof = 0
	if (missing(dir.pro)){
		fornow  = lapply(Y,function(x){ paste0("(", paste0(range(x),collapse="-"),")")})
		dir.sub =  paste0(c("Prof",paste0(substring(names(fornow),1,1),fornow)),collapse=".")
		#dir.sub =  paste0("pro.M(",paste0(range(M),collapse="-"),").A(",paste0(A,collapse=","),")")
		dir.pro = file.path(dir,dir.sub)
	} else {
		if (basename(dir.pro)==dir.pro || !grepl("^\\./|^[[:alpha:]]:/",dir.pro))  ## accept proper relative paths or those starting from a root
			dir.pro = file.path(dir,sub("^/","",dir.pro))
	}

	if (!file.exists(dir)) dir.create(dir)
	if (!file.exists(dir.pro)) dir.create(dir.pro)

	for (a in A){
		## Look for ss files
		afile = file.path(dir,paste0(prefix,"A",pad0(a,ncA),".ss"))
		sfile = file.path(dir,c("starter.base.ss", "forecast.ss", paste0(sub("Like","data",basename(dir)),".ss")))
		for (bfile in c(afile,sfile)) {
			if(!file.exists(bfile)) {
				if (file.exists(sub("/Like\\.","/MPD.",bfile)))
					file.copy(sub("/Like\\.","/MPD.",bfile), dir, copy.date=T)
				else
					stop(paste0("Cannot find file\n\t'", bfile, "'"))
			}
		}
		aline = readLines(afile)
		#grep(paste0(Yvar,collapse="|"),aline)

		for (i in 1:Nprof){
			Iprof = Iprof + 1
			yval = sapply(Y,function(x){x[i]})
			ychr = sapply(Ydec,function(x){x[i]})
			ychr = sapply(1:length(ychr),function(f){
				if(ipY[f]==0) ochr = show0(x=round(ychr[f],dpY[f]), n=dpY[f]) 
				else          ochr = pad0(x=ychr[f], n=ipY[f]+dpY[f], f=dpY[f])
			return(ochr)
#browser() ;return()
			})

## ipref need serious debugging
			#ipref = paste0(sub(Ypref,paste0(sub("0\\.","",show0(ychr,dpY)),"."),prefix),"A",pad0(a,ncA))
			ipref = paste0("control.",substring(gsub("control|\\.","",prefix),1,1), paste0(sub("^0\\.","",ychr),collapse=""),".A",pad0(a,ncA),collapse="")
#browser() ;return()
			.flush.cat(paste0("Processing run '", ipref, "' ..."), "\n")
			ifile = paste0(ipref,".ss")
			## for some reason, certain combos of M & A do not converge unless M is nudged
			if (strSpp=="WWR" && M==0.03 && A==45) smudge = 0.00001 ## necessary for WWR
			else smudge = 0 
			iline = aline
			for (l in 1:length(yval))
				iline = gsub(Yvar[l], yval[l]+smudge, iline)
			writeLines(iline, con=file.path(dir,ifile))

			## Modify starter file to reflect new run
			starter = readLines(file.path(dir,"starter.base.ss"))
			#dline   = grep(paste("data",pre.run.rwt,sep="."),starter)
			cline   = grep("^control",starter)
			#starter[dline] = sub(paste0("data\\.",l.run,"\\.",l.rwt),paste0("data.",new.run.rwt),starter[dline])
			starter[cline] = sub("^control.+\\.ss",paste0(ipref,".ss"),starter[cline])
			writeLines(starter, con=file.path(dir,"starter.ss"))
#browser();return()

			#expr = paste("mess = shell(cmd=\"awatea -ind ",mfile,argsMPD,"\", wait=TRUE, intern=TRUE)",sep="")
			setwd(dir)
			if (clean) cleanup()
			expr = paste("mpd = shell(cmd=\"ss\", wait=TRUE, intern=TRUE)",sep="")
			.flush.cat("   ", expr, "\n")
			eval(parse(text=expr))
			writeLines(mpd ,con=file.path(dir,"ss.mpd"))

			for (jfile in keep){
				if (file.exists(file.path(dir,jfile))){
					if (jfile=="ss.par"){
						## r4ss' SS_profile seems to rename ss.par; repeat here but still copy ss.par below
						kfile = paste0(jfile,"_",Iprof,".sso")
						file.copy(from=file.path(dir,jfile), to=file.path(dir.pro,kfile), overwrite=TRUE, copy.date=TRUE)
					}
					#else {
					jbits = strsplit(jfile,"\\.")[[1]]
					#jbits = c(jbits,"sumtingwong")  ## just for testing
					jlen  = length(jbits)-1
					jbits[jlen] = paste0(jbits[jlen],Iprof)
					kfile = paste0(jbits,collapse=".")
					#}
					file.copy(from=file.path(dir,jfile), to=file.path(dir.pro,kfile), overwrite=TRUE, copy.date=TRUE)
				}
			}
			if (clean) cleanup()
			rubbish = gc(verbose=FALSE)
		} ## end m loop (natural mortality)
	} ## end a loop (maximum age for plus class)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~repeatMPD


## runADnuts----------------------------2022-06-07
##  Based on Chris Grandis' R code that was used to manually run adnuts.
##  Originally see Cole Monnahan:
##    https://rdrr.io/github/Cole-Monnahan-NOAA/adnuts/src/R/samplers.R
##  Note: adnuts needs serious path intervention (RH 201021)
## ------------------------------------------CG|RH
runADnuts.off = function(path=getwd(), run="01.01", model="base", 
   rda_name, add2rda=FALSE, tag="REBSN", 
   parallel=TRUE, cores=6, syscall=c(FALSE,rep(TRUE,5)),
   iters=list(test=15,pilot=100,mle=500,update=500,pcburn=0.25),
   ssexe="ss", ssdir, ssfiles=c("control.ss","data.ss"), stupid_shiny=FALSE)
{
	## Check to see if file exists and keep a copy if it does before removing it.
	sweep.files = function(x) {
		rubbish = 0
		for (i in x) {
			if (file.exists(i)) {
				fstamp = paste0(sub("[[:space:]]","(",gsub("[[:punct:]]","",substring(file.info(i)$ctime,3,16))),")")
				ext    = tools::file_ext(i)
				pre    = gsub(paste0("\\.",ext,"$"),"",i)
				backup = paste0(pre,"-",fstamp,".",ext)
				if (!file.exists(backup))
					file.copy(from=i, to=backup, overwrite=TRUE, copy.date=TRUE)
				file.remove(i)
				rubbish = rubbish +1
			}
		}
		invisible(paste0(rubbish, " files removed"))
	}
	## Generate dispersed initial values from MLE estimates
	my_sample_inits = function(path, cv=0.025, reps) {
		fit   = .read_mle_fit("ss", path=path)  ##
		mle   = fit$est[1:fit$nopar]##; names(mle) = fit$par.names
		inits = lapply(1:reps, function(i) {
			ival = sapply(mle,function(x) { rnorm(1, mean=x, sd=abs(cv*x)) })
		})
		names(inits) = 1:reps
		return(inits)
	}

	start_time <- Sys.time()
	cwd  = getwd()
	poop = function(){ setwd(cwd); gdump = gc(verbose=FALSE) }
	on.exit(poop())
	setwd(path); d.path=getwd(); setwd(d.path)  ## Just in case a relative path is specified

	d.run  = file.path(path, run)
	if (!file.exists(d.run)) 
		stop(paste0("Need to start with an existing Run\n\t",d.run))
	d.mpd    = file.path(d.run, paste0("MPD.",model))
	if (!file.exists(d.mpd)) 
		stop(paste0("Need to start with an existing MPD\n\t",d.mpd))
#browser();return()
	file.copy(from=paste0(d.mpd,"/starter.ss"), to=paste0(d.mpd,"/starter.adnuts.ss"), overwrite=TRUE, copy.date=TRUE)
	d.nuts   = file.path(d.run, paste0("NUTS.",model))
	if (!file.exists(d.nuts)) 
		dir.create(d.nuts)
	d.mcmc   = file.path(d.run, paste0("MCMC.",model))
	if (!file.exists(d.mcmc)) 
		dir.create(d.mcmc)
	poop = function(){setwd(d.mcmc)}  ## set to MCMC directory
	if(missing(rda_name))
		rda_name = paste0(run,".nuts.",model,".rda")
		#rda_name = paste0(run,".",model,".rda")
	d.rda = file.path(d.run, rda_name)  ## store rda files in the run directory

	ssfiles = unique(c("starter.ss","forecast.ss",ssfiles))
	if (missing(ssdir)) ssdir = d.mpd
	if ( all(file.exists(paste0(ssdir,"/",ssfiles)))) {
		file.copy(from=paste0(ssdir,"/",ssfiles), to=paste0(d.nuts,"/",ssfiles), overwrite=TRUE, copy.date=TRUE)
		file.copy(from=paste0(ssdir,"/",ssfiles), to=paste0(d.mcmc,"/",ssfiles),  overwrite=TRUE, copy.date=TRUE)
		starter = readLines(file.path(d.mpd,"starter.adnuts.ss"))
		starter[grep("#_mcmc_burn",starter)] = sub("^[[:digit:]] #_mcmc_burn", paste0(0," #_mcmc_burn"), grep("#_mcmc_burn",starter,value=T))
		for (icon in file.path(c(d.nuts,d.mcmc),"starter.ss"))
			writeLines(starter, con=icon)
	} else
		stop("SS files: '", paste0(ssfiles,collapse="', '"),"'\n\tnot found in ", ssdir)
	
	mess  = "Hessian does not appear to be positive definite"
	dash  = paste0(rep("-",62),collapse="")
#browser();return()

	setwd(d.mpd)
	## ----------------------------------------------------------------
	if (syscall[1]) {
		## This step is needed to create files like 'covar.sso' and 'ss_summary.sso'
		## Should be done already but user might want to run again
#browser();return()
		if (ssdir!=d.mpd)
			file.copy(from=paste0(ssdir,"/",ssfiles), to=paste0(d.mpd,"/",ssfiles),   overwrite=TRUE, copy.date=TRUE)
		.flush.cat("SysCall 1 -- Running MPD for best fits...\n")
		#toast = setdiff(list.files("."),c(ssfiles))
		toast = setdiff(list.files("."),c(ssfile, basename(list.dirs("."))))  ## list files always sees subdirectories as files
		if (length(toast)>0)
			file.remove(toast)
		sys1 = system(paste0(ssexe), intern=TRUE)
		writeLines(sys1,"./sys1.txt")
		if (any(grepl(mess,sys1)))
			stop(paste0("\n",dash,"\nSystem call 1: ",mess,"\n",dash))
	}
	#file.remove(list.files(d.mpd, pattern="\\.log$"))  ## for now
	replist = SS_output(dir=d.mpd, verbose=FALSE)
#browser();return()
	ttput(replist)
	vals.save = c("replist")

	setwd(d.nuts)

	## Chains to run in parallel
	reps  <- min(cores, parallel::detectCores() - 1)
	set.seed(2022)
	seeds <- sample(1:1e4, size = reps)
	sweep.files(d.rda)
#browser();return()

	## ----------------------------------------------------------------
	if (syscall[2]) {
		## This step is needed before syscall 3
		.flush.cat("SysCall 2 -- Running precursor mcmc test...\n")
		fdump = file.remove(setdiff(list.files(),ssfiles))  ## Start fresh if syscall 2
		gdump = gc(verbose=FALSE)
		sys2 = system(paste0(ssexe, " -nox -iprint 200 -mcmc ", iters$test), intern=TRUE)
		writeLines(sys2,"./sys2.txt")
		if (any(grepl(mess,sys2)))
			stop(paste0("\n",dash,"\nSystem call 2: ",mess,"\n",dash))
	}
	## Then run parallel RWM chains as a first test to ensure
	## mcmc itself is working properly, or that model is converging in mcmc space
	thin   <- 10
	iter   <- iters$pilot * thin ## iter is per core !!!
	#warmup <- ceiling(iter/4)
	warmup <- ceiling(0.25 * iter)
	inits  <- NULL               ## start chains from MLE
	inits  <- my_sample_inits(path=d.mpd, reps=reps)
	vals.save = c(vals.save,c("seeds","iters"))

	## ----------------------------------------------------------------
	if (syscall[3]){
		.flush.cat("SysCall 3 -- Running pilot RWM chains...\n")
		gdump = gc(verbose=FALSE)
		pilot  <- sample_admb(model=basename(ssexe), path=d.nuts, iter=iter, init=inits, chains=reps, warmup=warmup, seeds=seeds, thin=thin, mceval=FALSE, duration=NULL, parallel=parallel, cores=reps, control=NULL, algorithm="RWM")
		ttput(pilot)
		save("pilot","inits","warmup","thin","seeds", file="pilot.rda")
	} else {
		ttget(pilot)
		if (!exists("pilot", inherits=F))
			load("pilot.rda")
	}

	## Check convergence and slow mixing parameters
#mon  <- rstan::monitor(pilot$samples, warmup=pilot$warmup, print=FALSE)
#ttput(mon)
	## max(mon[,'Rhat'])
	## min(mon[,'n_eff'])
	## Examine the slowest mixing parameters
#slow <- names(sort(mon[,"n_eff"]))[1:10]
#ttput(slow)
#vals.save = c(vals.save,c("pilot","mon","slow"))

#browser();return()
	#so("pairs_admb.r","synth")
	#pairs_admb(fit=pilot, pars=slow, label.cex=1)
	#pairs_admb(fit = pilot, pars = c("MGparm[1]", "SR_parm[1]", "SR_parm[2]")) ## must be specific to Hake
	#pairs_admb(fit = pilot, pars = c("SR_parm[1]", "Q_parm[1]", "selparm[1]"), label.cex=1)
	
	## After regularizing run NUTS chains. First reoptimize to get the
	## correct mass matrix for NUTS. Note the -hbf 1 argument. This is a
	## technical requirement b/c NUTS uses a different set of bounding
	## functions and thus the mass matrix will be different.
	## ----------------------------------------------------------------
	if (syscall[4]) {
		.flush.cat("SysCall 4 -- Reoptimize to get correct mass matrix...\n")
		#mess  = sys4 = "Hessian does not appear to be positive definite"
		#nsys4 = 0; isys4 = iters$test
		#while (any(grepl(mess,sys4)) && nsys4<5){ 
		#	nsys4 = nsys4 + 1
		#	iters$test = isys4 * nsys4
		#	.flush.cat(paste0("  System call 4, loop ", nsys4,": trying ", iters$test, " iterations\n"))
		gdump = gc(verbose=FALSE)
		sys4  = system(paste0(ssexe, " -hbf 1 -nox -iprint 100 -mcmc ", iters$test), intern=TRUE) 
		writeLines(sys4,"./sys4.txt")
		if (any(grepl(mess,sys4))) {
			message( paste0("\n",dash,"\nSystem call 4: ",mess,"\n",dash) )
			#replist <- SS_output(dir=d.nuts, verbose=F, printstats=F)
			#rep.names = row.names(replist$parameters[!is.na(replist$parameters[,"Active_Cnt"]),])
			fit <- .read_mle_fit("ss", path=d.nuts)
			par.names = fit$par.names; #names(par.names) = rep.names
			hes <- adnuts:::.getADMBHessian(path=d.nuts)
			ev  <- eigen(hes)
			hdiag = diag(hes); names(hdiag) = fit$par.names
			zbad  = hdiag < 0
			.flush.cat("Paramaters with concave (negative) 2nd derivatives :\n")
			.flush.cat( paste0(paste0("\t",sapply( grep(T,zbad), function(i) { paste0(names(hdiag)[i], " : ",hdiag[i], collapse="") }), "\n"),collapse="") )
			browser();return()
			#stop(paste0("\n",dash,"\nSystem call 4: ",mess,"\n",dash))
		}
		#}
#browser();return()
	}

	## Use default MLE covariance (mass matrix) and short parallel NUTS chains started from the MLE.
	## ----------------------------------------------------------------
	if (syscall[5]) {
		.flush.cat("SysCall 5 -- Using MLE covariance and short parallel NUTS chains...\n")
		sweep.files(c("ss.psv","unbounded.csv"))
		warmup <- ceiling(iters$pcburn * iters$mle)  ## Monnahan et al. (2019)
		gdump  <- gc(verbose=FALSE)
		inits  <- my_sample_inits(path=d.mpd, reps=reps)

##=== Break here to run on the command line ===
browser();return()
		nuts.mle <-  sample_admb(model=ssexe, path=d.nuts, iter=iters$mle, init=inits, chains=reps, warmup=warmup, seeds=seeds, thin=1, mceval=FALSE, duration=NULL, parallel=parallel, cores=reps, control=list(metric="mle", adapt_delta=0.9), algorithm="NUTS")
		ttput(nuts.mle)
		save("nuts.mle","inits","warmup","seeds", file="nuts.mle.rda")
	} else {
		ttget(nuts.mle)
		if (!exists("nuts.mle", inherits=F))
			load("nuts.mle.rda")
	}
	#vals.save = c(vals.save,c("nuts.mle"))
	## Check for issues like slow mixing, divergences, max treedepths with
	## ShinyStan and pairs_admb as above. Fix and rerun this part as needed.
	## launch_shinyadmb(nuts.mle)
#browser();return()
	
	## Once acceptable, run again for inference using updated mass matrix. Increase
	## adapt_delta toward 1 if you have divergences (runs will take longer).
	## Note this is in unbounded parameter space
	## The following, nuts.updated, was used for inferences in this appendix
	#vals.save = c(vals.save,c("nuts.mle"))
	mass  <- nuts.mle$covar.est
	inits <- sample_inits(nuts.mle, reps)
	## ----------------------------------------------------------------
	if (syscall[6]) {
		.flush.cat("SysCall 6 -- Run final NUTS MCMC using updated mass matrix...\n")
		setwd(d.mcmc)
		## Check to see if psv file exists and keep a copy if it does
		sweep.files(c("ss.psv","unbounded.csv"))
		#toast = setdiff(list.files("."),c(ssfiles,"mcmc",".RData"))
		toast = setdiff(list.files("."),c(ssfiles,"mcmc",".RData", basename(list.dirs("."))))  ## list files always sees subdirectories as files
#browser();return()
		if (length(toast)>0)
			file.remove(toast)
		warmup = ceiling(iters$pcburn * iters$update)  ## Monnahan et al. (2019)
		thin   = ifelse(is.null(iters$thin), 1, iters$thin)
		butter = c("admodel.cov","admodel.hes","ss.par","ss.cor")
		if ( all(file.exists(paste0(d.nuts,"/",butter))) )
			zb = file.copy(from=paste0(d.nuts,"/",butter), to=paste0(d.mcmc,"/",butter), overwrite=TRUE, copy.date=TRUE)
		else
			stop("MPD precursor directory is missin:\n\t'",paste0(butter[!xb],collapse="', '"),"'")
		gdump = gc(verbose=FALSE)
		nuts.updated <- sample_admb(model=ssexe, path=d.mcmc, iter=iters$update, init=inits, chains=reps, warmup=warmup, seeds=seeds, thin=thin, mceval=TRUE, duration=NULL, parallel=parallel, cores=reps, control=list(metric=mass, adapt_delta=0.9), algorithm="NUTS")
		ttput(nuts.updated)
		save("nuts.updated","inits","warmup","seeds","thin","mass", file="nuts.updated.rda")
	} else {
		ttget(nuts.updated)
		if (!exists("nuts.updated", inherits=F))
			load("nuts.updated.rda")
	}

	#vals.save = c(vals.save,c("nuts.updated","mass","inits"))
	tmpenv <- new.env()
	if (file.exists(d.rda) && add2rda) {
		load(d.rda, envir=tmpenv)
		oldtags = ls(envir=tmpenv)
	}
	mess = paste0(tag,"=list(",paste0(paste0(vals.save,"=",vals.save),collapse=","),"); ttput(", tag, "); tput(", tag, ",tenv=tmpenv); save(list=ls(tmpenv),file=\"",d.rda,"\", envir=tmpenv)")
	eval(parse(text=mess))
	end_time <- Sys.time()
	elapsed_time = end_time - start_time
	cat("Elapsed time: ", elapsed_time, " ", attributes(elapsed_time)$units, "\n")
	f.history = file.path(path, "runHistory.tex")
	if (file.exists(f.history))
		file.copy(from=f.history, to=file.path(d.run,"run.history.txt"), overwrite=TRUE, copy.date=TRUE)
	gdump = gc(verbose=FALSE)
	return(invisible(elapsed_time))

	## Have to get rid of column "lp__" to make shinyadmb work (bad programming in '.validate_sampler_params')
	#if (stupid_shiny) {
	#	fit = nuts.updated
	#	tmp_params = fit$sampler_params
	#	nms_params = lapply(tmp_params,function(x){x[,setdiff(colnames(x),"lp__")]})
	#	fit$sampler_params = nms_params
	#	ttput(fit)
	#	launch_shinyadmb(fit)
	#}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~runADnuts


## runSweave----------------------------2022-06-06
##  Run Sweave code to build pdfs for MPD and MCMC runs (not appendix)
## ---------------------------------------------RH
runSweave = function (d.model=getwd(), d.sweave, type="MPD", figs.only=FALSE)
{
	if (missing(d.sweave))
		d.sweave = "C:/Users/haighr/Files/Projects/R/Develop/PBSsynth/Authors/Rcode/develop"
	f.sweave = paste0(basename(d.model),ifelse(figs.only,".figs",""),".Rnw")
#browser();return()

	if (!file.copy(from=file.path(d.sweave, paste0("sweave",type,ifelse(figs.only,".figs",""),".Rnw")), to=file.path(d.model,f.sweave), overwrite=TRUE, copy.date=TRUE))
		stop("Failed to copy Sweave file to directory\n\t",d.model)

	## Other 'log' files (e.g., from Sweave) interfere with r4ss::SS_output
	if (file.exists(sub("\\.Rnw",".log",f.sweave)))
		file.remove(sub("\\.Rnw",".log",f.sweave))

	Sweave(file.path(d.model, f.sweave))
	mpd.sweave = shell(cmd=paste0("texify --pdf --synctex=1 --clean ",sub("\\.Rnw",".tex",f.sweave)), wait=TRUE, intern=TRUE)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~runSweave


## stab---------------------------------2021-04-19
## Quantile tabulation summary using significant digits
## --------------------------------------------AME
stab = function(xx.MCMC, dig=3, quants3=tcall(quants3), print=T) {  ## dig is number sig digits
	out = paste0( c( prettyNum(signif(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=options()$big.mark), 
		" & ", prettyNum(signif(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=options()$big.mark),
		" & ", prettyNum(signif(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=options()$big.mark)), collapse="")
	if (print) print(out)
	invisible(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~stab


## weightAF-----------------------------2022-02-24
##  Weight age frequencies using harmonic mean ratio method
## ---------------------------------------------RH
weightAF = function(replist, fleets, abase="agedbase", afld="Nsamp_adj",
   hfld="effN", rbase="recruit", rfld="pred_recr")
{
	agebase = replist[[abase]]
	fltbase = agebase[is.element(agebase$Fleet,fleets),]
	lstbase = split(fltbase, fltbase$Fleet)
	hmean  = sapply(lstbase,function(x){ 1/mean(1/x[,hfld]) })
	amean  = sapply(lstbase,function(x){ mean(x[,afld]) })
	wmean  = sapply(lstbase,function(x){ 1/mean(1/x[,hfld]) / mean(x[,afld]) })
#browser();return()

	## Recruitment predictions (flaky so do not use 'wadj')
	recbase = replist[[rbase]]
	Rmax = max(recbase[,rfld])
	#M = replist$parameters["NatM_p_1_Fem_GP_1","Value"]
	midx = grep("^Nat(.+)Fem",rownames(replist$parameters))
	M = replist$parameters[midx,"Value"]
	wsub = 1/((Rmax/10^floor(log10(Rmax)))*M^0.7)
	wadj = wmean / wsub

	w.mcallister = list(agedat=fltbase, hmean=hmean, amean=amean, wmean=wmean, Rmax=Rmax, M=M, wsub=wsub, wadj=wadj)
#browser();return()
	ttput(w.mcallister)
	return(w.mcallister)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~weightAF


## Flush the cat down the console (change to '.flash.cat' to avoid conflict with function in PBStools)
## Note: `.flush.cat' already in PBStools but this package is not assumed to be loaded.
.flash.cat = function(...) { cat(...); flush.console() }
