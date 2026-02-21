## ===============================================2024-08-27
## SS3 APPENDIX FUNCTIONS (Apps F & G)
## ------------------------------
## agileDT...............Produce agile decision tables for projections
## gatherMCMC............Gather MCMCs for functions below
## getSS.rdevs...........Get reruitment deviations from replist
## load_extra_mcmc.......Load extra MCMC information from Report files generated for every sample
## makeFSARfigs..........Make figures for the new FSAR
## plotSS.compo..........Make composite figures
## plotSS.senso..........Make sensitivity figures
## predictRec............Fit rockfish recruitment to environmental indices
## tabSS.compo...........Make base case tables
## tabSS.decision........Make decision tables
## tabSS.senso...........Make sensitivity tables
## ===============================================

## agileDT------------------------------2025-11-06
##  Produce agile decision tables for projections 
##  compared to reference points (using function 'findTarget')
##  Attempt to exapnd this to single-area models (RH 251024)
## ---------------------------------------------RH
agileDT <- function(compo, currYear=2026, projYear=2036, Ngen=3, gen1=25,
   area=c("5ABC","3CD","5DE"), cp=paste0("CC.", pad0(1:3,2)), 
   useRlow=FALSE, RPbase="B0", pdf=FALSE, 
   onepage=TRUE, tables.per.page=3, outnam, sigdig=2 )
{
	## Subfunctions ===============================
	## All B targets include area dimension
	subarray <- function(arr, ipos, inam, kpos) {  ## kpos = keep these dimension positions
		adim  = dim(arr)
		anam  = dimnames(arr)
		apos  = 1:length(adim)
		inam  = paste0("^",inam,"$")  ## to get exact match
		if (any(grepl(inam,anam[[ipos]]))) {
			indx = grep(inam,anam[[ipos]])
			if (length(indx)!=1) { message("Fix 'ipos'"); browser(); return() }
			if (missing(kpos)) kpos = apos[-ipos]
			nnam = anam[kpos]           ## new dimnames
			ndim = sapply(nnam, length) ## new dimension
			mess = rep("",length(adim))
			mess[ipos] = indx
			mess = paste0("arr[",paste0(mess,collapse=","),"]")
			mess = paste0("out = array(", mess, ", dim=ndim, dimnames=nnam)")
#browser();return()
			eval(parse(text=mess))
			return(out)
		} else {
			return(arr)
		}
	}  ## end function 'subarray'
	## ============================end subfunctions

	## Main function ==============================
	unpackList(compo, scope="L")

	## Get single-area information
	B.runs      = as.numeric(sapply(strsplit(rownames(avgPA),"\\."),function(x){x[1]}))
	B.run.rwt   = sapply(strsplit(rownames(avgPA),"\\."),function(x){paste0(x[1],".",x[2],".",x[3])}) ## vector of Runs and Rwts and ver
	B.run.ord   = unique(B.runs)
	B.run.nmc   = table(B.runs)[as.character(B.run.ord)]
	B.run.num   = rep(1:length(B.run.nmc), B.run.nmc)
	use.run.rwt = unique(B.run.rwt)

	## Take code from 'tabSS.compo.r'
	if (RPbase=="BMSY") {
		B.mcmc = data.frame (
			B0         = avgRP[,"B0"],
			Bcurr      = avgRP[,"Bcurr"],
			Bcurr.B0   = avgRP[,"Bcurr"] / avgRP[,"B0"],
			ucurr      = avgRP[,"ucurr"],
			umax       = apply(avgTS[,,"ut"],1,max,na.rm=T), ## for each mcmc sample across the time series
			MSY        = avgRP[,"MSY"],
			Bmsy       = avgRP[,"Bmsy"],
			LRP        = avgRP[,"LRP"],
			USR        = avgRP[,"USR"],
			Bcurr.Bmsy = avgRP[,"Bcurr"] / avgRP[,"Bmsy"],
			Bmsy.B0    = avgRP[,"Bmsy"] / avgRP[,"B0"],
			umsy       = avgRP[,"umsy"],
			ucurr.umsy = avgRP[,"ucurr"] / avgRP[,"umsy"]
		)
	}
	if (RPbase=="B0") {
		## eventually grab code from SweaveMCMC.Rnw
		B.mcmc = data.frame (
			B0         = avgRP[,"B0"],
			Bcurr      = avgRP[,"Bcurr"],
			Bcurr.B0   = avgRP[,"Bcurr"] / avgRP[,"B0"],
			ucurr      = avgRP[,"ucurr"],
			umax       = apply(avgTS[,,"ut"],1,max,na.rm=T), ## for each mcmc sample across the time series
			Ytrp       = avgRP[,"Ytgt"],
			#Btrp       = 0.40 * avgRP[,"B0"],  ## this should be OK when B0=SSB_1933 (VIRG) but not OK when B0=SSB_1935 (first year)
			Btrp       = avgRP[,"Btgt"],
			utrp       = avgRP[,"utgt"],
			#LRP        = 0.16 * avgRP[,"B0"],
			#USR        = 0.32 * avgRP[,"B0"],
			LRP        = 0.4 * avgRP[,"Btgt"],
			USR        = 0.8 * avgRP[,"Btgt"],
			Bcurr.Btrp = avgRP[,"Bcurr"] / avgRP[,"Btgt"],
			ucurr.utrp = avgRP[,"ucurr"] / avgRP[,"utgt"]
		)
	}
	tabQmed = t(apply(B.mcmc,2,quantile,tcall(quants5),na.rm=T))  ## not sure yet if we need this here

	A.collect = list('BC'=B.mcmc)

	## Extract area-based information if it exists
	if (exists("xavgRP") && length(xavgRP)>1) {
		## colnames(xavgRP) = "Bmsy" "Fmsy" "MSY"  "umsy" "Btgt"  "Ftgt"  "Ytgt"  "utgt"  "LRP"  "USR"  "Vmsy" "Vtgt"  "B0"   "V0"   "pVB"
		## dimnames(xavgTS)[3] =  "B"      "BtBtgt"  "BtBmsy" "D"      "fR"     "R"      "u"      "ututgt"  "utumsy" "V"
		if (RPbase=="BMSY") {
			for (a in dimnames(xavgRP)[[3]]) {
				A.mcmc = data.frame (
					## Needs revision (see below)
					B0         = avgRP[,"B0"],
					Bcurr      = avgRP[,"Bcurr"],
					Bcurr.B0   = avgRP[,"Bcurr"] / avgRP[,"B0"],
					ucurr      = avgRP[,"ucurr"],
					umax       = apply(avgTS[,,"ut"],1,max,na.rm=T), ## for each mcmc sample across the time series
					MSY        = avgRP[,"MSY"],
					Bmsy       = avgRP[,"Bmsy"],
					LRP        = avgRP[,"LRP"],
					USR        = avgRP[,"USR"],
					Bcurr.Bmsy = avgRP[,"Bcurr"] / avgRP[,"Bmsy"],
					Bmsy.B0    = avgRP[,"Bmsy"] / avgRP[,"B0"],
					umsy       = avgRP[,"umsy"],
					ucurr.umsy = avgRP[,"ucurr"] / avgRP[,"umsy"]
				)
				tabAmed = t(apply(A.mcmc,2,quantile,tcall(quants5),na.rm=T)) 
				A.collect[[a]] = A.mcmc
			} ## end area a loop
		} ## end RPbase 'BMSY'

		if (RPbase=="B0") {
			for (a in dimnames(xavgRP)[[3]]) {
				A.mcmc = data.frame (
					B0         = xavgRP[,"B0",a],
					Bcurr      = xavgTS[,as.character(currYear),"B",a],
					Bcurr.B0   = xavgTS[,as.character(currYear),"B",a] / xavgRP[,"B0",a],
					ucurr      = xavgTS[,as.character(currYear),"u",a],
					umax       = apply(xavgTS[,,"u",a],1,max,na.rm=T), ## for each mcmc sample across the time series
					Ytrp       = xavgRP[,"Ytgt",a],
					#Btrp       = 0.40 * xavgRP[,"B0",a],  ## this should be OK when B0=SSB_1933 (VIRG) but not OK when B0=SSB_1935 (first year)
					Btrp       = xavgRP[,"Btgt",a],
					utrp       = xavgRP[,"utgt",a],
					#LRP        = 0.16 * xavgRP[,"B0",a],
					#USR        = 0.32 * xavgRP[,"B0",a],
					LRP        = 0.4 * xavgRP[,"Btgt",a],
					USR        = 0.8 * xavgRP[,"Btgt",a],
					Bcurr.Btrp = xavgTS[,as.character(currYear),"B",a] / xavgRP[,"Btgt",a],
					ucurr.utrp = xavgTS[,as.character(currYear),"u",a] / xavgRP[,"utgt",a]
				)
				tabAmed = t(apply(A.mcmc,2,quantile,tcall(quants5),na.rm=T)) 
				A.collect[[a]] = A.mcmc
			} ## end area a loop
		} ## end RPbase 'B0'
	} ## end extra area collection

	areas = names(A.collect)

	B.collect = list()  ## colect various area-specific values
	## Now start geting decision tables from collected data frames
	for (i in 1:length(A.collect)) {
		ii = areas[i]
		B.mcmc = A.collect[[ii]]
		B.collect[["Nmcmc"]][[ii]] = Nmcmc  = nrow(B.mcmc)
		is.good = apply(B.mcmc,1,function(x){all(is.finite(x))})
#browser();return()
		if (sum(is.good)!=nrow(B.mcmc))
			.flush.cat(paste0(nrow(B.mcmc)-sum(is.good)," unsuitable records dropped."), "\n")
		B.mcmc  = B.mcmc[is.good,]
		B.collect[["Nmcmc.good"]][[ii]] = Nmcmc.good = nrow(B.mcmc)
		tput(Nmcmc.good) ## to use in 'tabSS.decision.r'

		## Determine subset vector for low recruitment events
		R.mcmc = avgTS[is.good,,"Rt"]
		if (useRlow) {
			if (spp.code=="BOR")
				Rlow   = R.mcmc[,"2017"] < quantile(R.mcmc[,"2017"], qRlow)
			else { 
				R.mean = apply(R.mcmc,1,mean)  ## mean recruitment 1935-2019
				Rlow   = R.mean < quantile(R.mean, qRlow)
			}
		} else {
			Rlow = rep(TRUE,nrow(R.mcmc))
		}
		## Subset vector of chosen runs to use in the base composite
		Ruse  = grepl(paste0(paste0("^",use.run.rwt),collapse="|"),rownames(B.mcmc))
		## Combined subset vector
		Rsub  = Rlow & Ruse
		if (length(Rsub)==nrow(B.mcmc)) {
			B.mcmc = B.mcmc[Rsub,]
		} else {
			message("sumtingwong with Rsub in 'agileDT'"); browser(); return()
		}
		B.collect[["Nmcmc.dtab"]][[ii]] = Nmcmc.dtab = nrow(B.mcmc)
		unpackList(B.mcmc)  ## to get the various RPs
		tput(Nmcmc.dtab)   ## to use (maybe) in 'tabSS.decision.r'
		
#browser();return()

		if (RPbase=="BSMY") {
			targets = c("Bmsy","B0","Bcurr","umsy","ucurr")
			Blst = list(
				list(ratio=0.4,target=Bmsy),
				list(ratio=0.8,target=Bmsy),
				list(ratio=1.0,target=Bmsy),
				list(ratio=1.0,target=Bcurr),
				list(ratio=0.2,target=B0),
				list(ratio=0.4,target=B0),
				##---COSEWIC---
				list(ratio=0.5,target=B0),
				list(ratio=0.7,target=B0),
				list(ratio=0.3,target=switch(ii, 'BC'=avgTS[,,"Bt"], xavgTS[,,"B",ii])),
				list(ratio=0.5,target=switch(ii, 'BC'=avgTS[,,"Bt"], xavgTS[,,"B",ii]))
			)
			## And have to do separately for u_t < u_MSY:
			Ulst = list(
				list(ratio=1, target=umsy),
				list(ratio=1, target=ucurr)
			)
			names(Blst)=c("0.4Bmsy","0.8Bmsy","Bmsy","Bcurr","0.2B0","0.4B0","0.5B0","0.7B0","0.3Gen","0.5Gen")
			names(Ulst) = c("umsy","ucurr")
		}
		if (RPbase=="B0") {
			targets = c("Btrp","B0","Bcurr","utrp","ucurr")
			Blst = list(
				list(ratio=0.16,target=B0),
				list(ratio=0.32,target=B0),
				list(ratio=0.40,target=B0),
				list(ratio=1.0, target=Bcurr),
				##---COSEWIC---
				list(ratio=0.5,target=B0),
				list(ratio=0.7,target=B0),
				list(ratio=0.3,target=switch(ii, 'BC'=avgTS[,,"Bt"], xavgTS[,,"B",ii])),
				list(ratio=0.5,target=switch(ii, 'BC'=avgTS[,,"Bt"], xavgTS[,,"B",ii]))
			)
			## And have to do separately for u_t < u_MSY:
			Ulst = list(
				list(ratio=1, target=utrp),
				list(ratio=1, target=ucurr)
			)
			names(Blst)=c("0.16B0","0.32B0","Btrp","Bcurr","0.5B0","0.7B0","0.3Gen","0.5Gen")
			names(Ulst) = c("utrp","ucurr")
		}
		if (ii %in% area) {
			Btab  = subarray(xavgPJ[,as.character(currYear:projYear),,ii,cp,drop=FALSE], ipos=3, inam="B", kpos=c(1,2,5))
			Utab  = subarray(xavgPJ[,as.character(currYear:projYear),,ii,cp,drop=FALSE], ipos=3, inam="u", kpos=c(1,2,5))
			cvec  = apply(xavgCP[,ii,],"proj",mean)
		} else {
			Btab  = subarray(avgPJ[,as.character(currYear:projYear),,cp,drop=FALSE], ipos=3, inam="Bt", kpos=c(1,2,4))
			Utab  = subarray(avgPJ[,as.character(currYear:projYear),,cp,drop=FALSE], ipos=3, inam="ut", kpos=c(1,2,4))
			cvec  = apply(apply(avgCP,c("yr","proj"),sum),"proj",mean)
		}
		if (length(cp)==1 && all(cp=="AC.00"))
			cvec  = cvec[1]
		else
			cvec  = cvec[-1]
		big.mark = ifelse(is.null(options()$big.mark), ",", options()$big.mark)
		B.collect[["cvec"]][[ii]] = cvec
		B.collect[["Cvec"]][[ii]] = Cvec  = formatC(cvec, digits=0, format="f", big.mark=big.mark)

		## Populate BRP -- P(B > RPs) and URP -- P(u < RPs) in same j loop (CPs)
		for (j in 1:length(cp)) {
			jj = cp[j]
			ijBRP = lapply(Blst,function(x,aa,pp) {  ## don't need argument 'aa' any more
				#targ = if (length(dim(x$target))==2) x$target[,aa] else x$target[,,aa]
				targ = x$target
				out = findTarget(Btab[,,pp], ratio=x$ratio, target=targ, retVal="p.hi", op=">", yrG=Ngen*gen1)
#.flush.cat(out , "\n")
				return(out)
				}, aa=ii, pp=jj)
			ijURP = lapply(Ulst,function(x,aa,pp){
				#targ = if (length(dim(x$target))==2) x$target[,aa] else x$target[,,aa]
				targ = x$target
				findTarget(Utab[,,pp], ratio=x$ratio, target=targ, retVal="p.hi", op="<", yrG=Ngen*gen1)
				}, aa=ii, pp=jj)
#browser();return()
			if (i==1 && j==1) {
				BRP = array(NA, dim=c(length(cp), length(names(ijBRP[[1]])), length(A.collect), length(ijBRP)), dimnames=list(cp=cp, year=names(ijBRP[[1]]), area=areas, refpt=names(ijBRP) ) )
				URP = array(NA, dim=c(length(cp), length(names(ijURP[[1]])), length(A.collect), length(ijURP)), dimnames=list(cp=cp, year=names(ijURP[[1]]), area=areas, refpt=names(ijURP) ) )
			}
			sapply(names(ijBRP), function(xx){ BRP[jj,,ii,xx] <<- ijBRP[[xx]] })
			sapply(names(ijURP), function(xx){ URP[jj,,ii,xx] <<- ijURP[[xx]] })
		} ## end j loop (catch policy)
	}  ## end i loop (area)
#browser();return()

	lsum = 0; lout = character(); agile=list()
	RPT = list(BRP=BRP, URP=URP)  ## list of reference point (decision) tables
	lapply (RPT, function (xtab) {
		for (l in 1:length(dimnames(xtab)$refpt)) {
			ll   = dimnames(xtab)$refpt[l]
			prb  = subarray(xtab, ipos=4, inam=ll, kpos=1:3)
			dns  = dimnames(prb); dns[[2]] = c("CC(t/y)",dns[[2]])
			prb2 = array(NA, dim=sapply(dns,length), dimnames=dns)
			prb2[,-1,] = prb
			prb2[,1,]  = as.matrix(do.call("cbind", lapply(B.collect$cvec, data.frame, stringsAsFactors=FALSE)))  ## insert catches for each policy by area
			lll = sub("Gen", paste0("B_{0,\\\\text{G}=",Ngen,"}"), 
				sub("curr", paste0("_{",currYear,"}"), 
				sub("B0", "B_0", 
				sub("trp","_\\\\text{TRP}", 
				sub("msy","_\\\\text{MSY}", 
				sub("USR", "0.8Bmsy", 
				sub("LRP", "0.4Bmsy",
			ll )))))))
			lll = paste0("$", lll, "$")
			lltab = paste0("tabDT.",ll)
			tabcap = paste0("Probability ", ifelse(any(grepl("^B", dimnames(xtab)$refpt)), "$B_t$ > ", "$u_t$ < "), lll)
			tex = texArray(prb2, sigdig=sigdig, use.round=T, zero="0", use.row.names=F, name.row.names="", add.header.column=T, new.header=areas, outnam=ll, table.caption=tabcap)
			agile[[lltab]] <<- tex
#browser();return()
			lout <<- c(lout, paste0(ll, c(".tex","+.tex")))
			## Aside: Geometry can be changed
			#texfile = readLines(paste0(ll, "+.tex"))
			#gline = grep("geometry",texfile)
			#texfile[gline] = sub("\\]\\{geometry\\}", ", paperheight=3in]{geometry}", texfile[gline])
			#writeLines(texfile, con="test.tex")
			#cmd1 = paste0("pdflatex.exe -quiet test.tex")
			#cmd2 = paste0("latexmk.exe  -c test.tex")
			if (onepage) {
				## Collect tex code to display decision tables in one pdf
				if (l==1) {
					onefile = tex$texfile
					onefile = onefile[-length(onefile)]  ## get rid of \\end{document}
				} else {
					onefile = c(onefile, ifelse(l%%tables.per.page==1, "\\clearpage", "\\bigskip"), tex$tabfile)
					#onefile = c(onefile, ifelse(l%%3==1, "\\clearpage", "\\bigskip"), tex$tabfile)
				}
			} else {
				cmd1 = paste0("pdflatex.exe -quiet ", ll, "+.tex")
				cmd2 = paste0("latexmk.exe  -c ", ll, "+.tex")
				if (pdf) {
					system(command=cmd1, intern=TRUE, wait=TRUE)
					system(command=cmd2, intern=TRUE, wait=TRUE)
				}
			}  ## end if onepage
			lsum =lsum + 1
		}  ## end l loop
	}) ## end lapply through RPT

#browser();return()

	tput(agile) ## save tex table results for use by 'tabSS.decision.r'
	save("agile", "A.collect", "B.collect", file=paste0(outnam, ".rda"))

	lout.remove = paste0(lout, collapse=" ")
	system(command=paste0("rm -f ", lout.remove), intern = TRUE)
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~agileDT


## gatherMCMC---------------------------2025-11-06
##  Gather MCMCs from base component runs or sensitivity runs
## ---------------------------------------------RH
gatherMCMC <- function( mcdir=".", type="compo", strSpp="SGR",
   basedir="C:/Users/haighr/Files/GFish/PSARC25/SGR/Data/SS3/SGR2025",
   vyrs=1933, ryrs=1935:2026, pyrs=2026:2036, valTS=c("SSB","F","Recr","RecrDev"), 
   valRP=c("SSB_MSY","annF_MSY","Dead_Catch_MSY"), debug=FALSE)
{
	## -----------------------------------------------
	## All values occur at start of year so B2022 is B at end of 2021 or start of 2022;
	##   u2022 must be u in 2021 after final catch was specified for 2021 in the data file.
	## Ian Taylor (Aug 9, 2021):
	##   Everything is start of year throughout.
	##   Setting 2021 as the final year of the model allows 2021 catch.
	##   The values for derived quantities like SSB_2021 and Bratio_2021 represent the beginning of the year values at the start of 2021 in the main time series.
	##   The SSB_2022 is one year later but will only depend on catches up through 2021, so changes in 2022 catch won't impact the 2022 quantities.
	## Chantel Wetzel (Jun 27, 2022)
	##   The F value in 2023 would be calculated from the total catch at the end of the model year.
	##   If the final year in your model with catches input into the data file is 2022, then the F_2023 would be based off the forecast catch from that year.
	## Chantel Wetzel (Jul 12, 2023)
	##   RH: If I specify year end=2023 in 'data.ss' and specify a 2023 catch in 'data.ss', while in 'forecast.ss' I specify a 2024 catch,
	##       are these the same catches? (because end-year 2023 catch = model start-year 2024 catch). 
	##   CW: Those catches would be removed from different years. 
	##       Any catches input in the forecast file will be removed from the specified year as long as it falls within the forecast period (year, season, fleet, catch (or F).
	##       Based on your description, the first forecast year would be 2024 (final data year + 1).
	## Rowan Haigh (Jul 13, 2023)
	##   ryrs = reconstructed years, needs the first projected year (endyr+1) to be the 'current' year (e.g., 2024 for POP 2023)
	##   pyrs = projected years, uses 2nd projected year onwards (e.g., 2025:2034 for POP 2023)
	## RH (230719) function 'pad0' automatically pads according to largest number in a vector (since 2011-12-12)
	## RH (230720) Need to run 'load_extra_mcmc.r' before collecting extra MCMC stuff (perhaps obvious but...)
	## RH (251027) ryrs and pyrs can have overlapping currYr: e.g., ryrs=1935:2026, pyrs=2026:2036 (SGR 2025)
	##  fixes single-area CP and PJ years; multi-area CP and PJ years already correct (calculated by 'prepCP()' and run in Linux)
	## -----------------------------------------------
	ayrs  = .su(c(vyrs, ryrs, pyrs)); nyrs = length(ayrs)
	currYr = rev(ryrs)[1]; byr = cyr = as.character(currYr)
	virgYr = vyrs        ; vyr = as.character(virgYr)
	## Add in virgin year, mostly for B0 (RH 251105)
	vryrs = c(vyrs, ryrs)

	mclst = Nmcmc = rowN = valPAs = valLLs = list()
	for (m in 1:length(mcdir)){
		mdir = file.path( sub("/$","",basedir), sub("/$","",mcdir[m]) )
		mm   = substring(basename(mdir),6) #,10)
		if(!file.exists(mdir)) {
			.flush.cat(paste0("   Skipping MCMC: run ", mm), "\n")
			next
		} else {
			.flush.cat(paste0("Getting MCMC: Run ", mm), "\n")
		}
		#if (penso) {
		#	## Special versions for PJS
		#	mm = substring(mdir,regexpr("Run",mdir)+3,regexpr("Run",mdir)+7)
		#}
		mclst[[mm]] = list()
		cp = "AC.00"

		## Gather MPD information
		d.mpd = sub("[a-z]$","",sub("\\.mh.+","",sub("\\.hm.+","",sub("\\.nuts.+","",sub("MCMC","MPD",mdir)))))
#browser();return()
		mpd = SS_output(dir=d.mpd, verbose=F, printstats=F)  ## 'Joining...' seems to occur with D-M MPDs
		Fmethod =  mpd$F_method

		## Grab likelihood components
		run = as.numeric(strsplit(mm,"\\.")[[1]][1]); names(run)="Run"
		LL.fleet  = mpd$likelihoods_by_fleet
		LL.fleet.ex = LL.fleet[is.element(LL.fleet$Label,"Surv_like"),grep("OTHER",mpd$FleetNames,invert=T,value=T)]
		names(LL.fleet.ex) = sapply(strsplit(names(LL.fleet.ex),"_"), function(x) {
			if (length(x)==3 && x[2]=="FISHERY") x[2]=x[3]
			#out = paste0(x[1],"_",substring(x[2],1,3))
			out = paste0(x[1],"_",x[2])
			return(out) })
		#names(LL.fleet.ex) = sub("TRAWL","CPUE",names(LL.fleet.ex))
		#names(LL.fleet.ex)[grep("TRAWL",names(LL.fleet.ex))]="CPUE_BT"
		LL.used   = mpd$likelihoods_used
		LL.used.ex  = LL.used[c("Survey","Age_comp","Recruitment","TOTAL"),"values"]
		names(LL.used.ex) = c("Index","AF","Recruit","Total")
		LL.compo = c(run, LL.fleet.ex, LL.used.ex)
		LL.compo = unlist(LL.compo)
#browser();return()
		mclst[[mm]][["LL"]] = LL.compo
		valLLs[[mm]] = names(LL.compo)  ## collect all names

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
		if (type %in% c("senso","penso")) CP = "AC"
		for (n in length(CP)){
			nn = CP[n]
			if (nn=="AC") next
			ndir = file.path(mdir,nn)
			#ncps = dir (ndir)[1]  ## name of catch policy subdirectories (debugging)
			ncps = dir (ndir)  ## name of catch policy subdirectories
			cp = c(cp, paste0(nn,".",sub("^[[:alpha:]]+","",ncps)))
			mdir = c(mdir, file.path(ndir,ncps))
		} ## end n loop for CP
		for (n in 1:length(mdir)) {
			mmm   = mdir[n]
			ccc   = cp[n]
			.flush.cat(paste0("...extracting catch policy: ", ccc), "\n")
			if (file.exists(paste0(mmm,"/posteriors.sso"))) {
				mmc   = SSgetMCMC(mmm, verbose=F)
			} else if (file.exists(paste0(mmm,"/sso/posteriors.sso"))) {
				mmc   = SSgetMCMC(paste0(mmm,"/sso"), verbose=F)
			} else {
				stop ("Cannot find MCMC sso files")
			}
			colnames(mmc) = convPN(colnames(mmc))

			## Sometimes need a subset of mmc
			if (file.exists(paste0(mmm,"/mcmc.ts.sub.rda"))) {
				load(paste0(mmm,"/mcmc.ts.sub.rda"))
				smcmc = attributes(mcmc.ts.sub)$samples  ## available samples from original reps
				if (!is.null(smcmc) && length(smcmc) != length(mmc$Iter)) {
					Amcmc = mmc$Iter ## all mcmc iterations
					amcmc = 1:length(Amcmc)
					names(amcmc) = pad0(Amcmc, floor(log10(max(Amcmc))) + 1)
					zmcmc = amcmc[smcmc]
					mmc = mmc[zmcmc,]
				}
			}
#browser();return()
#if (n==2) {browser();return()}

			msid  = mmc$Iter
			## ii    = as.character(msid)
			## pad0 automatically pads according to largest number in a vector
			## but if one subset is from early samples (low numbers), the code will break
			npad0 = floor(log10(max(msid))) + 1
			nmc   = length(msid)
			#run  = as.numeric(strsplit(mm,"\\.")[[1]][1])
			run   = mm
			ii = xii = paste0(run,".", pad0(msid,npad0))   ## xii = ii for extra mcmc
#browser(); return()
#print(length(ii))
			if (m==1){
				Nmcmc[[n]] = nmc
				rowN[[n]]  = ii
			} else {
				Nmcmc[[n]] = Nmcmc[[n]] + nmc
				rowN[[n]]  = c(rowN[[n]], ii)
			}

			if (n==1) { ## only need to look at the base catch policy for PA and RP
			## Parameters matrix
				matPA = array(NA, dim=c(nmc, length(valPA)), dimnames=list(sid=ii, val=valPA))
				for (k in 1:length(valPA)) {
					kk = valPA[k]
					vecPA = mmc[[kk]]
					if (is.null(vecPA)) next
					matPA[ii,kk] = vecPA
				}
				mclst[[mm]][["PA"]] = matPA

				## Reference points matrix
				matRP = array(NA, dim=c(nmc, length(valRP)), dimnames=list(sid=ii, val=valRP))
				for (k in 1:length(valRP)) {
					kk = valRP[k]
					vecRP = mmc[[kk]]
					if (is.null(vecRP)) next  ## just live with columns NA if specified valRPs are not in the data
					matRP[ii,kk] = vecRP
				}
				mclst[[mm]][["RP"]] = matRP

				## Time series matrix of population reconstruction
				matTS = array(NA, dim=c(nmc, length(vryrs), length(valTS)), dimnames=list(sid=ii, yr=vryrs, val=valTS))
				for ( k in 1:length(valTS) ) {
					kk = valTS[k]
					for (j in 1:length(vryrs)) {
						jj = as.character(vryrs[j])
						if (j==1)
							jj = "Virgin"
						vecTS = mmc[[paste0(kk,"_",jj)]]
						if (is.null(vecTS)) next
						jj = as.character(vryrs[j])  ## revert to year
						matTS[ii,jj,kk] = vecTS
#if(j==2025) {browser();return()}
					}
				}
				mclst[[mm]][["TS"]] = matTS
			}
#browser();return()

			## Projection series matrix (average catch and whatever catch policies are available)
			matPJ = array(NA, dim=c(nmc, length(pyrs), length(valTS)), dimnames=list(sid=ii, yr=pyrs, val=valTS))
			for ( k in 1:length(valTS) ) {
				kk = valTS[k]
				for (j in pyrs) {
					jj = as.character(j)
					vecPJ = mmc[[paste0(kk,"_",jj)]]
					if (is.null(vecPJ)) next
					matPJ[ii,jj,kk] = vecPJ
				}
			}
			mclst[[mm]][["PJ"]][[ccc]] = matPJ

			## Catch policy matrix
			matCP = array(NA, dim=c(nmc, length(pyrs)), dimnames=list(sid=ii, yr=pyrs))
			for (j in pyrs) {
				jj = as.character(j)
				vecCP = mmc[[paste0("ForeCatch_",jj)]]
				if (is.null(vecCP)) next
				matCP[ii,jj] = vecCP
			}
			mclst[[mm]][["CP"]][[ccc]] = apply(round(matCP),2,function(x){xx=as.numeric(names(rev(sort(table(x))))[1]); xx[is.na(xx) | !is.finite(xx)]=NA; return(xx)})
			## Note; 0-catch policy is sometimes set to a very small amount; sometimes other catch policies have small amounts added. WTF?
			## PJS thinks it may be a difference between v.2.30.16 and 3.30.17
			#mclst[[mm]][["CP"]][[ccc]] = apply(matCP,2,function(x){xx=mean(x,na.rm=T); xx[is.na(xx) | !is.finite(xx)]=NA; return(xx)})
			## CAUTION: Mode may not be sufficiently robust at higher CPs
#browser();return()

			if (file.exists(paste0(mmm,"/mcmc.posts.rda"))) {
				## Binary rda files were generated by 'load_extra_mcmc.r'
				load(paste0(mmm,"/mcmc.posts.rda"))
				load(paste0(mmm,"/mcmc.ts.sub.rda"))  ## for the catch-by-area table (cattab)

				## Collect extra reference points
				rpdata = c("MSY.mcmc","B0.mcmc","V0.mcmc","pVB.mcmc")
				rpnames = c(dimnames(MSY.mcmc[[1]])[[2]], gsub("\\.mcmc$", "", rpdata)[-1])  ## these may contain old Btgt names (e.g. B40)
				anames = names(MSY.mcmc)
#browser();return()
				#snumbs = as.numeric(sub("^s", "", dimnames(MSY.mcmc[[1]])[[1]]))     ## even though all except MSY have the right names
				#snames = paste0(run,".",pad0(snumbs, floor(log10(max(snumbs))) + 1)) ## standardise to mesh with xavgRP
				snames = xii ## run and original iteration number
				dnames = list(mcmc=snames, val=rpnames, area=anames)
				xmatRP = array(NA, dim=c(length(snames), length(rpnames), length(anames)), dimnames=dnames)
				for (k in 1:length(rpdata)) {
					kk  = rpdata[k]
					kdat = get(kk)
					for (a in anames) {
						if (kk %in% c("MSY.mcmc")) {
							adat = kdat[[a]]
							ii = colnames(adat)
						} else {
							adat = kdat[,a,drop=FALSE]
							ii = gsub("\\.mcmc$", "", kk)
						}
						xmatRP[,ii,a] = as.matrix(adat)
					}
				}
				dimnames(xmatRP)$val = sub("40$","tgt",dimnames(xmatRP)$val)  ## convert B40 to Btgt, etc. (if needed)
				mclst[[mm]][["xRP"]][[ccc]] = xmatRP

				## Collect extra time series (note: contains reconstructed and projected years)
				## Add in virgin year (RH 251105)
				tsdata  = setdiff(ls(pattern="\\.mcmc$"), rpdata)
				tsnames = gsub("\\.mcmc$", "", tsdata)
				#snames  = dimnames(get(tsdata[1])[[1]])[[1]]
				#snumbs = as.numeric(sub("^s", "", snames))
				#snames = paste0(run,".",pad0(snumbs, floor(log10(max(snumbs))) + 1))  ## standardise to mesh with xavgTS
				snames = xii ## run and original iteration number
				ynames  = dimnames(get(tsdata[1])[[1]])[[2]]
				dnames = list(mcmc=snames, year=ynames, val=tsnames, area=anames)
				xmatTS = array(NA, dim=c(length(snames),length(ynames),length(tsnames),length(anames)), dimnames=dnames)
				for (k in 1:length(tsnames)) {
					kk   = tsnames[k]
					kkk  = tsdata[k]
					kdat = get(kkk)
					for (a in anames) {
						adat = kdat[[a]]
#if (n==2) {browser();return()}
						jj = intersect(colnames(adat), colnames(xmatTS))
						xmatTS[,jj,kk,a] = as.matrix(adat[,jj])
					}
				}
				dimnames(xmatTS)$val = sub("40$","tgt",dimnames(xmatTS)$val)  ## convert B40 to Btgt, etc. (if needed)
				mclst[[mm]][["xTS"]][[ccc]] = xmatTS

				## Collect extra time series projections (just a subset of xmatTS)
				cpyrs = .su(c(ryrs[length(ryrs)],pyrs))  ## include the current year for convenience later  ## 251027 now already included
				xmatPJ = xmatTS[,as.character(cpyrs),,,drop=F]
				mclst[[mm]][["xPJ"]][[ccc]] = xmatPJ

				## Collect area-based catch policies
				xmatCP = cattab[as.character(cpyrs),]
				## To name the dimnames, must convert data.frame to matrix (je ne sais pas pourquoi)
				xmatCP = as.matrix(xmatCP)
				names(dimnames(xmatCP)) = c("year","area")
				mclst[[mm]][["xCP"]][[ccc]] = xmatCP
#if (n==1) {browser();return()}
			}  ## end collecting extras
		}  ## end n loop for mdir (catch policies)
		mclst[[mm]][["MPD"]] = P.mpd  ## for each base run
	}  ## end m loop for mcdir (base component runs)
	## ----- END initial data collection -----

	## Grab parameters across all runs
	if (strSpp %in% c("CAR","POP","YTR","SGR")) {
		valPA = unique(unlist(valPAs))
		valPA = c(grep("R0",valPA,value=T), grep("Rdist",valPA,value=T), grep("R0|Rdist|theta",valPA,invert=T,value=T), grep("theta",valPA,value=T))
		valLL = unique(unlist(valLLs))
		chunk = "Index|AF|Recruit|Total"
		valLL = c(grep(chunk,valLL,invert=T,value=T), grep(chunk,valLL,value=T))
	} else {
		.flush.cat("Include species '", strSpp, " in parameter grab", sep="", "\n")
		browser();return()
	}

	mpdPA = list()
	dim1  = Nmcmc[[1]] ## number of rows for base case (cumulative if more than one base run component)
	dnam1 = rowN[[1]]  ## base case row names
	iv    = max(unlist(Nmcmc))
	ivnam = .su(unlist(rowN))
#browser();return()

	avgPA = array(NA, dim=c(dim1, length(valPA)), dimnames=list(mcmc=dnam1, par=valPA))
	#vRP   = c("Bcurr","B0","20B0","40B0","Fcurr","ucurr","MSY","Bmsy","LRP","USR","Fmsy","umsy")
	## B40, Y40, etc are values at the TRP=0.4B0  -- assumes Btgt in control file is 40% of B0 (may not always be the case)
	## For now, use shorthand like 'B40' to mean 'Btgt' even if 'Btgt'=0.xxB0; probably better to call it 'Btgt'
	## Shorthand labels like '32B40' are problematic because it's actually '32B0' or '80B40' (80% of 0.4B0)
	vRP   = c("Bcurr","B0","20B0","40B0","MSY","Bmsy","40Bmsy","80Bmsy","Ytgt","Btgt","80Btgt","40Btgt","ucurr","Fcurr","umsy","Fmsy","utgt","Ftgt")
	avgRP = array(NA, dim=c(dim1, length(vRP)), dimnames=list(mcmc=dnam1, val=vRP))
	#vTS   = c("Bt","BtB0","BtBmsy","Ft","FtFmsy","ut","utumsy","Rt","Rtdev")
	vTS   = c("Bt","BtB0","BtBmsy","BtBtrp","Rt","Rtdev","Ft","ut","FtFmsy","utumsy","FtFtrp","ututrp")
	avgTS = array(NA, dim=c(dim1, length(vryrs), length(vTS)), dimnames=list(mcmc=dnam1, yr=vryrs, val=vTS))
	avgLL =array(NA, dim=c(length(valLL), length(mclst)), dimnames=list(ll=valLL, run=names(mclst)))
	## Average projections have to be dimensioned by the projection with the most MCMCs
	## (RH 251027) multi-area projections include current year (e.g., SGR 2026) but single-area projections don't; pyrs can now overla ryrs
	avgPJ = array(NA, dim=c(iv, length(pyrs), length(vTS), length(cp)), dimnames=list(mcmc=ivnam, yr=pyrs, val=vTS, proj=cp))
	avgCP = array(NA, dim=c(length(pyrs), length(mclst), length(cp)), dimnames=list(yr=pyrs, run=names(mclst), proj=cp))
	#avgPJ = array(NA, dim=c(iv, length(pyrs)+1, length(vTS), length(cp)), dimnames=list(mcmc=ivnam, yr=c(currYr,pyrs), val=vTS, proj=cp))
	#avgCP = array(NA, dim=c(length(pyrs)+1, length(mclst), length(cp)), dimnames=list(yr=c(currYr,pyrs), run=names(mclst), proj=cp))

	## Check for extra mcmc stuff (why is the last mm being used here? need to use the first, assuming it's the base case)
	aref = names(mclst)[1]
	if ("xRP" %in% names(mclst[[aref]])) {
		xavgRP = array(NA, dim=c(dim1, dim(mclst[[aref]][["xRP"]][[cp[1]]])[-1]), dimnames=c(list(mcmc=dnam1), dimnames(mclst[[aref]][["xRP"]][[cp[1]]])[-1]) )
	}
	if ("xTS" %in% names(mclst[[aref]])) { ## includes projections
		xavgTS = array(NA, dim=c(dim1, dim(mclst[[aref]][["xTS"]][[cp[1]]])[-1]), dimnames=c(list(mcmc=dnam1), dimnames(mclst[[aref]][["xTS"]][[cp[1]]])[-1]) )
	}
	## (RH 251027) multi-area projections include current year (e.g., SGR 2026) but single-area projections don't; will attempt to fix above in avgPJ
	if ("xPJ" %in% names(mclst[[aref]])) { ## includes projections
		xavgPJ = array(NA, dim=c(dim1, dim(mclst[[aref]][["xPJ"]][[cp[1]]])[-1], length(cp)), dimnames=c(list(mcmc=dnam1), dimnames(mclst[[aref]][["xPJ"]][[cp[1]]])[-1], list(proj=cp) ))
	}
	if ("xCP" %in% names(mclst[[aref]])) { ## includes projections
		xavgCP = array(NA, dim=c(dim(mclst[[aref]][["xCP"]][[cp[1]]]), length(cp)), dimnames=c(dimnames(mclst[[aref]][["xCP"]][[cp[1]]]), list(proj=cp) ))
	}

	## Build the composite run ('model average')
	for (m in 1:length(mclst)){
		mm    = names(mclst)[m]
.flush.cat(mm, "\n")
		if (debug) matrix(mclst[[mm]][["MPD"]][valPA],nrow=1)
		mpdPA = rbind(mpdPA, matrix(mclst[[mm]][["MPD"]][valPA], nrow=1))
		mcPA  = mclst[[mm]][["PA"]]
		mcRP  = mclst[[mm]][["RP"]]
		mcTS  = mclst[[mm]][["TS"]]
		mcPJ  = mclst[[mm]][["PJ"]]
		mcCP  = mclst[[mm]][["CP"]]
		#run   = as.numeric(strsplit(mm,"\\.")[[1]][1])
		run   = mm
		#npad0 = floor(log10(max(as.numeric(rownames(mcTS))))) +1
		iii = rownames(mcTS)
		#iii   = paste0(run,".",pad0(as.numeric(rownames(mcTS)), npad0)) ## but padding automatically detects largest number
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
		## Populate reference point array
		use.Btgt = FALSE
		avgRP[iii,"Bcurr"]  = mcTS[,cyr,"SSB"]
		avgRP[iii,"B0"]     = mcTS[,1,"SSB"]          ## first year should be 'vyr' (virgin year 1933)
#if (m==2) {browser();return()}
		avgRP[iii,"20B0"]   = 0.2 * mcTS[,1,"SSB"]
		avgRP[iii,"40B0"]   = 0.4 * mcTS[,1,"SSB"]
		avgRP[iii,"MSY"]    = mcRP[,"Dead_Catch_MSY"]
		avgRP[iii,"Bmsy"]   = mcRP[,"SSB_MSY"]
		avgRP[iii,"40Bmsy"] = 0.4 * mcRP[,"SSB_MSY"]
		avgRP[iii,"80Bmsy"] = 0.8 * mcRP[,"SSB_MSY"]
		if ("SSB_Btgt" %in% colnames(mcRP) && !all(is.na(mcRP[,"SSB_Btgt"]))) {
			use.Btgt = TRUE
			avgRP[iii,"Ytgt"] = mcRP[,"Dead_Catch_Btgt"]
			avgRP[iii,"Btgt"] = mcRP[,"SSB_Btgt"]
			avgRP[iii,"80Btgt"] = 0.8 * mcRP[,"SSB_Btgt"]
			avgRP[iii,"40Btgt"] = 0.4 * mcRP[,"SSB_Btgt"]
		}
#browser();return()
		if (Fmethod==1) {
			## Pope's approximation (discrete) ## Pope's method: check below (not tested)
			avgRP[iii,"ucurr"] = mcTS[,byr,"F"]
			avgRP[iii,"Fcurr"] = -log(1-mcTS[,byr,"F"])
			avgRP[iii,"umsy"]  = mcRP[,"annF_MSY"]
			avgRP[iii,"Fmsy"]  = -log(1-mcRP[,"annF_MSY"])
			if (use.Btgt) {
				avgRP[iii,"utgt"]  = mcRP[,"annF_Btgt"]
				avgRP[iii,"Ftgt"]  = -log(1-mcRP[,"annF_Btgt"])
			}
		} else {
#browser();return()
			## Baranov or Hybrid (continuous)
			avgRP[iii,"Fcurr"] = mcTS[,byr,"F"]
			avgRP[iii,"ucurr"] = 1 - exp(-mcTS[,byr,"F"])
			avgRP[iii,"Fmsy"]  = mcRP[,"annF_MSY"]
			avgRP[iii,"umsy"]  = 1-exp(-mcRP[,"annF_MSY"])
			if (use.Btgt) {
				avgRP[iii,"Ftgt"]  = mcRP[,"annF_Btgt"]
				avgRP[iii,"utgt"]  = 1-exp(-mcRP[,"annF_Btgt"])
			}
		}

		## Populate time series array
		avgTS[iii,dimnames(mcTS)$yr,"Bt"]     = mcTS[,,"SSB"]
		avgTS[iii,dimnames(mcTS)$yr,"BtB0"]   = mcTS[,,"SSB"] / avgRP[iii,"B0"]
		avgTS[iii,dimnames(mcTS)$yr,"BtBmsy"] = mcTS[,,"SSB"] / avgRP[iii,"Bmsy"]
		if (use.Btgt) {
			avgTS[iii,dimnames(mcTS)$yr,"BtBtrp"] = mcTS[,,"SSB"] / avgRP[iii,"Btgt"]
		}
		avgTS[iii,dimnames(mcTS)$yr,"Rt"]     = mcTS[,,"Recr"]
		avgTS[iii,dimnames(mcTS)$yr,"Rtdev"]  = mcTS[,,"RecrDev"]
		if (Fmethod==1) {
			## Pope's approximation (discrete) ## Pope's method: check below (not tested)
			avgTS[iii,dimnames(mcTS)$yr,"ut"]     = mcTS[,,"F"]
			avgTS[iii,dimnames(mcTS)$yr,"utumsy"] = mcTS[,,"F"] / avgRP[iii,"umsy"]
			avgTS[iii,dimnames(mcTS)$yr,"Ft"]     = -log(1 - mcTS[,,"F"])
			avgTS[iii,dimnames(mcTS)$yr,"FtFmsy"] = -log(1 - mcTS[,,"F"]) / avgRP[iii,"Fmsy"]
			if (use.Btgt) {
				avgTS[iii,dimnames(mcTS)$yr,"FtFtrp"] = -log(1 - mcTS[,,"F"]) / avgRP[iii,"Ftgt"]
				avgTS[iii,dimnames(mcTS)$yr,"ututrp"] = mcTS[,,"F"] / avgRP[iii,"utgt"]
			}
		} else {
			## Baranov or Hybrid (continuous)
			avgTS[iii,dimnames(mcTS)$yr,"Ft"]     = mcTS[,,"F"]
			avgTS[iii,dimnames(mcTS)$yr,"FtFmsy"] = mcTS[,,"F"] / avgRP[iii,"Fmsy"]
			avgTS[iii,dimnames(mcTS)$yr,"ut"]     = 1 - exp(-mcTS[,,"F"])
			avgTS[iii,dimnames(mcTS)$yr,"utumsy"] = (1 - exp(-mcTS[,,"F"])) / avgRP[iii,"umsy"]
			if (use.Btgt) {
				avgTS[iii,dimnames(mcTS)$yr,"FtFtrp"] = mcTS[,,"F"] / avgRP[iii,"Ftgt"]
				avgTS[iii,dimnames(mcTS)$yr,"ututrp"] = (1 - exp(-mcTS[,,"F"])) / avgRP[iii,"utgt"]
			}
		}

		## Loop through catch policies for each run
		for (n in 1:length(cp)){
			ccc = cp[n]
			avgCP[as.character(pyrs),mm,ccc] =  mcCP[[ccc]]
			ii  =  xii = dimnames(mcPJ[[ccc]])$sid
#browser();return()
			#iii = dimnames(avgRP)$mcmc
			kkk = dimnames(mcPJ[[ccc]])$yr
#if (n==2) {browser();return()}
			## Populate projection array
			avgPJ[ii,kkk,"Bt",ccc]      = mcPJ[[ccc]][ii,,"SSB"]
			avgPJ[ii,kkk,"BtB0",ccc]   = mcPJ[[ccc]][ii,kkk,"SSB"] / avgRP[ii,"B0"]
			avgPJ[ii,kkk,"BtBmsy",ccc] = mcPJ[[ccc]][ii,kkk,"SSB"] / avgRP[ii,"Bmsy"]
			if (use.Btgt) {
				avgPJ[ii,kkk,"BtBtrp",ccc] = mcPJ[[ccc]][ii,kkk,"SSB"] / avgRP[ii,"Btgt"]
			}
			avgPJ[ii,kkk,"Rt",ccc]      = mcPJ[[ccc]][ii,kkk,"Recr"]
			avgPJ[ii,kkk,"Rtdev",ccc]   = mcPJ[[ccc]][ii,kkk,"RecrDev"]
			if (Fmethod==1) {
				## Pope's approximation (discrete) ## Pope's method: check below (not tested)
				avgPJ[ii,kkk,"ut",ccc]      = mcPJ[[ccc]][ii,kkk,"F"]
				avgPJ[ii,kkk,"utumsy",ccc] = mcPJ[[ccc]][ii,kkk,"F"] / avgRP[ii,"umsy"]
				avgPJ[ii,kkk,"Ft",ccc]      = -log(1 - mcPJ[[ccc]][ii,kkk,"F"])
				avgPJ[ii,kkk,"FtFmsy",ccc] = -log(1 - mcPJ[[ccc]][ii,kkk,"F"]) / avgRP[ii,"Fmsy"]
				if (use.Btgt) {
					avgPJ[ii,kkk,"ututrp",ccc] = mcPJ[[ccc]][ii,kkk,"F"] / avgRP[ii,"utgt"]
					avgPJ[ii,kkk,"FtFtrp",ccc] = -log(1 - mcPJ[[ccc]][ii,kkk,"F"]) / avgRP[ii,"Ftgt"]
				}
			} else {
				## Baranov or Hybrid (continuous)
				avgPJ[ii,kkk,"Ft",ccc]      = mcPJ[[ccc]][ii,kkk,"F"]
				avgPJ[ii,kkk,"FtFmsy",ccc] = mcPJ[[ccc]][ii,kkk,"F"] / avgRP[iii,"Fmsy"]
				avgPJ[ii,kkk,"ut",ccc]      = 1 - exp(-mcPJ[[ccc]][ii,kkk,"F"])
				avgPJ[ii,kkk,"utumsy",ccc] = (1 - exp(-mcPJ[[ccc]][ii,kkk,"F"])) / avgRP[ii,"umsy"]
				if (use.Btgt) {
					avgPJ[ii,kkk,"FtFtrp",ccc] = mcPJ[[ccc]][ii,kkk,"F"] / avgRP[iii,"Ftgt"]
					avgPJ[ii,kkk,"ututrp",ccc] = (1 - exp(-mcPJ[[ccc]][ii,kkk,"F"])) / avgRP[ii,"utgt"]
				}
			}
		}  ## end n loop (catch policies)

		## ----- Collect extra mcmc data if they exist -----
		## Extra reference points
		if ("xRP" %in% names(mclst[[mm]])) {
			xmcRP = mclst[[mm]][["xRP"]][[cp[1]]]           ## grab the base runs' reference points (default catch policy)
			xtemp = xmcRP[xii,,,drop=FALSE]                 ## xii = original iteration names for the extra MCMCs
			jjj   = intersect(dimnames(xavgRP)[[2]],dimnames(xtemp)[[2]])  ## (RH 251012) in SGR 2025, we added B0 ref points for base run after we ran sensitivities
#browser();return()
			atemp = dimnames(xmcRP[xii,,,drop=FALSE])$area  ## need specific areas for certain runs (e.g., single-area runs)
			xavgRP[xii,jjj,atemp] = xtemp[xii,jjj,atemp]    ## need to match derived parameter values
		} else {
			xavgRP = NA
		}
		## Extra time series by area
		if ("xTS" %in% names(mclst[[1]])) {
			xmcTS = mclst[[mm]][["xTS"]][[cp[1]]]
			xtemp = xmcTS[xii,,,,drop=FALSE]                 ## xii = original iteration names for the extra MCMCs
			atemp = dimnames(xmcTS[xii,,,,drop=FALSE])$area  ## need specific areas for certain runs (e.g., single-area runs)
			jjj   = intersect(dimnames(xavgTS)[[2]],dimnames(xtemp)[[2]])  ## (RH 251106) only base run has the virgin year 1933 (SGR 2025)
			kkk   = intersect(dimnames(xavgTS)[[3]],dimnames(xtemp)[[3]])  ## (RH 251012) in SGR 2025, we added B0 ref points for base run after we ran sensitivities
			xavgTS[xii,jjj,kkk,atemp] = xtemp[xii,jjj,kkk,atemp]   ## need to match derived parameter values
#if(m==2) {browser();return()}
		} else {
			xavgTS = NA
		}
		## Extra projections
		if ("xPJ" %in% names(mclst[[1]])) {
			## Loop through catch policies for each run and collect projections
			for (n in 1:length(cp)){
				ccc   = cp[n]
				xmcPJ = mclst[[mm]][["xPJ"]][[ccc]]
				xtemp = xmcPJ[xii,,,,drop=FALSE]                 ## xii = original iteration names for the extra MCMCs
				atemp = dimnames(xmcPJ[xii,,,,drop=FALSE])$area  ## need specific areas for certain runs (e.g., single-area runs)
				jjj   = intersect(dimnames(xavgPJ)[[2]],dimnames(xtemp)[[2]])  ## (RH 251106) only base run has the virgin year 1933 (SGR 2025)
				kkk   = intersect(dimnames(xavgPJ)[[3]],dimnames(xtemp)[[3]])  ## (RH 251012) in SGR 2025, we added B0 ref points for base run after we ran sensitivities
#if(m==2) {browser();return()}
				xavgPJ[xii,jjj,kkk,atemp,ccc] = xtemp[xii,jjj,kkk,atemp] ## need to match derived parameter values
#browser();return()
				#xavgPJ[xii,,,,ccc] = xmcPJ[xii,,,]
			}
		} else {
			xavgPJ = NA
		}
#if (m==2) {browser();return()}
		## Extra catch policies
		if ("xCP" %in% names(mclst[[1]])) {
			if (m==1) {
				## Loop through catch policies for each run and collect area-based catch policies
				for (n in 1:length(cp)){
					ccc   = cp[n]
					xmcCP = mclst[[mm]][["xCP"]][[ccc]]
					xavgCP[,,ccc] = xmcCP
				}
			}
		} else {
			xavgCP = NA
		}
	}  ## end m loop (mclst, base components)
#browser();return()

	storage.mode(mpdPA)="double"
	rownames(mpdPA) = names(mclst); colnames(mpdPA) = valPA  ## because of occasional extra parameters
	if (type=="compo")
		out = list(ampdPA=mpdPA, avgLL=avgLL, avgPA=avgPA, avgRP=avgRP, avgTS=avgTS, avgPJ=avgPJ, avgCP=avgCP, xavgRP=xavgRP, xavgTS=xavgTS, xavgPJ=xavgPJ, xavgCP=xavgCP)
	else if (type=="senso")
		out = list(smpdPA=mpdPA, senLL=avgLL, senPA=avgPA, senRP=avgRP, senTS=avgTS, senPJ=avgPJ, senCP=avgCP, xsenRP=xavgRP, xsenTS=xavgTS, xsenPJ=xavgPJ, xsenCP=xavgCP)
	else if (type=="penso")
		out = list(pmpdPA=mpdPA, penLL=avgLL, penPA=avgPA, penRP=avgRP, penTS=avgTS, penPJ=avgPJ, penCP=avgCP, xpenRP=xavgRP, xpenTS=xavgTS, xpenPJ=xavgPJ, xpenCP=xavgCP)
	else if (type=="renso")
		out = list(rmpdPA=mpdPA, renLL=avgLL, renPA=avgPA, renRP=avgRP, renTS=avgTS, renPJ=avgPJ, renCP=avgCP, xrenRP=xavgRP, xrenTS=xavgTS, xrenPJ=xavgPJ, xrenCP=xavgCP)
	else
		out = "sumtingwong"
#browser();return()
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~gatherMCMC


## getSS.rdevs-------------------------2023-10-13
##  Get reruitment deviations from replist
## ----------------------------------------r4ss|RH
getSS.rdevs <- function (replist, forecast=FALSE, minyr=-Inf, maxyr=Inf) 
{
	parameters    <- replist$parameters
	recruit       <- replist$recruit
	startyr       <- replist$startyr
	endyr         <- replist$endyr
	sigma_R_in    <- replist$sigma_R_in
	recdevEarly   <- parameters[substring(parameters$Label, 1, 13) %in% c("Early_RecrDev"), ]
	early_initage <- parameters[substring(parameters$Label, 1, 13) %in% c("Early_InitAge"), ]
	main_initage  <- parameters[substring(parameters$Label, 1, 12) %in% c("Main_InitAge"), ]
	recdev        <- parameters[substring(parameters$Label, 1, 12) %in% c("Main_RecrDev"), ]
	recdevFore    <- parameters[substring(parameters$Label, 1, 8) == "ForeRecr", ]
	recdevLate    <- parameters[substring(parameters$Label, 1, 12) == "Late_RecrDev", ]
	ttput(recruit)
#browser();return()

	recdev$Yr <- as.numeric(substring(recdev$Label, 14))
	if (nrow(recdevEarly) > 0) {
		recdevEarly$Yr <- as.numeric(substring(recdevEarly$Label, 15))
	}
	else {
		recdevEarly$Yr <- integer(0)
	}
	if (nrow(early_initage) > 0) {
		early_initage$Yr <- startyr - as.numeric(substring(early_initage$Label, 15))
		recdevEarly <- rbind(early_initage, recdevEarly)
	}
	if (nrow(main_initage) > 0) {
		main_initage$Yr <- startyr - as.numeric(substring(main_initage$Label, 14))
		recdev <- rbind(main_initage, recdev)
	}
	if (nrow(recdevFore) > 0) {
		recdevFore$Yr <- as.numeric(substring(recdevFore$Label, 10))
	} else {
		recdevFore$Yr <- NULL
	}
	if (nrow(recdevLate) > 0) {
		recdevLate$Yr <- as.numeric(substring(recdevLate$Label, 14))
		recdevFore <- rbind(recdevLate, recdevFore)
	}
	Yr <- c(recdevEarly$Yr, recdev$Yr, recdevFore$Yr)
	if (forecast) {
		goodyrs <- ifelse(Yr >= minyr & Yr <= maxyr, TRUE, FALSE)
	} else {
		goodyrs <- Yr <= endyr + 1 & Yr >= minyr & Yr <= maxyr
	}
	alldevs <- rbind(recdevEarly, recdev, recdevFore)[goodyrs,]
	recdevs =  alldevs[,c("Yr","Value","Parm_StDev")]
#browser();return()
	ttput(recdevs)
	return(invisible(recdevs))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getSS.rdevs


## load_extra_mcmc----------------------2025-11-06
##  Load extra MCMC information from Report files 
##  generated for every sample.
##  Only use for area-based models at this point.
## 'run' is only used as a tag in this function.
## Modified for SGR 2025
## ------------------------------------------CG|RH
load_extra_mcmc <- function(dir.mcmc=".", dir.extra="./sso", 
   quants5=c(0.05,0.25,0.5,0.75,0.95),  RC=c(TRUE,FALSE), loadCP=FALSE,
   virginyr=1933, startyr=1935, run="28v2c", areas=c("5ABC","5DE","3CD"), Fmethod=3, benchmark=1,
   plot=TRUE, png=FALSE, pngres=400, PIN=c(9,9), lang="e",
   show.combo=TRUE, vertical=FALSE)
{
	cwd = getwd()
	on.exit(setwd(cwd))

	## RC = boolean to use Reports and/or CompReports (comps not really needed)
	# on.exit(gc(verbose=FALSE))  ## garbage collection is a painfully slow process

	## Subfunctions--------------------------------
	## Extract a table from the report file
	makeRepTab = function(dat, pat1, pat2, off1=2, off2=2, hup=1) {
		begin.idx = grep(pat1, dat) + off1
		end.idx   = grep(pat2, dat) - off2
		head.idx  = begin.idx - hup
#browser();return()
		head.dat   = gsub("^\\s+|\\s+$","",dat[head.idx])  ## get rid of leading and trailing whitespace
		col.nams  = strsplit(head.dat, split="\\s+")[[1]]
		dat.idx    = begin.idx:end.idx
		data.dat   = gsub("^\\s+|\\s+$","",dat[dat.idx])  ## get rid of leading and trailing whitespace
		data.tmp  = strsplit(data.dat, split="\\s+")
		data.tab  = do.call("rbind", lapply(data.tmp, function(x) {
			if(is.null(dim(x))) x=t(x)
			data.frame(x, stringsAsFactors=FALSE)
			}) )
		colnames(data.tab) = col.nams
		rownames(data.tab) = dat.idx
		return(data.tab)
	}
	orderRuns = function(rnames, get0pad.only=FALSE) {
		rdir   = dirname(rnames)
		rnames = basename(rnames)
		rlens  = regexpr("\\d+", rnames)
		oldnum = as.numeric(substring(rnames, rlens, rlens + attributes(rlens)$match.length - 1))
		npad0  = floor(log10(max(oldnum))) + 1
		if (get0pad.only)
			return(npad0)
		newnum = pad0(oldnum, n=npad0)
		outnam = paste0(rdir, "/", rnames[order(newnum)])
		attr(outnam,"npad0") = npad0
		return(outnam)
	}
	## Pad individula run based on names of all runs
	padRun = function(onefile, allfiles) {
		bname  = basename(onefile)
		z      = regexpr("\\d+", bname)
		npad0  = orderRuns(allfiles, get0pad.only=TRUE)
		outnam = pad0(substring(bname, z, z + attributes(z)$match.length - 1), npad0)  ## check that order of reps makes sense
		return(outnam)
	}
	make.ts <- function(reps) {
		## Make time series object (mcmc.ts)
		## ---------------------------------
		.flush.cat("Extracting 'mcmc.ts' object from 'reps' object\n")
		tic(); mcmc.ts = lapply(reps, function(x) { makeRepTab(dat=x, pat1="^TIME_SERIES", pat2="^SPR_SERIES") }); toc()  ## ~5 min for 2000 MCMC samples
		#mcmc.exp = lapply(reps, function(x) { makeRepTab(dat=x, pat1="^EXPLOITATION report:14", pat2="^CATCH report:15", off1=13, hup=5) })
		trun = run  ## time series run
		save("mcmc.ts", "trun", file="mcmc.ts.rda")
	}
	make.ts.sub <- function(mcmc.ts) {
		## Extract time series values (mcmc.ts.sub)
		.flush.cat("Deriving 'mcmc.ts.sub' from 'mcmc.ts'\n")
		mcmc.ts.sub = list()
		cpats.ts = c("^SpawnBio$", "Bio_smry", "Recruit", "dead\\(B", "Hrate|^F", "sel\\(B")  ## Appears to be no area-specific recdev (because R is coastwide)
		nmcmc = length(mcmc.ts)
		npad0 = floor(log10(nmcmc)) + 1   ## for actual number of MCMC samples

		for (i in 1:length(cpats.ts)) {
			ival = cpats.ts[i]
			ipos = grep(ival,colnames(mcmc.ts[[1]]))
			ii   = colnames(mcmc.ts[[1]])[ipos]  ## get all columns with label pattern
			iii  = sub(":_[1-9]", "", ii[1])     ## just getting one label for a later merge
			if (iii=="F" && Fmethod>1) {         ## convert to a harvest rate
				ilst = lapply(mcmc.ts, function(x) {
					xtab = apply(x[,ii,drop=FALSE], 2, function(xchr){
						xval = as.numeric(xchr); u = 1 - exp(-xval); return(u) })  ## F is a character at this point
					rownames(xtab) = rownames(x)
					colnames(xtab) = gsub("F_apical","Hrate",colnames(xtab))
					return(as.data.frame(xtab))
				})
				iii = "Hrate_apical"
			} else {
				## everything else
				ilst = lapply(mcmc.ts, function(x) { x[,ii,drop=FALSE] })
			}
			if (all(sapply(sapply(ilst,dim),is.null))) {
				ilst = lapply(ilst,as.numeric)
			} else {
				ilst = lapply(ilst, function(x) { #if (i==4) {browser();return()};
					apply(sapply(x,as.numeric),1,sum) 
				})
			}
			cols.init = mcmc.ts[[1]][,init.cols]
			cols.init[,c(1,2,4)] = sapply(cols.init[,c(1,2,4)], as.numeric)  ## leave hard-wired for now
			idat = data.frame(cols.init, do.call("cbind",lapply(ilst, data.frame, stringsAsFactors=FALSE)))
			colnames(idat)[-(init.cols)] = paste0("s", pad0(1:nmcmc, npad0))
			zbad = is.na(idat[,-(init.cols)])
			idat[,-(init.cols)][zbad] = NA
			mcmc.ts.sub[[iii]] = idat
		} ## end i loop (cpats.ts)

		## Vulnerable Biomass & Exploitation rate
		#VB  = data.frame(cols.init, mcmc.ts.sub[['sel(B)']][,-(init.cols)]/mcmc.ts.sub[['Hrate']][,-(init.cols)])
		VB = data.frame(cols.init, mcmc.ts.sub[['Bio_smry']][,-(init.cols)])
		u  = data.frame(cols.init, mcmc.ts.sub[['sel(B)']][,-(init.cols)] / mcmc.ts.sub[['Bio_smry']][,-(init.cols)])
		if (ncol(VB)==(init.len + 1))
			colnames(VB)[init.len+1] = colnames(u)[init.len+1] = "s1"
		zbad = is.na(VB[,-(init.cols)])
		VB[,-(init.cols)][zbad] = NA
		mcmc.ts.sub[["VB"]] = VB
		zbad = is.na(u[,-(init.cols)])
		u[,-(init.cols)][zbad] = NA
		mcmc.ts.sub[["u"]] = u
		F = u
		F[,-(init.cols)] = -log(1-F[,-(init.cols)])
#browser();return()

		## Gather catch for later (sample catches are over the place!)
		get_mode <- function(v) {  ## thanks Google AI
			uniqv <- unique(v)
			uniqv[which.max(tabulate(match(v, uniqv)))]
		}
		selB    = mcmc.ts.sub[['sel(B)']]
		zcat    = is.element(selB$Era,c("TIME","FORE"))
		catmode = apply(selB[zcat,-init.cols],1,get_mode)  ## better to get the modal catch rather than the first one (which is oftem compromised)
		catarea = cbind(selB[zcat,init.cols], catmode=catmode)
		catlst  = split(catarea[,"catmode"], catarea[,"Area"])
		## Find the intended forecast catch (RH 250513)
		#zFORE = is.element(selB$Era,"FORE")
		#cFORE = selB[zFORE,-init.cols]
		#mFORE = apply(cFORE, 1, max, na.rm=TRUE)
		#aFORE = selB[,c("Area",colnames(selB)[init.len + 1])]
		#aFORE[zFORE,2] = mFORE
		##catlst = split(selB[zcat,min(nmcmc,100)], selB[zcat,"Area"])  ## need to make sure that simulation catch is not reduced WTF?
		#catlst = split(aFORE[zcat,2], aFORE[zcat,"Area"])  ## need to make sure that simulation catch is not reduced by SS3
		cattab = do.call("cbind", lapply(catlst, data.frame, stringsAsFactors=FALSE))
		colnames(cattab) = areas
		rownames(cattab) = catarea[grep(1,catarea$Area),"Yr"] #.su(selB$Yr[zcat])
		cattab$total = apply(cattab,1,sum)

		## Depletion (B/B0, where B0 = B1933)
		SB = DB = mcmc.ts.sub[["SpawnBio"]]
		#B0 = SB[is.element(SB$Yr,startyr),]  ## use startyr (1935) not VIRG
		B0 = SB[is.element(SB$Yr,virginyr),]  ## switch to virgin year for SGR 2025 (RH 251105)
		for (i in unique(SB$Area)) {
			z  = is.element(SB$Area,i) & is.element(SB$Era,c("VIRG","INIT","TIME","FORE"))
			z0 = is.element(B0$Area,i) & is.element(B0$Era,"TIME")
			SBmat = as.matrix(DB[z,-(init.cols)])  ## matrix
			B0vec = unlist(B0[z0,-(init.cols)])    ## vector
			DBmat = sweep(SBmat, 2, B0vec, "/")
			DB[z,-(init.cols)] = DBmat
		}
		mcmc.ts.sub[["DB"]] = DB
	
		## Fraction recruitment (taken from r4ss::SSplotTimeseries)
		Recr   = mcmc.ts.sub[["Recruit_0"]]
		yvals  = Recr[,-(init.cols)]
		yvals2 = as.data.frame(array(NA, dim=c(length(Recr$Yr),nmcmc), dimnames=list(rownames(Recr),paste0("s",pad0(1:nmcmc,npad0)) ) ) )
		for (iyr in 1:nrow(yvals)) {
			y <- Recr$Yr[iyr]
			yvals2[iyr,] <- apply(yvals[Recr$Yr == y,],2,sum)
		}
		Frec <- data.frame(cols.init, yvals/yvals2)
		mcmc.ts.sub[["Frec"]] = Frec
		srun = run  ## sub ts run
		attr(mcmc.ts.sub,"samples") = names(mcmc.ts)
		save(list=c("mcmc.ts.sub","cattab","srun"), file="mcmc.ts.sub.rda")
	}
	## End subfunctions----------------------------

	## Extract relevant subset of time series values, if they exist, and bypass reconstruction
	dir.output = ifelse(loadCP, dir.extra, dir.mcmc)  ## if using CPs, change the output target directory
	setwd (dir.output)
	init.cols = 1:4; init.len = length(init.cols)  ## Initial columns commonly used by all MCMC samples; assumes they're sequential from 1 to n

	if (!any(RC) && file.exists(paste0(dir.output, "/mcmc.ts.sub.rda"))) {  ## always load it to be safe (unless any RC)
		.flush.cat("Loading mcmc.ts.sub object from saved binary\n")
		if (exists("mcmc.ts.sub", envir=.GlobalEnv))
			rm("mcmc.ts.sub", envir=.GlobalEnv)
		tic(); load("mcmc.ts.sub.rda", envir=.GlobalEnv); toc()
	} else if (!any(RC) && file.exists(paste0(dir.output, "/mcmc.ts.rda"))) {  ## load mcmc.ts if mcmc.ts.sub does not exist (unless any RC)
		.flush.cat("Loading mcmc.ts object from saved binary\n")
		if (exists("mcmc.ts", envir=.GlobalEnv))
			rm("mcmc.ts", envir=.GlobalEnv)
		if (exists("mcmc.ts.sub", envir=.GlobalEnv))
			rm("mcmc.ts.sub", envir=.GlobalEnv)
		tic(); load("mcmc.ts.rda", envir=.GlobalEnv); toc()
	} else {
		## Start the reconstruction from scratch
		## r4ss outputs names with non-sequential numbers if #iterations > 9999 (fixes padding of zeroes to 4)
		## Re-order extra files just in case (e.g., # iterations = 20,000):
		## Get the number of Report.sso files in the directory
		dir_list  = dir(dir.extra, full.names = TRUE)
		repfiles  = grep("/Report_mce_.*$", dir_list, value = TRUE)
		if (length(repfiles)==0)
			stop("Need the set of extra report files (on Linux)")
		repfiles  = orderRuns(grep("/Report_mce_.*$", dir_list, value = TRUE))
		## Get the number of CompReport.sso files in the directory
		compfiles  = grep("/CompReport_mce_.*$", dir_list, value = TRUE)
		if (RC[2]) {
			if (length(compfiles)==0)
				stop("Need the set of extra composition report files (on Linux)")
			else
				compfiles = orderRuns(grep("/CompReport_mce_.*$", dir_list, value = TRUE))
		}
		replist   = complist = list()
		#dump = gc(verbose=FALSE)

		## Report -- only do this once for gawd's sake
		if (!file.exists("reps.rda") && length(repfiles)>0)  ## check for need to start from scratch
			RC[1] = TRUE
		if (RC[1]) {  ## Reports files
			.flush.cat(paste0("Reading ", length(repfiles), " 'Reports' (ascii files)"), "\n")
			for (i in 1:length(repfiles)) {
				ii = padRun(repfiles[i], repfiles)
				replist[[ii]] = Kmisc::read(repfiles[i])  ## Kmisc::readlines crashes R
			}
			dump = gc(verbose=FALSE)
			tic(); reps = lapply(replist, function(x) { strsplit(x,split=ifelse(grepl("\\r\\n",x), "\\r\\n", "\\n"))[[1]] }); toc()
			#tic(); reps = mclapply(replist, function(x) { strsplit(x,split="\\r\\n")[[1]] }); toc()            ## 66.81 secs
			#tic(); reps = lapply(replist, function(x) { stringr::str_split(x, pattern="\\r\\n")[[1]] }); toc() ## 83.22 secs
			dump = gc(verbose=FALSE)
			rrun = run  ## reps run
			save("reps", "rrun", file="reps.rda")
		}
		## CompReport -- only do this once for gawd's sake (but really don't need it for POP|SGR)
		if (RC[2]) {  ## CompReports files
			.flush.cat(paste0("Reading ", length(compfiles), " 'CompReports' (ascii files)"), "\n")
			for (i in 1:length(compfiles)) {
				ii = padRun(compfiles[i], compfiles)
				complist[[ii]] = Kmisc::read(compfiles[i])  ## Kmisc::readlines crashes R
			}
			dump = gc(verbose=FALSE)
			tic(); comps = lapply(complist, function(x) { strsplit(x,split="\\r\\n")[[1]] }); toc() 
			dump = gc(verbose=FALSE)
			crun = run  ## comps run
			save("comps", "crun", file="comps.rda")
		}
		## Gather the time series elements
		## NOTE: parallel processing much slower than lapply
		## require(snow)
		## cl <- makeCluster(spec=detectCores()-1, type="SOCK")
		## tic(); out=clusterApply(cl, reps, function(x) {makeRepTab(dat=x, pat1="^TIME_SERIES", pat2="^SPR_SERIES")} ); toc() ## 635.52 sec = 10.6 min
		## stopCluster(cl)
	
		## Make time series object (mcmc.ts)
		make.ts(reps)  ## (RH 250922) now a subfunction
		tic(); load("mcmc.ts.rda", envir=.GlobalEnv); toc()
		## Make time series subset object (mcmc.ts.sub)
		make.ts.sub (mcmc.ts)  ## (RH 250922) now a subfunction
		tic(); load("mcmc.ts.rda", envir=.GlobalEnv); toc()
	}  ## end reconstruction of Reports into reps.rda, mcmc.ts.rda, and mcmc.ts.sub

	## ==========================================================================
	## Check to make sure reconstruction binaries mcmc.ts and mcmc.ts.sub were made
	if (!file.exists(paste0(dir.output, "/mcmc.ts.sub.rda")) ) {  ## recreate mcmc.ts.sub
		if (!file.exists(paste0(dir.output, "/mcmc.ts.rda")) ) {  ## recreate mcmc.ts.sub
			if (!exists("reps", envir=.GlobalEnv)) {
				if (!file.exists(paste0(dir.output, "/reps.rda")) ) {  ## recreate reps
					stop ("Either 'reps.rda' was not transferred to Windows (intentionally) or\n   reconstruction of Report files from 'sso' was not done in Linux")
				} else {
					tic(); load("reps.rda", envir=.GlobalEnv); toc()
				}
			}
			make.ts(reps)
			.flush.cat("Loading mcmc.ts object from saved binary\n")
			if (exists("mcmc.ts", envir=.GlobalEnv))
				rm("mcmc.ts", envir=.GlobalEnv)
			tic(); load("mcmc.ts.rda", envir=.GlobalEnv); toc()
		}
		make.ts.sub (mcmc.ts)
		.flush.cat("Loading mcmc.ts.sub object from saved binary\n")
		if (exists("mcmc.ts.sub", envir=.GlobalEnv))
			rm("mcmc.ts.sub", envir=.GlobalEnv)
		tic(); load("mcmc.ts.sub.rda", envir=.GlobalEnv); toc()
	} else {
		tic(); load("mcmc.ts.sub.rda", envir=.GlobalEnv); toc()
	}
		
#browser();return()
	#nmcmc = length(mcmc.ts) ## may not get loaded if code goes straight to 'mcmc.ts.sub'
	nmcmc = dim(mcmc.ts.sub[[1]])[2] - 4 ## because first 4 columns are c('Area,'Yr','Era','Seas') 
	npad0 = floor(log10(nmcmc)) + 1   ## for actual number of MCMC samples

	## Start the processing for figures etc.
	mcmc = mcmc.ts.sub ## just to save on typing
	options(scipen=10) ## stop displaying scientific notation (at least for the first 10 significant digits)
	##  "SpawnBio" "Recruit_0" "dead(B)" "Hrate" "sel(B)" "VB" "DB" "Frec" 
	## "dead(B)" and "sel(B)" don't change in MCMC samples
	plist = qlist = character()
	mcmc.names = c("SpawnBio", "Recruit_0", "u", "VB", "DB", "Frec")
	for (i in mcmc.names) {
		ii = switch(i, 'SpawnBio'="spawning", 'Recruit_0'="recruits", 'u'="exploitation", 'VB'="vulnerable", 'DB'="depletion", 'Frec'="frecruit")
		iii = switch(i, 'SpawnBio'="Spawning Biomass (kt)", 'Recruit_0'="Recruits (millions age-0 fish)", 'u'="Exploitation Rate (/y)", 'VB'="Vulnerable Biomass (kt)", 'DB'="Depletion (Bt/B0)", 'Frec'="Fraction Recruits")
		iv  = switch(i, 'SpawnBio'="B", 'Recruit_0'="R", 'u'="u", 'VB'="V", 'DB'="D", 'Frec'="fR")
		imcmc = mcmc[[i]]
		if (i == "SpawnBio") {
			## Extract SS3's estimation of VIRG biomass (V0, not vulnerable biomass) for allocation of MSY by area (suggested by PJS 230413)
			## Also calculate ORF's version of B0
			#B0.mcmc  = t(imcmc[is.element(imcmc$Era, c("TIME")) & is.element(imcmc$Yr, startyr),-c(init.cols)])
			B0.mcmc  = t(imcmc[is.element(imcmc$Era, c("VIRG")) & is.element(imcmc$Yr, virginyr),-c(init.cols)])  ## (RH 251105) switch to virgin year
			colnames(B0.mcmc) = areas
			B0.qmcmc = apply(B0.mcmc, 2, quantile, quants5, na.rm=TRUE)
			V0.mcmc  = t(imcmc[is.element(imcmc$Era, c("VIRG")),-c(init.cols)])
			colnames(V0.mcmc) = areas
			pVB.mcmc = t(apply(V0.mcmc, 1, function(x) { x / sum(x) }))    ## proportion Virgin Biomass used for allocating MSY
			if (length(areas)==1)
				pVB.mcmc = t(pVB.mcmc)
			colnames(pVB.mcmc) = areas
			V0.qmcmc  = apply(V0.mcmc, 2, quantile, quants5, na.rm=TRUE)
			pVB.med = V0.qmcmc['50%',] / sum(V0.qmcmc['50%',])  ## proportion SS3 VIRG B by area using median
			V0.mean = apply(V0.mcmc[,-c(init.cols)], 1, mean, na.rm=TRUE)
			pVB.mn = V0.mean / sum(V0.mean)                   ## proportion SS3 VIRG B by area using mean
			plist = c(qlist, c("B0.mcmc","V0.mcmc","pVB.mcmc"))
			qlist = c(qlist, c("B0.qmcmc","V0.qmcmc","pVB.med","pVB.mn"))
		}
		imcmc = imcmc[is.element(imcmc$Era, c("VIRG","TIME","FORE")),]  ## (RH 251105) include virgin values (1933 for SGR 2025)
		if (i=="DB") {
			## This only needs to be done if DB calculated in make.ts.sub() was using first year instead of virgin year
			## Just repeat the calculation to be sure (RH 251106)
			## Need to do this by area
			ssb = mcmc[["SpawnBio"]]
			ssb = ssb[is.element(ssb$Era, c("VIRG","TIME","FORE")),]
			for (a in sort(unique(ssb$Area))) {
				za = is.element(ssb$Area, a)
				ssb_vals = ssb[za,-c(init.cols)]
				ssb_virg = ssb_vals[1,]
				ssb_depl = sweep(as.matrix(ssb_vals),2,as.matrix(ssb_virg),FUN="/")
				if (sum(za) == nrow(ssb_depl)) {
					imcmc[za,-c(init.cols)] = ssb_depl
				} else {
					message("sumtingwong : nrows not same for 'imcmc' and 'ssb_depl'")
					browser(); return()
				}
			} ## end area loop
#browser();return()
		} ## end recalc of DB
		amcmc = split(imcmc, imcmc$Area)
		iyrs  = amcmc[[1]]$Yr
		endyr = rev(amcmc[[1]]$Yr[is.element(amcmc[[1]]$Era, c("TIME"))])[1]  ## first forecast year is actually the endyr of the model
		foryr = rev(amcmc[[1]]$Yr[is.element(amcmc[[1]]$Era, c("FORE"))])[1]

		pmcmc = lapply(amcmc, function(x){
			xx = t(x[,-c(init.cols)])
			colnames(xx) = iyrs
			return(xx)
		})
		qmcmc = lapply(amcmc, function(x){
			xx = apply(x[,-c(init.cols)], 1, quantile, quants5, na.rm=TRUE)
			colnames(xx) = iyrs
			return(xx)
		})
		names(pmcmc) = names(qmcmc) = areas
#if (iv=="D") {browser();return()}
		assign( paste0(iv, ".mcmc"), pmcmc)
		assign( paste0(iv, ".qmcmc"), qmcmc)
		plist = c(plist, paste0(iv, ".mcmc"))  ## posteriors added to a list for saving later
		qlist = c(qlist, paste0(iv, ".qmcmc"))
	}

	## Apply PJS allocation to MSY objects
	## -----------------------------------
	if (file.exists(paste0(dir.extra,"/derived_posteriors.sso")))
		dmcmc = SSgetMCMC (dir=dir.extra)
	else if (file.exists(paste0(dir.mcmc,"/derived_posteriors.sso")))
		dmcmc = SSgetMCMC (dir=dir.mcmc)
	else
		stop("Cannot find file 'derived_posteriors.sso'")

	## Sometimes need a subset of dmcmc
	smcmc = attributes(mcmc.ts.sub)$samples  ## available samples from original reps
	if (length(smcmc) != length(dmcmc$Iter)) {
		Amcmc = dmcmc$Iter #seq(thin, iters, thin)  ## all mcmc iterations
		Nmcmc = 1:length(Amcmc)
		pad.thai = floor(log10(max(Amcmc))) + 1  ## get common padding (RH 230905)
		names(Nmcmc) = pad0(Amcmc, pad.thai)
		smcmc = pad0(as.numeric(smcmc), pad.thai)
		Smcmc = Nmcmc[smcmc]
		dmcmc = dmcmc[Smcmc,]
	}
	#msy.mcmc = dmcmc[,grep("(^SSB|^annF|^Dead).*_MSY$",colnames(dmcmc),value=TRUE)]  ## I believe that F method 1 (Pope's approximation) give discrete F
	msy.mcmc = switch(benchmark,
		{ dmcmc[,grep("(^SSB|^annF|^Dead).*(_MSY$|_Btgt$)",colnames(dmcmc),value=TRUE)] },
		{ dmcmc[,grep("(^SSB|^annF|^Dead).*(_MSY$|_F01$)",colnames(dmcmc),value=TRUE)] }
	)
	colnames(msy.mcmc) = sub("Dead_Catch_MSY","MSY", sub("annF_MSY","Fmsy", sub("SSB_MSY","Bmsy",colnames(msy.mcmc))))
	colnames(msy.mcmc) = switch(benchmark,
		{ sub("Dead_Catch_Btgt","Ytgt", sub("annF_Btgt","Ftgt", sub("SSB_Btgt","Btgt",colnames(msy.mcmc)))) },
		{ sub("Dead_Catch_F01","Y01", sub("annF_F01","F01", sub("SSB_F01","B01",colnames(msy.mcmc)))) }
	)
	msy.mcmc[is.na(msy.mcmc)] = NA  ## change NaN to NA

	if (Fmethod==1) {
		colnames(msy.mcmc)[grep("Fmsy",colnames(msy.mcmc))] = "umsy"
		msy.mcmc$Fmsy = -log(1 - msy.mcmc$umsy)
		collect = c("Bmsy","Fmsy","MSY","umsy")
		if ("Ftgt" %in% colnames(msy.mcmc)) {
			colnames(msy.mcmc)[grep("Ftgt",colnames(msy.mcmc))] = "utgt"
			msy.mcmc$Ftgt = -log(1 - msy.mcmc$utgt)
			collect = c(collect, "Btgt","Ftgt","Ytgt","utgt")
		}
		if ("F01" %in% colnames(msy.mcmc)) {
			colnames(msy.mcmc)[grep("F01",colnames(msy.mcmc))] = "u01"
			msy.mcmc$F01 = -log(1 - msy.mcmc$u01)
			collect = c(collect, "B01","F01","Y01","u01")
		}
		msy.mcmc = msy.mcmc[,collect]
	} else {
		msy.mcmc$umsy = 1 - exp(-msy.mcmc$Fmsy)
		collect = c("Bmsy","Fmsy","MSY","umsy")
		if ("Ftgt" %in% colnames(msy.mcmc)) {
			msy.mcmc$utgt = 1 - exp(-msy.mcmc$Ftgt)
			collect = c(collect, "Btgt","Ftgt","Ytgt","utgt")
		}
		if ("F01" %in% colnames(msy.mcmc)) {
			msy.mcmc$u01 = 1 - exp(-msy.mcmc$F01)
			collect = c(collect, "B01","F01","Y01","u01")
		}
		msy.mcmc = msy.mcmc[,collect]
	}
	msy.mcmc$LRP = 0.4 * msy.mcmc$Bmsy
	msy.mcmc$USR = 0.8 * msy.mcmc$Bmsy
#tget(cattab)  ## temporary for debuggng : remove before send code to Linux
	endcat = cattab[as.character(endyr),"total"]  ## end-year catch for the BC coast
	msy.mcmc$Vmsy = endcat/msy.mcmc$umsy
	if ("Ftgt" %in% colnames(msy.mcmc))
		msy.mcmc$Vtgt = endcat/msy.mcmc$utgt
	if ("F01" %in% colnames(msy.mcmc))
		msy.mcmc$V01 = endcat/msy.mcmc$u01

	MSY.mcmc = MSY.qmcmc = list()
	## Allocate based on proportion Virgin Biomass (=VB ***** NOT Vulnerable Biomss)
	for (a in areas) {
		paVB = pVB.mcmc[,a]  ## proportion Virgin Biomass by area for allocating MSY (PJS 230413)
		amsy.mcmc  = sweep(msy.mcmc,1,paVB,"*")
		## Cannot do this to umsy so have to recalculate from Vmsy
		endcat = cattab[as.character(endyr), a]   ## end-year catch for each area
		amsy.mcmc$umsy = endcat/amsy.mcmc$Vmsy
		amsy.mcmc$umsy = pmin(amsy.mcmc$umsy, 0.99)  ## sometimes ad hoc calculation exceeds 1 (RH 251105)
		amsy.mcmc$Fmsy = -log(1 - amsy.mcmc$umsy) ## adjust area-specific Fmsy also, calculating from adjusted umsy
		if ("Ftgt" %in% colnames(amsy.mcmc)) {
			amsy.mcmc$utgt = pmin(endcat/amsy.mcmc$Vtgt, 0.99)
			amsy.mcmc$Ftgt = -log(1 - amsy.mcmc$utgt)
		}
		if ("F01" %in% colnames(amsy.mcmc)) {
			amsy.mcmc$u01 = pmin(endcat/amsy.mcmc$V01, 0.99)
			amsy.mcmc$F01 = -log(1 - amsy.mcmc$u01)
		}
		MSY.mcmc[[a]]  = amsy.mcmc
		MSY.qmcmc[[a]] = sapply(MSY.mcmc[[a]],quantile,quants5, na.rm=TRUE)
	}
#browser();return()
	#msy.names = c("BtLRP","BtUSR","BtBmsy","utumsy")
	msy.names = c("BtBmsy","utumsy")
	if ("Ftgt" %in% colnames(msy.mcmc))  msy.names = c(msy.names, c("BtBtgt","ututgt"))
	if ("F01" %in% colnames(msy.mcmc))   msy.names = c(msy.names, c("BtB01","utu01"))

	for (i in msy.names) {
		ii = paste0(i, ".mcmc")
		iii = paste0(i, ".qmcmc")
		assign(ii, list())
		assign(iii, list())
		for (a in areas) {
			aa = substring(i,1,1); bb = sub("^Bt|^ut","",i)
			pmess = paste0("num=",aa,".mcmc[[\"",a,"\"]]; den=MSY.mcmc[[\"",a,"\"]][,\"",bb,"\"]; ",ii,"[[\"",a,"\"]]=sweep(num,1,den,\"/\")")
#if (i==msy.names[3]) {browser();return()}
			eval(parse(text=pmess))
			qmess = paste0(iii,"[[\"",a,"\"]] = apply(",ii,"[[\"",a,"\"]],2,quantile,quants5,na.rm=T)")
			eval(parse(text=qmess))
		}
	}
	plist = c(plist, c("MSY.mcmc", paste0(msy.names,".mcmc")))   ## posteriors names
	qlist = c(qlist, c("MSY.qmcmc", paste0(msy.names,".qmcmc"))) ## quantiles names
#browser();return()

	save(list=plist, file="mcmc.posts.rda")  ## saves posteriors like B.mcmc
	save(list=qlist, file="mcmc.quants.rda") ## save quantiles of posteriors

	if (plot) {
		plotArea = function(qlist, area, hline=NULL, cp=NULL) {
			scale = ifelse (i %in% c("SpawnBio","Recruit_0", "VB"), 1000., 1.)
			ylim  =  range(unlist(qlist[area])/scale, na.rm=TRUE); ylim[1] = 0
#browser();return()
			acols = character()
			plot(0,0, xlim=xlim, ylim=ylim, type="n", xaxs="i", xaxt="n", xlab="", ylab="", cex.axis=1.2, cex.lab=1.5)
			if (!is.null(hline)) {
				hcol = rep("slategray",length(hline))
				if (i %in% c("BtBmsy","DB")) hcol = c("red","green4")
				hcol[match(1,hline)] = "purple"
				abline(h=hline, col="gainsboro", lty=1)
				abline(h=hline, col=hcol, lty="37")
			}
			axis (1, at=seq(1900,2100,5), tcl=-0.2, labels=FALSE)
			axis (1, at=seq(1900,2100,10), tcl=-0.5, labels=ifelse(par()$mfg[1]==par()$mfg[3], TRUE, FALSE), cex.axis=1.4)
			for (a in 1:length(area)){
				aa    = area[a]
				amcmc = qlist[[aa]] / scale
				acol  = switch(aa, '5ABC'="blue", '3CD'="green3", '5DE'="red")
				acols = c(acols,acol)
				env.pyrs.x = c(pyrs, rev(pyrs))
				env.pyrs.y = c(amcmc[1,as.character(pyrs)], rev(amcmc[5,as.character(pyrs)]))
				polygon(env.pyrs.x, env.pyrs.y, col=lucent(acol,0.2), border=FALSE)
				lines(pyrs, amcmc[1,as.character(pyrs)], col=acol, lty=2, lwd=1)
				lines(pyrs, amcmc[5,as.character(pyrs)], col=acol, lty=2, lwd=1)
				lines(pyrs, amcmc[2,as.character(pyrs)], col=acol, lty=3, lwd=1)
				lines(pyrs, amcmc[4,as.character(pyrs)], col=acol, lty=3, lwd=1)
				lines(pyrs, amcmc[3,as.character(pyrs)], col=acol, lty=1, lwd=3)
				env.fyrs.x = c(fyrs, rev(fyrs))
				env.fyrs.y = c(amcmc[1,as.character(fyrs)], rev(amcmc[5,as.character(fyrs)]))
				polygon(env.fyrs.x, env.fyrs.y, col=lucent("gainsboro",0.2), border=FALSE)
				lines(fyrs, amcmc[1,as.character(fyrs)], col=acol, lty=2, lwd=1)
				lines(fyrs, amcmc[5,as.character(fyrs)], col=acol, lty=2, lwd=1)
				lines(fyrs, amcmc[2,as.character(fyrs)], col=acol, lty=3, lwd=1)
				lines(fyrs, amcmc[4,as.character(fyrs)], col=acol, lty=3, lwd=1)
				lines(fyrs, amcmc[3,as.character(fyrs)], col=acol, lty=1, lwd=3)
			}
			mtext("Year", side=1, outer=T, line=2.75, cex=1.4)
			mtext(iii, side=2, outer=T, line=0.5, cex=1.4)
			big.mark = ifelse(is.null(options()$big.mark), ",", options()$big.mark)
#browser();return()
			legtxt   = if (is.null(cp)) area else paste0(area, " (proj.catch = ", formatC(as.numeric(ifelse(length(areas)==1,cp,cp[area])),digits=0,format="f",big.mark=big.mark), " t/y)")
			addLegend(0.99, 0.98, legend=legtxt, xjust=1, yjust=1, lty=1, col=acols, seg.len=3, bty="n")
		} ## end function 'plotArea'

		scroll.names = c("utumsy", "BtBmsy", mcmc.names)
		if ("F01" %in% colnames(msy.mcmc))  scroll.names = c("utu01", "BtB01", scroll.names)
		if ("Ftgt" %in% colnames(msy.mcmc))  scroll.names = c("ututgt", "BtBtgt", scroll.names)
		for ( i in scroll.names ) {
			ii = switch(i, 'SpawnBio'="spawning", 'Recruit_0'="recruits", 'u'="exploitation", 'VB'="vulnerable", 'DB'="depletion", 'Frec'="frecruit", 'BtBmsy'="BtBmsy", 'utumsy'="utumsy", 'BtBtgt'="BtBtgt", 'ututgt'="ututgt", 'BtB01'="BtB01", 'utu01'="utu01")
			iii = switch(i, 'SpawnBio'="Spawning Biomass (kt)", 'Recruit_0'="Recruits (millions age-0 fish)", 'u'="Exploitation Rate (/y)", 'VB'="Vulnerable Biomass (kt)", 'DB'="Depletion (Bt/B0)", 'Frec'="Fraction Recruits",
			'BtBmsy'=expression(italic(B)[italic(t)]/italic(B)[MSY]), 'utumsy'=expression(italic(u)[italic(t)]/italic(u)[MSY]),
			'BtBtgt'=expression(italic(B)[italic(t)]/0.4~italic(B)[0]), 'ututgt'=expression(italic(u)[italic(t)]/italic(u)[0.4~italic(B)~0]),
			'BtB01'=expression(italic(B)[italic(t)]/italic(B)[F01]), 'utu01'=expression(italic(u)[italic(t)]/italic(u)[F01]) )
			iv  = switch(i, 'SpawnBio'="B", 'Recruit_0'="R", 'u'="u", 'VB'="V", 'DB'="D", 'Frec'="fR", 'BtBmsy'="BtBmsy", 'utumsy'="utumsy", 'BtBtgt'="BtBtgt", 'ututgt'="ututgt", 'BtB01'="BtB01", 'utu01'="utu01")

			## Plot quantiles for each factor by area
			xlim  = range(iyrs)
			#ylim  =  range(unlist(qmcmc))
			pyrs  = startyr:endyr
			fyrs  = endyr:foryr
			qmcmc = get(paste0(iv,".qmcmc"))
			cp    = cattab[as.character(endyr+1),areas]  ## catch policies

			if (vertical) {
				rc = if(length(areas)==1) c(1,1) else c(length(areas)+ifelse(show.combo,1,0),1) ## 230822
			} else {
				rc = if(length(areas)==1) c(1,1) else .findSquare(length(areas)+1) ## 230822
			}
			fout.e = paste0("extra.mcmc(", ii, ")")
			if (prod(rc)>1){
				fout.e = paste0(fout.e, "-", prod(rc), "panels")
				if (!vertical)  fout.e = paste0(fout.e, "(square)")
			}
#browser();return()
			for (l in lang) {
				changeLangOpts(L=l)
				fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				if (png) png(filename=paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				#expandGraph(mfrow=c(ifelse(length(areas)==1,1,length(areas)+1),1), mar=c(0,2,0,1), oma=c(4,2.5,1,0), mgp=c(1.6,0.75,0))
				expandGraph(mfrow=rc, mar=c(0,2,0,1), oma=c(4,2.5,1,0), mgp=c(1.6,0.75,0))
				hline = if (i=="BtBmsy") rep(list(c(0.4,0.8)),length(areas))
					else if (i=="utumsy") rep(list(c(1)),      length(areas))
					else if (i=="DB")     rep(list(c(0.2,0.4)),length(areas))
					else if (i %in% c("ututgt","utu01")) lapply(qmcmc, function(x) { x["95%", as.character(endyr)] })
					else NULL
				for (a in 1:length(areas))
					plotArea(qlist=qmcmc, area=areas[a], hline=hline[[a]], cp=cp)
				if (length(areas)>1 && show.combo)
					plotArea(qlist=qmcmc, area=areas)
				if (png) dev.off()
			}; eop()
#browser();return()
		} ## end i loop for plots
	} ## end if plot
	## check L436
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~load_extra_mcmc


## makeFSARfigs-------------------------2025-12-01
##  Make figures for the new FSAR
## ---------------------------------------------RH
#makeFSARfigs <- function (xTS, xRP, xPJ, years=1935:2024, TAC,
makeFSARfigs <- function (envo, years=1935:2025, TAC, RPbase="B0",
   png=F, pngres=400, PIN=c(10,7), lang="e",
   fig.p4=TRUE, fig.snail=FALSE, fig.catch=FALSE, fig.hbar=FALSE, fig.proj=FALSE)
{
	vomit <- function() { gc(verbose=FALSE); while ("ee" %in% search()) detach(ee) }; #resetGraph() }
	on.exit( vomit() )
	#par.old = par(no.readonly=TRUE)
	#on.exit(par(par.old))
	for (e in 1:length(envo)) {
		keep.pars = ls()  ## need to clear environment of objects on subsequent loops through envo
		ee = envo[[e]]
		attach(ee)
		rm(list=intersect(ls(.GlobalEnv), ls(ee)), envir=.GlobalEnv)  ## remove objects in global environment masking objects in e1
		unpackList(compo)

		xTS=xavgTS; xRP=xavgRP; xPJ=xavgPJ
		if (is.na(xTS)) {  ## the single-area model
			xTS=avgTS; xRP=avgRP; xPJ=avgPJ
		}

		narea    = length(area.names)  ## from e2
		col.area = c("blue","darkgreen","red","black")
		bg.area  = c("cyan","green","pink","gainsboro")
		lty.area = rep(1,4) #c(2, 3, 5, 1)
		pch.area = c(22:24, 21)# c(15,18,17,16) #
		col.refs = c(.colBlind[c("redpurple","bluegreen","skyblue")]) ## LRP, USR, RR
		lty.refs = c(4,5,2)
	
		## Get catch time series from the base run
		replist = SS_output(dir=central.mpd.dir, verbose=F, printstats=F)
		repcats = split(replist$catch, replist$catch$Fleet)
		fleetcat = lapply(repcats, function(x){x[,c("Obs")]})
		catch    = do.call("cbind", lapply(fleetcat, data.frame, stringsAsFactors=FALSE))
		rownames(catch) = .su(replist$catch$Yr)
		#if (all(is.numeric(replist$catch$Area))) { ## not using extra area info; assume single-area model
		if (length(unique(replist$catch$Area))==1) { ## not using extra area info; assume single-area model
			if (ncol(catch)>1) { ## could be multiple fisheries in single area
				catch = as.data.frame(matrix(apply(catch,1,sum), ncol=1))
				rownames(catch) = .su(replist$catch$Yr)
			}
			areas = "BC" ##.su(sapply(strsplit(replist$catch$Fleet_Name,split="_"),function(x){rev(x)[1]})) ## can have numerous fisheries in one area
			## Need to expand arrays with area dimension for code later on
			xRP.new = array(NA, dim=c(dim(xRP),length(areas)), dimnames=c(dimnames(xRP),list(area=areas)))
	#browser();return()
			xRP.new[1:nrow(xRP),1:ncol(xRP),areas] = xRP
			xRP.old = xRP; xRP = xRP.new
			xTS.new = array(NA, dim=c(dim(xTS),length(areas)), dimnames=c(dimnames(xTS),list(area=areas)))
			xTS.new[1:nrow(xTS),1:ncol(xTS),1:dim(xTS)[3],areas] = xTS
			xTS.old = xTS; xTS = xTS.new
	#browser();return()
		} else {
			areas = .su(replist$catch$Area)  ## from multi-area model (only valid when there is one fleet per area)
		}
		colnames(catch) = areas
		nyrs   = length(years)
		Cscale = 1.; catch = round(catch/Cscale,3)
		if (!missing(TAC))
			TAC[,-1] = TAC[,-1]/Cscale
		## Default figures (for SGR 2025) -- 4-panel
		## -----------------------------------------
		if (fig.p4) {
			for (a in 1:ncol(catch)) {
				aa = colnames(TAC)[-1][a]
				fout.e = paste0("FSAR.4panel.", spp.code, ".", aa)
				for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
					changeLangOpts(L=l)
					#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
					fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
					if (png) { 
						clearFiles(paste0(fout,".png"))
						png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
					}
					expandGraph(mfrow=c(2,2), mar=c(2.75,3.75,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
	
					## Plot 1 : catch (top left)
					plot(0,0, xlim=range(years), ylim=c(0,max(catch[,a],na.rm=T)), xlab=linguaFranca("Year",l), ylab="", cex.axis=1.0, cex.lab=1.2, las=1)
					axis(1, at=seq(1935,2025,5), labels=F, tcl=-.3)
					mtext(linguaFranca(paste0("Catch (", ifelse(Cscale==1000,"k",""), "t)"),l), side=2, line=2.5, cex=1.1)
					za = is.element(rownames(catch), years)
					points(TAC$Year, TAC[,aa], pch= convUTF("\\u{2584}"), col="orange")
					lines(years, catch[za,a], col=col.area[a], lwd=2)
					addLegend(ifelse(strSpp %in% c("405"),0.5,0.95), 0.975, col=c(col.area[a],"orange"), lty=c(1,2), legend=linguaFranca(c(paste0("Catch (", ifelse(Cscale==1000,"kilo",""), "tonnes)"), paste0("TAC (", ifelse(Cscale==1000,"kilo",""), "tonnes)")),l), bty="n", xjust=1, lwd=3, seg.len=3)
					addLabel(0.05, 0.95, "(A)", cex=1.5, col="black", adj=c(0,1))
					addLabel(0.05, 0.85, linguaFranca(aa,l), cex=1.5, col=col.area[a], adj=c(0,1))
#browser();return()
	
	#				## Plot 2 : biomass (top right)
	#				yy = as.character(years)
	#				Bscale = 1000.
	#				B.qts   = apply(xTS[,yy,"B",aa],2, quantile, probs=ttcall(quants3), na.rm=T) / Bscale
	#				BRP.qts = apply(xRP[,c("LRP","USR"),aa], 2, quantile, probs=ttcall(quants3), na.rm=T) / Bscale
	#				plot(0,0, xlim=range(years), ylim=c(0,max(B.qts,na.rm=T)), xlab="Year", ylab="", cex.axis=1.0, cex.lab=1.2, las=1)
	#				axis(1, at=seq(1935,2025,5), labels=F, tcl=-.3)
	#				mtext(paste0("Female Spawning Biomass (", ifelse(Bscale==1000,"k",""), "t)"), side=2, line=2.2, cex=1.1)
	#				polygon(c(years[c(1,nyrs)], years[c(nyrs,1)]), c(BRP.qts[c(1,1),"LRP"],BRP.qts[c(3,3),"LRP"]), col=lucent(col.refs[1],0.35), border=FALSE)
	#				lines(years, rep(BRP.qts[2,"LRP"],length(years)), col="black", lwd=1, lty=4)
	#				polygon(c(years[c(1,nyrs)], years[c(nyrs,1)]), c(BRP.qts[c(1,1),"USR"],BRP.qts[c(3,3),"USR"]), col=lucent(col.refs[2],0.35), border=FALSE)
	#				lines(years, rep(BRP.qts[2,"USR"],length(years)), col="black", lwd=1, lty=5)
	#				lines(years, B.qts[1,], col="gainsboro", lwd=1, lty=1) ## just to add a bit more emphasis
	#				lines(years, B.qts[1,], col=col.area[a], lwd=1, lty=3)
	#				lines(years, B.qts[3,], col="gainsboro", lwd=1, lty=1) ## just to add a bit more emphasis
	#				lines(years, B.qts[3,], col=col.area[a], lwd=1, lty=3)
	#				lines(years, B.qts[2,], col=col.area[a], lwd=2)
	#				addLegend(0.95, 0.975, col=c(rep(col.area[a],2),rep("black",2)), lty=c(1,3,lty.refs[2:1]), legend=c("Median biomass","90% credibility envelope", "Median USR","Median LRP"), bty="n", xjust=1)
	#				addLabel(0.05, 0.95, "(B)", cex=1.5, col="black", adj=c(0,1))
	
					## Plot 2 : biomass relative to Bmsy (top right)
					yy = as.character(years)
					B.qts   = apply(xTS[,yy,switch(RPbase,'BMSY'="BtBmsy",'B0'="BtB0",'BTRP'="BtBtrp"),aa],2, quantile, probs=ttcall(quants3), na.rm=T) 
					plot(0,0, xlim=range(years), ylim=c(0,max(B.qts,na.rm=T)), xlab=linguaFranca("Year",l), ylab="", cex.axis=1.0, cex.lab=1.2, las=1)
					axis(1, at=seq(1935,2025,5), labels=F, tcl=-.3)
					ylab = linguaFranca(paste0("Spawning Biomass relative to ", RPbase), l)
					ylab = sub("B(MSY|RMD|0|TRP|PRC)","italic(B)[\\1]", gsub("\\s+","~",ylab))
	#browser();return()
					mtext(eval(parse(text=paste0("expression(", ylab, ")"))), side=2, line=2.3, cex=1.1)
					#mtext(linguaFranca(expression(Spawning~Biomass~relative~to~italic(B)[MSY]),l), side=2, line=2.3, cex=1.1)
					BURP = c(0.4, 0.8)
					if (RPbase %in% c("B0"))
						BURP = 0.4 * BURP
					lines(years, rep(BURP[1],length(years)), col=col.refs[1], lwd=1.5, lty=4)
					lines(years, rep(BURP[2],length(years)), col=col.refs[2], lwd=1.5, lty=5)
					lines(years, B.qts[1,], col="gainsboro", lwd=1, lty=1) ## just to add a bit more emphasis
					lines(years, B.qts[1,], col=col.area[a], lwd=1, lty=3)
					lines(years, B.qts[3,], col="gainsboro", lwd=1, lty=1) ## just to add a bit more emphasis
					lines(years, B.qts[3,], col=col.area[a], lwd=1, lty=3)
					lines(years, B.qts[2,], col=col.area[a], lwd=2)
					addLegend(0.95, 0.975, col=c(rep(col.area[a],2),col.refs[2:1]), lty=c(1,3,lty.refs[2:1]), legend=linguaFranca(c("Median relative biomass","90% credibility envelope", "USR","LRP"),l), bty="n", xjust=1, lwd=1.5, seg.len=3)
					addLabel(0.05, 0.95, "(B)", cex=1.5, col="black", adj=c(0,1))
#browser();return()
	
					## Plot 3 : exploitation (bottom left)
					yy = as.character(years)
	#browser();return()
					u.qts   = apply(xTS[,yy, grep("^u$|^ut$",dimnames(xTS)[[3]]), aa],2, quantile, probs=ttcall(quants3), na.rm=T)
					U.qts   = apply(xTS[,yy,switch(RPbase,'BMSY'="utumsy",'B0'="ututrp",'BTRP'="ututrp"),aa],2, quantile, probs=ttcall(quants3), na.rm=T)
					useRR = TRUE
					if(useRR) {
						u.qts.raw = u.qts ## retain just in case (not used further)
						u.qts     = U.qts ## use ut/umsy (PJS suggestion)
					}
					plot(0,0, xlim=range(years), ylim=c(0,max(u.qts,na.rm=T)), xlab=linguaFranca("Year",l), ylab="", cex.axis=1.0, cex.lab=1.2, las=1)
					axis(1, at=seq(1935,2025,5), labels=F, tcl=-.3)
					if (useRR) {
						ylab = linguaFranca(paste0("Exploitation rate relative to ", switch(RPbase,'BMSY'="uMSY",'B0'="uTRP",'BTRP'="uTRP")), l)
						ylab = sub("u(MSY|RMD|TRP|PRC)", "italic(u)[\\1]", gsub("\\s+","~",ylab))
						mtext(eval(parse(text=paste0("expression(", ylab, ")"))), side=2, line=2.3, cex=1.1)
						lines(years, rep(1,length(years)), col=col.refs[2], lwd=1.5, lty=5)
					} else {
						mtext(linguaFranca("Exploitation rate (per year)",l), side=2, line=2.5, cex=1.1)
					}
					lines(years, u.qts[1,], col="gainsboro", lwd=1, lty=1) ## just to add a bit more emphasis
					lines(years, u.qts[1,], col=col.area[a], lwd=1, lty=3)
					lines(years, u.qts[3,], col="gainsboro", lwd=1, lty=1) ## just to add a bit more emphasis
					lines(years, u.qts[3,], col=col.area[a], lwd=1, lty=3)
					lines(years, u.qts[2,], col=col.area[a], lwd=2)
					if (useRR) {
						if (strSpp %in% c("405"))
							addLegend(0.025, 0.9, col=c(rep(col.area[a],2),col.refs[2]), lty=c(1,3,lty.refs[2]), legend=linguaFranca(c("Median relative exploitation","90% credibility envelope", "RRR"),l), bty="n", xjust=0, lwd=1.5, seg.len=3)
						else
							addLegend(0.95, 0.975, col=c(rep(col.area[a],2),col.refs[2]), lty=c(1,3,lty.refs[2]), legend=linguaFranca(c("Median relative exploitation","90% credibility envelope", "RRR"),l), bty="n", xjust=1, lwd=1.5, seg.len=3)
					} else {
						addLegend(0.95, 0.975, col=c(rep(col.area[a],2)), lty=c(1,3), legend=linguaFranca(c("Median exploitation","90% credibility envelope"),l), bty="n", xjust=1, lwd=1.5, seg.len=3)
					}
					addLabel(0.05, 0.95, "(C)", cex=1.5, col="black", adj=c(0,1))
#browser();return()
	
					## Plot 4 : recruitment (bottom right)
					yy = as.character(years)
					Rscale = 1000.
					R.qts   = apply(xTS[,yy, grep("^R$|^Rt$",dimnames(xTS)[[3]]), aa],2, quantile, probs=ttcall(quants3), na.rm=T) / Rscale
					plot(0,0, xlim=range(years), ylim=c(0,max(R.qts,na.rm=T)), xlab=linguaFranca("Year",l), ylab="", cex.axis=1.0, cex.lab=1.2, las=1)
					axis(1, at=seq(1935,2025,5), labels=F, tcl=-.3)
					mtext(linguaFranca(paste0("Recruitment (", ifelse(Rscale==1000,"millions","1000s"), " age-0 fish)"),l), side=2, line=2.3, cex=1.1)
					lines(years, R.qts[1,], col="gainsboro", lwd=1, lty=1) ## just to add a bit more emphasis
					lines(years, R.qts[1,], col=col.area[a], lwd=1, lty=3)
					lines(years, R.qts[3,], col="gainsboro", lwd=1, lty=1) ## just to add a bit more emphasis
					lines(years, R.qts[3,], col=col.area[a], lwd=1, lty=3)
					lines(years, R.qts[2,], col=col.area[a], lwd=2)
					addLegend(ifelse(strSpp %in% c("405"),0.65,0.95), 0.975, col=c(rep(col.area[a],2)), lty=c(1,3), legend=linguaFranca(c("Median recruitment","90% credibility envelope"),l), bty="n", xjust=1, lwd=1.5, seg.len=3)
					addLabel(0.05, 0.95, "(D)", cex=1.5, col="black", adj=c(0,1))
					if (png) dev.off()
				}; eop()
	#browser();return()
			} ## end a (catch)
		} ## end include 4 panel
browser();return()
	
		## Stock status 4-panel plot for GMU
		## ---------------------------------
		if(fig.snail) {
			fout.e = paste0("FSAR.4GMU.subareas")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				if (png) {
					clearFiles(paste0(fout,".png"))
					png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				}
				expandGraph(mfrow=c(2,2), mar=c(2.75,3.75,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		
				## Panel 1 -- stock status (taken from 'plotSS.compo.r')
				load("Bbase.for.FSAR.rda")
				subs = 1:length(Mnams)
				rm.cst = T
				if (rm.cst) {
					subs  = 2:4
					Bbase = lapply(Bbase,function(x){x[subs]})
				}
				out = compBmsy(Bspp=Bbase, spp=L1, boxwidth=0.4, medcol=col.area[subs-1], boxfill=boxfill[subs], boxlim=boxlim, whisklwd=2, staplelwd=3, Mnams=gsub("\\s+","\n",Mnams[subs]), figgy=list(win=TRUE), pngres=pngres, spplabs=F, t.yr=currYear, left.space=20, top.space=1, fout="base.subarea.status", calcRat=F, lang=l, add=T, cex.lab=0.8, rlty=c(4,2))
				addLabel(0.975, 0.975, "(A)", cex=1.5, col="black", adj=c(1,1)) 
		#browser();return()
		
				## Panel 2-4 -- phase plots for 5ABC, 3CD, 5DE
				yy = as.character(years)
				## Need to recreate object found in 'mcmc.quants.rda' from 'load_extra_mcmc.r")
				BtBmsy = apply(xTS[,yy,"BtBmsy",,drop=FALSE], c(2,4), quantile, probs=ttcall(quants5), na.rm=T)
				utumsy = apply(xTS[,yy,"utumsy",,drop=FALSE], c(2,4), quantile, probs=ttcall(quants5), na.rm=T)
				BtBmsy.qmcmc = utumsy.qmcmc = list()
				for (a in 1:dim(BtBmsy)[3]) {
					aa = dimnames(BtBmsy)$area[a]
					BtBmsy.qmcmc[[aa]] = BtBmsy[,,aa]
					utumsy.qmcmc[[aa]] = utumsy[,,aa]
				}
		#browser();return()
				for (i in 1:dim(BtBmsy)[3]) {
					#ii = switch(i, "5ABC", "3CD", "5DE")
					ii = dimnames(BtBmsy)$area[i]
					plotSnail(BtBmsy.qmcmc, utumsy.qmcmc, yrs=1935:2023, p=c(0.05,0.95), xLim=NULL, yLim=NULL, ngear=NULL, assYrs=assYrs, outs=F, Cnames=NULL, Lwd=2, ptypes="win", outnam=paste0("snail.",ii), lang=l, labYrs=c(1950,seq(1960,1975,5), assYrs), onepanel=F, subarea=i, add=T)
					addLabel(0.975, 0.975, paste0("(",LETTERS[i+1],")"), cex=1.5, col="black", adj=c(1,1))
				}
		#browser();return()
				if (png) dev.off()
			}; eop()
		}
	
		## Catch all in one for catch section (proposal)
		## ---------------------------------------------
		if (fig.catch) {
			catch$'4' = apply(catch,1,sum,na.rm=T)
			aa = c(colnames(TAC)[-1],"Coast")
			fout.e = "FSAR.Catch.All"
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				if (png) {
					clearFiles(paste0(fout,".png"))
					png(paste0(fout,".png"), units="in", res=pngres, width=12, height=6)
				}
				expandGraph(mfrow=c(1,1), mar=c(2.75,3.75,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
				plot(0,0, xlim=range(years), ylim=c(0,max(catch,na.rm=T)), xlab="Year", ylab="", cex.axis=1.0, cex.lab=1.2, las=1)
				axis(1, at=seq(1935,2025,5), labels=F, tcl=-.3)
				axis(2, at=seq(1,50,1), labels=F, tcl=-.2)
				mtext(paste0("Catch (", ifelse(Cscale==1000,"k",""), "t)"), side=2, line=2.2, cex=1.2)
				for (a in 1:ncol(catch)) {
					za = is.element(rownames(catch), years)
					lines(years, catch[za,a], col=col.area[a], lty=lty.area[a], lwd=2)
					points(years, catch[za,a], col=col.area[a], bg=bg.area[a], pch=pch.area[a], lwd=0.75, cex=0.8)
					#points(years, catch[za,a], col="gainsboro", bg=col.area[a], pch=pch.area[a], lwd=0.5, cex=0.9) #col.area[a]
					#points(years, catch[za,a], col="gainsboro", bg=bg.area[a], pch=pch.area[a], lwd=0.5, cex=1) #col.area[a]
					#points(years, catch[za,a], col=bg.area[a], pch=pch.area[a], lwd=0.5, cex=1) #col.area[a]
				}
				addLegend(0.95, 0.975, pch=pch.area, col=col.area, pt.bg=bg.area, lty=lty.area, legend=aa, bty="n", xjust=1, seg.len=3)
				if (png) dev.off()
			}; eop()
		}
	
		## Merged stock status for PJS
		## ---------------------------
		if (fig.hbar) {
			for (ss.type in c(1,2,3)[3]) {
				Mnams.all = character()
				## Load stock status choice
				if (ss.type %in% c(1,3)) {
					load("Bbase.for.FSAR.rda")
					Bss = Bbase
					Mnams.all = c(Mnams.all, Mnams)
				}
				## Splice together another stock status set
				if (ss.type %in% c(2,3)) {
					if (!"Bsens.for.PJS.rda" %in% list.files())
						stop("Need 'Bsens.for.PJS.rda' from AppF")
					load("Bsens.for.PJS.rda")
					Bss = Bsens
					Mnams.all = c(Mnams.all, Mnams)
				}
				if (ss.type %in% c(3)) {
					Bss = Bbase
					Bss[[1]] = c(Bss[[1]], Bsens[[1]])  ## merge
					boxfill  = rep(boxfill,2)
					medcol   = rep(medcol,2)
				}
				fout.e = paste0("PJS.SS.", switch(ss.type, "subareas", "single.areas", "both.areas"))
				for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
					changeLangOpts(L=l)
					#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
					fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
					if (png) {
						clearFiles(paste0(fout,".png"))
						png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
					}
					expandGraph(mfrow=c(1,1), mar=c(3.5,6,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
	
					subs = 1:length(Mnams)
					rm.cst = T
					if (rm.cst) {
						subs  = grep("CST|Base",names(Bss[[1]]),invert=T)
						Bss = lapply(Bss, function(x){ x[subs] })
					}
	#browser();return()
					Mnams = gsub("\\s+", "\n", sub("A[0-9]\\s+\\(R\\d+\\)\\ss+", "S", Mnams.all[subs]))
					out = compBmsy(Bspp=Bss, spp=L1, boxwidth=0.4, medcol=medcol[subs], boxfill=boxfill[subs], boxlim=c(0,6.25), whisklwd=2, staplelwd=3, Mnams=Mnams, figgy=list(win=TRUE), pngres=pngres, spplabs=F, t.yr=currYear, left.space=7, top.space=1, fout="stock.status", calcRat=F, lang=l, add=T, cex.axis=1.2, cex.lab=1.5, rlty=c(4,2))
					#addLabel(0.975, 0.975, "(A)", cex=1.5, col="black", adj=c(1,1)) 
					if (png) dev.off()
				}; eop()
	#browser();return()
			}
		}
	
		## Spawning biomass with projections
		## --------------------
		if (fig.proj) {
			tsfld   = "B"  #"BtBmsy" ## BtBmsy for TSC report 2023
			is.biomass = tsfld %in% c("B","V")
			catpol3 = c("CC.01","CC.03","CC.07")  ## POP 2023
			cpcol   = c("green3","orange","red")
			yr.proj = dimnames(xPJ)$year
			yr.late = 2015:2023
			yr.main = setdiff(dimnames(xTS)$year, c(yr.late, yr.proj)) 
			## Adjust periods by adding a year for connection in plots
			yr.main = c(yr.main, yr.late[1])
			yr.late = c(yr.late, yr.proj[1])
	#browser();return()
			Bscale = ifelse (is.biomass, 1000., 1.)
			Qmain   = apply(xTS[,yr.main,tsfld,], c("year","area"), quantile, probs=ttcall(quants3), na.rm=T) / Bscale
			Qlate   = apply(xTS[,yr.late,tsfld,], c("year","area"), quantile, probs=ttcall(quants3), na.rm=T) / Bscale
			Qproj   = apply(xPJ[,yr.proj,tsfld,,catpol3], c("year","area","proj"), quantile, probs=ttcall(quants3), na.rm=T) / Bscale
			#BRP.qts = apply(xRP[,c("LRP","USR"),], c("area"), quantile, probs=ttcall(quants3), na.rm=T) / Bscale   ## LRP|USR not used
			areas   = dimnames(Qmain)$area
			narea   = length(areas)
			myrs    = as.numeric(yr.main)
			pyrs    = as.numeric(yr.proj)
			ayrs    = .su(c(myrs,pyrs))
			nyrs    = length(ayrs)
			## ----------------------------------------------
			## Subroutine to plot envelopes
			## Added argument elty (250805) for Karen Dwyer
			## ----------------------------------------------
			addEnvelope <- function(mat, ecol, elty=c(1,3)) {
				xx   = rownames(mat)
				yy   = colnames(mat)
				yrs  = as.numeric(yy)
				nyrs = length(yrs)
				polygon(c(yrs,rev(yrs)), c(mat["5%",yy],rev(mat["95%",yy])), col=ifelse(ecol=="black","ghostwhite",lucent(ecol,0.10)), border=FALSE)
				lines(yrs, mat["5%",yy], col="gainsboro", lwd=1, lty=elty[1]) ## just to add a bit more emphasis
				lines(yrs, mat["5%",yy], col=ecol, lwd=1, lty=elty[2])
				lines(yrs, mat["95%",yy], col="gainsboro", lwd=1, lty=elty[1]) ## just to add a bit more emphasis
				lines(yrs, mat["95%",yy], col=ecol, lwd=1, lty=elty[2])
				lines(yrs, mat["50%",yy], col=ecol, lwd=2, lty=elty[1])
				#lines(yrs, mat["50%",yy], col=ecol, lwd=ifelse(elty==1,2,1), lty=elty[1])
	#browser();return()
			}
			fout.e = paste0("FSAR.", tsfld, ".proj.subarea")
			if (diff(PIN)<0)      { orient="landscape"; rc=c(1,narea) }
			else if(diff(PIN)>=0) { orient="portrait" ; rc=c(narea,1) }
			fout.e = paste0(fout.e, ".", orient)
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				if (png) {
					clearFiles(paste0(fout,".png"))
					png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				}
				expandGraph(mfrow=rc, mar=c(2.75,3.75,0.5,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
	
				for (a in 1:narea) {
					aa = areas[a]
					plot(0,0, xlim=range(ayrs), ylim=c(0,max(Qmain[,,aa],na.rm=T)), xlab=linguaFranca("Year",l), ylab="", cex.axis=1.2, cex.lab=1.5, las=1)
					axis(1, at=seq(1935,max(years),5), labels=F, tcl=-.3)
					if ((orient=="portrait" && a==2) || (orient=="landscape")){
						if (is.biomass)
							mtext(linguaFranca(paste0("Spawning Biomass (", ifelse(Bscale==1000,"k",""), "t)"), l), side=2, line=2.2, cex=1.1)
						else {
							ylab = "Spawning Biomass italic(B)[italic(t)] relative to italic(B)[MSY]"
							ylab = linguaFranca(ylab,l)
							mtext(eval(parse(text=paste0("expression(",gsub(" ","~",ylab),")"))), side=2, line=1.75, cex=1.1)
						}
					}
					## Skip LRP|USR envelopes
					#polygon(c(years[c(1,nyrs)], years[c(nyrs,1)]), c(BRP.qts[c(1,1),"LRP"],BRP.qts[c(3,3),"LRP"]), col=lucent(col.refs[1],0.35), border=FALSE)
					#lines(years, rep(BRP.qts[2,"LRP"],length(years)), col="black", lwd=1, lty=4)
					#polygon(c(years[c(1,nyrs)], years[c(nyrs,1)]), c(BRP.qts[c(1,1),"USR"],BRP.qts[c(3,3),"USR"]), col=lucent(col.refs[2],0.35), border=FALSE)
					#lines(years, rep(BRP.qts[2,"USR"],length(years)), col="black", lwd=1, lty=5)
					addEnvelope(Qmain[,,aa], ecol="black")
					addEnvelope(Qlate[,,aa], ecol="black")
					cp = numeric()
					for (p in 1:length(catpol3)) {
						plty = switch(p, c(2,2), c(1,3), c(4,4)) ## (RH 250805) added for Karen Dwyer
						pp = catpol3[p]
						cp = c(cp, xavgCP[,aa,pp][1])
						addEnvelope(Qproj[,,aa,pp], ecol=cpcol[p], elty=plty)
					}
					## Add reference points if appropriate
					if (tsfld %in% "BtBmsy")
						abline(h=c(0.4,0.8), lty=c(4,5), col=.colBlind[c("redpurple","bluegreen")])
					cpleg = paste0(cp," t")
					#if (a==1)
						#addLegend(0.01, ifelse(orient=="portrait",0.975,0.2), col=c("black","blue"), lty=1, lwd=2, legend=paste0(c("Main","Late")," recruitment"), title=linguaFranca("Reconstruction",l), bty="n", xjust=0, cex=1.1)
					addLegend(0.75, 0.975, col=cpcol, lty=c(2,1,4), lwd=2, legend=cpleg, title=linguaFranca("Projected catch",l), bty="n", xjust=1, cex=1.2, seg.len=3.5)
	#browser();return()
					addLabel(0.95, 0.95, aa, cex=1.5, col=col.area[a], adj=c(1,1))
				}
	#browser();return()
				if (png) dev.off()
			}; eop()
	#browser();return()
		} ## end fig.proj
		rm(list=setdiff(ls(),c(keep.pars,"keep.pars")))  ## need to clear the objects in first loop so they don't mask those in second environment
		detach(ee)
	} ## end envo loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~makeFSARfigs


## plotSS.compo-------------------------2025-12-03
## Make Composite Figures
## ---------------------------------------------RH
#plotSS.compo <- function(compo, spp.code="SGR", istock="3area", 
plotSS.compo <- function(envo, 
   subset=NULL, redo.figs=FALSE, redo.panels=FALSE, 
   ptypes, lang, pngres=400, PIN=c(9,9), gmu=TRUE)
{
	vomit <- function() { gc(verbose=FALSE); while ("ee" %in% search()) detach(ee) }; #resetGraph() }
	on.exit( vomit() )

	## Subfunctions-------------------------
	## Swap areas and parameters in list
	## Adapted from https://stackoverflow.com/questions/26890615/rearrange-hierarchy-in-a-list-of-lists-in-r
	swapList <- function(x) {
		i <- 1:length(x); ii = names(x)
		j <- 1:length(x[[1]]); jj = names(x[[1]])
		swap <- lapply(j, function(jj) lapply(i, function(ii) {x[[ii]][[jj]]}))
		names(swap) = jj
		for (jjj in jj)
			names(swap[[jjj]]) = ii
#browser();return()
		return(swap)
	}
##--------------------------subfunctions

	## unpack environemntal objects : https://stackoverflow.com/questions/76057103/unpack-contents-of-an-r-environment-object-to-current-working-environment
	#list2env(as.list(envo), environment())
	for (e in 1:length(envo)) {
		keep.pars = ls()  ## need to clear environment of objects on subsequent loops through envo
		ee = envo[[e]]
		attach(ee)
		rm(list=intersect(ls(.GlobalEnv), ls(ee)), envir=.GlobalEnv)  ## remove objects in global environment masking objects in e1
		#unpackList(stock[[istock]][["Controls"]])
		#unpackList(compo)
		mess = paste0("unpackList(tcall(data.compo.figs.", istock, "))")  ## created in 'AppG_Rcode.R'
		createFdir(lang)

		modYrsChar  = as.character(modYrs)
		## Diagnostics for select parameters
		#P.names = colnames(avgPA)
		P.names = colnames(avgPA)  = gsub("\\_gp\\([0-9]\\)", "", colnames(avgPA))  ## just for M
		P.runs      = as.numeric(sapply(strsplit(rownames(avgPA),"\\."),function(x){x[1]}))
		P.run.ord   = unique(P.runs)
		P.run.nmc   = table(P.runs)[as.character(P.run.ord)]
		P.run.num   = rep(1:length(P.run.nmc), P.run.nmc)
		use.run.rwt = is.element(P.runs, P.run.ord) ## default use all runs but a subset might be used for management advice
	
		catpol3 = NULL #c("AC.00")  ## SGR 202%: 1500t; use this when CC policies have not yet been run.
		#catpol3 = c("CC.01","CC.03","CC.09")  ## SGR 2025: 0, 1500, 3000t
		#catpol3 = NULL ## comment this out when you have catch policies
		#proYrs = .su(c(currYear, proYrs))
sumting=T
if (sumting) {
		## Envelope plots
		## (RH 251001) add utumsy as an envelope plot rather than a quantile boxplot
		for (epar in c("Bt","BtB0","ututrp")[2]) {#[1:3]) {
			.flush.cat(epar, " envelope","\n", sep="")
			elab = switch(epar, 
				'Bt'     = expression(paste("Spawning biomass   ", group("(",italic(B)[italic(t)],")"))),
				'BtBmsy' = expression(paste("Stock status   ", group("(",italic(B)[italic(t)] / italic(B)[MSY],")"))),
				'BtB0'   = expression(paste("Depletion   ", group("(",italic(B)[italic(t)] / italic(B)[0],")"))),
				'ut'     = expression(paste("Harvest rate   ", group("(",italic(u)[italic(t)],")"))),
				'utumsy' = expression(paste("Exploitation status   ", group("(",italic(u)[italic(t)] / italic(u)[MSY],")"))),
				'ututrp' = expression(paste("Exploitation status   ", group("(",italic(u)[italic(t)] / italic(u)[TRP],")")))
			)
			refpts = switch(epar,
				'Bt'     = list(LRP=NULL,USR=NULL,TRP=NULL,RRR=NULL),
				'BtBmsy' = list(LRP=0.4, USR=0.8),
				'BtB0'   = list(LRP=0.16, USR=0.32, TRP=0.40),
				'ut'     = list(LRP=NULL,USR=NULL,TRP=NULL,RRR=NULL),
				'utumsy' = list(RRR=1),
				'ututrp' = list(RRR=1)
			)
			if (!exists("areaTS")) { ##|| length(areaTS)==10) {  ## should be able to handle a single-area model (RH 251028)
				catpol = avgPJ[,,epar,catpol3,drop=T]  ## YTR: see projections using catches of (0,min,max)
				plotSS.pmcmc(get(paste0(epar,".mcmc")), yrs=modYrs, pyrs=proYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=elab, outnam=paste0(prefix,"compo.", epar), xyType="envelope", ptypes=ptypes, catpol=catpol, refpts=refpts, gmu=(epar!="Bt"))  ## (RH 240912) label ref pts for GMU
			} else {
				areas = names(areaTS)
				areaPJ.swap = swapList(areaPJ)
				areaTS.swap = swapList(areaTS)
				if (is.null(catpol3) || all(catpol3=="AC.00")) {
					catpol  = NULL
					ncatpol = 1
				} else {
					catpol  = areaPJ.swap[[epar]]
					catpol  = lapply(catpol, function(x){x[,,catpol3]})  ## NA if just AC.01
					ncatpol = length(catpol3)
				}
				outnam = paste0(prefix, "compo.", epar, ".", ifelse(ncatpol==1,"AC","CC") )
				plotSS.pmcmc(obj=areaTS.swap[[paste0(epar,".mcmc")]], pqs=tcall(quants5), yrs=modYrs, pyrs=proYrs, lang=lang, cex.axis=ifelse(length(areaTS)==1,1.2,1), cex.lab=ifelse(length(areaTS)==1,1.5,1.2), yLab=elab, outnam=outnam, xyType="envelope", ptypes=ptypes, catpol=catpol, refpts=refpts, gmu==(epar!="Bt"))
			}
		} ## end epar (envelope plots)
browser();return()

		## Quantile box plots
		for (epar in c("ut","Rt")) { #,"utumsy"
			.flush.cat(epar, " quantile plot","\n", sep="")
			elab = switch(epar, 
				'ut'     = expression(paste("Exploitation rate   ", group("(",italic(u)[italic(t)],")"))),
			#	'utumsy' = expression(paste("Exploitation status   ", group("(",italic(u)[italic(t)] / italic(u)[MSY],")"))),
				'Rt'     = expression(paste("Recruitment   ", group("(",italic(R)[italic(t)],")"))) )
			refpts = switch(epar,
				'ut'     = list(LRP=NULL,USR=NULL,TRP=NULL,RRR=NULL),
			#	'utumsy' = list(RRR=1),
				'Rt'     = list(LRP=NULL,USR=NULL,TRP=NULL,RRR=NULL)
			)
			if (!exists("areaTS")) { ##|| length(areaTS)==1) {
				plotSS.pmcmc(get(paste0(epar,".mcmc")), yrs=modYrs, pyrs=proYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=elab, outnam=paste0(prefix,"compo.",epar), xyType="quantBox", ptypes=ptypes, refpts=refpts)
			} else {
				areas = sub("CST","BC",names(areaTS))
				areaPJ.swap = swapList(areaPJ)
				areaTS.swap = swapList(areaTS)
#browser();return()
				plotSS.pmcmc(areaTS.swap[[paste0(epar,".mcmc")]], yrs=modYrs, pyrs=proYrs, lang=lang, cex.axis=ifelse(length(areaTS)==1,1.2,1), cex.lab=ifelse(length(areaTS)==1,1.5,1.2), yLab=elab, outnam=paste0(prefix,"compo.",epar), xyType="quantBox", ptypes=ptypes, refpts=refpts)
			}
		}
	
		## Recruitment deviations (quantile boxes)
		.flush.cat("Rtdev boxes","\n")
		if (!is.null(dim(Rtdev.mcmc)))
			Rtdev.mcmc = list(BC=Rtdev.mcmc)
		#area.names = "BC"; narea = 1  ## can only be coastwide if multi-area model (at least for now)
#browser();return()
		plotSS.pmcmc(Rtdev.mcmc, yrs=modYrs, pyrs=proYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=expression(paste("Recruitment deviations ", group("(",delta~italic(R)[italic(t)],")"))), outnam=paste0(prefix,"compo.Rtdev"), ptypes=ptypes, refpts=list(USR=0), gmu=FALSE, yLim=c(-2.75,3.5))

		## Plot MCMC diagnostics select parameters for each component run
		## --------------------------------------------------------------
		for (i in c(1)){
			ii   = P.names[i]
			iii  = gsub("[_|[:space:]]","",ii)
			P.i  = split(avgPA[,i], P.runs)[as.character(P.run.ord)]
			  ## splits by 1:length(run.num) in 'gather.compo.case.r' so retains correct order
			P.ii = data.frame(P.i)
			#colnames(P.ii) = paste(iii,paste0("R",names(P.i)),sep="_")
			#colnames(P.ii) = paste(ii,paste0("B",names(run.num)),sep="_")  ## RH 200508 (for subsets of B)
	
			## Trace plots
			.flush.cat("Parameter trace plots","\n")
			fout.e = paste0(prefix,"compo.", iii, ".traces")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=7, horizontal=FALSE,  paper="special")
					else if (p=="png") {
						clearFiles(paste0(fout,".png"))
						png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
					}
					par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
					#colnames(P.ii) = paste(iii,paste0(switch(l, 'e'="R",'f'="E"),names(P.i)),sep="_")
					colnames(P.ii) = paste(iii,base.lab,sep=" ")
					panelTraces(mcmc=P.ii, mpd=ampdPA[,i], xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, cex.strip=1.5, same.limits=ifelse(i%in%c(1),F,T), lang=l)
					if (p %in% c("eps","png")) dev.off()
				} ## end p (ptypes) loop
			}; eop()
#browser();return()
	
			## Split chains (specify number of chains used in NUTS MCMCs) -- function found in PBSsynth
			.flush.cat("Parameter split chains","\n")
			nchains   = 8
			#col.trace = rep(rev(c("red", "blue", "black", "green", "orange", "purple", "cyan", "gold")), nchains)[1:nchains]
			col.trace = rep(rev(c("red", "blue", "black", "green", "orange", "purple", "brown", "pink")), nchains)[1:nchains]
			fout.e = paste0(prefix,"compo.", ii, ".chains")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=7, horizontal=FALSE,  paper="special")
					else if (p=="png") {
						clearFiles(paste0(fout,".png"))
						png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
					}
					par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))  ## mar and oma ignored, fixed in call to `mochaLatte'
					panelChains(mcmc=P.ii, nchains=nchains, axes=TRUE, pdisc=0, between=list(x=0, y=0), col.trace=col.trace, xlab="Parameter Value", ylab="Cumulative Frequency", cex.axis=1.2, cex.lab=1.4, cex.strip=1.5, yaxt="n", lang=l, lwd.trace=2)
					addLegend(0.975, 0.8, legend=paste0("chain ",1:nchains), lwd=2, seg.len=3, col=col.trace, bty="n", xjust=1, yjust=1)
					if (p %in% c("eps","png")) dev.off()
				} ## end p (ptypes) loop
			}; eop()
#browser();return()
	
			## ACF plots
			.flush.cat("Parameter ACFs","\n")
			fout.e = paste0(prefix,"compo.", ii, ".acfs")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(fout,".eps"), width=8, height=8, horizontal=FALSE,  paper="special")
					else if (p=="png") {
						clearFiles(paste0(fout,".png"))
						png(paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
					}
					plotACFs(P.ii, lag.max=60, lang=l, lwd=2)
					if (p %in% c("eps","png")) dev.off()
				} ## end p (ptypes) loop
			}; eop()
#browser();return()
		}


		## BtBmsy UtUmsy snail trail plot
		## ------------------------------
		.flush.cat("Stock status snail trail","\n")
		z = use.run.rwt
		if (RPbase=="BMSY") {
			BoverBrp = avgTS[z,,"BtBmsy"]
			UoverUrp = avgTS[z,,"utumsy"]
		} else if (RPbase=="B0") {
			BoverBrp = avgTS[z,,"BtB0"]
			UoverUrp = avgTS[z,,"ututrp"]
		}
		outnam    = paste0(prefix,"compo.snail")
		plotSnail(BoverBrp, UoverUrp, yrs=modYrs, p=tcall(quants3)[c(1,3)], xLim=NULL, yLim=NULL, ngear=length(gseries), assYrs=assYrs, outs=F, Cnames=gseries, ptypes=ptypes, outnam=outnam, lang=lang, labYrs=c(1950,seq(1960,1975,5), assYrs), refpts=RPbase)
	
		## Prepare base composite and components for function 'compBmsy'
		## ------------------------------------------------------------
		.flush.cat("Stock status biomass bars","\n")
		Bbase = list()
		L1    = toupper(ifelse(spp.code=="REBS", istock, spp.code))
		Bbase[[L1]] = list()
		ibase  = toupper(istock)  ## for SGR 2025, departs from previous use
		if (L1==ibase)
			ibase  = toupper(area.name)
	
		z  = use.run.rwt
		zz = sapply(split(z,P.runs),unique)[as.character(P.run.ord)]
		Ubase = sum(zz)
		runlab = base.lab[zz]
	
#browser();return()
	
		## (RH 251002) New nonsense to deal with alternate reference points for SGR
		if (!exists("RPbase")) RPbase = "BMSY"
		if (RPbase=="BMSY") {
			Btemp = data.frame(run=P.runs, BcurrBmsy=avgTS[,as.character(currYear),"BtBmsy"])[z,]
			Bcurr.Brp = Btemp[,"BcurrBmsy"]
		} else if (RPbase=="B0") {
			Btemp = data.frame(run=P.runs, BcurrB0=avgTS[,as.character(currYear),"BtB0"])[z,]
			Bcurr.Brp = Btemp[,"BcurrB0"]
		} else {
			stop ("'RPbase' specified is not supported")
		}
		Bbase[[L1]][[ibase]] = Bcurr.Brp ## composite from multiple bases
	
		if (length(P.run.ord)>1){
			for (i in 1:length(P.run.ord[zz])) {
				ii = P.run.ord[i]
				if (!is.element(ii, Btemp$run)) next
				iBtemp = Btemp[is.element(Btemp$run,ii),]
				iii = runlab[i]
				Bbase[[L1]][[paste0("R",ii)]] = iBtemp$BcurrBrp
			}
		}
		mam = FALSE
		if (exists("areaTS") && length(areaTS) > 1){  ## list object
			## Multi-area model
			mam   = TRUE
			areas = setdiff(names(areaTS), ibase) ## may be touchy
			areaTS.swap = swapList(areaTS)
			Bt.area = areaTS.swap[["Bt.mcmc"]]
			for (a in areas) {
				if (!(a %in% area.names))
					next
				else {
					aBcurr = Bt.area[[a]][,as.character(currYr)]
					if (a %in% area.names)
						aBrp  = switch( RPbase, 'BMSY'=xavgRP[,"Bmsy",a], 'B0'=xavgRP[,"B0",a] ) ## 3-D array
					#else 
					#	aBmsy  = avgRP[,"Bmsy"]  ## 2-D array
					aBcurrBrp = aBcurr/aBrp
					Bbase[[L1]][[a]] = aBcurrBrp
				}
			}
		}
		## the section above needs some serious cleaning, but works to this point
	
		L1nam = ifelse(L1=="BSR","REBS North", ifelse(L1=="RER","REBS South", L1))
		#Mnams = c(paste0(L1nam," Composite"), gsub("_| AE=3","",base.runs.lab))
		if(length(P.run.ord)==1)
			Mnams = paste0(L1nam, " Base Run")
		else
			Mnams = c(paste0(L1nam, " Composite"), runlab)
		## Add in projections -- CAUTION: for now, the code is hard-wired to have both CC and HR at year 2
		if (spp.code=="BOR") {
			pY = 2
			Bbase[[L1]][["CC-proj"]] = B.proj[["CC"]][,as.character(currYear+pY)]/Bmsy
			Bbase[[L1]][["HR-proj"]] = B.proj[["HR"]][,as.character(currYear+pY)]/Bmsy
			Mnams = c(Mnams, paste0("Comp.",c("CC ","HR "), pY,"y"))
		}
		N = length(Mnams)
	
		if (spp.code=="BOR") {
			bord=c(1:N); medcol=c(ifelse(istock==names(stock)[1],"blue","red"),rep("grey30",N-Ubase), rep("gold",2)); 
			boxfill=c(ifelse(istock==names(stock)[1],"aliceblue","mistyrose"), rep("grey95",N-Ubase), rep("lightyellow",2))
		} else {
			bord=c(1:N); medcol=c(ifelse(istock==names(stock)[1],"blue","red"),rep("grey30",ifelse(N==1,0,Ubase))); 
			boxfill=c(ifelse(istock==names(stock)[1],"aliceblue","mistyrose"), rep("grey95",ifelse(N==1,0,Ubase)))
		}
		if (mam) {  ## multi-area model
			#Mnams = c(Mnams, paste0("Subarea ",areas))  ## areas now appears to include the whole coast
			Mnams = c(Mnams, paste0("Subarea ",area.names))
			N = length(Mnams)
			bord=c(1:N)
			medcol=c("orange","blue","red","green3")  ## SGR 2025
			boxfill=c("lemonchiffon","aliceblue","mistyrose","honeydew")
		}
		
		boxlim = c(0, max(sapply(Bbase[[L1]],quantile,tcall(quants5)[5],na.rm=T)) )
		boxlim[2] = boxlim[2] + diff(boxlim) * 0.04
		#ptypes="png"; scenario=0; so("compBmsy.r")
	
		if ("win" %in% ptypes) resetGraph()
		par(mfrow=c(1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		mess =  paste0("list(",paste0(paste0(ptypes,"=TRUE"),collapse=","),")")

		out = compBmsy(Bspp=Bbase, spp=L1, boxwidth=0.5, medcol=medcol, boxfill=boxfill, boxlim=boxlim, whisklwd=2, staplelwd=3, Mnams=Mnams[bord], width=9, height=6, figgy=eval(parse(text=mess)), pngres=pngres, spplabs=F, t.yr=currYear, left.space=8, top.space=0.5, fout=paste0(prefix,"compo", ifelse(Ubase>=1,".stock.","."), "status"), calcRat=F, lang=lang, ratios=switch(RPbase, 'BMSY'=c(0.4,0.8), 'B0'=c(0.16,0.32)), refpt=switch(RPbase, 'BMSY'="MSY", 'B0'=0)  )

		save("Bbase","L1","medcol","boxfill","boxlim","Mnams","bord","mess","RPbase", file=paste0("Bbase.for.FSAR.", istock, ".rda"))

#browser();return()
} ## end sumting

		## Make quantile plots of component-run parameters and quantities
		## --------------------------------------------------------------
		if (redo.panels) {
			#P.collect = c(1:16) ## YTR: collect all parameters
			#P.collect = c(1:6,seq(16,ncol(avgPA),3)) ## SGR: collect major parameters + mu  ## only works for area3
			P.collect = grep("R0|Rdist|M\\_|BH|mu", colnames(avgPA))  ## SGR 2025
			P.pars    = data.frame( avgPA[,P.collect] )
			colnames(P.pars) = colnames(avgPA)[P.collect]

			Q.pars   = switch(RPbase,
				'BMSY' = data.frame(
				Bcurr      = avgRP[,"Bcurr"],
				B0         = avgRP[,"B0"],
				Bcurr.B0   = avgRP[,"Bcurr"]/avgRP[,"B0"],
				MSY        = avgRP[,"MSY"],
				Bmsy       = avgRP[,"Bmsy"],
				Bmsy.B0    = avgRP[,"Bmsy"]/avgRP[,"B0"],
				Ucurr      = avgRP[,"ucurr"],
				Umsy       = avgRP[,"umsy"],
				Umax       = apply(avgTS[,,"ut"],1,max,na.rm=T) ## for each mcmc sample across the time series
				#Ucurr.Umsy = avgRP[,"ucurr"]/avgRP[,"umsy"]
				),
				'B0' = data.frame(
				Bcurr      = avgRP[,"Bcurr"],
				B0         = avgRP[,"B0"],
				Bcurr.B0   = avgRP[,"Bcurr"]/avgRP[,"B0"],
				Ytrp       = avgRP[,"Ytgt"],
				TRP        = avgRP[,"Btgt"],
				#Btrp.B0    = avgRP[,"Btgt"]/avgRP[,"B0"],  ## this just evaluates to 0.4 (or it should)
				LRP       = avgRP[,"40Btgt"],  ## or some other value to fill in blank panel
				Ucurr      = avgRP[,"ucurr"],
				Utrp       = avgRP[,"utgt"],
				Umax       = apply(avgTS[,,"ut"],1,max,na.rm=T) ## for each mcmc sample across the time series
				#Ucurr.Utrp = avgRP[,"ucurr"]/avgRP[,"utgt"]
				) )
			nchains = length(P.run.ord)
	
			names(Q.pars) = sub("\\.", " / ", gsub("U","u", 
				gsub("Ucurr",paste0("U",currYear-1), 
				gsub("Bcurr",paste0("B",currYear), names(Q.pars)))))
	
			if (spp.code %in% c("YTR","SGR")){
				## Just split the one base run into 8 chains
				nchains=8
				boxfill = colorRampPalette(c("lightblue1","lightblue4"))(nchains)
			}
			if (spp.code %in% c("YMR")){
				## nchains determined above
				nchains = nchains; ngroups = 1  ## no secondary axis of uncertainty
				boxfill = c("cyan","green","coral","dodgerblue","yellow")
			} else if (spp.code %in% c("WWR")){
				nchains = 9; ngroups = nchains/3
				boxfill = paste0(rep(c("cyan","green","coral"),each=ngroups), rep(1:ngroups,3)) 
			} else if (spp.code %in% c("BOR")){
				nchains = 3; ngroups = nchains/3
				boxfill = paste0(rep(c("cyan","green","salmon"),each=ngroups), rep(1:ngroups,3)) 
			} else if (spp.code %in% c("REBS")){
				if (istock %in% c("BSR","RER")) {
					#nchains = 9; ngroups = nchains/3
					#boxfill = paste0(rep(c("cyan","green","coral"),each=ngroups), rep(1:ngroups,3)) 
					nchains = Ubase; ngroups = ceiling(as.numeric(names(run.num)[use.num])/3)    ## RH 200508 (for subsets of B)
					colrnum = as.vector(unlist(sapply(split(ngroups,ngroups), function(x){1:length(x)})))
					boxfill = paste0(c("cyan","green","coral")[ngroups], colrnum ) 
				}
			}
			xlim = range(1:nchains)+c(-0.75,0.75)
	#browser();return()
	
			fout = fout.e = paste0(prefix,"compo.pars.qbox")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				for (p in ptypes) {
					if (p=="png") {
						clearFiles(paste0(fout,".png"))
						png(paste0(fout,".png"), units="in", res=pngres, width=8, height=9)
					}
					pflds  = setdiff(colnames(P.pars), "sumting")
					#pflds  = c("q_6",pflds)
					if (spp.code %in% c("YTR","SGR"))
						xfac = paste0("c", 1:nchains)  ## for chains of YTR and SGR
					else
						xfac = paste0("B",names(run.num[use.num]))
					if (length(xfac) != nchains) {
						message("sumtingwong : length of xfac != number chains"); browser(); return()
					}
					panelBoxes(P.pars[,pflds], nchains=nchains, xlab=linguaFranca("Base Run",l), ylab=linguaFranca("Parameter estimates",l), cex.strip=1.2, cex.axis=1.1, cex.lab=1.5, outline=FALSE, xlim=xlim, boxfill=boxfill, xfac=xfac, mar=c(0,3.5,0,0), oma=c(4,2,0.5,1))  ## RH 200508 (for subsets of B)
					if (p %in% c("png","eps")) dev.off()
				}
			}; eop()
#browser();return()
	
			fout = fout.e = paste0(prefix,"compo.rfpt.qbox")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				for (p in ptypes) {
					if (p=="png") {
						clearFiles(paste0(fout,".png"))
						png(paste0(fout,".png"), units="in", res=pngres, width=8, height=8)
					}
					panelBoxes(Q.pars, nchains=nchains, xlab=linguaFranca("Base Run",l), ylab=linguaFranca("Derived Quantities",l), cex.strip=1.2, cex.axis=1.1, cex.lab=1.5, outline=FALSE, xlim=xlim, boxfill=boxfill, xfac=xfac, mar=c(0,3.8,0,0), oma=c(4,2,0.5,1))  ## RH 200508 (for subsets of B)
					if (p %in% c("png","eps")) dev.off()
				}
			}; eop()
#browser();return()
		} ## end redo.panels

		rm(list=setdiff(ls(),c(keep.pars,"keep.pars")))  ## need to clear the objects in first loop so they don't mask those in second environment
		detach(ee)
	} ## end envo loop
#browser();return()

	## ==============================================================
	## Prepare multiple stock base composites for function 'compBmsy'
	## If comparing two stocks or models (e.g., RSR,POP or SGR1 and SGR2)
	if (length(stock)>1) {
		for (i in 1:length(stolab)) {
			istock = stolab[i]
			load(paste0("Bbase.for.FSAR.", istock, ".rda"))
			if (i==1) {
				Bcompo = Bbase
				Mnamso = sub(" Run", paste0(": ",istock), Mnams)
				mdcolo = medcol
				bxfilo = boxfill
			} else {
				Bcompo[[1]][[names(Bbase[[1]])]] = Bbase[[1]][[names(Bbase[[1]])]]
				Mnamso = c(Mnamso, sub(" Run", paste0(": ",istock), Mnams))
				mdcolo = c(mdcolo, medcol)
				bxfilo = c(bxfilo, boxfill)
			}
		}
		## Very much geared to SGR 2025 (will need to change for future species)
		z = c(grep("COAST",names(Bcompo[[1]])), grep("COAST",names(Bcompo[[1]]),invert=TRUE))  ## Put coast first
		Bcompo[[1]] = Bcompo[[1]][z]
		Mnamso = Mnamso[z]
		mdcolo = mdcolo[z]; mdcolo[1] = "purple"
		bxfilo = bxfilo[z]; bxfilo[1] = "thistle1"

		#outnam = paste0("status.",paste0(paste0("(",paste(tolower(stospp),stolab,sep="."),")"),collapse=".vs."))
		outnam = paste0(stospp[1],".status.",paste0(paste0("(",stolab,")"),collapse=".vs."))

#browser();return()
		spp.code = unique(stospp)[1]  ## needs to be one species
		out    = compBmsy(Bspp=Bcompo, spp=spp.code, boxwidth=0.5, medcol=mdcolo, boxfill=bxfilo, whisklwd=2, staplelwd=2, Mnams=Mnamso, width=9, height=4, figgy=eval(parse(text=mess)), pngres=pngres, spplabs=F, t.yr=e1$currYear, left.space=c(9,9), top.space=1, fout=outnam, calcRat=F, lang=lang, ratios=switch(RPbase, 'BMSY'=c(0.4,0.8), 'B0'=c(0.16,0.32)), refpt=switch(RPbase, 'BMSY'="MSY", 'B0'=0))
	} ## end more than 1 stock
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.compo


## plotSS.senso-------------------------2025-12-04
## Make Sensitivity Figures
## ---------------------------------------------RH
plotSS.senso <- function(envo, #senso, spp.code="SGR", istock="SGR",
	subset=NULL, redo.figs=FALSE, redo.panels=FALSE,
	ptypes, lang, pngres=400, PIN=c(9,9),
	useRlow=FALSE)
{
	vomit <- function() { gc(verbose=FALSE); while ("ee" %in% search()) detach(ee) }; #resetGraph() }
	on.exit( vomit() )

	for (e in 1:length(envo)) {
		keep.pars = ls()  ## need to clear environment of objects on subsequent loops through envo
#browser();return()
		ee = envo[[e]]
		attach(ee)
		rm(list=intersect(ls(.GlobalEnv), ls(ee)), envir=.GlobalEnv)  ## remove objects in global environment masking objects in e1
		#unpackList(stock[[istock]][["Controls"]])
		#unpackList(senso)

		sen.fig     = "senso."
		Nsenso      = length(sen.run.num)
		all.run.num = c(exp.run.num, sen.run.num)
		## Provide 15 colours and line types
		col.senso   = c("green4", "blue", "red", "purple", "orange", "skyblue", "gold3", "darkolivegreen", "lightseagreen", "hotpink", "brown", "cyan", "tomato", "chartreuse3", "darkorchid1") # 
		lty.senso   = c("22", "44", "66", "13", "73", "1551", "1343", "2262", "3573", "654321", "12345678", "88", "17", "5937", "9876")
		col.senso   = rep(col.senso, Nsenso)[1:Nsenso]
		lty.senso   = rep(lty.senso, Nsenso)[1:Nsenso]
		names(col.senso) = names(lty.senso) = sen.run.num

		## Need to agglomerate parameters due to fleet offsets
		if (spp.code=="POP") {
			fleets = c("(TRAWL_5ABC)","QCS","WCVI","WCHG","NMFS")
			good   = c("^(M_)Female", "^(M_)Male", paste0("mu.+",fleets), paste0("varL.+",fleets), paste0("delta1.+",fleets))
			bad    = c("DM\\_theta")
			idx    = 2
		}
		if (spp.code=="YTR") {
			fleets = c("(TRAWL|BT_BC)","QCS","WCVI","NMFS")
			good   = c("^(M_|M2_)Female", "^(M_|M2_)Male", paste0("mu.+",fleets), paste0("varL.+",fleets), paste0("delta1.+",fleets))
			bad    = c("^M1\\_","sigmaR","MW\\_BC","beta6","varR","delta[2-9]","DM\\_theta")
			idx    = 3
		}
		if (spp.code=="SGR") {
			fleets = c("(TRAWL_5ABC)", "(TRAWL_5DE)", "(TRAWL_3CD)", "QCS", "WCHG", "WCVI")
			good   = c("^(M_|M1_|M2_)Female", "^(M_|M1_|M2_)Male", "BH\\_h", "LN\\(R0\\)", "Rdist", paste0("mu.+",fleets), paste0("varL.+",fleets), paste0("Delta.+",fleets), paste0("LnQ.+",fleets))
			bad    = c("HBLL","DM\\_theta")
			idx    = 3   ## number of strings separated by "." in rownames of PA
		}
		## This is a mess I will have to sort later
#print(istock)
#if (istock=="3area") {browser();return()}
		PA.mpd  = mergePA(smpdPA, good=good, bad=bad, idx=idx)
#browser();return()
		PA.mcmc = mergePA(senPA, good=good, bad=bad, idx=idx)
		#unpackList(PA.list, scope="L")
		S.mpd = PA.mpd[["PA"]]
		senPA = PA.mcmc[["PA"]]

		## Check to see if sensitivities are alternative runs (e.g., single-are models)
		is.penso = as.character(substitute(senso)) == "penso"
		if (is.penso) {
			sen.run.num = pen.run.num
			sen.rwt.num = pen.rwt.num
			sen.run.rwt = pen.run.rwt
			sen.lab     = pen.lab
			sen.long    = pen.long
			sen.fig     = "penso."
			col.senso   = c("blue", "green4", "red")  ## may need more when penso changes
		}
		#unpackList(tcall(data.compo.figs))  ## created in Rnw
		createFdir(lang)
	
		#modYrs      = startYear:currYear
		modYrsChar  = as.character(modYrs)
		## Diagnostics for select parameters
		P.names     = colnames(senPA) ## parameter names
		#S.mpd       = senPA #smpdPA; storage.mode(S.mpd)="double"     ## For some reason, this matrix is stored as a list (maybe due to NA values?)
		S.runs      = as.numeric(sapply(strsplit(rownames(senPA),"\\."),function(x){x[1]})) ## vector of Runs
		#S.run.rwt   = sapply(strsplit(rownames(senPA),"\\."),function(x){paste0(x[1],".",x[2],".",x[3])}) ## vector of Runs and Rwts and ver
		S.run.rwt   = sapply(strsplit(rownames(senPA),"\\."),function(x,i){
			paste0(x[1:i], collapse=".")
		}, i=ifelse(spp.code=="POP",2,3)) ## vector of Runs and Rwts and ver ## but POP does not have ver
		S.run.ord   = unique(S.runs)                           ## unique Run numbers (1st is central run)
		S.run.nmc   = table(S.runs)[as.character(S.run.ord)]   ## vector of sensitivity runs ordered
	
		S.num = match(S.run.ord[-1],sen.run.num)  ## Sensitivity numbers according to control file
		S.run.num = rep(c(0,S.num), S.run.nmc)  ## 1st run is Central Run (include in gatherMCMC senso)
		#S.run.num   = rep(1:length(S.run.nmc), S.run.nmc) - 1  ## 1st run is Central Run (include in gatherMCMC senso)
		CRS.num      = unique(S.run.num)
	#browser();return()
		## Look into subsetting later
		use.run.rwt = unique(S.run.rwt)
		#is.element(P.runs, P.run.ord) ## default use all runs but a subset might be used for management advice
	
		## Create Sensitivity labels
		S.prefix    = paste0("S",pad0(CRS.num,2)," (R", pad0(S.run.ord,2) ,") ")
		if (is.penso) {
			S.prefix    = paste0("A",pad0(CRS.num,1)," (R", pad0(S.run.ord,2) ,") ")
		}
		iCR         = grep( sub("\\.v[0-9]+[a-z]", "", exp.run.rwt), unique(B.index) )
		S.prefix[1] = paste0("B",iCR, " (R", strsplit(exp.run.rwt,"\\.")[[1]][1], ") ")
		S.labels    = paste0(S.prefix, gsub("\\_"," ", c("Base Run",sen.lab[S.num])))  ## Central Run if using a composite
		S.labels    = sub("\\\\pc", "%", S.labels)  ## just in case  wtf?  maybe reverse the subsitution (leave for now)
	#browser();return()
	
		## Function 'calcQs' now available in 'util.Funs.r'
		P.qnts = calcQs(senPA, ivec=S.run.rwt, ovec=S.run.ord)
		B.qnts = calcQs(senTS[,,"Bt"], ivec=S.run.rwt, ovec=S.run.ord)
		D.qnts = calcQs(senTS[,,"BtB0"], ivec=S.run.rwt, ovec=S.run.ord)
		U.qnts = calcQs(senTS[,,"ut"], ivec=S.run.rwt, ovec=S.run.ord)
		R.qnts = calcQs(senTS[,,"Rt"], ivec=S.run.rwt, ovec=S.run.ord)
		RD.qnts = calcQs(senTS[,,"Rtdev"], ivec=S.run.rwt, ovec=S.run.ord)
	#browser();return()
	
		L1 = if (spp.code %in% c("REBS")) istock else spp.code  ## RH 200416

sumting=F
if (sumting) {

		## Diagnostics for select parameters
		## ---------------------------------
		if (redo.figs) {
			for (i in grep("LN\\(R0\\)", P.names)) {
				ii   = P.names[i]
				iii  = gsub("[_|[:space:]]","",ii)
				P.i  = split(senPA[,i], S.runs)[as.character(S.run.ord)]
				## splits by 1:length(run.num) in 'gather.compo.case.r' so retains correct order
				P.ii = data.frame(P.i)
	#browser();return()
	
				#colnames.e = paste0(ii, ": ", S.labels)
				colnames.e =  S.labels
				colnames.f = linguaFranca(colnames.e, "f")
				fout.e = paste0(prefix, sen.fig, ii, ".traces")
				for (l in lang) {
					changeLangOpts(L=l)
					#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
					fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
					colnames(P.ii) = switch(l, 'e'=colnames.e, 'f'=colnames.f)
					for (p in ptypes) {
						if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=7, horizontal=FALSE,  paper="special")
						else if (p=="png") {
							clearFiles(paste0(fout,".png"))
							png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
						}
						panelTraces(mcmc=P.ii, mpd=S.mpd[,i], xlab="Samples", ylab=paste0("Parameter value ",ii), cex.axis=1.2, cex.lab=1.5, cex.strip=1.0, same.limits=ifelse(ii%in%c("LN(R0)"),F,T), lang=l, mar=c(0,2.5,0,0), oma=c(3.2,2,0.5,0.5), rc=.findSquare(ncol(P.ii)))
						if (p %in% c("eps","png")) dev.off()
					} ## end p (ptypes) loop
				}; eop()
#browser();return()
	
				nchains   = 8
				col.trace = rep(rev(c("red", "blue", "black", "green", "orange", "purple", "brown", "pink")), nchains)[1:nchains]
				fout.e = paste0(prefix, sen.fig, ii, ".chains")
				for (l in lang) {  
					changeLangOpts(L=l)
					fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
					#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
					colnames(P.ii) = switch(l, 'e'=colnames.e, 'f'=colnames.f)
					for (p in ptypes) {
						if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=7, horizontal=FALSE,  paper="special")
						else if (p=="png") {
							clearFiles(paste0(fout,".png"))
							png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
						}
						panelChains(mcmc=P.ii, nchains=nchains, axes=TRUE, pdisc=0, between=list(x=0, y=0), col.trace=col.trace, xlab=paste0("Parameter value ",ii), ylab="Cumulative Frequency", cex.axis=1.2, cex.lab=1.4, cex.strip=1.2, lang=l, mar=c(2,0,0,0), oma=c(3,4,0.5,0.5), rc=.findSquare(ncol(P.ii))) #, yaxt="n"
						if (p %in% c("eps","png")) dev.off()
					} ## end p (ptypes) loop
				}; eop()
#browser();return()
	
				fout.e = paste0(prefix, sen.fig, ii, ".acfs")
				for (l in lang) {  
					changeLangOpts(L=l)
					#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
					fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
					colnames(P.ii) = switch(l, 'e'=colnames.e, 'f'=colnames.f)
					for (p in ptypes) {
						if (p=="eps") postscript(paste0(fout,".eps"), width=8, height=8, horizontal=FALSE,  paper="special")
						else if (p=="png") {
							clearFiles(paste0(fout,".png"))
							png(paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
						}
						plotACFs(P.ii, lag.max=60, lang=l, rc=.findSquare(ncol(P.ii)))
						if (p %in% c("eps","png")) dev.off()
					} ## end p (ptypes) loop
				}; eop()
			}
		}
#browser();return()
} ## sumting

sumtingmore=T
if (sumtingmore) {
	
		## Create Rhat figures for the worst sensitivities
		##"C:/Users/haighr/Files/GFish/PSARC24/YTR/Data/SS3/YTR2024/Run02/MCMC.02.01.v1c"
		cwd = getwd()
		on.exit(setwd(cwd), add = TRUE)
		if (spp.code=="YTR") {
			use.base = FALSE
			sens.dir = paste0(base.dir, "/Run", pad0(unique(S.runs),2), "/MCMC.", unique(S.run.rwt))[-1]  ## remove base run
			#sens.use = c(1,2,5,6) ## YTR: split-M, dome sel, est sigmaR, use D-M
			sens.use = c(1,2,5,7) ## YTR: split-M, dome sel, est sigmaR, remove AE
			sens.run = S.run.ord[-1]
		}
		if (spp.code=="SGR") {
			use.base = TRUE
			sens.dir = paste0(base.dir, "/Run", pad0(unique(S.runs),2), "/MCMC.", unique(S.run.rwt))  ## keep base run
			sens.use = switch(istock,
				'3area'=c(8,11,13,4) + 1,  ## SGR: base, rm catch, rm GIG/NMFS, fix 5DE Rdist (251204)
				'coast'=c(10,9,5,1) + 1 ) ## SGR: base, rm catch, rm GIG/NMFS, fix 5DE Rdist (251204)
			exlax = paste0("S",pad0(sens.use-1,2), " :  ")
			exlax = gsub("S00", "B1", exlax)
			sens.run = S.run.ord
		}
		for ( i in 1:length(sens.use)) {
			ii = sens.dir[sens.use[i]]
			iii = paste0("s", pad0(sens.use[i]-use.base,2),"r",pad0(sens.run[sens.use[i]],2))
			fout.e = paste0(prefix, sen.fig, "rhat.", iii)
.flush.cat(iii, "\n")
			for (l in lang) {
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(fout,".eps"), width=10, height=6, horizontal=FALSE,  paper="special")
					else if (p=="png"){
						clearFiles(paste0(fout,".png"))
						png(paste0(fout,".png"), units="in", res=pngres, width=10, height=8)
					}
					calcRhat(dir=ii, recdevs=F, rhat.only=T, lang=l, xlim=c(0.99,1.05), exlax=exlax[i])#, barcols=barcols)
					if (p %in% c("png","eps")) dev.off()
				}
			}
#browser();return()
		}
#browser();return()
} ## sumtingmore

#} ## sumting

		## Make quantile plots of component-run parameters and quantities
		## --------------------------------------------------------------
		## Note: When outliers are plotted, these quantile plots take forever when sent to a png file.
		if (redo.panels) {
			P.collect   = switch(spp.code,
				'YMR'=c(1:5),
				'CAR'=1:6,
				'POP'=c(1,5:8,11,14,17,20),
				'YTR'=1:8,
				'SGR'=switch(istock, '3area'=c(1:6,8:13), 'coast'=c(1:4,5:9))
			) 
#browser();return()
			if (is.penso){
				## Fiddly amalgamation of parameters from single-area models for PJS
				if (spp.code=="POP") {
					P.collect = grep("^M_|BH_h|^mu", colnames(senPA)) ## collect M and mu
					P.cent.sens = data.frame( senPA[,P.collect] )
					colnames(P.cent.sens) = colnames(senPA)[P.collect]
					#P.exclude = grep("3CD$|5DE$",colnames(P.cent.sens),invert=TRUE)  ## get rid of 3CD and 5DE
					#P.cent.sens = P.cent.sens[,P.exclude]
					P.penso = as.data.frame(array(NA, dim=c(nrow(P.cent.sens),9), dimnames=list(sample=rownames(P.cent.sens), par=c("M_Female","M_Male","BH_h","mu_5ABC|3CD|5DE","mu_QCS","mu_WCVI","mu_WCHG","mu_GIG","mu_NMFS"))))
					P.penso[,1:3] = P.cent.sens[,1:3] ## M, and hC in same place across models
					for (i in colnames(P.penso)[4:9]) {
						ii  = sub("^mu_","",i)
						iii = grep(ii, colnames(P.cent.sens))
	#browser();return()
						P.penso[,i] = apply(P.cent.sens[,iii,drop=F], 1, sum, na.rm=T)
					}
					## Get rid of artifical zeroes created by appply 
					is.zero = P.penso==0 & !is.na(P.penso)
					P.penso[is.zero] = NA
	
				}
			} else {
				P.cent.sens = data.frame( senPA[,P.collect] )
				colnames(P.cent.sens) = colnames(senPA)[P.collect]
			}
			## dimnames(senRP)[2]
			## $val
 			##  [1] "Bcurr"  "B0"     "20B0"   "40B0"   "MSY"    "Bmsy"   "40Bmsy" "80Bmsy" "Ytgt"    "Btgt"    "80Btgt"  "40Btgt" 
			## [13] "ucurr"  "Fcurr"  "umsy"   "Fmsy"   "utgt"    "Ftgt"
#browser();return()
			useRP = apply(senRP,2,function(x){!all(is.na(x))})
			Q.cent.sens   = switch(RPbase,
				'BMSY' = data.frame(
					Bcurr      = senRP[,"Bcurr"],
					B0         = senRP[,"B0"],
					Bcurr.B0   = senRP[,"Bcurr"]/senRP[,"B0"],
					MSY        = senRP[,"MSY"],
					Bmsy       = senRP[,"Bmsy"],
					Bmsy.B0    = senRP[,"Bmsy"]/senRP[,"B0"],
					Ucurr      = senRP[,"ucurr"],
					Umsy       = senRP[,"umsy"],
					Umax       = apply(senTS[,,"ut"],1,max,na.rm=T) ## for each mcmc sample across the time series
					#Ucurr.Umsy = senRP[,"ucurr"]/senRP[,"umsy"]
				),
				'B0' = data.frame(
					Bcurr      = senRP[,"Bcurr"],
					B0         = senRP[,"B0"],
					Bcurr.B0   = senRP[,"Bcurr"]/senRP[,"B0"],
					Btrp       =
						if (useRP["Btgt"]) {
							senRP[,"Btgt"]
						} else {
							0.4 * senRP[,"B0"]
						},
					Ucurr      = senRP[,"ucurr"],
					Umax       = apply(senTS[,,"ut"],1,max,na.rm=T)     ## for each mcmc sample across the time series
				)
			)  ## end switch RPbase
			if (RPbase=="B0") {
				if (useRP["Ytgt"]) Q.cent.sens = data.frame(Q.cent.sens, Ytrp = senRP[,"Ytgt"])   ## SGR 2025: not common to all runs
				if (useRP["utgt"]) Q.cent.sens = data.frame(Q.cent.sens, Utrp = senRP[,"utgt"], Ucurr.Utrp = senRP[,"ucurr"]/senRP[,"utgt"])   ## SGR 2025: not common to all runs
			}
#browser();return()

			nchains = length(S.run.ord)
	
			names(Q.cent.sens) = sub("\\.", " / ", gsub("U","u", 
				gsub("Ucurr",paste0("U",currYear-1), 
				gsub("Bcurr",paste0("B",currYear), names(Q.cent.sens)))))
	
			if (spp.code %in% c("BOR","WWR")) {  ## this code seems to be old so will change without guarantee of compatibility
				P.pars = senPA[,c("R_0","h","q_1","mu_1","q_7","mu_7")]
			} else if (spp.code %in% c("REBS")) {
				if (istock=="BSR")
					P.pars = senPA[,c("R_0","mu_3","q_1","mu_1","q_2","mu_2")]
				if (istock=="RER")
					P.pars = senPA[,c("R_0","q_1","q_2","q_3","mu_4","mu_5")]
			} else if (is.penso) {
				P.pars = P.penso
			} else {
				P.pars  = P.cent.sens
			}
#browser();return()
	
			## Sometimes want to exclude sensitivities  ## RH 200416
			verboten = NULL
			if (spp.code=="BOR") {
				verboten = "no_CPUE"
				zap.runs = match(verboten,sen.lab)
				#zap.runs = grep("9",c(sen.run.nums))
			}
			if (spp.code=="YMR") {
				verboten = c("estimate_M", "start_Rdevs_in_1970")#, "steepness_h=0.5")
				zap.runs = match(verboten,sen.lab) ## + 1  ## to account for central run
			}
			if (!is.null(verboten)) {
				## Not really robust to mcsub choices
				nmcmc    = sen.mcsubs[zap.runs]
				bad.rows = lapply(1:length(zap.runs),function(i){ (nmcmc[[i]][1]:nmcmc[[i]][2]) + (zap.runs[i] * nmcmc[[i]][2]) })
				zap.rows = unlist(bad.rows)
				#P.pars[zap.rows,] = NA
				fn.ylim <- function(x, yzero=T){
					xr = extendrange(sapply(split(x,names(x)), quantile,quants5[c(1,5)], na.rm=TRUE))
					if (yzero) xr[1] = 0
					return(xr)
				}
			} else {
				## Sometimes the 95pc limits are outrageously high so use this function:
				##   hardwire index i for now
				##   hardwire yzero because passing yzero to yzero is problematic (promise already under evaluation)
				fn.ylim <- function(x, i=rep(1:(length(S.num)+1),each=2000), yzero=T){
					xr = range(x, na.rm=TRUE)
					xx = split(x,i);
#browser();return()
					if (yzero)
						xr[1] = 0
					else
						xr[1] = quantile(x,0.025,na.rm=T)
					xr[2] = quantile(x,0.975,na.rm=T)
					return(xr)
				}
			}
#browser();return()

			fout.e = paste0(prefix, sen.fig, "pars.qbox")
			for (l in lang) {  
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				LP.pars = P.pars
				names(LP.pars) = linguaFranca(names(LP.pars),l)
				for (p in ptypes) {
					if (p=="png") {
						clearFiles(paste0(fout,".png"))
						png(filename=paste0(fout,".png"), units="in", res=pngres, width=8, height=8)
					}
					#boxfill = c("gainsboro","green4","blue","red","purple","orange","skyblue","gold3","salmon3","red2","hotpink","darkred","dodgerblue") ## colours used in the trajectory plots (REBS)
					boxfill = c("gainsboro",col.senso[S.num]) ## colours used in the trajectory plots (YMR)
					#panelBoxes(P.pars, xlim=c(0.25,nchains+0.75), xlab=linguaFranca("Sensitivity Runs",l), ylab=linguaFranca("Parameter estimates",l), nchains=nchains, xfac=c(ifelse(NrefM>1,"CR","B1"), paste0("S",pad0(1:(nchains-1),2))), boxfill=boxfill, cex.strip=1.2, cex.axis=1.2, cex.lab=1.5, outline=FALSE)
					panelBoxes(LP.pars, xlim=c(0.25,nchains+0.75), xlab="", ylab=linguaFranca("Parameter Estimates",l), nchains=nchains, xfac=c(ifelse(NrefM>1,"CR","B1"), paste0(ifelse(is.penso,"A","S"), pad0(S.num,ifelse(is.penso,1,2))) ), boxfill=boxfill, cex.strip=1.2, cex.axis=1.2, cex.lab=1.5, outline=FALSE, rc=.findSquare(ncol(LP.pars)), fn.ylim=fn.ylim)
					mtext(linguaFranca(ifelse(is.penso,"Multi-area (B) vs. Single-area (A) Runs","Sensitivity Runs"),l), side=1, outer=T, line=2.5, cex=1.2, adj=0.55)
					if (p %in% c("png","eps")) dev.off()
				}
			}; eop()

			Q.pars = Q.cent.sens
			if (!is.null(verboten))
				Q.pars[zap.rows,] = NA
			fout.e = paste0(prefix, sen.fig, "rfpt.qbox")
			for (l in lang) {  
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				LQ.pars = Q.pars
				names(LQ.pars) = linguaFranca(names(LQ.pars),l)
				for (p in ptypes) {
					if (p=="png") {
						clearFiles(paste0(fout,".png"))
						png(filename=paste0(fout,".png"), units="in", res=pngres, width=8, height=8)
					}
					#boxfill = c("gainsboro","green4","blue","red","purple","orange","skyblue","gold3","salmon3","red2","hotpink","darkred","dodgerblue") ## colours used in the trajectory plots (REBS)
					boxfill = c("gainsboro", col.senso[S.num]) ## colours used in the trajectory plots (YMR)
					panelBoxes(LQ.pars, xlim=c(0.25,nchains+0.75), xlab="", ylab=linguaFranca("Derived Quantities",l), nchains=nchains, xfac=c(ifelse(NrefM>1,"CR","B1"), paste0(ifelse(is.penso,"A","S"), pad0(S.num,ifelse(is.penso,1,2))) ), boxfill=boxfill, cex.strip=1.2, cex.axis=1.1, cex.lab=1.5, outline=FALSE, fn.ylim=fn.ylim)
					mtext(linguaFranca(ifelse(is.penso,"Multi-area (B) vs. Single-area (A) Runs","Sensitivity Runs"),l), side=1, outer=T, line=2.5, cex=1.2, adj=0.55)
					if (p %in% c("png","eps")) dev.off()
				}
			}; eop()
		}
#browser();return()

		## Make plots of median trajectories
		## ---------------------------------
		ii   = as.character(startYear:currYear)
		bb   = list('Bt'=B.qnts, 'BtB0'=D.qnts, 'U'=U.qnts, 'R'=R.qnts, 'RD'=RD.qnts)
		bdat = lapply(bb,function(x){
			sapply (x, function(y){
				#if (!any(grepl("50", rownames(y)))) {browser();return()}
				yyy = rep(NA,length(ii)); names(yyy)=ii
				yy  = y[grep("50%", rownames(y)),]
				yyy[names(yy)] = yy
	#browser();return()
				return(yyy)
			})
		})
		bdat$label = S.labels
	
		tlty = c("solid", lty.senso[S.num])
		tcol = c("black",col.senso[S.num])
	
		#so("plotTraj.r","awatea")
		if (redo.figs) {
			traj.names = c("Bt","BtB0","U","R","RD")##[1:2] for les figs de french
			Ntraj = length(traj.names)
	
			for (k in 1:Ntraj){
				kk = traj.names[k]
				logR = FALSE
				traj.meds = bdat[[kk]]
				if (is.null(traj.meds)) next
				fout.e = paste0(prefix, sen.fig, "traj.",kk)
				traj = traj.meds[intersect(rownames(traj.meds),ii),]
				Nruns = ncol(traj)
				Nyrs  = nrow(traj)
				x     = as.numeric(rownames(traj))
				xlim  = range(x)
				ylim = range(traj,na.rm=TRUE)
				if (spp.code %in% c("SGR","YTR","POP")) {
					y0=ifelse(kk=="RD",F,T);   if (y0) ylim[1] = 0; if (kk=="RD") ylim[1]=-2; if (kk=="BtB0") ylim[1]=0;# if (kk=="Bt") ylim[2]=21000 #1.1*ylim[2]
					logR=F; if (kk=="R" && logR) ylim[1]=1
				} else if (spp.code=="CAR") {
					y0=ifelse(kk=="RD",F,T);   if (y0) ylim[1] = 0; if (kk=="RD") ylim[1]=-0.5; if (kk=="BtB0") ylim[1]=0; if (kk=="Bt") ylim[2]=21000 #1.1*ylim[2]
				} else if (spp.code=="YMR") {
					y0=ifelse(kk=="RD",F,T);   if (y0) ylim[1] = 0; if (kk=="RD") ylim[1]=-3; if (kk=="BtB0") ylim[1]=-0.2; if (kk=="Bt") ylim[2]=50000 #1.1*ylim[2]
				}
				## tcol and legtxt need to be reset inside language loop because they are altered when bdat is supplied
				tcol   = rep(tcol,Nruns)[1:Nruns]
				tlty   = rep(tlty,Nruns)[1:Nruns]
				legtxt = bdat$label
				for (l in lang) {  
					changeLangOpts(L=l)
					#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
					fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
					for (p in ptypes) {
						if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
						else if (p=="png") {
							clearFiles(paste0(fout,".png"))
							png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
						}
						par(mfrow=c(1,1), mar=c(3,3.5,1,0), oma=c(0,0,0,1), mgp=c(2,0.5,0))
						plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="",log=ifelse(kk=="R"&&logR,"y",""))
						if (kk %in% c("BtB0"))
							abline(h=c(0.16,0.32,0.4), lty=5, col="grey20") #col=c("salmon","darkorchid","navy"))
						if (kk %in% c("RD"))
							abline(h=c(0), lty=5, col="grey20") #col=c("salmon","darkorchid","navy"))
						for (j in Nruns:1) {
	#.flush.cat(j," ",jj," ", tlty[j], "\n")
							jj = colnames(traj)[j]
							y  = traj[,jj]
							lines(x,y, col=tcol[j], lwd=ifelse(j==1,3,2), lty=tlty[j])
						}
						#lines(x=as.numeric(names(bline)), y=if (jj=="R" && logR) log10(bline) else bline, col="black", lty=1, lwd=3)
						mtext(linguaFranca("Year",l), side=1, line=1.75, cex=ifelse(Nruns==1,1.5,1.5)) ## wtf?
						ylab = switch(kk, 'Bt'="Spawning Biomass", 'BtB0'="Spawning Biomass Depletion", 'B'="Spawning Biomass", 'VB'="Vulnerable Biomass", 'R'="Recruitment", 'RD'="Recruitment Deviation", 'U'="Exploitation Rate", "Unknown")
						#mtext(linguaFranca(ylabs[[jj]],l), side=2, line=1.8, cex=ifelse(Ntraj==1,1.5,1.2))
						mtext(linguaFranca(ylab,l), side=2, line=1.8, cex=ifelse(Ntraj==1,1.5,1.5))
						legtxt = gsub("_"," ",legtxt)
						if (spp.code %in% c("YTR","SGR")) {
							if (kk %in% c("sumting")) {
								xleg=0.99; yleg=0.99; xjust=1; yjust=1
							} else if (kk %in% c("BtB0")){
								xleg=0.04; yleg=0.04; xjust=0; yjust=0
							} else if (kk %in% c("U","R","RD")){
								xleg=0.04; yleg=0.98; xjust=0; yjust=1
							} else if (kk %in% c("Bt")){
								xleg=0.04; yleg=switch(spp.code, 'YTR'=0.855, 'SGR'=switch(istock, '3area'=0.95, 'coast'=0.4)); xjust=0; yjust=1
							}
						} else if (spp.code %in% c("POP")) {
							if (kk %in% c("Bt","BtB0", "U","R","RD")) {
								xleg=0.99; yleg=0.99; xjust=1; yjust=1
							} else {
								xleg=0.02; yleg=0.01; xjust=0; yjust=0
							}
						}
						addLegend(xleg, yleg, col=tcol, seg.len=5, legend=linguaFranca(legtxt,l), bty="o", box.col="grey", bg=ifelse(kk%in%c("BtB0"),"white","transparent"), xjust=xjust, yjust=yjust, lwd=2, lty=tlty, cex=0.8)
						#if (Ntraj==1) {
							axis(1,at=intersect(seq(1900,2500,5),x),labels=FALSE,tcl=-0.2)
							axis(1,at=intersect(seq(1900,2500,10),x),labels=FALSE)
						#}
						box()
						if (p %in% c("png","eps")) dev.off()
					}
				}; eop()
#browser();return()
			}  ## end k loop (Ntraj)
		}  ## end if redo.figs

		## Prepare sensitivity runs for function 'compBmsy' (stock status)
		## ---------------------------------------------------------------
		if (redo.figs) {
			Bsens = list()
			L1    = toupper(ifelse(spp.code=="REBS", istock, spp.code))
			Bsens[[L1]] = list()
			isens  = toupper(istock)
	
			z   = use.run.rwt  ## look into subsetting later -- perhaps just subset z
			zz  = is.element(use.run.rwt, z)
			zzz = is.element(S.run.rwt, z)
			#zz = sapply(split(z,S.runs),unique)[as.character(S.run.ord)]
			Ubase = sum(zz)
			runlab = S.labels[zz] ## created on L33

			## (RH 251002) New nonsense to deal with alternate reference points for SGR
			if (!exists("RPbase")) RPbase = "BMSY"
			if (RPbase=="BMSY") {
				Stemp = data.frame(run=S.runs[zzz], BcurrBrp=senTS[,as.character(currYear),"BtBmsy"])[zzz,]
			} else if (RPbase=="B0") {
				Stemp = data.frame(run=S.runs[zzz], BcurrBrp=senTS[,as.character(currYear),"BtB0"])[zzz,]
			} else {
				stop ("'RPbase' specified is not supported")
			}
			Bcurr.Brp = Stemp[is.element(Stemp$run, central.run),"BcurrBrp"]
			medCR = median(Bcurr.Brp, na.rm=T)
			Bsens[[L1]][[runlab[1]]] = Bcurr.Brp ## central run from composite for snesitivity comparison
			for (i in 2:length(S.run.ord[zz])) {
				ii = S.run.ord[i]
				if (!is.element(ii, Stemp$run)) next
				iStemp = Stemp[is.element(Stemp$run,ii),]
				iii = runlab[i]
				Bsens[[L1]][[iii]] = iStemp$BcurrBrp
			}
			Mnams = runlab
			## Check for restricting MCMCs to low recruits
			if (useRlow) {
				R.cent = stock[[istock]]$Base[[exp.run.rwt]]$currentMCMC$R
				if (spp.code=="BOR")
					Rlow   = R.cent[,"2017"] < quantile(R.cent[,"2017"], qRlow)
				else {
					R.mean = apply(R.cent,1,mean)  ## mean recruitment 1935-2020
					Rlow   = R.mean < quantile(R.mean, qRlow)
				}
				Mnams = c(Mnams, "low recruitment")
			}
			N = length(Mnams)
			bord=c(1:N); medcol=c(ifelse(istock==names(stock)[1],"blue","red"),rep("grey30",N-1)); 
			boxfill = c(ifelse(istock==names(stock)[1],"aliceblue","mistyrose"), rep("grey95",N-1))
			if (is.penso) {
				medcol = c("slategray", col.senso)
				boxfill = c("grey95", "aliceblue","honeydew","mistyrose")
			}
	
	#ptypes="png"; scenario=0; so("compBmsy.r")
			if ("win" %in% ptypes) resetGraph()
			mess =  paste0("list(",paste0(paste0(ptypes,"=TRUE"),collapse=","),")")
	
			Bspp = Bsens
			#Bspp = Bsens["BtBmsy"]; names(Bspp)=L1
			#if (spp.code=="BOR")
			#	Bspp[[L1]] = Bspp[[L1]][c(0,iisens)+1]  ## BOR 2021: remove no_CPUE
			#if (useRlow)
			#	Bspp[[L1]][[paste0(names(Bspp$BOR)[1],"_low")]] = Bspp[[L1]][[1]][Rlow]
			boxlim    = c(0, max(sapply(Bspp[[L1]],quantile,quants5[5],na.rm=T)) )
			boxlim[2] = boxlim[2] + diff(boxlim) * 0.04
			left.space = 
				if (istock %in% c("BCN")) c(7.5,10)
				else if (istock %in% c("YMR","CAR","POP")) c(12,16)
				else if (istock %in% c("YTR")) c(13,16)
				else if (spp.code %in% c("SGR")) c(11,14)
				else c(10,12)
			rcol =c(.colBlind[c("redpurple","bluegreen")], "blue")  ## for LRP, USR, median
			out = compBmsy(Bspp=Bsens, spp=L1, boxwidth=0.5, medcol=medcol, boxfill=boxfill, boxlim=boxlim, whisklwd=2, staplelwd=2, Mnams=Mnams[bord], width=9, height=6, figgy=eval(parse(text=mess)), pngres=pngres, spplabs=F, t.yr=currYear, left.space=left.space, top.space=1.25, fout=paste0(prefix, sub("\\.$","",sen.fig), ifelse(Ubase>1,".stock.","."), "status",ifelse(useRlow,"+","")), calcRat=F, lang=lang, rlty=c(4,5,3), rlwd=c(2,2,1), cex.axis=1, cex.lab=1.4, rcol=rcol, ratios=switch(RPbase, 'BMSY'=c(0.4,0.8,medCR), 'B0'=c(0.16,0.32,medCR)), refpt=switch(RPbase, 'BMSY'="MSY", 'B0'=0) )
#browser();return()
			save("Bsens","L1","medcol","boxfill","boxlim","Mnams","bord","mess", file=paste0("Bsens.for.", istock, ".rda"))
		}
		rm(list=setdiff(ls(),c(keep.pars,"keep.pars")))  ## need to clear the objects in first loop so they don't mask those in second environment
		detach(ee)
	} ## end envo loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.senso


## predictRec --------------------------2025-10-23
##  Fit rockfish recruitment to an environmental index.
## ---------------------------------------------RH
predictRec <- function(rec, indices, mos=1:12,
   ryrs=NULL, rfun=median, ifun=mean,
   nmcmc=2000, refit=TRUE, polyno=2, outnam,
   png=FALSE, pngres=400, PIN=c(10,8), lang=c("f","e"))
{
	## Subfunctions----------------------
	## Polynomial regression (Michal Cukrowski)
	## https://towardsdatascience.com/polynomial-autoregression-improve-your-forecasts-in-2-minutes-746d8b57d896
	polyreg <- function(x, y, p, ar){
		eval(parse(text="require(dplyr)"))
		xdata=sapply(1:p, function(i){x^i})
		colnames(xdata)=paste0("x.",1:p)
		data = data.frame(y=y,xdata)
		data1 = data.frame(y=y)
		for (i in 1:ar){ 
			data[paste0("L", i)] <- Hmisc::Lag(data1$y, i) #here we are adding lags
			for (j in 2:ar) {
				data[paste0("L", i, ".", j)] <- Hmisc::Lag(data1$y, i)^j  ## from here on we are adding polynomial lags
			}
		}
		data = data %>% dplyr::mutate_at(names(data), ~(scale(.) %>% as.vector))
		data <- na.omit(data)  #get rid of the NAs
		#boxplot(data)
		model <- lm(y ~ ., data)
		model2 <- MASS::stepAIC(model, k=2) # k=log(nrow(train)) - if you want to use BIC
#browser();return()
		print(summary(model2))
		#refit <- forecast::Arima(data$y, model=ar)
		#test = data
		fcast <- predict(model2, data)
		plot(data$y, type='l')
		#lines(refit$fitted, col='purple')
		lines(fcast, col='red'); points(fcast, col="red")
	}
	merde <- function(dat, iter=10000, chains=3) {  ## [not use] the function crashes because rstan is crap
		modres = brm(
			bf(response ~ poly(value, 2) + ar(time = time)),
			data=dat, iter=iter, chains=chains,
			prior = c(
				set_prior("normal(0, 1)", class = "ar"),
				set_prior("normal(0, 10)", class = "b"),
				set_prior("student_t(3, 0, 2)", class = "sigma"),
				set_prior("normal(0, 10)", class = "Intercept")
			),
			#backend = "cmdstan"
		)
	}
	## Penn State : STAt 510 Applied time Series Analysis
	## https://online.stat.psu.edu/stat510/lesson/8/8.1
	goPenn <- function(regmodel, x, y) {  ## [not used]
		eval(parse(text="require(astsa)"))
		trend = time(y)
		#regmodel=lm(y~trend+x) # first ordinary regression.
		dtx = residuals(lm (x~time(x)))  ## change in x
#browser();return()
		regmodel= lm(y~ trend + dtx) # first ordinary regression with detrended x, dtx.
		summary(regmodel) # This gives us the regression results
		astsa::acf2(resid(regmodel)) # ACF and PACF of the residuals
		adjreg = astsa::sarima (y, 0,0,1, xreg=cbind(trend, dtx)) # This is the adjustment regression with MA(1) residuals
		adjreg # Results of adjustment regression. White noise should be suggested.
		astsa::acf2(resid(adjreg$fit))
		eval(parse(text="require(nlme)"))
		glsfit = nlme::gls(y ~ dtx + trend, correlation=corARMA(form = ~ 1, p=0, q=1))
		summary(glsfit) #For comparison
		return(glsfit$fitted)
	}
	## end Subfuctions-------------------

	fnam = as.character(substitute(rec))
	colnames(rec) = gsub("Recr_|RecrDev_","", colnames(rec))
	if (!is.null(ryrs))
		rec = rec[,as.character(ryrs)]
	if (!is.null(rfun) && is.function(rfun)) {
		tmp.yr = lapply(rec, rfun, na.rm=T)
		rec.yr = t(do.call("rbind", lapply(tmp.yr, data.frame, stringsAsFactors=FALSE)))
		rownames(rec.yr) = as.character(substitute(rfun))
	} else {
		rec.yr = rec[1:nmcmc,]
	}
	#rc = PBSmodelling:::.findSquare(length(indices))
	rc=c(1,1)  ## language loop makes multi-panel plots difficult

	for (i in 1:length(indices)) {
		ii = indices[i]
		mess = c(
			paste0("load(\"", ii, ".rda\")"),
			paste0("index <- ", ii)
			#paste0("index <- index[is.element(index$month, ", deparse(mos), "),]")
		)
		.flush.cat(mess, "\n")
		eval(parse(text=paste0(mess, collapse="; ")))
		.flush.cat("   index years : ", min(index$year), " to ", max(index$year), "\n", sep="")
#browser();return()
		if ("month" %in% colnames(index)) {
			index <- index[is.element(index$month, mos),]
			## Deal with months spanning year changes (e.g., 11,12,1,2)
			if (any(diff(mos)<1)) {
				zlate = mos > rev(mos)[1]
				zearly = !zlate
				if (sum(zlate) > sum(zearly)) { ## change early-month years to previous year
					zmos = is.element(index$month, mos[zearly])
					index$year[zmos] = index$year[zmos] - 1
				} else 
				if (sum(zearly) >= sum(zlate)) { ## change late-month years to next year
					zmos = is.element(index$month, mos[zlate])
					index$year[zmos] = index$year[zmos] + 1
				}
			}
		} else {
			mess = paste0("Index '", ii, "' does not have a 'month' field")
			message("\n", mess, "\n")
		} ## end if 'month' is a field
		idx.grp = split(index$anomaly, index$year)
		idx.yr  = sapply(idx.grp, ifun)
		iyrs    = intersect(names(idx.yr), colnames(rec.yr))
		idx.yr  = idx.yr[iyrs]
		rec.yr  = rec.yr[,iyrs,drop=FALSE]
		yrs     = as.numeric(iyrs)
		logged  = FALSE
		if (grepl("R.mcmc", fnam)) {
			logged = TRUE
			rec.yr = log10(rec.yr)  ## for statistical reasons (Michal Cukrowski)
		}

		## Test Philina's brms model  (cannot install rstan using Rtools44, officially giving up RH:240610)
		#require(brms); require(cmdstanr)
		#testdat = cbind(response=log(apply(rec.yr,2,median)), value=idx.yr, time=yrs)
		#guano   = merde(testdat, iter=nmcmc, chains=3)

		if (!refit && !is.null(tcall(ifit)) && length(tcall(ifit))==nmcmc && isnull(attributes(tcall(ifit))$index,"nada")==ii) {
			tget(ifit)
		} else {
			ifit = list()
			for (j in 1:nrow(rec.yr)) {
				idx   = idx.yr
				jj    = rownames(rec.yr)[j]
				kk    = as.ts(yrs)
				jrec  = unlist(rec.yr[j,])
				## Various attempts discarded
				#jfit  = lm(jrec ~ poly(idx, polyno, raw=T))
				#fit.arima = auto.arima(y=jrec.raw)
				#fit.ar = ar(x=jrec.raw)
				#if (length(fit.ar$ar) > 1000)
				#	jrec = Arima(y=jrec.raw, order=c(length(fit.ar$ar),0,0))$fitted
				#else
				#	jrec = jrec.raw
				#jfit  = lm(jrec ~ idx)
				#kfit  = polyreg(x=idx, y=jrec, p=polyno, ar=ar)  ##testing polynomial regression (not what I want)
				#goPenn(jfit, x=idx, y=jrec)
				
				## Solution using nlme::gls (Sean Anderson 240612)
				data = data.frame(time=yrs-min(yrs)+1, y=jrec, x=idx)
				#jfit <- gls(y ~ poly(x,polyno), correlation=corAR1(form = ~ time), data=data)
				jmod <- substitute(gls(y ~ poly(x,polyno), correlation=corAR1(form = ~ time), data=data))
				jfit = eval(jmod)
				xnew  = seq(min(idx), max(idx), length=length(idx))
				newdata = data[,setdiff(colnames(data), "y")]
				newdata$x = xnew
#browser();return()
				jcoef = coef(jfit); jcoef=c(jcoef,Phi=intervals(jfit)$corStruct["Phi","est."])
				#ynew.lm  = jcoef[1] + jcoef[2] * xnew  ## only good for simple llinear model
				ynew = predict(jfit, newdata=newdata)
				#plot(idx, jrec); lines(xnew, ynew.lm, col="blue",lw=3); lines(xnew, ynew, col="red") ## diagnostic for linear fit
				ifit[[jj]] = list(idx=idx, jrec=jrec, jfit=jfit, xnew=xnew, jcoef=jcoef, ynew=ynew)
			}
			attr(ifit, "index") = ii
			tput(ifit)
		}  ## end if fit  or refit
		nfit = length(ifit)
		## Collect quantiles
		q5      = c(0.05, 0.25, 0.5, 0.75, 0.95)
		yobs.q  = apply(rec.yr, 2, quantile, q5)
#browser();return()
		#yobs.q  = apply(yobs.p, 2, quantile, q5)
		xnew.mc = lapply(ifit, function(x) { xx=as.data.frame(matrix(x$xnew,nrow=1)); return(xx) })
		xnew    = do.call("rbind", lapply(xnew.mc, data.frame, stringsAsFactors=FALSE, check.names=F))
		xnew.o  = xnew[1,]  ## they are all the same
		ynew.mc = lapply(ifit, function(x) { xx=as.data.frame(matrix(x$ynew,nrow=1)); return(xx) })
		ynew    = do.call("rbind", lapply(ynew.mc, data.frame, stringsAsFactors=FALSE, check.names=F))
		ynew.q  = apply(ynew, 2, quantile, q5)
		coef.mc = lapply(ifit, function(x) { xx=as.data.frame(matrix(x$jcoef,nrow=1)); return(xx) })
		coef    = do.call("rbind", lapply(coef.mc, data.frame, stringsAsFactors=FALSE, check.names=F))
		coef.q  = apply(coef, 2, quantile, q5)
		#cnams = c("a", paste0("b^", 1:polyno),"phi"); cnams=sub("b\\^1","b",cnams)
		cnams = c(paste0("b", 0:polyno),"phi")
		cnams = c(paste0(convUTF("\u{03B2}"), 0:polyno),convUTF("\u{03D5}")) ## test
		colnames(coef.q) = cnams  ## will need to change for other models
		xCL = as.vector(rbind(idx.yr,idx.yr,NA))
		yCL = as.vector(rbind(yobs.q[c(1,5),],NA))

		if (missing(outnam))
			fout.e = paste0("pred", sub("\\.mcmc","",fnam), ifelse(logged,"(log10)",""), "-", ii)
		else
			fout.e = outnam
#browser();return()
		for (l in lang) {
			changeLangOpts(L=l)
			## Need to process legtxt withing language loop
			legtxt = apply(coef.q, 2, function(x) { paste0( formatC(x[3],digits=3,format="g",flag="#"), " (", formatC(x[1],digits=3,format="g",flag="#"), ", ", formatC(x[5],digits=3,format="g",flag="#"), ")") })
			#legtxt = paste0("expression(", names(legtxt), "==", gsub("\\s+", "~", legtxt), ")")
			legtxt = paste0(names(legtxt), " = ", legtxt) #gsub("\\s+", "~~", legtxt))
			if (nrow(ynew)==1) 
				legtxt = gsub(r"{\s*\([^\)]+\)}","",legtxt)  ## remove credibility limits

			fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
			if (png) png(filename=paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
			if (i==1)
				expandGraph(mfrow=rc, mar=c(3,3.5,1,1), mgp=c(1.9,0.6,0))
			xlim = range(xnew.o)
			#ylim = quantile(unlist(ynew), probs = if(nfit==1) c(0,1) else c(0.01,0.99)) #range(rec.yr)
			ylim = range(c(yobs.q, unlist(ynew.q)), na.rm=T)
			#ylim = range(c(unlist(ynew.q)), na.rm=T)
			ylab = switch (fnam, 'R.mcmc'="Recruitment", 'Rdev.mcmc'="Recruitment Deviation", "Something")
			if (logged) ylab = paste0("log10 ", ylab)
			plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=linguaFranca(toupper(ii),l), ylab=linguaFranca(ylab,l), cex.axis=1.2, cex.lab=1.5)
			axis(1, at=pretty(xlim, n=10L), labels=FALSE, tcl=-0.3)
			polygon(x=c(0,0,par()$usr[c(2,2)]), y=par()$usr[c(3,4,4,3)], border=F, col=lucent(.colBlind["vermillion"],0.2)) #"darkseagreen1")
			polygon(x=c(par()$usr[c(1,1)],0,0), y=par()$usr[c(3,4,4,3)], border=F, col=lucent(.colBlind["blue"],0.2)) #"mistyrose")
			abline(v=0, col="gainsboro", lwd=2)
#browser();return()
			if (ii %in% c("ui51","ao"))
				addLegend(0.025, 0.975, legend=legtxt, bty="n", xjust=0)  ## forget using expressions -- too finicky and drops significant zeroes
			else
				addLegend(0.975, 0.975, legend=legtxt, bty="n", xjust=1)  ## forget using expressions -- too finicky and drops significant zeroes
			#addLegend(0.95, 0.95, legend=eval(parse(text=paste0("expression(", legtxt[1], ")"))), bty="n", xjust=1)
			#addLegend(0.95, 0.925, legend=eval(parse(text=paste0("expression(", legtxt[2], ")"))), bty="n", xjust=1)
	
			if (nfit==1) {
				unpackList(ifit[[1]], scope="L")
				points(idx, jrec, pch=21, cex=2.5, col="blue", bg="aliceblue",lwd=1)
				text(idx, jrec, labels=substr(names(idx),3,4), cex=0.8, col="black")
				lines(xnew, ynew, col="red", lwd=2)
			}
	#		for (j in 1:nfit) {
	#			jj = names(ifit)[j]
	#			unpackList(ifit[[jj]], scope="L")
	#			points(idx,jrec, pch=20, col="grey", cex=0.8)
	#		}
			if (nfit>1) {
				lcol = "black"
				lines(xCL, yCL, col=lucent("black",0.1), lwd=2) ## credibility limits by year
				lines(xnew.o, ynew.q[3,], lty=1,lwd=3, col=lcol)
				lines(xnew.o, ynew.q[2,], lty=3,lwd=2, col=lcol)
				lines(xnew.o, ynew.q[4,], lty=3,lwd=2, col=lcol)
				lines(xnew.o, ynew.q[1,], lty=2,lwd=2, col=lcol)
				lines(xnew.o, ynew.q[5,], lty=2,lwd=2, col=lcol)
				points(idx.yr, yobs.q[3,], pch=21, cex=2.5, col="black", bg="gainsboro",lwd=1)  ## median recruitment
				text(idx.yr, yobs.q[3,], labels=substr(names(idx.yr),3,4), cex=0.8, col="black")
##browser();return()
				#points(idx.yr, yobs.q[3,], pch=21, col="balck", bg="yellow
#browser();return()
			}
			box()
			if (png) dev.off()
		}; eop()  ## end l loop for lang
	} ## end i loop indices
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~predictRec


## tabSS.compo--------------------------2025-12-04
## Make Base Case Tables
## Note: u2023=u2022 (see 'gatherMCMC.r') so change 
##       labels here in rfpt tables to use 'prevYear'
## ---------------------------------------------RH
tabSS.compo <- function(envo, #istock="YTR", prefix="ytr.", compo,
  useRlow=FALSE, qRlow=0.25, sigdig=4)
{
	vomit <- function() { gc(verbose=FALSE); while ("ee" %in% search()) detach(ee) }; #resetGraph() }
	on.exit( vomit() )

	for (e in 1:length(envo)) {
		keep.pars = ls()  ## need to clear environment of objects on subsequent loops through envo
		## Note that objects generated below go into the function environment, not the model environments (e1 ad e2)
		ee = envo[[e]]
		attach(ee)
		rm(list=intersect(ls(.GlobalEnv), ls(ee)), envir=.GlobalEnv)  ## remove objects in global environment masking objects in e1
		#unpackList(stock[[istock]][["Controls"]])
		#unpackList(compo)

		modYrs      = startYear:currYear
		modYrsChar  = as.character(modYrs)
		## Diagnostics for select parameters
		P.mpd       = ampdPA; storage.mode(P.mpd)="double"     ## For some reason, this matrix is stored as a list (maybe due to NA values?)
		P.names     = colnames(avgPA) = sub("\\_gp\\([0-9]\\)", "", colnames(avgPA))
		P.runs      = as.numeric(sapply(strsplit(rownames(avgPA),"\\."),function(x){x[1]}))
		P.run.rwt   = sapply(strsplit(rownames(avgPA),"\\."),function(x){paste0(x[1],".",x[2],".",x[3])}) ## vector of Runs and Rwts
		P.run.ord   = unique(P.runs)
		P.run.nmc   = table(P.runs)[as.character(P.run.ord)]
		P.run.num   = rep(1:length(P.run.nmc), P.run.nmc)
		use.run.rwt = is.element(P.runs, P.run.ord) ## default use all runs but a subset might be used for management advice
		B.labs      = paste0(rep("B",NrefM),1:NrefM," (R",P.run.ord,")")
		if (spp.code %in% c("YMR"))
			P.ord=c(1,seq(2,10,2),seq(3,11,2))
		else
			P.ord = 1:ncol(ampdPA)

		##----------Table 1-----------
		## MCMC model parameters table
		##----------------------------
		tabPmed = t(apply(avgPA,2,quantile,tcall(quants5),na.rm=T))  ## arrays so cannot use sapply
		tabPmed = tabPmed[P.ord,]
		tabPmed = formatCatch(tabPmed,N=sigdig)
	
		names.pars = rownames(tabPmed)
		## (RH 251010) This substitution sequence was overhauled for SGR 2025

		## Deal with stupid catchabilities first
		zq = grep("LnQ", names.pars)
		if (any(zq)) {
			names.pars[zq] = sub("LnQ\\_base","\\\\log q", names.pars[zq])
			names.pars[zq] = sub("TRAWL_(BC|5ABC|5DE|3CD)", "TRAWL~\\1_", names.pars[zq ])
			names.pars[zq] = sub("\\(([0-9])\\)", "{\\1}", names.pars[zq ])
			names.pars[zq] = sapply(strsplit(names.pars[zq], split="\\_"), function(x){ paste0(x[1],"_",x[3],"~(\\text{",x[2],"})")})
		}
		## Fix up Rdist
		zd = grep("Rdist", names.pars)
		if (any(zd)) {
			oRdist = rownames(tabPmed)[zd]
			oRdist = sub("3",areas[3],sub("2",areas[2],sub("1",areas[1],oRdist)))
			names.pars[zd] = paste0("\\text{", gsub("\\_","~",oRdist), "}")
			names.pars[zd] = sub("area", "", names.pars[zd])
		}
		names.pars = sub("LN\\(R0)", "\\\\log R_{0}", names.pars)
		names.pars = sub("LN\\(R0)", "\\\\log R_{0}", names.pars)
		names.pars = sub("\\_GP1", "", names.pars)
		names.pars = sub("M\\_Female", "M_{1}~(\\\\text{Female})", names.pars)
		names.pars = sub("M\\_Male", "M_{2}~(\\\\text{Male})", names.pars)
		names.pars = sub("BH\\_(h)", "h~(\\\\text{BH steepness})", names.pars)
		## Fix selectivity
		zs = grep("mu|varL|[Dd]elta", names.pars)
		if (any(zs)) {
			names.pars[zs] = sub("\\)","}",sub("\\(","_{",gsub("\\_","~", names.pars[zs])))
			names.pars[zs] = sub("~(TRAWL|QCS|WCHG|WCVI|HS|GIG|NMFS)", "~(\\\\text{\\1", names.pars[zs])
			names.pars[zs] = paste0("\\", names.pars[zs], "})")
			names.pars[zs] = sub("\\\\varL\\_\\{([0-9])\\}", "\\\\log v_{\\\\text{L}\\1}", names.pars[zs])
		}
		rownames(tabPmed) =  paste(rep("$",nrow(tabPmed)),names.pars,rep("$",nrow(tabPmed)),sep="")
		#colnames(tabPmed) =  gsub("\\%","\\\\%",colnames(tabPmed))
	
		xtab.compo.pars = xtable(tabPmed, align="lrrrrr",
			label   = paste0("tab:",prefix,"base.pars"), digits = if (exists("formatCatch")) NULL else sigdig,
			caption = paste0("Base run (", istock, ")~: the ", texThatVec(tcall(quants5)), " quantiles for ", ifelse(NrefM>1,"pooled",""),
			" model parameters (defined in \\AppEqn) from MCMC estimation of \\numberstringnum{", NrefM, "} ",
			ifelse(NrefM>1, "component model runs of \\Nmcmc{} samples each.", "base run of \\Nbase{} samples.") ) )
		xtab.compo.pars.out = capture.output(print(xtab.compo.pars,
			caption.placement="top", 
			sanitize.rownames.function=function(x){x}, 
			add.to.row =list(pos=list(-1), command=c("\\\\[-1.0ex]")),
			size = "\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\footnotesize" ) )
		#tput(xtab.compo.pars.out)  ## no reason to save this in an environment

	# \\hline\\\\[-2.2ex] CC & 2023 & 2024 & 2025 & 2026 & 2027 & 2028 & 2029 & 2030 & 2031 \\\\[0.2ex]\\hline\\\\[-1.5ex]
		##----------Table 2-----------
		## MCMC Derived Parameters
		##----------------------------
		yrsUse = modYrs[c(1,length(modYrs)+c(-1,0))]
		yrsChr = as.character(yrsUse)

		if (RPbase=="BMSY") {
			B.mcmc = data.frame (
				B0         = avgRP[,"B0"],
				Bcurr      = avgRP[,"Bcurr"],
				Bcurr.B0   = avgRP[,"Bcurr"] / avgRP[,"B0"],
				ucurr      = avgRP[,"ucurr"],
				umax       = apply(avgTS[,,"ut"],1,max,na.rm=T),  ## for each mcmc sample across the time series (note 1933 ut=NA)
				MSY        = avgRP[,"MSY"],
				Bmsy       = avgRP[,"Bmsy"],
				LRP        = avgRP[,"LRP"],
				USR        = avgRP[,"USR"],
				Bcurr.Bmsy = avgRP[,"Bcurr"] / avgRP[,"Bmsy"],
				Bmsy.B0    = avgRP[,"Bmsy"] / avgRP[,"B0"],
				umsy       = avgRP[,"umsy"],
				ucurr.umsy = avgRP[,"ucurr"] / avgRP[,"umsy"]
			)
		}
		if (RPbase=="B0") {
			## eventually grab code from SweaveMCMC.Rnw
			B.mcmc = data.frame (
				B0         = avgRP[,"B0"],
				Bcurr      = avgRP[,"Bcurr"],
				Bcurr.B0   = avgRP[,"Bcurr"] / avgRP[,"B0"],
				ucurr      = avgRP[,"ucurr"],
				umax       = apply(avgTS[,,"ut"],1,max,na.rm=T), ## for each mcmc sample across the time series (note 1933 ut=NA)
				Ytrp       = avgRP[,"Ytgt"],
				#Btrp       = 0.40 * avgRP[,"B0"],  ## this is OK as long as B0=SSB_1933 (VIRG) and not SSB_1935
				Btrp       = avgRP[,"Btgt"],
				utrp       = avgRP[,"utgt"],
				#LRP        = 0.16 * avgRP[,"B0"],
				#USR        = 0.32 * avgRP[,"B0"],
				LRP        = 0.4 * avgRP[,"Btgt"],
				USR        = 0.8 * avgRP[,"Btgt"],
				Bcurr.Btrp = avgRP[,"Bcurr"] / avgRP[,"Btgt"],
				ucurr.utrp = avgRP[,"ucurr"] / avgRP[,"utgt"]
			)
		}
		tabQmed = t(apply(B.mcmc,2,quantile,tcall(quants5),na.rm=T)) 
#browser();return()

		## Extract area-based information if it exists
		if (exists("xavgRP") && length(xavgRP)>1) {
			## colnames(xavgRP) = "Bmsy" "Fmsy" "MSY"  "umsy" "Btgt"  "Ftgt"  "Ytgt"  "utgt"  "LRP"  "USR"  "Vmsy" "Vtgt"  "B0"   "V0"   "pVB"
			## dimnames(xavgTS)[3] =  "B"      "BtBtgt"  "BtBmsy" "D"      "fR"     "R"      "u"      "ututgt"  "utumsy" "V"
			if (RPbase=="BMSY") {
				A.collect = list('BC'=tabQmed)
				for (a in dimnames(xavgRP)[[3]]) {
					A.mcmc = data.frame (
						## Needs revision (see below)
						B0         = avgRP[,"B0"],
						Bcurr      = avgRP[,"Bcurr"],
						Bcurr.B0   = avgRP[,"Bcurr"] / avgRP[,"B0"],
						ucurr      = avgRP[,"ucurr"],
						umax       = apply(avgTS[,,"ut"],1,max,na.rm=T), ## for each mcmc sample across the time series
						MSY        = avgRP[,"MSY"],
						Bmsy       = avgRP[,"Bmsy"],
						LRP        = avgRP[,"LRP"],  ## LRP and USR automatically set to 0.4Bmsy and 0.8Bmsy in 'load_extra_mcmc.R' for multi-area models
						USR        = avgRP[,"USR"],
						Bcurr.Bmsy = avgRP[,"Bcurr"] / avgRP[,"Bmsy"],
						Bmsy.B0    = avgRP[,"Bmsy"] / avgRP[,"B0"],
						umsy       = avgRP[,"umsy"],
						ucurr.umsy = avgRP[,"ucurr"] / avgRP[,"umsy"]
					)
					tabAmed = t(apply(A.mcmc,2,quantile,tcall(quants5),na.rm=T)) 
					A.collect[[a]] = tabAmed
				} ## end area a loop
			} ## end RPbase 'BMSY'

			if (RPbase=="B0") {
				A.collect = list('BC'=tabQmed)
				for (a in dimnames(xavgRP)[[3]]) {
					A.mcmc = data.frame (
						B0         = xavgRP[,"B0",a],
						Bcurr      = xavgTS[,as.character(currYear),"B",a],
						Bcurr.B0   = xavgTS[,as.character(currYear),"B",a] / xavgRP[,"B0",a],
						ucurr      = xavgTS[,as.character(currYear),"u",a],
						umax       = apply(xavgTS[,,"u",a],1,max,na.rm=T), ## for each mcmc sample across the time series
						Ytrp       = xavgRP[,"Ytgt",a],
						#Btrp       = 0.40 * xavgRP[,"B0",a],  ## this is OK as long as B0=SSB_1933 (VIRG) and not SSB_1935
						Btrp       = xavgRP[,"Btgt",a],
						utrp       = xavgRP[,"utgt",a],
						#LRP        = 0.16 * xavgRP[,"B0",a],
						#USR        = 0.32 * xavgRP[,"B0",a],
						LRP        = 0.4 * xavgRP[,"Btgt",a],
						USR        = 0.8 * xavgRP[,"Btgt",a],
						Bcurr.Btrp = xavgTS[,as.character(currYear),"B",a] / xavgRP[,"Btgt",a],
						ucurr.utrp = xavgTS[,as.character(currYear),"u",a] / xavgRP[,"utgt",a]
					)
					tabAmed = t(apply(A.mcmc,2,quantile,tcall(quants5),na.rm=T)) 
					A.collect[[a]] = tabAmed
#browser();return()
				} ## end area a loop
			} ## end RPbase 'B0'
			## How to interleave a data frame (Google AI 251016)
			## Add an identifier and row number to each data frame
			interleaved_list <- lapply(seq_along(A.collect), function(i) {
				df <- as.data.frame(A.collect[[i]])
				df$df_id <- rep(i,nrow(df)) # Identifier for the original data frame
				df$row_num <- 1:nrow(df) # Original row number within each data frame
				return(df)
			})
			## Combine all data frames
			combined_df <- do.call(rbind, interleaved_list)
			## Sort by row_num and then by df_id to interleave
			interleaved_df <- combined_df[order(combined_df$row_num, combined_df$df_id), ]
			## Remove the helper columns if not needed
			interleaved_df$df_id <- NULL
			interleaved_df$row_num <- NULL
			tabQmed.original = tabQmed
			tabQmed = interleaved_df
			xareas  = dimnames(xavgRP)[[3]]
			rownames(tabQmed) = sub("3",paste0(".",xareas[3]), sub("2",paste0(".",xareas[2]), sub("1", paste0(".",xareas[1]), rownames(tabQmed) ) ) )
#browser();return()
		} ## end extra area collection

		tab.base.rfpt = formatCatch(tabQmed,N=sigdig-1)  ## use 3 instead of 4
		names.rfpt = rownames(tab.base.rfpt)

		names.rfpt =
			gsub("\\.(B|u)", "/\\1",
			gsub("_Trawl", "~(\\\\text{trawl})",
			gsub("_Other", "~(\\\\text{other})",
			gsub("trp", "_\\\\text{TRP}",
			gsub("msy", "_\\\\text{MSY}",
			gsub("LRP",  paste0(switch(RPbase, 'BMSY'="0.4B_{\\\\text{MSY}}", 'B0'="0.16B_0")),
			gsub("USR",  paste0(switch(RPbase, 'BMSY'="0.8B_{\\\\text{MSY}}", 'B0'="0.32B_0")),
			gsub("MSY",  paste0("\\\\text{MSY}"),
			gsub("VB",  "V",
			gsub("[Uu]max", "u_\\\\text{max}",
			gsub("ucurr",  paste0("u_{",prevYear,"}"),
			gsub("Bcurr",  paste0("B_{",currYear,"}"),
			gsub("B0", "B_{0}",
			names.rfpt)))))))))))))
#browser();return()
		if (exists("areas")) {
			za = grep(paste0(areas,collapse="|"), names.rfpt)
			if (any(za)) {
				apat = paste0("(",paste0(areas,collapse="|"),")")
				names.rfpt[za] = gsub(apat, "~~~\\1", names.rfpt[za])
				names.rfpt[za] = gsub(paste0("(.*)\\.~~~", apat), "~~~\\\\text{\\2}~~\\1", names.rfpt[za])
			}
		}
		rownames(tab.base.rfpt) =  paste(rep("$",nrow(tab.base.rfpt)),names.rfpt,rep("$",nrow(tab.base.rfpt)),sep="")

		#remove.rows = grep("~~0\\.16|~~0\\.32|~~B\\_\\\\text\\{TRP|~~u\\_\\\\text\\{max",rownames(tab.base.rfpt))
		remove.rows = grep("~~0\\.16|~~0\\.32|~~u\\_\\\\text\\{max",rownames(tab.base.rfpt))
#browser();return()
		if (length(remove.rows) > 0)
			tab.base.rfpt = tab.base.rfpt[-remove.rows,]   ## otherwise it removes all rows


		#remove.rows = grep("~~0\\.16|~~0\\.32|~~B\\_\\\\text\\{TRP|~~u\\_\\\\text\\{max",rownames(tab.base.rfpt))
		#if (length(remove.rows) > 0)
		#	tab.base.rfpt = tab.base.rfpt[-remove.rows,]   ## otherwise it removes all rows

		## Caption for Ref Point table
		cap.rfpt = paste0(
			"Base run (", istock, ")~: the ", texThatVec(tcall(quants5)), " quantiles of MCMC-derived quantities from \\Nbase{} samples ", 
			ifelse(NrefM>1,"pooled","")," from ", ifelse(NrefM>1, "component runs.", "a single base run."),
			" Definitions are: ",
			"$B_0$ -- unfished equilibrium spawning biomass (mature females), ",
			"$B_{", currYear, "}$ -- spawning biomass at the beginning of ", currYear, ", ",
			"$u_{", prevYear, "}$ -- exploitation rate (ratio of total catch to vulnerable biomass) in the middle of ", prevYear, ", ",
			"$u_\\text{max}$ -- maximum exploitation rate (calculated for each sample as the maximum exploitation rate from ",
			modYrs[1], "-", prevYear, "), ",
			switch(RPbase, 'BMSY'=paste0(c(
				"MSY -- maximum sustainable yield at equilibrium, ",
				"$B_\\text{MSY}$ -- equilibrium spawning biomass at MSY, ",
				"$u_\\text{MSY}$ -- equilibrium exploitation rate at MSY. ",
			),collapse=""), 'B0'=paste0(c(
				"$Y_\\text{TRP}$ -- equilibrium yield at target reference point (0.4$B_0$), ",
				"$B_\\text{TRP}$ -- equilibrium spawning biomass at target reference point (0.4$B_0$), ",
				"$u_\\text{TRP}$ -- equilibrium exploitation rate at target reference point (0.4$B_0$), "
			), collapse="") ),
			"All biomass values (and Yield) are in tonnes. ", refCC.sentence
		)
		xtab.compo.rfpt = xtable(tab.base.rfpt, align="lrrrrr",
			label=paste0("tab:",prefix,"base.rfpt"), digits=if (exists("formatCatch")) NULL else sigdig, caption=cap.rfpt )
		xtab.compo.rfpt.out = capture.output( 
			print(xtab.compo.rfpt,
			caption.placement="top", sanitize.rownames.function=function(x){x}, 
			add.to.row=list(pos=list(-1), command=c("\\\\[-1.0ex]" )),
			#add.to.row=list(pos=list(-1,3), command=c("\\\\[-1.0ex]", "\\hdashline \\\\[-1.75ex]")), hline.after=c(-1,0,5,nrow(xtab.compo.rfpt)) 
			hline.after =  c(-1,0,nrow(tab.base.rfpt)),
			size = "\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\footnotesize"
		)	)
#browser();return()
	
		## Component run likelihoods
		## -------------------------
		tabll = formatCatch(avgLL)
		rownames(tabll) =
			sub("^Index$", "Abundance Index",
			sub("^Recruit$", "Recruitment",
			sub("^AF$", "Age Frequency",
			sub("BT$",  "Bottom Trawl",
			sub("TRI$", "Triennial",
			sub("HIS$", "Historical",
			sub("SYN$", "Synoptic",
			gsub("_", " ",
		rownames(tabll) ))))))))

		xtab.cruns.ll = xtable(tabll, align=paste0("l", paste0(rep("r",dim(avgLL)[2]),collapse="")),
			label   = paste0("tab:",prefix,"log.likes"), digits=NULL, 
			caption = paste0("MPD log likelihood (LL) values reported by ", ifelse(NrefM>1, "component base runs", "the single base run"), " for survey indices, age composition (AF), recruitment, and total (not all LL components reported here)") )
		xtab.cruns.ll.out = capture.output(print(xtab.cruns.ll,  include.rownames=TRUE,
			caption.placement="top", sanitize.rownames.function=function(x){x}, add.to.row =list(pos = list(-1), command = c("\\\\[-1.0ex]")) ) )
		rowhead = grep("^\\s&",xtab.cruns.ll.out)
		if (length(rowhead)>0)
			xtab.cruns.ll.out[rowhead[1]] =  paste0(c("LL value", grep("^\\s&",xtab.cruns.ll.out,value=T)[1]), collapse="")
#browser();return()
		
		## Only process these if there is more than one base component run
		if (length(run.num) > 1) {
			## Parmeter comparisons among component runs
			##  (not used unless base case has more than one base run)
			## -----------------------------------------
			## Calculate quantiles for later (eventually make this function global)
			calcQs <- function (dat, ivec, ovec, qval=quants3) {
				F.data = as.data.frame(dat)
				F.list = split(F.data, ivec)
				rr.ord = match(ovec, as.numeric(substring(names(F.list),1,2)) )  ## order not so important for sensitivities
				F.list = F.list[rr.ord]
		#browser();return()
				F.qnts = lapply(F.list,function(x){
					z = apply(x,2,function(xx){!all(is.na(xx))})
					out = apply(x[z],2,quantile,qval,na.rm=T)
					return(out)
				}) ## lapply does not sort like split does
				return(F.qnts)
			}
			P.qnts = calcQs(avgPA, ivec=P.run.rwt, ovec=P.run.ord)
			P.med  = lapply(P.qnts,function(x){apply(x,2,function(xx){c(xx[2],xx[1],xx[3])})})
			P.mpd.mcmc = as.list(rownames(ampdPA)); names(P.mpd.mcmc) = rownames(ampdPA)
			for (i in names(P.mpd.mcmc)){
				P.mpd.mcmc[[i]] = rbind(ampdPA[i,P.ord,drop=FALSE], P.med[[i]][,P.ord])
			}
			P.vals = lapply(lapply(P.mpd.mcmc,t),formatCatch)
			P.tabs = lapply(P.vals, function(mat){
				out = t(t(apply(mat, 1, function(x) {
					paste0(c("|", c(x[1], "| ", x[2], " (", x[3],",",x[4],")")), collapse="")
					#paste0(c(x[1], " , ", x[2], " (", x[3],",",x[4],")"), collapse="")
				})))
				return(out)
			})
			P.tab = array(NA,dim=rev(dim(ampdPA)), dimnames=list(names.pars,names(P.qnts)))
			for (i in names(P.tabs)){
				P.tab[,i] = P.tabs[[i]]
			}
			## Strip out fleet names in parentheses:
			i.fleet = grep("TRAWL|QCS|WCVI|NMFS|HS|WCHG|HBLL|GIG", rownames(P.tab))
			rownames(P.tab)[i.fleet] =  gsub("\\s*~\\(.*?\\)$","",rownames(P.tab)[i.fleet])
			rownames(P.tab) = paste(rep("$",nrow(P.tab)),  rownames(P.tab), rep("$",nrow(P.tab)), sep="")
			colnames(P.tab) = B.labs
		
			xtab.cruns.pars = xtable(P.tab, align=paste0("l", paste0(rep("c",dim(P.tab)[2]),collapse="")),
				label   = paste0("tab:",prefix,"runs.pars"), digits = if (exists("formatCatch")) NULL else sigdig,
				caption = paste0("Base run (", istock, ")~: model parameter MPDs (delimited by `|') and MCMC medians (with 0.05 and 0.95 quantile limits) for ",
				ifelse(NrefM>1, paste0("each of the \\numberstringnum{", NrefM, "} component model runs of \\Nmcmc{} samples each."),
				"the base model run of \\Nbase{} samples.") ) )
			xtab.cruns.pars.out = capture.output(print(xtab.cruns.pars,
				caption.placement="top", sanitize.rownames.function=function(x){x}, add.to.row =list(pos = list(-1), command = c("\\\\[-1.0ex]")) ) )
	#browser();return()
		
			## Comparisons among derived quantities for component runs
			##  (not used unless base case has more than one base run)
			## -------------------------------------------------------
			Q.qnts = calcQs(B.mcmc, ivec=P.run.rwt, ovec=P.run.ord)
			Q.med  = lapply(Q.qnts,function(x){apply(x,2,function(xx){c(xx[2],xx[1],xx[3])})})
		
			#Q.mpd.mcmc = as.list(colnames(B.mcmc)); names(Q.mpd.mcmc) = colnames(B.mcmc)  ## have not collected MPD estimates of 
			#for (i in names(Q.mpd.mcmc)){
			#	Q.mpd.mcmc[[i]] = rbind(ampdPA[i,,drop=FALSE], Q.med[[i]])
			#}
			#Q.vals = lapply(lapply(Q.mpd.mcmc,t),formatCatch)
			Q.vals = lapply(lapply(Q.med,t), formatCatch, N=2)
			Q.tabs = lapply(Q.vals, function(mat){
				t(t(apply(mat, 1, function(x) {
					paste0(c(x[1], " (", x[2],",",x[3],")"), collapse="")
					#paste0(c("|", c(x[1], "| ", x[2], " (", x[3],",",x[4],")")), collapse="")
					#paste0(c(x[1], " , ", x[2], " (", x[3],",",x[4],")"), collapse="")
				})))
			})
			Q.tab = array(NA,dim=c(ncol(B.mcmc),NrefM), dimnames=list( dimnames(B.mcmc)[[2]],names(Q.qnts)))
			for (i in names(Q.tabs)){
				Q.tab[,i] = Q.tabs[[i]]
			}
			## Note that this table has not yet been expanded to include subareas (remove them for now)
			rownames(Q.tab) = names.rfpt[grep("5ABC|5DE|3CD", names.rfpt, invert=T)]
			i.fleet = grep("TRAWL|QCS|WCVI|NMFS|HS|WCHG|HBLL|GIG", rownames(Q.tab))
			if (length(i.fleet)>0)
				rownames(Q.tab)[i.fleet] =  gsub("\\s*~\\(.*?\\)$","",rownames(Q.tab)[i.fleet])
			rownames(Q.tab) = paste(rep("$",nrow(Q.tab)),  rownames(Q.tab), rep("$",nrow(Q.tab)), sep="")
			colnames(Q.tab) = B.labs
		
			xtab.cruns.rfpt = xtable(Q.tab, align=paste0("l", paste0(rep("r",dim(Q.tab)[2]),collapse="")),
				label   = paste0("tab:",prefix,"runs.rfpt"), digits = if (exists("formatCatch")) NULL else sigdig,
				caption = paste0("Base run (", istock, ")~: MCMC median (with 0.05 and 0.95 quantile limits) for derived model quantities for ",
				ifelse(NrefM>1, paste0("each of the \\numberstringnum{", NrefM, "} component model runs of \\Nmcmc{} samples each."),
				"the base model run of \\Nbase{} samples.") ) )
			xtab.cruns.rfpt.out = capture.output(print(xtab.cruns.rfpt, caption.placement="top", 
				sanitize.rownames.function = function(x){x},
				add.to.row = list(pos=list(-1,3), command=c("\\\\[-1.0ex]", "\\hdashline \n")),
				hline.after=c(-1,0,5,nrow(xtab.cruns.rfpt)) ) )
		}  ## end if more than 1 base run
#browser();return()

		outwithit = c("xtab.compo.pars.out", "xtab.compo.rfpt.out", "xtab.cruns.ll.out")
		if (length(run.num) > 1)
			outwithit = c( outwithit, "xtab.cruns.pars.out", "xtab.cruns.rfpt.out")
		save(list=outwithit, file=paste0(prefix,"compo.tabs", ifelse(useRlow, paste0(".Rlow(q=",pad0(qRlow*100,2),")"),""), ".rda") )
		## Tables are loaded when 'AppG_Rcode.R' is run (usually once), 
		##  but the tables may be regenerated and no longer sync with e1|e2
		for (i in outwithit) {
			mess = paste0("ee$", i, " <- ", i)
			eval(parse(text=mess))
		}
		#ee$xtab.compo.pars.out <- xtab.compo.pars.out
		#ee$xtab.compo.rfpt.out <- xtab.compo.rfpt.out
		#ee$xtab.cruns.ll.out   <- xtab.cruns.ll.out
		#ee$xtab.cruns.pars.out <- xtab.cruns.pars.out
		#ee$xtab.cruns.rfpt.out <- xtab.cruns.rfpt.out

		rm(list=setdiff(ls(),c(keep.pars,"keep.pars")))  ## need to clear the objects in first loop so they don't mask those in second environment
		detach(ee)
	} ## end envo loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tabSS.compo


## tabSS.decision-----------------------2025-12-10
## Make Decision Tables (Probabilities)
## ---------------------------------------------RH
tabSS.decision <- function(envo, #istock="YTR", prefix="ytr.", compo,
  useRlow=FALSE, qRlow=0.25, decdig=2, cp.type="CC",  ## HR not used in SS3 yet (although it could be)
  cp.num=c(1,3,6,9:15), use.agile=FALSE)
{
	vomit <- function() { gc(verbose=FALSE); while ("ee" %in% search()) detach(ee) }; #resetGraph() }
	on.exit( vomit() )
	require(xtable)

	for (e in 1:length(envo)) {
		keep.pars = ls()  ## need to clear environment of objects on subsequent loops through envo
		ee = envo[[e]]
		attach(ee)
		rm(list=intersect(ls(.GlobalEnv), ls(ee)), envir=.GlobalEnv)  ## remove objects in global environment masking objects in e1
		unpackList(stock[[istock]][["Controls"]])
		big.mark   = ifelse(is.null(options()$big.mark), ",", options()$big.mark)
	
		## Check to see if projections are alternative runs (e.g., low-recruitment)
		is.renso = as.character(substitute(compo)) == "renso"
		if (is.renso) {
			names(compo) = sub("rmpd","ampd",sub("ren","avg",names(compo)))
		}
		unpackList(compo)
	
		proYears    = (currYear+1):projYear
		refYears    = as.character(currYear:projYear)  ## include current year in decision tables
	#browser();return()
	
		#lab.recentCatchYears = texThatVec(rev(rev(rownames(catch))[1:num.recentCatchYears]))
		recentCatchMean = lapply(recentCatch,mean);
		num.recentCatchYears = lapply(recentCatch,length)
		recentHarvMean  = 0.01 ## dummy for now
	
		lab.recentCatchYears = lapply(recentCatch,function(x){texThatVec(names(x))})
		refCC.sentence = paste0(" For reference, the average catch over the last ", num.recentCatchYears[[1]], " years (", lab.recentCatchYears[[1]], ") was ", 
			if (length(recentCatchMean)==1)
				formatC(recentCatchMean[[1]], digits=0, format="f", big.mark=big.mark)
			else
				paste0(paste0(names(recentCatchMean),"=",sapply(recentCatchMean, round, dig=0)), collapse=", ")
		, "~t. ")
		#refHR.sentence = paste0(" For reference, the average harvest rate over the last ", num.recentCatchYears, " years (", lab.recentCatchYears, ") was ", round(recentHarvMean, dig=3), ". ")
	#browser();return()
	
		## Use this routine developed for multi-area models (now adapted to deal with single-area models)
		if (use.agile) {
			agileDT(compo=compo, onepage=F, cp=paste0("CC.", pad0(cp.num,2)), tables.per.page=1, useRlow=useRlow, RPbase=RPbase, outnam=paste0("agileDT-",istock), pdf=F)
			unpackList(tcall(agile))
			Nmcmc.good = tcall(Nmcmc.good)
			#Nmcmc = format(tcall(Nmcmc),big.mark=",")
		} else {  ## deprecated
			## Use old routine for single-stock models
			## Diagnostics for select parameters
			P.mpd       = ampdPA; storage.mode(P.mpd)="double"     ## For some reason, this matrix is stored as a list (maybe due to NA values?)
			P.names     = colnames(avgPA)
			P.runs      = as.numeric(sapply(strsplit(rownames(avgPA),"\\."),function(x){x[1]}))
			P.run.rwt   = sapply(strsplit(rownames(avgPA),"\\."),function(x){paste0(x[1],".",x[2],".",x[3])}) ## vector of Runs and Rwts and ver
			P.run.ord   = unique(P.runs)
			P.run.nmc   = table(P.runs)[as.character(P.run.ord)]
			P.run.num   = rep(1:length(P.run.nmc), P.run.nmc)
			#use.run.rwt = is.element(P.runs, P.run.ord) ## default use all runs but a subset might be used for management advice
			use.run.rwt = unique(P.run.rwt)
		
			B.mcmc = data.frame (
				B0         = avgRP[,"B0"],
				Bcurr      = avgRP[,"Bcurr"],
				B.currB0   = avgRP[,"Bcurr"] / avgRP[,"B0"],
				Ucurr      = avgRP[,"ucurr"],
				Umax       = apply(avgTS[,,"ut"],1,max), ## for each mcmc sample across the time series
				MSY        = avgRP[,"MSY"],
				Bmsy       = avgRP[,"Bmsy"],
				LRP        = avgRP[,"LRP"],
				USR        = avgRP[,"USR"],
				Bcurr.Bmsy = avgRP[,"Bcurr"] / avgRP[,"Bmsy"],
				Bmsy.B0    = avgRP[,"Bmsy"] / avgRP[,"B0"],
				Umsy       = avgRP[,"umsy"],
				Ucurr.Umsy = avgRP[,"ucurr"] / avgRP[,"umsy"]
			)
			Nmcmc = nrow(B.mcmc)
			## Sometimes Umsy=0 for some reason
			is.good = apply(B.mcmc,1,function(x){all(is.finite(x))})
			if (sum(is.good)!=nrow(B.mcmc))
				.flush.cat(paste0(nrow(B.mcmc)-sum(is.good)," unsuitable records dropped."), "\n")
			B.mcmc  = B.mcmc[is.good,]
			Nmcmc.good = nrow(B.mcmc)
		
			## Determine subset vector for low recruitment events
			R.mcmc = avgTS[is.good,,"Rt"]
			if (useRlow) {
				if (spp.code=="BOR")
					Rlow   = R.mcmc[,"2017"] < quantile(R.mcmc[,"2017"], qRlow)
				else { 
					R.mean = apply(R.mcmc,1,mean)  ## mean recruitment 1935-2019
					Rlow   = R.mean < quantile(R.mean, qRlow)
				}
			} else {
				Rlow = rep(TRUE,nrow(R.mcmc))
			}
			## Subset vector of chosen runs to use in the base composite
			Ruse  = grepl(paste0(paste0("^",use.run.rwt),collapse="|"),rownames(B.mcmc))
			## Combined subset vector
			Rsub  = Rlow & Ruse
		
			Ubase = sum(use.num)  ## may be different than Ubase in 'gather.compo.case.r' if Rlow subsets low-recruitment samples from posterior
	
			## -----------------------------------------------
			## 'refProbs3GenList' contains probabilities for the main
			##   decision tables (usually only use MSY-based RefPts).
			## Will include Bmsy and B0 RefPts already calculated above,
			##   e.g., refProbs3GenList$'0.4Bmsy' - refProbs$LRP = matrix of 0's
			## -----------------------------------------------
			so("findTarget.r","synth")
	
			Bmsy  = B.mcmc$Bmsy
			Bcurr = B.mcmc$Bcurr
			B0    = B.mcmc$B0
			Umsy  = B.mcmc$Umsy
			Ucurr = B.mcmc$Ucurr
		
			## List of targets:
			Tlst = list(
			  list(ratio=0.4,target=Bmsy),
			  list(ratio=0.8,target=Bmsy),
			  list(ratio=1.0,target=Bmsy),
			  list(ratio=1.0,target=Bcurr),
			  list(ratio=0.2,target=B0),
			  list(ratio=0.4,target=B0),
			  ##---COSEWIC---
			  list(ratio=0.5,target=B0),
			  list(ratio=0.7,target=B0),
			  list(ratio=0.3,target=avgTS[,,"Bt"]),
			  list(ratio=0.5,target=avgTS[,,"Bt"])
			)
			names(Tlst)=c("0.4Bmsy","0.8Bmsy","Bmsy","Bcurr","0.2B0","0.4B0","0.5B0","0.7B0","0.3Gen","0.5Gen")
		
			## And have to do separately for u_t < u_MSY:
			Ulst = list(
			  list(ratio=1, target=Umsy),
			  #list(ratio=1, target=Ucurr[[1]]*Usplit[1] + Ucurr[[2]]*Usplit[2])
			  list(ratio=1, target=Ucurr)
			)
			names(Ulst) = c("umsy","ucurr")
			
			## Testing COSEWIC criterion
			## findTarget(Vmat=B.pols[["CC"]][["600"]], target=Tlst[["0.5Gen"]]$target, ratio=Tlst[["0.5Gen"]]$ratio, retVal="p.hi", op=">", yrG=Ngen*gen1)
			## findTarget(Vmat=B.pols[["CC"]][["600"]], target=Tlst[["Bmsy"]]$target, ratio=Tlst[["Bmsy"]]$ratio, retVal="p.hi", op=">", yrG=Ngen*gen1)
			## findTarget(Vmat=B.pols[["CC"]][["600"]], target=Tlst[["Bcurr"]]$target, ratio=Tlst[["Bcurr"]]$ratio, retVal="p.hi", op=">", yrG=Ngen*gen1)
			## findTarget(Vmat=B.pols[["HR"]][["0.08"]], target=Bmsy, ratio=0.4, retVal="N", op=">", yrG=Ngen*gen1, plotit=T)
		
			BRPpList = URPpList = Ttab.conf = Utab.conf = list()
		
			cp = setdiff(dimnames(avgPJ)[[4]], "AC.00")
			cp = cp[cp.num]  ## subset catch policies for decision tables
			CP = round(sapply(avgCP[1,1,],mean,na.rm=T)[cp])  ## should already have been rounded in 'gatherMCMC'
	
			## create table for B > Bmsy (or whatever)
			## ---------------------------------------
			for (k in cp.type){
				BRPpList[[k]] = list()
				for(i in 1:length(cp)) {
					ii = cp[i]
					kk = ifelse (grepl("CC",ii), "CC", "HR")
					Bmat = avgPJ[is.good,as.character(proYears),"Bt",ii]  ## subset years until/if we can project further
					Btab = data.frame(Bcurr,Bmat)
					colnames(Btab)= currYear:projYear
					BRPp.temp = sapply(Tlst,function(x){findTarget(Btab, ratio=x$ratio, target=x$target, retVal="p.hi", op=">", yrG=Ngen*gen1)}, simplify=FALSE)
	#browser();return()
					if (i==1)
						BRPp.collect = BRPp.temp
					else {
						for (j in 1:length(BRPp.temp)){
							jj = names(BRPp.temp)[j]
							BRPp.collect[[jj]] = rbind(BRPp.collect[[jj]], BRPp.temp[[jj]])
						}
					}
				}
				BRPp.collect = lapply(BRPp.collect, function(x,rnam){ rownames(x) = rnam; return(x) }, rnam=CP)
				BRPpList[[kk]] = BRPp.collect
			}
	#browser();return()
	
			## create table for u < umsy (or ucurr)
			## ------------------------------------
			for (k in cp.type){
				URPpList[[k]] = list()
				for(i in 1:length(cp)) {
					ii = cp[i]
					kk = ifelse (grepl("CC",ii), "CC", "HR")
					Umat = avgPJ[is.good,as.character(proYears),"ut",ii]  ## subset years until/if we can project further
					Utab = data.frame(Ucurr,Umat)
					colnames(Utab)= currYear:projYear
					URPp.temp = sapply(Ulst,function(x){findTarget(Utab, ratio=x$ratio, target=x$target, retVal="p.hi", op="<", yrG=Ngen*gen1)}, simplify=FALSE)
					if (i==1) {
						URPp.collect = URPp.temp
					} else {
						for (j in 1:length(URPp.temp)){
							jj = names(URPp.temp)[j]
							URPp.collect[[jj]] = rbind(URPp.collect[[jj]], URPp.temp[[jj]])
						}
					}
				}
				URPp.collect = lapply(URPp.collect, function(x,rnam){ rownames(x) = rnam; return(x) }, rnam=CP)
				URPpList[[kk]] = URPp.collect
			}
	#browser();return()
	
			### Also calculate number of years to reach reference target points with a specified confidence
			### (keep code in case we need it in future)
			#Ttab.temp = Utab.temp = list()
			#for (j in c(0.50,0.65,0.80,0.95)) {
			#	jj = formatC(j,digits=2,format="f")
			#	Ttab.temp[[jj]]  = pmin(sapply(Tlst, function(x){sapply(B.pols[[i]], findTarget, ratio=x$ratio, target=x$target, conf=j,  retVal="N", op=">", yrG=Ngen*gen1)}), Ngen*gen1)
			#	Utab.temp[[jj]]  = pmin(sapply(Ulst, function(x){sapply(U.pols[[i]], findTarget, ratio=x$ratio, target=x$target, conf=j,  retVal="N", op="<", yrG=Ngen*gen1)}), Ngen*gen1)
			#}
			#Ttab.conf[[i]] = Ttab.temp
			#Utab.conf[[i]] = Utab.temp
		} ## end not agile
	
		##=====TABLE: DEFAULT DFO MSY REFERENCE POINTS=====
		## Just use and adapt code from runSweaveMCMC.Rnw
		projYearsNum = length(proYears)
		for.area = paste0(" (",gsub(")","]",gsub("\\(","[",name)),")")
		maxCatSentence = ""  ## only ever used for the ROL assessment (see run-masterMCMC.Snw)
		#Nmcmc = sum(Rsub)
	
		##=====GMU -- Guidance for Setting TAC=====
	
		adjUyrs <- function(tab){
			## adjust the exploitation years to reflect catch before start of current year
			colnames(tab) = as.character(as.numeric(colnames(tab))-1)
			return(tab)
		}
		## Change the label and caption of an agile tab
		adjAgile <- function(atab, tablab, tabcap) {
			apos  = grep("caption", atab); aline = atab[apos]
			aline = sub("label\\{(.+)\\}", paste0("label{", tablab,"}"), aline)
			aline = gsub("caption\\{(.*)\\\\label", paste0("caption{", tabcap,"} \\\\label"), aline)
			atab[apos] = aline
			return(atab)
		}
		## Change the rownames (catch policies) for non-agile decision tables to include big.mark delimiters
		adjRnam <- function(dtab) {
			dtab.adj = dtab
			rownames(dtab.adj) = formatC(as.numeric(rownames(dtab)), digits=0, format="f", big.mark=big.mark)
			return(dtab.adj)
		}
		prefix     = ifelse(is.renso, "low.", prefix)
		prefix.cap = ifelse(is.renso, paste0("Base run", ifelse(Narea>1," subareas","")," (0.5$R$)"), paste0("Base run", ifelse(Narea>1," subareas","")))
		Nmcmc.cap  = format(Nmcmc.good, big.mark=big.mark)
		
		## (RH 251027) Try new approach to collecting and printing numerous tables that are available (not what might be available).
		tabitha = ls(pattern="tabDT")
		tabitha = grep("Gen",tabitha,invert=T,value=T) ## drop the COSEWIC running historical points
		##  "tabDT.0.16B0" "tabDT.0.32B0" "tabDT.0.3Gen" "tabDT.0.5B0"  "tabDT.0.5Gen" "tabDT.0.7B0"  "tabDT.Bcurr" "tabDT.Btrp"   "tabDT.ucurr"  "tabDT.utrp" 
		tabDT.out = list()
		for (i in 1:length(tabitha)) {
			ii   = tabitha[i]
			iii  = sub("tabDT.","",ii)
			ilab  = switch(iii,  ## try to simplify acronyms (ilab[3])
				'0.16B0' = c("limit reference point (LRP) 0.16$B_0$", "P$(B_t > 0.16B_0)$", "LRP"),
				'0.32B0' = c("upper stock reference (USR) 0.32$B_0$", "P$(B_t > 0.32B_0)$", "USR"),
				'Btrp'   = c("target reference point (TRP) 0.4$B_0$", "P$(B_t > 0.4B_0)$",  "TRP"),
				'Bcurr'  = c(paste0("current spawning biomass (BCU) $B_{",currYear,"}$"), paste0("P$(B_t > B_{",currYear,"})$"),  "BCU"),
				'utrp'   = c("reference removal rate (RRR): exploitation rate at TRP 0.4$B_0$", "P$(u_t < u_\\\\text{TRP})$", "RRR"),
				'ucurr'  = c(paste0("current exploitation rate (UCU) $u_{",currYear,"}$"), paste0("P$(u_t < u_{",currYear,"})$"), "UCU"),
				## COSEWIC RPs
				'0.5B0' = c("COSEWIC reference criterion A2 Endangered (A2E) 0.5$B_0$", "P$(B_t > 0.5B_0)$", "A2E"),
				'0.7B0' = c("COSEWIC reference criterion A2 Threatened (A2T) 0.7$B_0$", "P$(B_t > 0.7B_0)$", "A2T"),
				c("sumting", "wong", "XXX") ## can add more if MSY RPs come back into vogue
			)
			eval(parse(text=paste0("itab <- ", ii)))
			client = ifelse(grepl("Gen|0\\.5|0\\.7",iii), "cosewic.", "gmu.")
			tablab = paste0("tab:", prefix, client, ilab[3],".CC")
			tabcap =
			if (ii %in% grep("tabDT\\.[Uu]",tabitha,invert=T,value=T)) {
				eval(paste0(prefix.cap, ": decision table for the ", ilab[1], " featuring current- and ", projYearsNum, "-year projections for a range of constant catch (CC) strategies (in tonnes). Values are ", ilab[2], ", i.e.~the probability of the spawning biomass (mature females) at the start of year $t$ being greater than ", ilab[3], ". ", refCC.sentence, maxCatSentence))
			} else if (ii %in% grep("tabDT\\.[Uu]",tabitha,value=T)) {
				eval(paste0(prefix.cap, ": decision table for the ", ilab[1], " featuring current- and ", projYearsNum, "-year projections for a range of constant catch (CC) strategies (in tonnes). Values are ", ilab[2], ", i.e.~the probability of the annual exploitation rate prior to the start of year $t$ being less than ", ilab[3], ". ", refCC.sentence, maxCatSentence))
			} else {
				"Sumtingwong"
			}
			iout = paste0(client, ilab[3], ".CC.out")
			itab.out = adjAgile(itab$tabfile, tablab=tablab, tabcap=tabcap)
			#tput(gmu.LRP.CCs.out)  ## no need to store it
			tabDT.out[[ilab[3]]][[iout]] <- itab.out
			tabDT.out[[ilab[3]]][["label"]] <- tablab
			tabDT.out[[ilab[3]]][["caption"]] <- tabcap
		}
#browser();return()
		save("tabDT.out", file=paste0(prefix,"decision.tabs", ifelse(useRlow, paste0(".Rlow(q=",pad0(qRlow*100,2),")"),""), ".rda") )
		.flush.cat("Finished building decision tables.","\n")
	
		rm(list=setdiff(ls(),c(keep.pars,"keep.pars")))  ## need to clear the objects in first loop so they don't mask those in second environment
		detach(ee)
	} ## end envo loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tabSS.decision


## tabSS.senso--------------------------2025-12-04
## Make Sensitivity Tables
## Note: u2023=u2022 (see 'gatherMCMC.r') so change 
##       labels here in rfpt tables to use 'prevYear'
## ---------------------------------------------RH
#tabSS.senso <- function(istock="YTR", prefix="ytr.", senso, sigdig=4)
tabSS.senso <- function(envo, sigdig=4)
{
	vomit <- function() { gc(verbose=FALSE); while ("ee" %in% search()) detach(ee) }; #resetGraph() }
	on.exit( vomit() )

	for (e in 1:length(envo)) {
		keep.pars = ls()  ## need to clear environment of objects on subsequent loops through envo
		ee = envo[[e]]
		attach(ee)
		rm(list=intersect(ls(.GlobalEnv), ls(ee)), envir=.GlobalEnv)  ## remove objects in global environment masking objects in e1
		#unpackList(stock[[istock]][["Controls"]])
		#unpackList(senso)
		fleets.all = fleets

		## Need to agglomerate parameters due to fleet offsets
		if (spp.code=="YTR") {
			fleets = c("(TRAWL|BT_BC)","QCS","WCVI","NMFS")
			good = c("^(M_|M2_)Female", "^(M_|M2_)Male", paste0("mu.+",fleets), paste0("varL.+",fleets), paste0("delta1.+",fleets))
			bad  = c("^M1\\_","sigmaR","MW\\_BC","beta6","varR","delta[2-9]","DM\\_theta")
		}
		if (spp.code=="SGR") {
			fleets = c("(TRAWL_5ABC)", "(TRAWL_5DE)", "(TRAWL_3CD)", "QCS", "WCHG", "WCVI")
			good   = c("^(M_|M1_|M2_)Female", "^(M_|M1_|M2_)Male", "BH\\_h", "LN\\(R0\\)", "Rdist\\_area\\([12]", paste0("mu.+",fleets), paste0("varL.+",fleets), paste0("Delta.+",fleets), paste0("LnQ.+",fleets))
			bad    = c("HBLL","DM\\_theta","Rdist\\_area\\(3")
			idx    = 3   ## number of strings separated by "." in rownames of PA
		}

		PA.mpd = mergePA(smpdPA, good=good, bad=bad)
#browser();return()
		PA.mcmc = mergePA(senPA, good=good, bad=bad)
		#unpackList(PA.list, scope="L")
		S.mpd = PA.mpd[["PA"]]
		senPA = PA.mcmc[["PA"]]
	
		## Use same preliminary code from 'plotSS.senso'
		## ---------------------------------------------
		modYrs      = startYear:currYear
		modYrsChar  = as.character(modYrs)
		## Diagnostics for select parameters
		P.names     = colnames(senPA) ## parameter names
		#S.mpd       = smpdPA; storage.mode(S.mpd)="double"     ## For some reason, this matrix is stored as a list (maybe due to NA values?)
		S.runs      = as.numeric(sapply(strsplit(rownames(senPA),"\\."),function(x){x[1]})) ## vector of Runs
		S.run.rwt   = sapply(strsplit(rownames(senPA),"\\."),function(x){paste0(x[1],".",x[2],".",x[3])}) ## vector of Runs and Rwts
		S.run.ord   = unique(S.runs)                           ## unique Run numbers (1st is central run)
		S.run.nmc   = table(S.runs)[as.character(S.run.ord)]   ## vector of sensitivity runs ordered
		S.run.num   = rep(1:length(S.run.nmc), S.run.nmc) - 1  ## 1st run is Central Run
		S.num       = unique(S.run.num)
		Nsens       = length(S.num) - 1
		## Look into subsetting later
		use.run.rwt = unique(S.run.rwt)
		#is.element(P.runs, P.run.ord) ## default use all runs but a subset might be used for management advice
		if (spp.code %in% c("YMR"))
			P.ord=c(1,seq(2,10,2),seq(3,11,2),12,13)
		else
			P.ord = 1:ncol(S.mpd)

		## Create Sensitivity labels
		S.prefix    = paste0("S",pad0(S.num,2)," (R", pad0(S.run.ord,2) ,") ")
		iCR         = grep(paste0(strsplit(exp.run.rwt,split="\\.")[[1]][1:2],collapse="."), unique(B.index))
		S.prefix[1] = paste0("B",iCR, " (R", strsplit(exp.run.rwt,"\\.")[[1]][1], ") ")
		sen.lab     = sen.lab[1:(length(S.num)-1)]
		S.labels    = paste0(S.prefix, gsub("\\_"," ", c("Central Run",sen.lab)))
		S.labels    = sub("\\\\pc", "%", S.labels)  ## just in case
#browser();return()

		P.qnts = calcQs(senPA, ivec=S.run.rwt, ovec=S.run.ord)
	
		##----------Table 1-----------
		## MCMC model parameters table
		##----------------------------
		tabPmed = array(NA, dim=c(length(P.names),length(use.run.rwt)), dimnames=list(P.names,use.run.rwt))
		for (k in 1:length(use.run.rwt)){
			kk = use.run.rwt[k]
			qq = P.qnts[[kk]]
			ii = qq[grep("50%", rownames(qq)),]
			tabPmed[names(ii), kk] = ii
		}
		tabPmed = tabPmed[P.ord,]
		tab.sens.pars = formatCatch(tabPmed,N=sigdig, na="-")
		colnames(tab.sens.pars) = gsub("\\s+","",S.prefix)

		test = names.pars = rownames(tab.sens.pars)
		## Deal with stupid catchabilities first
		zq = grep("LnQ", names.pars)
		if (any(zq)) {
			names.pars[zq] = sub("LnQ\\_","\\\\log q~", names.pars[zq])
			names.pars[zq] = sub("TRAWL_(BC|5ABC|5DE|3CD)", "TRAWL_\\1_", names.pars[zq ])
			names.pars[zq] = sub("\\_$", "", names.pars[zq ])
			names.pars[zq] = sapply(strsplit(names.pars[zq], split="~"),  function(x){ paste0(x[1],"~(\\text{",gsub("\\_","~",x[2]),"})")})
		}
		## Fix up Rdist
		zd = grep("Rdist", names.pars)
		if (any(zd)) {
			oRdist = rownames(tabPmed)[zd]
			oRdist = sub("3",areas[3],sub("2",areas[2],sub("1",areas[1],oRdist)))
			names.pars[zd] = paste0("\\text{", gsub("\\_","~",oRdist), "}")
			names.pars[zd] = sub("area", "", names.pars[zd])
		}
		names.pars = sub("LN\\(R0)", "\\\\log R_{0}", names.pars)
		names.pars = sub("\\_GP1", "", names.pars)
		names.pars = sub("M(1)?\\_Female", "M_{1}~(\\\\text{Female})", names.pars)
		names.pars = sub("M(1)?\\_Male", "M_{2}~(\\\\text{Male})", names.pars)
		names.pars = sub("BH\\_(h)", "h~(\\\\text{BH steepness})", names.pars)
		## Fix selectivity
		zs = grep("mu|varL|[Dd]elta", names.pars)
		if (any(zs)) {
			names.pars[zs] = sub("\\)","}",sub("\\(","_{",gsub("\\_","~", names.pars[zs])))
			names.pars[zs] = sub("~(TRAWL|QCS|WCHG|WCVI|HS|GIG|NMFS)", "~(\\\\text{\\1", names.pars[zs])
			names.pars[zs] = paste0("\\", names.pars[zs], "})")
#browser();return()
			names.pars[zs] = sub("\\\\varL", "\\\\log v_{\\\\text{L}}", names.pars[zs])
		}
		rownames(tab.sens.pars) =  paste(rep("$", nrow(tab.sens.pars)), names.pars,rep("$",nrow(tab.sens.pars)), sep="")
#browser();return()

		sen.leg = paste0("Sensitivity runs: ", paste0("S", formatC(S.num[-1], width=0, format="d", flag="0"),"~= ", gsub("\\&","\\\\&", gsub("\\%","\\\\%", gsub("_"," ",sen.lab))), collapse=", "))
		if (spp.code=="REBS") {  ## need to add Other fishery because it has ages but no CPUE index
			nseries = length(iseries)
			iseries[nseries] = paste0(iseries[nseries],"/Trawl fishery")
			iseries = c(iseries, "Other fishery")
		}
		cap.par = paste0(
			name, "~: median values of MCMC samples for the primary estimated parameters, ",
			"comparing the ", ifelse(NrefM>1,"central","base"), " run to ", Nsens, " sensitivity runs (\\Nmcmc{} samples each). R~= Run, S~= Sensitivity. ",
			"Numeric subscripts other than those for $R_0$ and $M$ indicate the following gear types $g$: ",
			texThatVec(paste0(c(1, match(iseries,fleets)[-1]),"~= ",iseries),simplify=F), ". ", sen.leg
		)
		xtab.sens.pars = xtable(tab.sens.pars, align=paste0("c",paste0(rep("r",Nsens+1),collapse="")),
			label   = paste0("tab:",prefix,"sens.pars"), digits = if (exists("formatCatch")) NULL else sigdig,
			caption = cap.par )
		xtab.sens.pars.out = capture.output( print(xtab.sens.pars, caption.placement="top",
			sanitize.rownames.function=function(x){x},
			add.to.row =list(pos = list(-1), command = c("\\\\[-1.0ex]")),
			size = "\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\footnotesize"
		) )
		#tput(xtab.sens.pars.out)
#browser();return()

sumting = T
if (sumting) {

		## Deal with extra sensitivity parameters
		if (!is.null(PA.mcmc$exPA)) {
			exPA = PA.mcmc$exPA
			Pex.qnts = calcQs(exPA, ivec=substring(rownames(exPA),1,9), ovec=as.numeric(unique(substring(rownames(exPA),1,2))) )
			#test = sapply(Pex.qnts,function(x){x[grep("50\\%",rownames(x)),]})
			qnts3 = sapply(Pex.qnts,function(x){
				z  = grep("^5\\%$|^50\\%$|^95\\%$",rownames(x))
				xx = x[z,,drop=FALSE]
				return(xx)
			})
			ragged = lapply(qnts3, function(x){
				#if(is.null(names(x))) names(x)="sigmaR"  ## hack for YTR
				xx = formatCatch(as.data.frame(x))
				vv = apply(xx,2,function(xxx){
					paste0(xxx[2], " (", xxx[1], ", ", xxx[3], ")")
				})
				return(vv)
			})
			## Start constructing table
			jseries = match(names(ragged), rownames(S.mpd))  ## remove base run
			cap.par2 = paste0(
				name, "~: median values of MCMC samples for remaining estimated parameters from ",
				length(ragged), " sensitivity runs (\\Nmcmc{} samples each). ",
				"R~= Run, S~= Sensitivity, M1~= natural mortality for young fish (ages 1 to 9)."
			)
			ragtab = as.character()
			ragtab = c(
				#"\\setlength{\\tabcolsep}{2pt}",  ## added in Sweave
				"\\begin{table}[!h]",
				"\\centering",
				paste0("\\caption{", cap.par2, "}"),
				paste0("\\label{tab:", prefix, "sens.pars2}"),
				#"\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\small",
				"\\begin{tabular}{ll}",
				"\\hline \\\\ [-1.5ex]",
				"{\\bf Parameter} & {\\bf median (90\\% CI)} \\\\ [1ex]",
				"\\hline \\\\ [-1.5ex]"
			)
			for (i in 1:length(ragged)) {
				ii = names(ragged)[i]
				iii = grep(paste0("R",substring(ii,1,2)), S.labels, value=T)
				iii = gsub("\\&","\\\\&", gsub("\\%","\\\\%", iii))
				ragtab = c(ragtab, paste0(iii, " & \\\\"))
				for (j in 1:length(ragged[[i]])) {
					jj = ragged[[i]][j]
					jjj = gsub("\\_","~",names(jj))
					jjjj = ragged[[i]][[j]]
					ragtab = c(ragtab, paste0(jjj, " & ", jjjj, " \\\\"))
				}
			}
	#browser();return()
			ragtab = c(ragtab, c(
				"\\hline",
				"\\end{tabular}",
				#"\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\normalsize",
				"\\end{table}"
			))
			xtab.sens.pars2.out = ragtab
			tput(xtab.sens.pars2.out)  ## may need fixing
	#browser();return()
		} else {
			xtab.sens.pars2.out = NULL
		}
} ## end sumting

#browser();return()

		## MCMC MSY-based quantities  (RH 251009) we don't have TRP info from the sensitivities
		## -------------------------
		ls2df <- function(x) {
			nn = as.character(substitute(x))
			xx = as.data.frame(x)
			colnames(xx) = paste(nn,colnames(xx),sep="_")
			return(xx)
		}
		#D.qnts = calcQs(senTS[,,"BtB0"], ivec=S.run.rwt, ovec=S.run.ord)
		#U.qnts = calcQs(senTS[,,"ut"], ivec=S.run.rwt, ovec=S.run.ord)
		#R.qnts = calcQs(senTS[,,"Rt"], ivec=S.run.rwt, ovec=S.run.ord)
		#RD.qnts = calcQs(senTS[,,"Rtdev"], ivec=S.run.rwt, ovec=S.run.ord)
	
		if (RPbase=="BMSY") {
			Q.mcmc = data.frame (
				B0         = senRP[,"B0"],
				Bcurr      = senRP[,"Bcurr"],
				Bcurr.B0   = senRP[,"Bcurr"] / senRP[,"B0"],
				ucurr      = senRP[,"ucurr"],
				umax       = apply(senTS[,,"ut"],1,max,na.rm=T), ## for each mcmc sample across the time series
				MSY        = senRP[,"MSY"],
				Bmsy       = senRP[,"Bmsy"],
				LRP        = senRP[,"LRP"],
				USR        = senRP[,"USR"],
				Bcurr.Bmsy = senRP[,"Bcurr"] / senRP[,"Bmsy"],
				Bmsy.B0    = senRP[,"Bmsy"] / senRP[,"B0"],
				umsy       = senRP[,"umsy"],
				ucurr.umsy = senRP[,"ucurr"] / senRP[,"umsy"]
			)
		}
		if (RPbase=="B0") {
			Q.mcmc = data.frame (
				B0         = senRP[,"B0"],
				Bcurr      = senRP[,"Bcurr"],
				Bcurr.B0   = senRP[,"Bcurr"] / senRP[,"B0"],
				ucurr      = senRP[,"ucurr"],
				umax       = apply(senTS[,,"ut"],1,max,na.rm=T), ## for each mcmc sample across the time series
				#Ytrp      = senRP[,"Ytgt"],
				Btrp       = 0.40 * senRP[,"B0"],  ## this should be OK when B0=SSB_1933 (VIRG) but not OK when B0=SSB_1935 (first year)
				LRP        = 0.16 * senRP[,"B0"],
				USR        = 0.32 * senRP[,"B0"],
				Bcurr.Btrp = senRP[,"Bcurr"] / (0.40 * senRP[,"B0"])
				#Btrp.B0    = (0.40 * senRP[,"B0"]) / senRP[,"B0"],  ## just 0.4
				#umsy       = senRP[,"umsy"],
				#ucurr.umsy = senRP[,"ucurr"] / senRP[,"umsy"]
			)
		} 
		Q.mcmc.sens = split(Q.mcmc, S.run.rwt)  ## split rearranges order of Runs alphabetically
		Q.mcmc.sens = Q.mcmc.sens[use.run.rwt]  ## force back to original order
		tabQmed = sapply(Q.mcmc.sens, function(Q){
			sapply(Q, median, na.rm=T)
		})
#browser();return()

		## Extract area-based information if it exists
		if (exists("xsenRP") && length(xsenRP)>1) {
			## colnames(xsenRP) = "Bmsy" "Fmsy" "MSY"  "umsy" "Btgt"  "Ftgt"  "Ytgt"  "utgt"  "LRP"  "USR"  "Vmsy" "Vtgt"  "B0"   "V0"   "pVB"
			## dimnames(xsensTS)[3] =  "B"      "BtBtgt"  "BtBmsy" "D"      "fR"     "R"      "u"      "ututgt"  "utumsy" "V"
			if (RPbase=="BMSY") {
				A.mcmc = data.frame (
					## Needs revision (see below)
					B0         = senRP[,"B0"],
					Bcurr      = senRP[,"Bcurr"],
					Bcurr.B0   = senRP[,"Bcurr"] / senRP[,"B0"],
					ucurr      = senRP[,"ucurr"],
					umax       = apply(senTS[,,"ut"],1,max,na.rm=T), ## for each mcmc sample across the time series
					MSY        = senRP[,"MSY"],
					Bmsy       = senRP[,"Bmsy"],
					LRP        = senRP[,"LRP"],
					USR        = senRP[,"USR"],
					Bcurr.Bmsy = senRP[,"Bcurr"] / senRP[,"Bmsy"],
					Bmsy.B0    = senRP[,"Bmsy"] / senRP[,"B0"],
					umsy       = senRP[,"umsy"],
					ucurr.umsy = senRP[,"ucurr"] / senRP[,"umsy"]
				)
			}
			if (RPbase=="B0") {
				A.collect = list('BC'=tabQmed)
				
				for (a in dimnames(xsenRP)[[3]]) {
#browser();return()
					A.mcmc = data.frame (
						B0         = xsenRP[,"B0",a],
						Bcurr      = xsenTS[,as.character(currYear),"B",a],
						Bcurr.B0   = xsenTS[,as.character(currYear),"B",a] / xsenRP[,"B0",a],
						ucurr      = xsenTS[,as.character(currYear),"u",a],
						umax       = apply(xsenTS[,,"u",a],1,max,na.rm=T), ## for each mcmc sample across the time series
						#Ytrp       = xsenRP[,"Ytgt",a],
						Btrp       = 0.40 * xsenRP[,"B0",a],  ## this should be OK when B0=SSB_1933 (VIRG) but not OK when B0=SSB_1935 (first year)
						LRP        = 0.16 * xsenRP[,"B0",a],
						USR        = 0.32 * xsenRP[,"B0",a],
						Bcurr.Btrp = xsenTS[,as.character(currYear),"B",a] / (0.40 * xsenRP[,"B0",a])
						#Btrp.B0    = (0.40 * senRP[,"B0"]) / senRP[,"B0"],  ## just 0.4
						#umsy       = senRP[,"umsy"],
						#ucurr.umsy = senRP[,"ucurr"] / senRP[,"umsy"]
					)
					A.mcmc.sens = split(A.mcmc, S.run.rwt)  ## split rearranges order of Runs alphabetically
					A.mcmc.sens = A.mcmc.sens[use.run.rwt]  ## force back to original order
					tabAmed = sapply(A.mcmc.sens, function(Q){
						sapply(Q, median, na.rm=T)
					})
					A.collect[[a]] = tabAmed
				} ## end area a loop
			} ## end RPbase 'B0'
			## How to interleave a data frame (Google AI 251016)
			## Add an identifier and row number to each data frame
			interleaved_list <- lapply(seq_along(A.collect), function(i) {
				df <- as.data.frame(A.collect[[i]])
				df$df_id <- rep(i,nrow(df)) # Identifier for the original data frame
				df$row_num <- 1:nrow(df) # Original row number within each data frame
				return(df)
			})
			## Combine all data frames
			combined_df <- do.call(rbind, interleaved_list)
			## Sort by row_num and then by df_id to interleave
			interleaved_df <- combined_df[order(combined_df$row_num, combined_df$df_id), ]
			## Remove the helper columns if not needed
			interleaved_df$df_id <- NULL
			interleaved_df$row_num <- NULL
			tabQmed.original = tabQmed
			tabQmed = interleaved_df
			xareas  = dimnames(xsenRP)[[3]]
			rownames(tabQmed) = sub("3",paste0(".",xareas[3]), sub("2",paste0(".",xareas[2]), sub("1", paste0(".",xareas[1]), rownames(tabQmed) ) ) )
#browser();return()
		} ## end extra area collection

		tab.sens.rfpt = formatCatch(tabQmed,N=sigdig-1)  ## use 3 instead of 4
		colnames(tab.sens.rfpt) = gsub(" +","",S.prefix)
#browser();return()

		names.rfpt = rownames(tab.sens.rfpt)
		## Use routine from 'make.base.tabs.r':
		names.rfpt =
			gsub("\\.(B|u)", "/\\1",
			gsub("_Trawl", "~(\\\\text{trawl})",
			gsub("_Other", "~(\\\\text{other})",
			gsub("trp", "_\\\\text{TRP}",
			gsub("msy", "_\\\\text{MSY}",
			gsub("LRP",  paste0(switch(RPbase, 'BMSY'="0.4B_{\\\\text{MSY}}", 'B0'="0.16B_0")),
			gsub("USR",  paste0(switch(RPbase, 'BMSY'="0.8B_{\\\\text{MSY}}", 'B0'="0.32B_0")),
			gsub("MSY",  paste0("\\\\text{MSY}"),
			gsub("VB",  "V",
			gsub("[Uu]max", "u_\\\\text{max}",
			gsub("ucurr",  paste0("u_{",prevYear,"}"),
			gsub("Bcurr",  paste0("B_{",currYear,"}"),
			gsub("B0", "B_{0}",
			names.rfpt)))))))))))))
		if (exists("areas")) {
			za = grep(paste0(areas,collapse="|"), names.rfpt)
			if (any(za)) {
				apat = paste0("(",paste0(areas,collapse="|"),")")
				names.rfpt[za] = gsub(apat, "~~~\\1", names.rfpt[za])
				names.rfpt[za] = gsub(paste0("(.*)\\.~~~", apat), "~~~\\\\text{\\2}~~\\1", names.rfpt[za])
			}
		}
		rownames(tab.sens.rfpt) =  paste(rep("$",nrow(tab.sens.rfpt)),names.rfpt,rep("$",nrow(tab.sens.rfpt)),sep="")
		#remove.rows = grep("max|0\\.16|0\\.32|~~B\\_\\\\text\\{TRP",rownames(tab.sens.rfpt))
		remove.rows = grep("~~0\\.16|~~0\\.32|~~B\\_\\\\text\\{TRP|~~u\\_\\\\text\\{max",rownames(tab.sens.rfpt))
		if (length(remove.rows) > 0)
			tab.sens.rfpt = tab.sens.rfpt[-remove.rows,]   ## otherwise it removes all rows

		cap.rfpt = paste0(
			name, "~: medians of MCMC-derived quantities from the ", ifelse(NrefM>1,"central","base"), " run and ", Nsens,
			" sensitivity runs (\\Nmcmc{} samples each) from their respective MCMC posteriors. Definitions are: ",
			"$B_0$                = unfished equilibrium spawning biomass (mature females), ",
			#"$V_0$ = unfished equilibrium vulnerable biomass (males and females), ",
			"$B_{", currYear, "}$ = spawning biomass at the end of ",currYear, ", ",
			"$u_{", currYear, "}$ = exploitation rate (ratio of total catch to vulnerable biomass) in the middle of ", currYear, ", ",
			"$u_\\text{max}$      = maximum exploitation rate (calculated for each sample as the maximum exploitation rate from ",
			startYear , " to ", currYear, "), ",
			switch(RPbase, 'BMSY'=c(
				"MSY -- maximum sustainable yield at equilibrium, ",
				"$B_\\text{MSY}$ -- equilibrium spawning biomass at MSY, ",
				"$u_\\text{MSY}$ -- equilibrium exploitation rate at MSY. ",
			), 'B0'=c(
				"$B_\\text{TRP}$ -- equilibrium spawning biomass at target reference point (0.4$B_0$), "
			)),
			"All biomass values are in tonnes. ", sen.leg
		)
#browser();return()
	
		xtab.sens.rfpt = xtable(tab.sens.rfpt, align=paste0("l",paste0(rep("r",Nsens+1),collapse="")),
			label   = paste0("tab:",prefix,"sens.rfpt"), digits = if (exists("formatCatch")) NULL else sigdig,
			caption = cap.rfpt )
		xtab.sens.rfpt.out = capture.output( 
			#print(xtab.sens.rfpt, tabular.environment="longtable", floating=FALSE, caption.placement="top",
			print(xtab.sens.rfpt,  caption.placement="top",
			sanitize.rownames.function=function(x){x},
			#add.to.row=list(pos=list(-1,3,switch(RPbase,'BMSY'=11,'B0'=5)), command=c("\\\\[-1.0ex]", "\\hdashline \\\\[-1.75ex]", "\\hdashline \\\\[-1.75ex]")),
			add.to.row=list(pos=list(-1), command=c("\\\\[-1.0ex]" )),
			#add.to.row = list(pos = list(0), command = 
			#	paste0("\\hline \\endhead", "\\hline",   "\\multicolumn{", ncol(tab.sens.rfpt) + 1, "}{l}",
			#	"{\\footnotesize Continued on next page} ", "\\endfoot ", "\\endlastfoot ")), # Custom longtable commands
			hline.after =  c(-1,0,nrow(xtab.sens.rfpt)),
			size = "\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\footnotesize"
			)
		)
#browser();return()

		## Sensitivity run likelihoods
		## ---------------------------
		cutout = switch(spp.code, 'YTR'="HBLL|BT_|MW_|CPUE_", 'SGR'="Index|AF|Recruit|Total")
#browser();return()
		senLL = senLL[grep(cutout,rownames(senLL),value=T,invert=T),]  ## YTR 2024 : get rid of these rows
		tab.sens.ll = formatCatch(senLL,N=sigdig,na="-",zero="-")
		names.ll = rownames(senLL)
		names.ll = sub("^Run","Sen.Run",names.ll)
		#names.ll = gsub("\\_"," ",names.ll)
		rownames(tab.sens.ll) = names.ll
		LL.senso = t(tab.sens.ll)
		LL.senso[,"Sen.Run"] = sub("\\s+$", "", S.prefix)
		LL.senso = data.frame(Sen.Run=LL.senso[,"Sen.Run"], Label = c("base run", gsub("_"," ",sen.lab)), LL.senso[,-c(1)])
#browser();return()

		xtab.sens.ll = xtable(LL.senso, align=paste0("l", "p{0.6in}p{1.0in}", paste0(rep("p{0.6in}",dim(LL.senso)[2]-2),collapse="")),
			label   = paste0("tab:",prefix,"sens.ll"), digits=NULL, 
			caption = paste0("MPD log likelihood (LL) values reported by ", 
			ifelse(NrefM>1,"central","base"), " and sensitivity runs for survey indices, age composition (AF), recruitment, and total (not all LL components reported here)") )
		xtab.sens.ll.out = capture.output(print(xtab.sens.ll,  include.rownames=FALSE,
			caption.placement="top", sanitize.rownames.function=function(x){x}, 
			add.to.row = list(pos = list(-1,1), command = c("\\\\[-0.5ex]", "\\hdashline \\\\[-1.75ex]")),
			size = "\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\footnotesize"
		) )
		#tput(xtab.sens.ll.out)
#browser();return()
	
		outwithit = c("xtab.sens.pars.out", "xtab.sens.rfpt.out", "xtab.sens.ll.out")
	#browser();return()
		if (!is.null(xtab.sens.pars2.out))
			outwithit = c(outwithit, "xtab.sens.pars2.out")
		save(list=outwithit, file=paste0(prefix,"senso.tabs.rda"))
		## Tables are loaded when 'AppG_Rcode.R' is run (usually once), 
		##  but the tables may be regenerated and no longer sync with e1|e2
		for (i in outwithit) {
			mess = paste0("ee$", i, " <- ", i)
			eval(parse(text=mess))
		}
		collect=c("tab.sens.pars","tab.sens.rfpt")
		for (i in collect) 
			eval(parse(text=paste0("stock[[istock]][[\"Sens\"]][[\"",i,"\"]] = ",i)))

		rm(list=setdiff(ls(),c(keep.pars,"keep.pars")))  ## need to clear the objects in first loop so they don't mask those in second environment
		detach(ee)
	} ## end envo loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tabSS.senso
