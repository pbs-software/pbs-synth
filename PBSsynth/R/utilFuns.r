##==========================================================
## PBS Stock Synthesis utility functions:
## --------------------------------------
## calcEQyield...........Calculate MSY by area/region
## calcMA................Calculate indicators for fits to mean age data.
## calcMSY...............Calculate MSY (very rough equivalent to Awatea function)
## calcStdRes............Calculate standardised residuals for robustified normal likelihood
## convPN................Convert parameter names from SS to Awatea
## cquantile.vec.........Calculate cumulative quantile as a vector
## extract.between.......Extract character strings between two delimiters.
## findTarget............Derive decision tables for reference points, etc.
## gatherMCMC............Gather and combine MCMC (and MPD) values from compo and/or senso runs
## getNpan...............Get panel number when inside a multi-panel plot.
## getSS.rdevs...........Get reruitment deviations from replist.
## getSS.control.........Get select values out of the control file.
## importCor.............Import SS parameter correlations (mod.PBSawatea)
## importEva.............Import SS Hessian eigenvlaues (mod.PBSawatea)
## importPar.............Import all SS parameters (mod.PBSawatea).
## is.numStr.............Check if strings can be converted to numerics.
## load_extra_mcmc.......Load extra MCMC information from Report files generated for every sample.
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

## calcEQyield--------------------------2023-04-04
##  Ian Taylor's suggestion for MSY by area (230404)
##  Works for MPD but not MCMC samples because 
##    'equil_yield' table missing from sample replists.
##  Conversation with PJS 230413: use proportions from
##    area-specific B0 to allocate MSY, etc.
## ------------------------------------------IT|RH
calcEQyield = function(replist, areas=c("5ABC","3CD","5DE"))
{
	equil_yield = replist$equil_yield
	parameters  = replist$parameters

	## Exclude SPRloop number 3 as it has different fleet allocations
	EQY_table   = equil_yield[!is.element(equil_yield$SPRloop,3),]
	msy.idx     = which.max(EQY_table$Tot_Catch)
	MSY_area    = EQY_table[msy.idx, grep("Dead",colnames(EQY_table))]
	MSY_prop    = MSY_area / sum(MSY_area)

	Rdist  = parameters$Value[grep("^RecrDist.+month_1$",parameters$Label)]
	Rdexp  = exp(Rdist)
	Rprop  = Rdexp / sum(Rdexp)

	Rdevnam = grep("^RecrDist.+DEVadd",parameters$Label,value=T)
	yrs     = as.numeric(substring(Rdevnam,nchar(Rdevnam)-3,nchar(Rdevnam)))
	ayrs    = split(yrs, cumsum(c(1, diff(yrs) != 1)))
	uyrs    = .su(yrs)
	Rdevadd = array(0, dim=c(length(uyrs),length(areas)), dimnames=list(year=uyrs,area=areas))
	Rdevval = parameters$Value[grep("^RecrDist.+DEVadd",parameters$Label)]
	names(Rdevval) = yrs
	Adevval = split(Rdevval, cumsum(c(1, diff(yrs) != 1)))
	for (i in 1:length(Adevval)) {
		ival = Adevval[[i]]
		Rdevadd[names(ival),i] = ival
	}
	Rdevse   = c(parameters$Value[grep("^RecrDist.+dev_se$",parameters$Label)], 1)
	Rdevadj  = sweep(Rdevadd,2,Rdevse,"*")
	Rdistyr  = sweep(Rdevadj,2,Rdist,"+")
	Rdexpyr  = exp(Rdistyr)
	Rpropyr  = t(apply(Rdexpyr,1,function(x){x/sum(x)}))
#browser();return()
	return(EQY_area)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcEQyield


## calcMA-------------------------------2023-02-17
## Calculate indicators for fits to mean age data.
## -----------------------------------------PJS|RH
calcMA = function(runs=1, rwts=0:1, overwrite=FALSE, fleets.lab, fleets.af,
   cwd="C:/Users/haighr/Files/GFish/PSARC23/POP/Data/SS/POP2023")
{
	so("plotSS.francis.r","synth")
	require(r4ss)
	mean.age.ind = list()
	for (i in runs){
		ii = pad0(i,2)
		for (j in rwts){
			jj = pad0(j,2)
			setwd(paste0(cwd,"/Run",ii,"/MPD.",ii,".",jj))
			replist=SS_output(dir=getwd(), verbose=F, printstats=F)
			francis = plotSS.francis(replist, "age", fleet=fleets.af, printit=F, plotit=F, png=F, outnam="rubbish", lang="e")  ##bt_cpue, qcs_syn, wcvi_syn, nmfs_tri

			if (overwrite)
				write.csv(francis$agedat, file=paste0("mean.ages.",ii,".",jj,".csv"))
			sum.res.flt = sapply(split(francis$agedat$Std.res,francis$agedat$Fleet),sum)
			sum.res.tot = sum(sum.res.flt)
#browser();return()
			pjs.res.flt = sapply(split(francis$agedat,francis$agedat$Fleet),function(x){sum(abs(x$Obsmn-x$Expmn))})
			pjs.res.tot = sum(pjs.res.flt)
			age.ind = data.frame(Run_Rwt=paste0("R.",ii,".",jj), Fleet=c(fleets.lab[fleets.af],"Total"), Sum_Std_Res=c(sum.res.flt,sum.res.tot), Sum_PJS_Res=c(pjs.res.flt,pjs.res.tot))
			mean.age.ind = rbind(mean.age.ind, age.ind)
		}
	}
	return(mean.age.ind)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcMA


## calcMSY------------------------------2023-03-29
## method = method of fishing (e.g., trawl, other, hook and line, midwater trawl, etc.)
## ------------------------------------------AH|RH
calcMSY = function(replist, strategy, method=1, proj_gears=TRUE)
{
	## Subfunctions -----------------------------------------

	MSY = function(MSY.in)
	{
		##****************************************************
		## Calculates MSY (MSYtab stored in .PBStoolEnv)
		## Modified from Allan Hicks' Awatea code in ModelParams.cpp
		##****************************************************
		dIC.out  = determInitCond(dIC.in=MSY.in); ##//deterministic initial conditions
		unpackList(dIC.out)
		Nmat  = Ninit
		#SBold = SBnew = SB0;
		#Nold  = Nnew  = sum(Nmat);
		theCatch = 0;
		oYDP.start = c(MSY.in, list(Catch=theCatch, Nmat=Nmat, VB=VB0, SSB=SB0, alpha=alpha, beta=beta))

		for (theCatch in seq(StartStrategy, EndStrategy, StepStrategy)) {
			if(StrategyType == 2) {
				harvRate = theCatch; ## because working over harvest rates and theCatch is changed below, save for output later
			}
			oYDP.out = oneYearDetermProj(oYDP.start);
			unpackList(oYDP.out, scope="L")
			Nold     = sum(N2)
			oYDP.in  = c(MSY.in, list(Catch=theCatch, Nmat=N2, VB=VB, SSB=SSB, alpha=alpha, beta=beta))
			oYDP.out = oneYearDetermProj(oYDP.in);
			unpackList(oYDP.out, scope="L")
			Nnew     = sum(N2) #Nnew = sum(Nmat);
			diffN    = abs(Nold-Nnew);
			nproj=1;
#browser();return()
			while(diffN > MSYtol && nproj < MSYmaxIter) {
				Nold = Nnew;
				if(StrategyType == 2) {
					theCatch = harvRate; ## //because working over harvest rates and theCatch is changed by reference
				}
				unpackList(oYDP.out, scope="L")
				oYDP.in  = c(MSY.in, list(Catch=theCatch, Nmat=N2, VB=VB, SSB=SSB, alpha=alpha, beta=beta))
				oYDP.out = oneYearDetermProj(oYDP.in);
				unpackList(oYDP.out, scope="L")
				Nnew     = sum(N2) 
				diffN    = abs(Nold-Nnew);
				nproj = nproj + 1
			}
			ttget(MSYtab)
			iii = MSY.in$nmcmc
			kkk = ifelse(harvRate==0, "0", show0(harvRate,3))  ## theCatch = harvest rate
			MSYtab[iii,,kkk] = c(harvRate*VB,VB,SSB,nproj) ## c("Yield","VB","SB","Niter")
			ttput(MSYtab)
		}
	} ## end subfun MSY

	## Strategy selectivity (uses EndYear selectivity only)
	strategySelect = function(va, u, method, year, Nsexes, Nages)  # void model_parameters::strategySelect() {
	{
		##*******************************************************
		## finds the selectivity for projections and yields
		## fills in strategyva
		## Modified from Allan Hicks' Awatea code in ModelParams.cpp
		##*******************************************************
		unpackList(sS.in)
		strategyva = array(0, dim=c(Nsexes,Nages), dimnames=list(sex=1:Nsexes, age=1:Nages));
		endyear    = as.character(year)
		Nmethods   = length(method)

		sumu = 0;
		## // selectivity at age for the projections
		for (midx in 1:Nmethods) {
			meth = as.character(method[midx])
			for (sex in 1:Nsexes) {
				for (age in 1:Nages) {
					if (proj_gears) {
						## Use end-year harvest rates by method
						strategyva[sex,age] = strategyva[sex,age] + u[meth,endyear] * va[meth,endyear,sex,age];
					} else {
						## Use user-supplied harvest rate proportions
						strategyva[sex,age] = strategyva[sex,age] + strategyuproportion[meth] * va[meth,endyear,sex,age];
					}
				} ## end age loop
			} ## end sex loop
			sumu = sumu + u[meth,endyear];
		} ## end method loop
		if (proj_gears)
			strategyva = strategyva/sumu
		return(strategyva)
	} ## end subfun strategySelect

	determInitCond = function(dIC.in)  ## dmatrix model_parameters::determInitCond(double& VB0, double& SSB0)
	{
		## Modified from Allan Hicks' Awatea code in ModelParams.cpp
		unpackList(dIC.in, scope="L")
		## // Weight and Fecundity for the projection --> EndYear +1 (from setup_variable_parameters)
		wProj    = w[as.character(endyr+1),,]
		fecProj  = fec[as.character(endyr+1),]
		survival = aperm(survVec, c(2,1));
		vuln     = strategyva;
		h        = SR[grep("steep",names(SR))]
		R0       = exp(SR[grep("R0",names(SR))])
		if ( any(grepl("Rprop",names(SR))) )
			R0    = R0 * SR[grep("Rprop",names(SR))]  ## Allocation of recruitment by area
#browser();return()

		Ninit = array(0, dim=c(Nsexes, Nages), dimnames=list(sex=1:Nsexes, age=0:(Nages-1)))  ## Ninit.initialize();
		VB0 = TB0 = 0;

		## From the Appendix E:
		A3 = 3 * A - 1
		for (sex in 1:Nsexes) {
			Na0 = RecFraction[sex] * R0 * exp(-(0:A3) * M[sex])  ## start using App E equations
			NA0 = sum( Na0[(A+1):A3] + (Na0[A3] * exp(-M[sex])) / (1 - exp(-M[sex])) )
			N0  = c( Na0[1:A], NA0)  ## Adjusted for 0-age 
			Ninit[sex,] = N0
			## Difference between reported fecundity and calculated fecundity f=w*m*eggs. Why?
			if (sex==1) { ## females
				SB0 = sum( fec[as.character(startyr),] * N0 )  ## closer to SS3 calculation but doesn't make sense; however f=w*m*eggs) so perhaps is good
			}
			TB0 = TB0 + sum( w[as.character(startyr),sex,] * maturity[,1] * N0 )  ## use same maturity?
			VB0 = VB0 + sum( w[as.character(startyr),sex,] * maturity[,1] * N0 * vuln[sex,] )  ## use same maturity?
		} ## end sex loop

		## //Compute Recruitment Parameters
		alpha = (SB0/R0) * (1-((h-0.2) / (0.8*h)));
		beta  = (5*h-1) / (4*h*R0);

		## Update MSYtab with virgin VB and SSB
		ttget(MSYtab)
		MSYtab[nmcmc,,"Init"] = c(0,VB0,SB0,0)
		ttput(MSYtab)
		return( list(Ninit=Ninit, SB0=SB0, VB0=VB0, TB0=TB0, alpha=alpha, beta=beta) );
	} ## end subfun determInitCond

	oneYearDetermProj = function(oYDP.in) 
	{
		##*******************************************************
		## projects the population one year for a certain catch or harvest rate
		##   if StrategyType is 2, then it is a harvest rate
		## uses same selectivity as projections (strategySelect function)
		## this function was written by Allan Hicks on 3/DEC/07
		## Modified from Allan Hicks' Awatea code in ModelParams.cpp
		##*********************************************************
		unpackList(oYDP.in, scope="L")

		## Weight and Fecundity for the projection --> EndYear +1 (from setup_variable_parameters)
		wProj   = w[as.character(endyr+1),,] #w[EndYear+1];
		fecProj = fec[as.character(endyr+1),]  #fec[EndYear+1];

		## Survival curve
		survVec   = aperm(survVec, c(2,1));  ## to make consistent with code below
		survVec05 = aperm(survVec05, c(2,1));

		## Recruitment
		rec = SSB / (alpha + beta*SSB); # value(SSB/(alpha+beta*SSB));

		if (StrategyType == 2) {  ## hr = harvest rate
			hr = Catch;
			Catch = VB*hr;
		}
		else {
			if(VB>0) { hr=Catch/VB; }
			else     { hr=0; }
			if(hr>1) { hr=1; }
		}

		SSB  = VB = 0;
		Z    = M + hr ## there must be natural mortality as time progresses
		Z    = rep(Z,length(Nsexes))[1:Nsexes] 
		N2   = Nmat
		aa   = 2:(Nages-1)
		aaa  = 2:Nages
		for (sex in 1:Nsexes) {
			N2[sex, aa] = N2[sex, aa-1] * exp(-Z[sex])  ## assumes no time-varying F or M
			N2[sex, Nages] = N2[sex, Nages-1] * exp(-Z[sex]) + N2[sex, Nages] * exp(-Z[sex])
			N2[sex,aaa][N2[sex,aaa]<0] = 0
			VB = VB + sum( N2[sex,aaa] * wProj[sex,aaa] * strategyva[sex,aaa] ) ## * survVec05[sex,Nages];
			if (sex==1) ## females
				SSB = SSB + sum( N2[sex,aaa] * fecProj[aaa] )
			N2[sex,1] = rec * RecFraction[sex];
			N2[sex,1][N2[sex,1]<0] = 0
			VB = VB + sum( N2[sex,1] * wProj[sex,1] * strategyva[sex,1] ) ## * survVec05[sex,Nages];
			if (sex==1) ## females
				SSB = SSB + sum( N2[sex,1] * fecProj[1] )
		} ## end sex loop
		return( list(N2=N2, VB=VB, SSB=SSB) )
	} ## end subfun oneYearDetermProj

	## End subfunctions--------------------------------------

	## Main function calcMSY
	## ======================================================
	if (missing(strategy)) {
		strategy = list(
			StrategyType  = 2,     ## Strategy Type (1=constant catch; 2=harvest rate)
			MSYmaxIter    = 15000, ## Maximum number of iterations
			MSYtol        = 0.01,  ## Tolerance for convergence
			StartStrategy = 0,     ## Start Strategy (lower bound of catch for projections)
			EndStrategy   = 0.200, ## End Strategy (upper bound of catch for projections)
			StepStrategy  = 0.001  ## Step Strategy (interval to use for catch projections)
		)
	}
	unpackList(strategy)
	StartYear  = 1935
	EndYear    = 2023   ## -99 for MSY : if(ProjectYear>EndYear) MSYsetup()  in 'ModelParams.cpp'
	A          = 60     ## Age of plus class

	parameters = replist$parameters
	## Determine recruitment distribution by area
	Rdist      = parameters[grep("RecrDist.+month_1$",parameters$Label),"Value"]
	Rprop      = if (length(Rdist)==0) 1 else exp(Rdist)/sum(exp(Rdist))
	## For now, set U strategy to Rprop
	strategyuproportion = Rprop
	## Grab stock Recruitment pars R0 and h (and perhaps others)
	SRpars     = parameters[grep("R0|steep",parameters$Label),c("Label","Value")]
	SR         = apply(SRpars, 1, function(x) { par=x["Value"]; names(par)=x["Label"]; return(as.numeric(par)) })
	SR         = c(SR, 'Rprop'=sum(Rprop[method]))

	## Exploitation -- limit methods to commercial fleets
	fmeths   = grep("FISHERY",replist$FleetNames)
	fnames   = grep("FISHERY",replist$FleetNames,value=TRUE)
	names(fmeths) = fnames
	Nmethods = length(fmeths)
	exploit  = replist$exploitation
	fyears   = .su(exploit$Yr)  # .su(x) = sort(unique(x))
	u        = array(0, dim=c(length(fmeths),length(fyears)), dimnames=list(method=fmeths,year=fyears));
	for (i in fmeths) {
		ii  = as.character(i)
		iii = names(fmeths)[i]
		for (j in fyears) {
			jj = as.character(j)
			u[ii,jj] = unlist(exploit[exploit$Yr%in%j, iii])
		}
	}

	agesel = replist$ageselex  ## contains various Factors
	agesel.byfac = split(agesel, agesel$Factor)

	## Selectivity
	Asel   = agesel.byfac[["Asel"]]
	Asel   = Asel[is.element(Asel$Fleet, fmeths),]
	ameths = .su(Asel$Fleet)
	ayears = .su(Asel$Yr)   ; Nyears   = length(ayears)
	asexes = .su(Asel$Sex) ; Nsexes   = length(asexes)
	fage   = colnames(Asel)[-(1:grep("Label",colnames(Asel)))]
	ages   = as.numeric(fage); Nages = length(ages)
	va     = array(0, dim=c(Nmethods,Nyears,Nsexes,Nages), dimnames=list(method=ameths, year=ayears, sex=asexes, age=ages));
	for (i in ameths) {
		ii = as.character(i)
		for (j in ayears) {
			jj = as.character(j)
			for (k in asexes) {
				kk = as.character(k)
				va[ii,jj,kk,fage] = unlist(Asel[Asel$Fleet%in%i & Asel$Yr%in%j & Asel$Sex%in%k, fage])
			}
		}
	}

	## Body Weight
	bodywt = agesel.byfac[["bodywt"]]
	bodywt = bodywt[is.element(bodywt$Fleet, fmeths),]
	bmeths = .su(bodywt$Fleet)
	byears = .su(bodywt$Yr)  ; Nyears = length(byears)  ## overwrite previous value
	bsexes = .su(bodywt$Sex) ; Nsexes = length(bsexes)
	fage   = colnames(bodywt)[-(1:grep("Label",colnames(bodywt)))]
	ages   = as.numeric(fage); Nages = length(ages)
	w      = array(0, dim=c(Nmethods,Nyears,Nsexes,Nages), dimnames=list(method=bmeths, year=byears, sex=bsexes, age=ages));
	for (i in bmeths) {
		ii = as.character(i)
		for (j in byears) {
			jj = as.character(j)
			for (k in bsexes) {
				kk = as.character(k)
				w[ii,jj,kk,fage] = unlist(bodywt[bodywt$Fleet%in%i & bodywt$Yr%in%j & bodywt$Sex%in%k, fage])
			}
		}
	}

	## Fecundity (females, no specific fleet)
	fecund = agesel.byfac[["Fecund"]]
	fyears = .su(fecund$Yr)
	fage   = colnames(fecund)[-(1:grep("Label",colnames(fecund)))]
	ages   = as.numeric(fage)
	fec    = array(0, dim=c(length(fyears),length(ages)), dimnames=list(year=fyears, age=ages));
	for (i in fyears) {
		ii = as.character(i)
		fec[ii,fage] = unlist(fecund[fecund$Yr%in%i, fage])
	}

	## Natural Mortality and Survival (assume not time-varying)
	natM   = replist$Natural_Mortality  ## contains M by area and growth pattern (which can be area-specific)
	mareas = .su(natM$Area);        Nareas = length(mareas)
	mgpats = .su(natM$Bio_Pattern); Ngpats = length(mgpats)  ## growth patterns
	msexes = .su(natM$Sex);         Nsexes = length(msexes)
	mage   = colnames(natM)[-(1:grep("Era",colnames(natM)))]
	mages  = as.numeric(mage);      Nages = length(mages)
	Mvec   = array(0, dim=c(Nareas,Nsexes,Nages), dimnames=list(area=mareas, sex=msexes, age=mages));  ## assume one growth pattern for now
	for (i in mareas) {
		ii = as.character(i)
		for (j in msexes) {
			jj = as.character(j)
			Mvec[ii,jj,mage] = unlist(natM[natM$Area%in%i & natM$Sex%in%j & natM$Yr%in%EndYear, mage])
		}
	}
	sumting = apply(Mvec,1:2, function(y) {
		x = as.numeric(names(y)); #browser();return(); 
		S = if (all(diff(y)==0)) exp(-y * (x-x[1])) else exp(-y) ## not sure if this is right if M is calculated by age (or even if SS3 can calculate M by age)
		return(S)
	})
	sumtinghalf = apply(Mvec,1:2, function(y) {
		x = as.numeric(names(y)); #browser();return(); 
		S = if (all(diff(y)==0)) exp(-0.5*y * (x-x[1])) else exp(-0.5*y) ## not sure if this is right if M is calculated by age (or even if SS3 can calculate M by age)
		return(S)
	})
	## I prefer dim order to be age, sex, area 
	## https://stackoverflow.com/questions/10679131/how-to-change-order-of-array-dimensions
	Svec   <- aperm(sumting, c(1,3,2))
	Svec05 <- aperm(sumtinghalf, c(1,3,2))
	Mvec   <- aperm(Mvec, c(3,2,1))

	## Maturity
	control  = list.files(".", pattern="control.+ss$")
	cfile    = r4ss::SS_readctl(control)
	maturity = t(cfile$Age_Maturity)
	rownames(maturity) = sub("Age_","",rownames(maturity))
	colnames(maturity) = "Mat"

	## Objects for functions:
	sS.in      = list(va=va, u=u, method=method, year=EndYear, Nsexes=Nsexes, Nages=Nages)
	strategyva = strategySelect(sS.in)

	MSY.in = list(startyr=StartYear, endyr=EndYear, Nsexes=Nsexes, Nages=Nages, w=w[method,,,], fec=fec, survVec=Svec[,,method], survVec05=Svec05[,,method], strategyva=strategyva, SR=SR, RecFraction = c(0.5,0.5), maturity=maturity, M=Mvec[1,,method])

	## Setup MSY collection array
	U       = seq(StartStrategy, EndStrategy, StepStrategy)
	nProj   = length(U)
	collect = c("Yield","VB","SB","Niter")
	Nmcmc   = 10
	Npad0   = floor(log10(Nmcmc))+1
	MSYtab = array(NA, dim=c(Nmcmc,length(collect),nProj+1), dimnames=list(mcmc=1:Nmcmc, val=collect, U=c("Init",show0(U,3))) )
	ttput(MSYtab) ## store in .PBStoolEnv

	MSY.in = c(MSY.in, list(nmcmc=0))
	
	for (i in 1:1){ #Nmcmc) {  ## replaced by list.files for Report_mce_xxxx.sso
		MSY.in$nmcmc = i
		sumting = MSY(MSY.in)
	}
	ttget(MSYtab)
	msy.idx = which.max(MSYtab[1,"Yield",])
	MSY = MSYtab[1,,msy.idx]

#browser();return()

	return(list(MSYtab=MSYtab, MSY=MSY, msy.idx=msy.idx))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcMSY


## calcStdRes---------------------------2022-11-01
##  This implements standardised residuals for the
##  Awatea implementation of the Fournier robustified
##  normal likelihood for proportions at length. 
##  Based on PJS summary of CASAL doc and ACH change to length.
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

		## For Dirichlet-Multinomial, SS code "SS_write_report.tpl" shows:
		##   show_Pearson = value((ocomp - ecomp) / sqrt(ecomp * (1.0 - ecomp) / nsamp * (nsamp + dirichlet_Parm) / (1. + dirichlet_Parm))); // Pearson for Dirichlet-multinomial using negative-exponential parameterization
		## But this is too complicated to replicated here (for now)
	}
	## Pearson residuals truncated
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


## convPN-------------------------------2023-02-14
##  Convert SS3 parameter names to Awatea (or simpler) names.
##  2023 modification for 3-area POP.
## ---------------------------------------------RH
convPN = function(pnams) {
	cnams =
		sub("RecrDist_GP_(\\d+)_area_(\\d+)_month_(\\d+)","Rdist_area(\\2)",
		sub("\\)_Age_P(\\d+)","_\\1)",
		sub("_steep","_h",
		sub("NatM_break_(\\d+)?_Fem_GP_(\\d+)","M\\1_Female_GP\\2",
		sub("NatM_break_(\\d+)?_Mal_GP_(\\d+)","M\\1_Male_GP\\2",
		sub("NatM_p_1_Fem_GP_(\\d+)|NatM_uniform_Fem_GP_(\\d+)", "M_Female_gp(\\1\\2)",
		sub("NatM_p_1_Mal_GP_(\\d+)|NatM_uniform_Mal_GP_(\\d+)", "M_Male_gp(\\1\\2)",
		sub("Early_RecrDev|Main_RecrDev|Late_RecrDev|ForeRecr","RecrDev",
		sub("(_TRAWL|_OTHER)?_FISHERY(_5ABC|_3CD|_5DE)?","\\1\\2", 
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
		pnams)))))))))))))))))))))))))
		onams   = sapply(strsplit(cnams,"_"), function(x){
		if(length(x)==3 && x[1] %in% c("mu","beta2","varL","varR","beta6")) {
			if (grepl("5ABC|3CD|5DE", x[3])) {
				paste0(x[1], sub("5ABC|3CD|5DE","",x[3]), "_", x[2], "_", sub("\\(\\d+\\)","",x[3]))
#browser();return()
			} else {
				paste0(x[1],x[3],"_",x[2])
			}
		} else {
			paste0(x,collapse="_")
		}
	})
#browser();return()
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
## Chantel Wetzel (Jun 27, 2022)
##   The F value in 2023 would be calculated from the total catch at the end of the model year.
##   If the final year in your model with catches input into the data file is 2022, then the F_2023 would be based off the forecast catch from that year.
## ---------------------------------------------RH
gatherMCMC = function( mcdir=".", type="compo", strSpp="CAR",
   basedir="C:/Users/haighr/Files/GFish/PSARC22/CAR/Data/SS/YMR2022",
   ryrs=1935:2023, pyrs=2024:2033, valTS=c("SSB","F","Recr","RecrDev"), valRP=c("SSB_MSY","annF_MSY","Dead_Catch_MSY"))
{
	ayrs  = c(ryrs, pyrs); nyrs = length(ayrs)
	currYr = rev(ryrs)[1]; byr = cyr = as.character(currYr)
	#prevYr = currYr-1;     byr = as.character(prevYr) ## see notes above
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


## getSS.rdevs-------------------------2023-05-02
##  Get reruitment deviations from replist
## ----------------------------------------r4ss|RH
getSS.rdevs = function (replist, forecast=FALSE, minyr=-Inf, maxyr=Inf) 
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


## load_extra_mcmc----------------------2023-04-14
##  Load extra MCMC information from Report files 
##  generated for every sample.
##  Only use for area-based models at this point.
## ------------------------------------------CG|RH
load_extra_mcmc <- function(dir.mcmc=".", dir.extra="./sso", 
   quants5=c(0.05,0.25,0.5,0.75,0.95),  RC=c(TRUE,FALSE),
   startyr=1935, areas=c("5ABC","3CD","5DE"), Fmethod=1, 
   png=F, pngres=400, PIN=c(9,6), lang="e" ) # posts_file_name="posteriors.sso", ...)
{
	## RC = boolean to use Reports and/or CompReports
	#on.exit(gc(verbose=FALSE))  ## garbage collection is a painfully slow process

	## Subfunctions--------------------------------
	## Extract a table from the report file
	makeRepTab = function(dat, pat1, pat2)
	{
		begin.idx = grep(pat1, dat) + 2
		end.idx   = grep(pat2, dat) - 2
		head.idx  = begin.idx - 1
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
#browser();return()

	# Get the number of Report.sso files in the directory
	dir_list  = dir(dir.extra, full.names = TRUE)
	repfiles  = grep("/Report_mce_.*$", dir_list, value = TRUE)
	compfiles = grep("/CompReport_mce_.*$", dir_list, value = TRUE)
	replist = complist = list()
	#dump = gc(verbose=FALSE)

	## Only do this once for gawd's sake
	if (RC[1]) {  ## Reports files
		if (!file.exists("reps.rda")) {
			.flush.cat(paste0("Reading ", length(repfiles), " 'Reports' (ascii files)"), "\n")
			for (i in 1:length(repfiles))
				replist[[i]] = Kmisc::read(repfiles[i])  ## Kmisc::readlines crashes R
			dump = gc(verbose=FALSE)
			tic(); reps = lapply(replist, function(x) { strsplit(x,split=ifelse(grepl("\\r\\n",x), "\\r\\n", "\\n"))[[1]] }); toc()              ## 66.82 secs
			#tic(); reps = mclapply(replist, function(x) { strsplit(x,split="\\r\\n")[[1]] }); toc()            ## 66.81 secs
			#tic(); reps = lapply(replist, function(x) { stringr::str_split(x, pattern="\\r\\n")[[1]] }); toc() ## 83.22 secs
			dump = gc(verbose=FALSE)
			save("reps", file="reps.rda")
		} else {
			if (!exists("reps", envir=.GlobalEnv)) {
				.flush.cat("Loading Reports object from saved binary\n")
				tic(); load("reps.rda", envir=.GlobalEnv); toc()
			}
		}
	}
	## Only do this once for gawd's sake
	if (RC[2]) {  ## CompReports files
		if (!file.exists("comps.rda")) {
			.flush.cat(paste0("Reading ", length(compfiles), " 'CompReports' (ascii files)"), "\n")
			for (i in 1:length(compfiles))
				complist[[i]] = Kmisc::read(compfiles[i])  ## Kmisc::readlines crashes R
			dump = gc(verbose=FALSE)
			tic(); comps = lapply(complist, function(x) { strsplit(x,split="\\r\\n")[[1]] }); toc() 
			dump = gc(verbose=FALSE)
			save("comps", file="comps.rda")
		} else {
			if (!exists("comps", envir=.GlobalEnv)) {
				.flush.cat("Loading CompReports object from saved binary\n")
				tic(); load("comps.rda", envir=.GlobalEnv); toc()
			}
		}
	}
	## Gather the time series elements
	## Aside: parallel processing much slower than lapply
	# require(snow)
	# cl <- makeCluster(spec=detectCores()-1, type="SOCK")
	# tic(); out=clusterApply(cl, reps, function(x) {makeRepTab(dat=x, pat1="^TIME_SERIES", pat2="^SPR_SERIES")} ); toc() ## 635.52 sec = 10.6 min
	# stopCluster(cl)

	## Make time series object
	if (!file.exists("mcmc.ts.rda")) {
		.flush.cat("Extracting 'mcmc.ts' object from 'reps' object\n")
		tic(); mcmc.ts = lapply(reps, function(x) { makeRepTab(dat=x, pat1="^TIME_SERIES", pat2="^SPR_SERIES") }); toc()  ## ~5 min for 2000 MCMC samples
		save("mcmc.ts", file="mcmc.ts.rda")
	} else {
		if (!exists("mcmc.ts", envir=.GlobalEnv)) {
			.flush.cat("Loading mcmc.ts object from saved binary\n")
			tic(); load("mcmc.ts.rda", envir=.GlobalEnv); toc()
		}
	}
	nmcmc = length(mcmc.ts)
	npad0 = floor(log10(nmcmc))+1

#	## Make MSY object (report 54 missing from Report_mcs*.sso files!)
#	if (!file.exists("mcmc.msy.rda")) {
#		.flush.cat("Extracting 'mcmc.msy' object from 'reps' object\n")
#		tic(); mcmc.msy = lapply(reps, function(x) { makeRepTab(dat=x, pat1="^TIME_SERIES", pat2="^SPR_SERIES") }); toc()  ## ~5 min for 2000 MCMC samples
#		save("mcmc.msy", file="mcmc.msy.rda")
#	} else {
#		if (!exists("mcmc.msy", envir=.GlobalEnv)) {
#			.flush.cat("Loading mcmc.msy object from saved binary\n")
#			tic(); load("mcmc.msy.rda", envir=.GlobalEnv); toc()
#		}
#	}

	## Extract time series values
	if (file.exists("mcmc.ts.sub.rda")) {
		if (!exists("mcmc.ts.sub", envir=.GlobalEnv)) {
			.flush.cat("Loading mcmc.ts.sub object from saved binary\n")
			tic(); load("mcmc.ts.sub.rda", envir=.GlobalEnv); toc()
		}
	} else {
		.flush.cat("Deriving 'mcmc.ts.sub' from 'mcmc.ts'\n")
		mcmc.ts.sub = list()
		#cnams.ts = rev(sort(grep("^SpawnBio$|Recruit|dead\\(B|Hrate",colnames(mcmc.ts[[1]]),value=TRUE)))
		cpats.ts = c("^SpawnBio$", "Recruit", "dead\\(B", "Hrate", "sel\\(B")
		for (i in cpats.ts) {
			ii = grep(i,colnames(mcmc.ts[[1]]),value=TRUE)
			iii = sub(":_[1-9]", "", ii[1])
			ilst = lapply(mcmc.ts, function(x) {x[,ii]})
			if (all(sapply(sapply(ilst,dim),is.null)))
				ilst = lapply(ilst,as.numeric)
			else
				ilst = lapply(ilst, function(x) { apply(sapply(x,as.numeric),1,sum) } )
	
			cols.init = mcmc.ts[[1]][,1:4]
			cols.init[,c(1,2,4)] = sapply(cols.init[,c(1,2,4)], as.numeric)
	
			idat = data.frame(cols.init, do.call("cbind",lapply(ilst, data.frame, stringsAsFactors=F)))
			colnames(idat)[-(1:4)] = paste0("s", pad0(1:nmcmc, npad0))
			zbad = is.na(idat[,-(1:4)])
			idat[,-(1:4)][zbad] = NA
			mcmc.ts.sub[[iii]] = idat
		}
		## Vulnerable Biomass
		VB  = data.frame(cols.init, mcmc.ts.sub[['sel(B)']][,-(1:4)]/mcmc.ts.sub[['Hrate']][,-(1:4)])
		zbad = is.na(VB[,-(1:4)])
		VB[,-(1:4)][zbad] = NA
		mcmc.ts.sub[["VB"]] = VB

		## Gather catch for later
		selB   = mcmc.ts.sub[['sel(B)']]
		zcat   = is.element(selB$Era,c("TIME","FORE"))
		catlst = split(selB[zcat,100], selB[zcat,"Area"])  ## need to make sure that simulation catch is not reduced
		cattab = do.call("cbind", lapply(catlst, data.frame, stringsAsFactors=FALSE))
		colnames(cattab) = areas
		rownames(cattab) = .su(selB$Yr[zcat])
		cattab$total = apply(cattab,1,sum)
#browser();return()

		## Depletion (B/B0, where B0=B1935)
		SB = DB = mcmc.ts.sub[["SpawnBio"]]
		B0 = SB[is.element(SB$Yr,startyr),]  ## use startyr (1935) not VIRG
		for (i in unique(SB$Area)) {
			z  = is.element(SB$Area,i) & is.element(SB$Era,c("VIRG","INIT","TIME","FORE"))
			z0 = is.element(B0$Area,i) & is.element(B0$Era,"TIME")
			SBmat = as.matrix(DB[z,-(1:4)])  ## matrix
			B0vec = unlist(B0[z0,-(1:4)])    ## vector
			DBmat = sweep(SBmat, 2, B0vec, "/")
			DB[z,-(1:4)] = DBmat
		}
		mcmc.ts.sub[["DB"]] = DB
	
		## Fraction recruitment (taken from r4ss::SSplotTimeseries)
		Recr   = mcmc.ts.sub[["Recruit_0"]]
		yvals  = Recr[,-(1:4)]
		yvals2 = as.data.frame(array(NA, dim=c(length(Recr$Yr),nmcmc), dimnames=list(rownames(Recr),paste0("s",pad0(1:nmcmc,npad0)) ) ) )
		for (iyr in 1:nrow(yvals)) {
			y <- Recr$Yr[iyr]
			yvals2[iyr,] <- apply(yvals[Recr$Yr == y,],2,sum)
		}
		Frec <- data.frame(cols.init, yvals/yvals2)
		mcmc.ts.sub[["Frec"]] = Frec
		save(list=c("mcmc.ts.sub","cattab"), file="mcmc.ts.sub.rda")
	}

#	## Extract MSY values
#	if (file.exists("mcmc.msy.sub.rda")) {
#		if (!exists("mcmc.msy.sub", envir=.GlobalEnv)) {
#			.flush.cat("Loading mcmc.msy.sub object from saved binary\n")
#			tic(); load("mcmc.msy.sub.rda", envir=.GlobalEnv); toc()
#		}
#	} else {
#		.flush.cat("Deriving 'mcmc.msy.sub' from 'mcmc.ts'\n")
#		mcmc.ts.sub = list()
#	}

	mcmc = mcmc.ts.sub ## just to save on typing
	options(scipen=10) ## stop displaying scientific notation (at least for the first 10 significant digits)
	##  "SpawnBio"  "Recruit_0" "dead(B)"   "Hrate"     "sel(B)"    "VB"        "DB"        "Frec" 
	## "dead(B) and "sel(B)" don't change in MCMC samples
	plist = qlist = character()
	mcmc.names = c("SpawnBio", "Recruit_0", "Hrate", "VB", "DB", "Frec")
	for (i in mcmc.names) {
		ii = switch(i, 'SpawnBio'="spawning", 'Recruit_0'="recruits", 'Hrate'="harvest", 'VB'="vulnerable", 'DB'="depletion", 'Frec'="frecruit")
		iii = switch(i, 'SpawnBio'="Spawning Biomass (kt)", 'Recruit_0'="Recruits (millions age-0 fish)", 'Hrate'="Harvest Rate (/y)", 'VB'="Vulnerable Biomass (kt)", 'DB'="Depletion (Bt/B0)", 'Frec'="Fraction Recruits")
		iv  = switch(i, 'SpawnBio'="B", 'Recruit_0'="R", 'Hrate'="u", 'VB'="V", 'DB'="D", 'Frec'="fR")
		imcmc = mcmc[[i]]
		if (i == "SpawnBio") {
			## Extract SS3's estimation of VIRG biomass (V0, not vulnerable biomass) for allocation of MSY by area (suggested by PJS 230413)
			## Also calculate ORF's version of B0
			B0.mcmc  = t(imcmc[is.element(imcmc$Era, c("TIME")) & is.element(imcmc$Yr, startyr),-c(1:4)])
			colnames(B0.mcmc) = areas
			B0.qmcmc = apply(B0.mcmc, 2, quantile, quants5, na.rm=TRUE)
			V0.mcmc  = t(imcmc[is.element(imcmc$Era, c("VIRG")),-c(1:4)])
			colnames(V0.mcmc) = areas
			pVB.mcmc = t(apply(V0.mcmc, 1, function(x) { x / sum(x) }))  ## used for allocating MSY
			V0.qmcmc  = apply(V0.mcmc, 2, quantile, quants5, na.rm=TRUE)
			pVB.med = V0.qmcmc['50%',] / sum(V0.qmcmc['50%',])  ## proportion SS3 VIRG B by area using median
			V0.mean = apply(V0.mcmc[,-c(1:4)], 1, mean, na.rm=TRUE)
			pVB.mn = V0.mean / sum(V0.mean)                   ## proportion SS3 VIRG B by area using mean
			plist = c(qlist, c("B0.mcmc","V0.mcmc","pVB.mcmc"))
			qlist = c(qlist, c("B0.qmcmc","V0.qmcmc","pVB.med","pVB.mn"))
		}
		imcmc = imcmc[is.element(imcmc$Era, c("TIME","FORE")),]
		amcmc = split(imcmc, imcmc$Area)
		iyrs  = amcmc[[1]]$Yr
		endyr = rev(amcmc[[1]]$Yr[is.element(amcmc[[1]]$Era, c("TIME"))])[1]  ## first forecast year is actually the endyr of the model
		foryr = rev(amcmc[[1]]$Yr[is.element(amcmc[[1]]$Era, c("FORE"))])[1]
		pmcmc = lapply(amcmc, function(x){
			xx = t(x[,-c(1:4)])
#browser();return()
			colnames(xx) = iyrs
			return(xx)
		})
		qmcmc = lapply(amcmc, function(x){
			xx = apply(x[,-c(1:4)], 1, quantile, quants5, na.rm=TRUE)
			colnames(xx) = iyrs
			return(xx)
		})
		names(pmcmc) = names(qmcmc) = areas
		assign( paste0(iv, ".mcmc"), pmcmc)
		assign( paste0(iv, ".qmcmc"), qmcmc)
		plist = c(plist, paste0(iv, ".mcmc"))
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
	msy.mcmc = dmcmc[,grep("(^SSB|^annF|^Dead).*_MSY$",colnames(dmcmc),value=TRUE)]  ## I believe that F method 1 (Pope's approximation) give discrete F
	colnames(msy.mcmc) = sub("Dead_Catch_MSY","MSY", sub("annF_MSY","Fmsy", sub("SSB_MSY","Bmsy",colnames(msy.mcmc))))
	msy.mcmc[is.na(msy.mcmc)] = NA  ## change NaN to NA

	if (Fmethod==1)
		colnames(msy.mcmc)[grep("Fmsy",colnames(msy.mcmc))] = "umsy"
	else
		msy.mcmc$umsy = 1 - exp(-msy.mcmc$Fmsy)
	msy.mcmc$LRP = 0.4 * msy.mcmc$Bmsy
	msy.mcmc$USR = 0.8 * msy.mcmc$Bmsy
	endcat = cattab[as.character(endyr),"total"]  ## end-year catch for the BC coast
	msy.mcmc$Vmsy = endcat/msy.mcmc$umsy

	MSY.mcmc = MSY.qmcmc = list()
	for (a in areas) {
		paVB = pVB.mcmc[,a]
		amsy.mcmc  = sweep(msy.mcmc,1,paVB,"*")
		## Cannot do this to umsy so have to recalculate from Vmsy
		endcat = cattab[as.character(endyr), a]  ## end-year catch for each area
		amsy.mcmc$umsy = endcat/amsy.mcmc$Vmsy
		MSY.mcmc[[a]]  = amsy.mcmc
		MSY.qmcmc[[a]] = sapply(MSY.mcmc[[a]],quantile,quants5, na.rm=T)
	}
#browser();return()
	msy.names = c("BtLRP","BtUSR","BtBmsy","utumsy")
	for (i in msy.names) {
		ii = paste0(i, ".mcmc")
		iii = paste0(i, ".qmcmc")
		assign(ii, list())
		assign(iii, list())
		for (a in areas) {
			aa = strsplit(i,split="t")[[1]]
			pmess = paste0("num=",aa[1],".mcmc[[\"",a,"\"]]; den=MSY.mcmc[[\"",a,"\"]][,\"",aa[2],"\"]; ",ii,"[[\"",a,"\"]]=sweep(num,1,den,\"/\")")
			eval(parse(text=pmess))
			qmess = paste0(iii,"[[\"",a,"\"]] = apply(",ii,"[[\"",a,"\"]],2,quantile,quants5,na.rm=T)")
			eval(parse(text=qmess))
		}
	}
	plist = c(plist, c("MSY.mcmc", paste0(msy.names,".mcmc")))   ## posteriors
	qlist = c(qlist, c("MSY.qmcmc", paste0(msy.names,".qmcmc"))) ## quantiles

	save(list=plist, file="mcmc.posts.rda")
	save(list=qlist, file="mcmc.quants.rda")
#browser();return()

	plotArea = function(qlist, area, hline=NULL) {
		scale = ifelse (i %in% c("SpawnBio","Recruit_0", "VB"), 1000., 1.)
		ylim  =  range(unlist(qlist[area])/scale); ylim[1] = 0
		acols = character()
		plot(0,0, xlim=xlim, ylim=ylim, type="n", xaxs="i", xaxt="n", xlab="", ylab="", cex.axis=1.2, cex.lab=1.5)
		if (!is.null(hline)) {
			hcol = rep("slategray",length(hline))
			hcol[match(0.4,hline)] = "red"
			hcol[match(0.8,hline)] = "green4"
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
		addLegend(0.98, 0.98, legend=area, xjust=1, yjust=1, lty=1, col=acols, seg.len=3, bty="n")
	} ## end function 'plotArea'

	for (i in c("utumsy", "BtBmsy",mcmc.names)) {
		ii = switch(i, 'SpawnBio'="spawning", 'Recruit_0'="recruits", 'Hrate'="harvest", 'VB'="vulnerable", 'DB'="depletion", 'Frec'="frecruit", 'BtBmsy'="BtBmsy", 'utumsy'="utumsy")
		iii = switch(i, 'SpawnBio'="Spawning Biomass (kt)", 'Recruit_0'="Recruits (millions age-0 fish)", 'Hrate'="Harvest Rate (/y)", 'VB'="Vulnerable Biomass (kt)", 'DB'="Depletion (Bt/B0)", 'Frec'="Fraction Recruits", 'BtBmsy'=expression(italic(B)[italic(t)]/italic(B)[MSY]), 'utumsy'=expression(italic(u)[italic(t)]/italic(u)[MSY]))
		iv  = switch(i, 'SpawnBio'="B", 'Recruit_0'="R", 'Hrate'="u", 'VB'="V", 'DB'="D", 'Frec'="fR", 'BtBmsy'="BtBmsy", 'utumsy'="utumsy")

		## Plot quantiles for each factor by area
		xlim  = range(iyrs)
		#ylim  =  range(unlist(qmcmc))
		pyrs  = startyr:endyr
		fyrs  = endyr:foryr
		qmcmc = get(paste0(iv,".qmcmc"))

		fout.e = paste0("extra.mcmc(", ii, ")")
		for (l in lang) {
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			if (png) png(filename=paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			expandGraph(mfrow=c(4,1), mar=c(0,2,0,1), oma=c(4,2.5,1,0), mgp=c(1.6,0.75,0))
			hline = if (i=="BtBmsy") c(0.4,0.8)
				else if (i=="utumsy") 1
				else if (i=="depletion") c(0.2,0.4)
				else NULL
			plotArea(qlist=qmcmc, area=areas[1], hline=hline)
			plotArea(qlist=qmcmc, area=areas[2], hline=hline)
			plotArea(qlist=qmcmc, area=areas[3], hline=hline)
			plotArea(qlist=qmcmc, area=areas)
			if (png) dev.off()
		}; eop()
#browser();return()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~load_extra_mcmc


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


## prepMPD------------------------------2023-05-18
##  Prepare MPD runs for POP likelihood analysis.
## ---------------------------------------------RH
prepMPD = function(run.rwt, #n.run, n.rwt, l.run, l.rwt,
	d.base = "C:/Users/haighr/Files/GFish/PSARC23/POP/Data/SS/POP2023",
	w=NULL, cvpro=NULL, modify=TRUE)
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
	if (modify) {
		.flush.cat(paste0("***Edit the ss files in\n\t",d.mpd,"\nbefore proceeding with an MPD fit in SS."),"\n\n")
	} else {
		##file.copy sumting
	}
#	if (!is.null(w) || !is.null(cvpro)){
#	}
	setwd(d.mpd)
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
##  Switched to Chris Grandin's 'run_adnuts' function on his Linux
##    server for the 2023 POP assessment. (RH 230323)
## ------------------------------------------CG|RH
runADnuts=function(path=getwd(), run="01.01", model="base", 
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

