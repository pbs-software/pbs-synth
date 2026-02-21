##================================================2025-07-31
## PBS Stock Synthesis utility functions:
## --------------------------------------
## calcEQyield...........Calculate MSY by area/region
## calcMA................Calculate indicators for fits to mean age data.
## calcMSY...............Calculate MSY (very rough equivalent to Awatea function)
## calcQs................Calculate quantiles from MCMC posteriors
## calcStdRes............Calculate standardised residuals for robustified normal likelihood
## convPN................Convert parameter names from SS to Awatea
## doRetros..............Do retrospective analyses (wrapper for the r4ss' function 'retro')
## extract.between.......Extract character strings between two delimiters.
## findTarget............Derive decision tables for reference points, etc.
## fRcalc................Calculate fraction allocation of R0 by area (and year)
## importCor.............Import SS parameter correlations (mod.PBSawatea)
## importEva.............Import SS Hessian eigenvlaues (mod.PBSawatea)
## importPar.............Import all SS parameters (mod.PBSawatea).
## mergePA...............Agglomerate parameters due to fleet offsets
## prepCP................Prepare Catch Policies -- 'CC'=constant catch
## prepMPD...............Prepare MPD runs for likelihood analysis
## quickDT...............Produce quick decision tables for Bt/LRP, Bt/USR, etc.
## repeatMPD ............Repeat MPDs for axes of uncertainty to visualise likelihood density
## runSweave.............Run Sweave code to build pdfs for MPD and MCMC runs (not appendix)
## tabDQs................Tables of Derived Quantities
## weightAF..............Weight age frequencies using harmonic mean ratio method
##==========================================================

## calcEQyield--------------------------2023-04-04
##  Ian Taylor's suggestion for MSY by area (230404)
##  Works for MPD but not MCMC samples because 
##    'equil_yield' table missing from sample replists.
##  Conversation with PJS 230413: use proportions from
##    area-specific B0 to allocate MSY, etc.
## ------------------------------------------IT|RH
calcEQyield <- function(replist, areas=c("5ABC","3CD","5DE"))
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
	#return(EQY_area)  ## this has not been defined !?
	EQY_area = list(EQY_table=EQY_table, MSY_area=MSY_area, MSY_prop=MSY_prop, Rdexpyr=Rdexpyr, Rpropyr=Rpropyr)  ## (RH 230822) may need changing
	return(EQY_area)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcEQyield


## calcMA-------------------------------2024-07-16
## Calculate indicators for fits to mean age data.
## -----------------------------------------PJS|RH
calcMA <- function(runs=1, rwts=0:1, vers=1, overwrite=FALSE, fleets.lab, fleets.af,
   cwd="C:/Users/haighr/Files/GFish/PSARC24/YTR/Data/SS3/YTR2024", use.cwd.as.is=FALSE)
{
	so("plotSS.francis.r","synth")
	#require(r4ss)
	mean.age.ind = list()
	for (i in runs) {
		ii = pad0(i,2)
		for (j in rwts) {
			jj = pad0(j,2)
			for (k in vers) {
				kk =paste0("v",k)
				if (!use.cwd.as.is)
					setwd(paste0(cwd, "/Run", ii, "/MPD.", ii, ".", jj, ".", kk))
				replist=SS_output(dir=getwd(), verbose=F, printstats=F)
				francis = plotSS.francis(replist, "age", fleet=fleets.af, printit=F, plotit=F, png=F, outnam="rubbish", lang="e")  ##bt_cpue, qcs_syn, wcvi_syn, nmfs_tri

				if (overwrite)
					write.csv(francis$agedat, file=paste0("mean.ages.",ii,".",jj,".",kk,".csv"))
				sum.res.flt = sapply(split(francis$agedat$Std.res,francis$agedat$Fleet),sum)
				sum.res.tot = sum(sum.res.flt)
#browser();return()
				pjs.res.flt = sapply(split(francis$agedat,francis$agedat$Fleet),function(x){sum(abs(x$Obsmn-x$Expmn))})
				pjs.res.tot = sum(pjs.res.flt)
				age.ind = data.frame(Run_Rwt_Ver=paste0("R.",ii,".",jj,".",kk), Fleet=c(fleets.lab[fleets.af],"Total"), Sum_Std_Res=c(sum.res.flt,sum.res.tot), Sum_PJS_Res=c(pjs.res.flt,pjs.res.tot))
				mean.age.ind = rbind(mean.age.ind, age.ind)
			}  ## end i loop
		}  ## end j loop
	}  ## end k loop
	return(mean.age.ind)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcMA


## calcMSY------------------------------2023-03-29
## method = method of fishing (e.g., trawl, other, hook and line, midwater trawl, etc.)
## ------------------------------------------AH|RH
calcMSY <- function(replist, strategy, method=1, proj_gears=TRUE)
{
	## Subfunctions -----------------------------------------

	MSY <- function(MSY.in)
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
	strategySelect <- function(va, u, method, year, Nsexes, Nages)  # void model_parameters::strategySelect() {
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

	determInitCond <- function(dIC.in)  ## dmatrix model_parameters::determInitCond(double& VB0, double& SSB0)
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

	oneYearDetermProj <- function(oYDP.in) 
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


## calcQs-------------------------------2024-07-11
## Calculate quantiles from MCMC posteriors
##   ivec = vector 'run.rwt.ver' with elements for every record
##   ovec = vector model run, one value for each model
## ---------------------------------------------RH
calcQs <- function (dat, ivec, ovec)
{
	fnam = as.character(substitute(dat))
	F.data = as.data.frame(dat)
	F.list = split(F.data, ivec)
	rr.ord = match(ovec, as.numeric(substring(names(F.list),1,2)) )  ## order not so important for sensitivities
	F.list = F.list[rr.ord]
	## lapply does not sort like split does
	F.qnts = lapply(F.list,function(x){
		z = apply(x,2,function(xx){!all(is.na(xx))})
		out = apply(x[,z,drop=FALSE],2,quantile,quants5,na.rm=T)
		return(out)
	})
	return(F.qnts)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcQs


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
	Z <- function (x, r=0.001){
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


## convPN-------------------------------2025-07-17
##  Convert SS3 parameter names to Awatea (or simpler) names.
## ---------------------------------------------RH
convPN <- function(pnams)
{
	cnams =
		sub("\\.$|_$", "",                                                             ## get rid of trailing delimiters
		sub("(mu|varL|varR)(.+)(_)?\\.(\\d+)\\.","\\1(\\4)\\2",                        ## try to fix scan of posteriors.sso
		sub("RecrDist_GP_(\\d+)_area_(\\d+)_month_(\\d+)","Rdist_area(\\2)",
		sub("\\)_Age_P(\\d+)","_\\1)",
		sub("_steep","_h",
		sub("NatM_break_(\\d+)?_Fem_GP_(\\d+)","M\\1_Female_GP\\2",
		sub("NatM_break_(\\d+)?_Mal_GP_(\\d+)","M\\1_Male_GP\\2",
		sub("NatM_p_1_Fem_GP_(\\d+)|NatM_uniform_Fem_GP_(\\d+)", "M_Female_GP\\1\\2",
		sub("NatM_p_1_Mal_GP_(\\d+)|NatM_uniform_Mal_GP_(\\d+)", "M_Male_GP\\1\\2",    ##  keep GP in case multiple growth patterns used
		sub("Early_RecrDev|Main_RecrDev|Late_RecrDev|ForeRecr","RecrDev",
		sub("(_TRAWL|_OTHER)?_FISHERY(_BC|_5ABC|_3CD|_5DE)?","\\1\\2", 
		sub("SYNOPTIC","",
		sub("HISTORIC(AL)?","",
		sub("TRIENNIAL","",
		sub("HBLL_NORTH","HBLLN_",
		sub("HBLL_SOUTH","HBLLS_",
		sub("^AgeSel_(\\d+)MaleatDogleg","delta3(\\1)",       ## (RH 250717) SGR 2025 : selectivity offset
		sub("^AgeSel_(\\d+)MaleatMaxage","delta4(\\1)",       ## (RH 250717) SGR 2025 : selectivity offset
		sub("^AgeSel_(\\d+)?(Male|Fem)_Scale","Delta5(\\1)",  ## parameter offset Delta 5
		sub("^AgeSel_(\\d+)?(Male|Fem)_Final","Delta4(\\1)",  ## parameter offset Delta 4
		sub("^AgeSel_(\\d+)?(Male|Fem)_Descend","Delta3(\\1)",## parameter offset Delta 3
		sub("^AgeSel_(\\d+)?(Male|Fem)_Ascend","Delta2(\\1)", ## parameter offset Delta 2
		sub("^AgeSel_(\\d+)?(Male|Fem)_Peak","Delta(\\1)",    ## parameter offset Delta 1
		sub("^Age_DblN_end_logit","final",                    ## beta 6
		sub("^Age_DblN_top_logit","plateau",                  ## beta 2
		sub("^Age_DblN_descend_se","varR",                    ## beta 4
		sub("^Age_DblN_ascend_se","varL",                     ## beta 3
		sub("^Age_DblN_peak","mu",                            ## beta 1
		sub("^SR_","",
		pnams)))))))))))))))))))))))))))))
#print(pnams); cat("\n"); print(cnams)
#browser();return()
		onams   = sapply(strsplit(cnams,"_"), function(x){
		if(length(x)==3 && x[1] %in% c("mu","plateau","varL","varR","final")) {   ## note: no conversion yet for beta 5 (initial)
			if (grepl("^BC|^5ABC|^3CD|^5DE", x[3])) {
				paste0(x[1], sub("^BC|^5ABC|^3CD|^5DE","",x[3]), "_", x[2], "_", sub("\\(\\d+\\)","",x[3]))
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


## doRetros-----------------------------2025-10-14
##  Basically a wrapper for the r4ss' function 'retro'.
## ---------------------------------------------RH
doRetros <- function(strSpp="405", assyr=2025, stock="BC", newsubdir="retros",
   exe="C:/Users/haighr/Files/Archive/Bat/ss.exe")
{
	if (strSpp %in% c("437","CAR") && assyr %in% c(2022)){
		basedir = "C:/Users/haighr/Files/GFish/PSARC/PSARC_2020s/PSARC22/CAR/Data/SS/CAR2022"
		years=0:-10; run=24; rwt=1
	}
	if (strSpp %in% c("396","POP") && assyr %in% c(2023)){
		basedir = "C:/Users/haighr/Files/GFish/PSARC23/POP/Data/SS/POP2023"
		switch(stock,
			'5ABC' = { years=0:-13;  run=24; rwt=1; ver="1"},  ## POP 2023 -- go back to 2010 for all areas
			'3CD'  = { years=0:-13; run=25; rwt=1; ver="1" },  ## POP 2023
			'5DE'  = { years=0:-13; run=26; rwt=1; ver="1" }   ## POP 2023
		)
	}
	if (strSpp %in% c("418","YTR") && assyr %in% c(2024)){
		basedir = "C:/Users/haighr/Files/GFish/PSARC24/YTR/Data/SS3/YTR2024"
		switch(stock,
			#'BC'   = { years=0:-10;  run=2; rwt=1; ver="1"}  ## YTR 2024 RPR (with fecundity error)
			'BC'   = { years=0:-10;  run=2; rwt=1; ver="2"}  ## YTR 2024 revised post RPR
		)
	}
	if (strSpp %in% c("405","SGR") && assyr %in% c(2025)){
		basedir = "C:/Users/haighr/Files/GFish/PSARC25/SGR/Data/SS3/SGR2025"
		switch(stock,
			'3area'   = { years=0:-10;  run=28; rwt=1; ver="2"},  ## SGR 2025 multi-area model
			'coast'   = { years=0:-10;  run=29; rwt=1; ver="2"}   ## SGR 2025 single-area model
		)
	}
	## Generic for all species
	mpddir = paste0(basedir, "/Run", pad0(run,2), "/MPD.", pad0(run,2), ".", pad0(rwt,2), ".v", ver )
	retdir = paste0(basedir, "/Run", pad0(run,2), "/Retro.", pad0(run,2), ".", pad0(rwt,2), ".v", ver )
#browser();return()
	if (!dir.exists(retdir)) 
		dir.create(retdir)
	## Grab the input files:
	inputs = list.files(mpddir, pattern="(^control|^data|^forecast|^starter).*\\.ss$")
	file.copy(from=paste0(mpddir,"/",inputs), to=retdir, overwrite=T, copy.date=T)
#browser();return()
	retro(dir=retdir, years=years, newsubdir=newsubdir, exe=exe)
	retroModels = SSgetoutput( dirvec=file.path(retdir, newsubdir, paste("retro",years,sep="") ) )
	save("retroModels","strSpp","assyr","stock","years","run","rwt","retdir", file=paste0(retdir,"/retroModels.rda"))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~doRetros


## extract.between----------------------2020-09-17
##  Extract character strings between two delimiters; based on:
##  https://stackoverflow.com/questions/8613237/extract-info-inside-all-parenthesis-in-r
## ---------------------------------------------RH
extract.between <- function(x, open="{", close="}", first.only=TRUE)
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
findTarget <- function(Vmat, yrU=as.numeric(dimnames(Vmat)[[2]]), yrG=90, 
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


## fRcalc ------------------------------2025-07-29
##  Calculate fraction allocation of R0 by area (and year)
##  Note: need to run convPN() on colnames of d.mcmc
## ---------------------------------------------RH
fRcalc = function(d.mcmc, narea)
{
	calcParea = function(Rdist, Rdev, RdevSE) {
		Narea  = length(Rdist)
		Rdist.adj = Rdist + Rdev * RdevSE
		for (a in 1:Narea) {
			parea = exp(Rdist.adj[,a]) / t(t(apply(exp(Rdist.adj),1,sum)))
			if (a==1) Parea = parea
			else      Parea = cbind(Parea, parea)
		}
		colnames(Parea) = gsub("\\_DEVadd","",gsub("Rdist","Rprop",colnames(Rdev)))
		return(Parea)
	}
	fR.mcmc = NA
	if (narea>1) {
		fR = d.mcmc[,grep("Rdist\\_area\\(\\d)$", colnames(d.mcmc),value=T)]
		missRdist = setdiff (1:narea, gsub("[^0-9.-]", "", colnames(fR)) )
		fR[,gsub("\\d", missRdist, colnames(fR)[1])] = rep(0,nrow(fR))
		RdistDev = setdiff(grep("Rdist",colnames(d.mcmc),value=T), colnames(fR))
		if (any(grepl("dev_se",RdistDev))) {
			## SE for Rdist has been estimated
			RdistDevSE = grep("dev_se", RdistDev, value=TRUE)
			RdistDev   = setdiff(RdistDev, RdistDevSE)
			fR[,RdistDevSE] = d.mcmc[,RdistDevSE]
			fR[,gsub("\\d", missRdist, colnames(fR)[4])] = rep(1,nrow(fR))
		}
		RdistYrs = .su(as.numeric(revStr(substring(revStr(RdistDev),1,4))))
		fR.names = colnames(fR)[1:narea]
		## create a weird data.frame with every year's DEV as a column (wtf?)
		for (i in RdistYrs) {
			dev.name = paste0(fR.names,"_DEVadd_",i)
			fR[,dev.name] = 0
			for (j in dev.name) {
				if(!j %in% colnames(d.mcmc)) next
				fR[,j] = d.mcmc[,j]
			}
		}
		devSE.name = paste0(fR.names,"_dev_se")
		if (all(devSE.name %in% colnames(fR))) {
			iRdevSE    = fR[,devSE.name]
		} else {
			if (exists("replist")) {
				fixRdevSE = c(replist$parameters[grep("dev_se", rownames(replist$parameters)),"Value"], 1)
			} else {
				fixRdevSE = rep(1, narea)
			}
			iRdevSE   = as.data.frame(matrix(rep(fixRdevSE, nrow(fR)), nrow=nrow(fR), ncol=narea, dimnames=list(1:nrow(fR),devSE.name),byrow=T))
		}
		for (i in RdistYrs) {
			dev.name = paste0(fR.names,"_DEVadd_",i)
			iRdist   = fR[,fR.names]
			iRdev    = fR[,dev.name]
			parea    = calcParea(Rdist=iRdist, Rdev=iRdev, RdevSE=iRdevSE)
			fR       = cbind(fR, parea)
		}
		fR.mcmc = fR
	} ## end calculation when narea > 1
	return(fR.mcmc)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~fRcalc


## importCor----------------------------2020-09-28
## Import SS parameter correlations.
##  (Adapted from 'importCor' in PBSawatea.)
## ---------------------------------------------RH
importCor <- function(cor.file)
{
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


## importEva----------------------------2020-09-28
## Import SS Hessian eigenvlaues.
##  (Adapted from 'importEva' in PBSawatea.)
## ---------------------------------------------RH
importEva <- function(eva.file)
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


## importPar----------------------------2020-09-24
##  Import all SS parameters.
##  (Adapted from 'importPar' in PBSawatea.)
## ---------------------------------------------RH
importPar <- function(par.file) 
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


## importStd----------------------------2020-09-28
## Import SS table of standard deviations.
##  (Adapted from 'importStd' in PBSawatea.)
## ---------------------------------------------RH
importStd <- function(std.file, vnam="name") 
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


## mergePA------------------------------2025-10-30
## Agglomerate parameters due to fleet offsets
## Revised for SGR becaue previous algorith a hot mess
## ---------------------------------------------RH
## need to check if any bad entries (parameter stragglers)
## idx is the number of strings separated by "." in rownames of PA to 
## obtain a unique run identifier (e.g., idx=2 for POP [run.rwt], idx=3 for YTR [run.rwt.ver])
mergePA <- function(PA, good, bad, idx=3)
{
	colnames(PA) = gsub("\\_GP[0-9]$", "", gsub("M_(Female|Male)(\\_GP[0-9])$", "M1_\\1\\2", colnames(PA)))
	pnams = colnames(PA)
	exIdx = grep(paste0(bad,collapse="|"), pnams)
	exPA = PA[,exIdx]  ## put occasional parameters into separate table
	inPA = PA[,-exIdx] ## then remove them from the main PA table
#browser();return()
	for (i in 1:length(good)) {
		pnams = colnames(inPA)  ## this changes every interation
		ii = good[i]
		iii = grep(ii, pnams, value=TRUE)
		if (length(iii)==0) next
		pdat = inPA[,iii,drop=FALSE]
		if (i==1) {
			outPA = pdat
		} else {
			outPA = cbind(outPA, pdat)
		}
#if (i==3) {browser();return()}
		#cnam = unique(gsub("\\([0-9]\\)","",iii))  ## wtf?
		#if (length(cnam)>1) { .flush.cat("More than one P name; using the first: ", cnam[1], "\n"); cnam=cnam[1] }
		#psum = matrix(rowSums(pdat,na.rm=T), ncol=1, dimnames=list(rownames(pdat), cnam))
		#PA = PA[,setdiff(colnames(PA), iii)]
		#PA = cbind(PA,psum)
	}
#browser();return()
	## Clean up the main PA table
	if (spp.code=="YTR")
		colnames(PA) = gsub("\\_BC$", "", gsub("^delta1", "delta", colnames(PA)))
	if (spp.code %in% c("SGR")) {
		colnames(outPA) = gsub("LnQ\\_base", "LnQ", colnames(outPA))
		colnames(outPA) = gsub("Rdist\\_area\\(1\\)", "Rdist_5ABC", colnames(outPA))
		colnames(outPA) = gsub("Rdist\\_area\\(2\\)", "Rdist_5DE", colnames(outPA))
		colnames(outPA) = gsub("Rdist\\_area\\(3\\)", "Rdist_3CD", colnames(outPA))
		colnames(outPA) = gsub("\\([0-9]\\)", "", colnames(outPA))
	}

	## Clean up the excised (low occurrence) PA table  (still needs to be revised)
	exPA = exPA[rowSums(exPA,na.rm=T)>0,,drop=F]
	unam = unique(substring(rownames(exPA),1,9))
	unam = unique(sapply(strsplit(rownames(exPA),"\\."),function(x){paste0(x[1:idx],collapse=".")}))  ## seems unnecessary
	if (length(unam)==nrow(exPA)) {  ## MPD agglomeration could likely use the MCMC routing below
		idex  = apply(exPA,2,function(x){(1:length(x))[!is.na(x)][1]})
		exPA = exPA[,order(idex), drop=FALSE]
	} else {  ## manipulations for MCMCs
		rlab  = unam
		clab  = colnames(exPA)
		sexPA = array(NA, dim=c(length(rlab),length(clab)), dimnames=list(rlab,clab))
		for (i in 1:length(rlab)) {
			ii = rlab[i]
			iii = grep(ii,rownames(exPA))
			ival = apply(exPA[iii,,drop=FALSE],2,median,na.rm=T)
			sexPA[ii,] = ival
		}
		idex  = apply(sexPA,2,function(x){(1:length(x))[!is.na(x)][1]})
#browser();return(); 
		exPA = exPA[,order(idex)]
	}
	#colnames(exPA) = gsub("\\([0-9]\\)","",colnames(exPA))
#browser();return()
	return(list(PA=outPA,exPA=exPA))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~mergePA


## prepCP-------------------------------2025-10-22
##  Prepare Catch Policies
##    'CC'=constant catch, 'HR'=harvest rate (not implemented)
##  Note: 'fyrs' (forecast years) includes the current year 
##    (e.g., 2024=beginning of 2024, end of 2023), which is treated as a projection
##  Updated for:
##    Canary Rockfish 2022 (RH 220630)
##    Pacific Ocean Perach (RH 230713)
##    Yellowtail Rockfish  (RH 240704)
##    Silvergray Rockfish  (RH 251022)
## ---------------------------------------------RH
prepCP <- function(run.rwt, ver="", cp=list('BC'=4000), d.cp="CC", tag="", #tag=".nuts4K", 
	fyrs=2025:2035, season=1, fleet=1,
	d.base = "C:/Users/haighr/Files/GFish/PSARC24/YTR/Data/SS3/YTR2024",
	w=NULL, cvpro=NULL, linux=TRUE, onlinux=FALSE)
	## linux=T when MCMCs created on Linux and transferred to Windows
{
	padded = pad0(run.rwt,2)
	t.run=padded[1]; t.rwt=padded[2]
	tar.run.rwt = paste(t.run,t.rwt,sep=".")
	tar.run.rwt.ver = paste0(tar.run.rwt, ".v", ver)
	if (onlinux) {
		d.target = file.path(d.base, paste0("Run", t.run, "v", ver), "mcmc")
	} else {
		d.target = file.path(d.base, paste0("Run", t.run), paste0("MCMC.",tar.run.rwt.ver,tag))
	}
	poop <- function(){ setwd("/srv/haigh/SGR2025/R/"); gdump = gc(verbose=FALSE) }
	on.exit(poop())

	## Check for same number of catch policies by fleet
	##  (using sapply which reverts to lapply when uneven vectors in list)
	n.pol = sort(unique(sapply(cp,length)))
	if (length(n.pol)>1) stop("Revise catch policy input list to have equal numbers by fleet")

	need.files = c("starter.ss","forecast.ss", "ss.par","ss.cor","ss.psv", paste0("admodel.",c("hes","cov")))
	if (linux)  ## without Linux server
		need.files = c(need.files, "data.ss","control.ss")  ## CG's runADnuts requires simple names
	else
		need.files = c(need.files, paste0(c("control","data"),".",tar.run.rwt,".ss")) ## RH's runADnuts allows specific names

	for (i in 1:n.pol) {
		#if (i %in% c(1,2)) next  ## already done
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
			clearFiles( setdiff(list.files(d.catch.policy,full.names=TRUE),list.dirs(d.catch.policy)) )  ## need full names and explicit exclusion of directories
		}
		pee = file.copy(from=file.path(d.target,need.files), to=d.catch.policy, copy.date=TRUE, overwrite=TRUE)
		
		## Modify forecast file to reflect new catch policies
		forecast = readLines(file.path(d.catch.policy,"forecast.ss"))

		## NOTE: pattern 2 = Linux because MCMC's use MLE's '.ss_new' files (output)
		fline    = grep("#_fcast_nyrs|# N forecast years",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,(regexpr(" #_fcast_nyrs| # N forecast years",forecast[fline])-1)), length(fyrs), forecast[fline])
		fline    = grep("#_ctl_rule_ul|# Control rule inflection",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,(regexpr(" #_ctl_rule_ul| # Control rule inflection",forecast[fline])-1)), 0.002, forecast[fline])
		fline    = grep("#_ctl_rule_ll|# Control rule cutoff",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,(regexpr(" #_ctl_rule_ll| # Control rule cutoff",forecast[fline])-1)), 0.001, forecast[fline])
		fline    = grep("#_fcast_yr1|#(\\s+)?FirstYear",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,4), max(fyrs) + 100, forecast[fline])

		## Not really necessary when USWC rebuilder is set to 0, but Adam Langley resets these just in case
		fline    = grep("##_rebuild_catch_yr|# Rebuilder:  first",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,4), max(fyrs) + 101, forecast[fline])
		fline    = grep("#_rebuild_start_yr|# Rebuilder:  year",forecast)
		forecast[fline] = sub(substring(forecast[fline],1,4), max(fyrs) + 101, forecast[fline])

		f1     = grep("#_Yr Seas Fleet Catch",forecast)
		f2     = grep("#_end_catpols|-9999 [01] [01] [01]",forecast); f2 = rev(f2)[1]
		nyrs   = length(fyrs)
		#fleet  = names(cp)
		nfleet = length(fleet)
		newpol = data.frame(Yr=rep(fyrs,nfleet), Seas=rep(season,nyrs*length(fleet)), Fleet=rep(fleet,each=nyrs), Catch=rep(ii,each=nyrs))
		foopol = apply(newpol,1,paste0,collapse=" ")
		forecast = c(forecast[1:f1], foopol, forecast[f2:length(forecast)])
		writeLines(forecast, con=file.path(d.catch.policy,"forecast.ss"))
		.flush.cat(paste0("CP = ", i, " -- run mceval to generate catch policies.\n"))
		setwd(d.catch.policy)
		gdump = gc(verbose=FALSE)
		cmd = ifelse(onlinux,paste0(ss_executable, " -mceval"),"ss -mceval")
#browser();return()
		syssy = system(cmd, intern=TRUE)
		syssy = syssy[grep("Error -- base = 0", syssy, invert=T)]  ## get rid of stupid ADMB errors
		writeLines(syssy,paste0("./sys", pad0(i,4), ".txt"))
	}
	.flush.cat("C'est toute, folks\n")
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~prepCP


## prepMPD------------------------------2025-08-11
##  Prepare MPD runs from previous runs for analysis.
## ---------------------------------------------RH
prepMPD <- function(run.rwt.ver,
	d.base = "C:/Users/haighr/Files/GFish/PSARC25/SGR/Data/SS3/SGR2025",
	w=NULL, cvpro=NULL, modify=TRUE)
{
	padded = pad0(run.rwt.ver[c(1,2,4,5)],2)
	vers   = paste0("v",run.rwt.ver[c(3,6)])
	n.run=padded[1]; n.rwt=padded[2]; l.run=padded[3]; l.rwt=padded[4]
	pre.run.rwt = paste(l.run, l.rwt, vers[2], sep=".")
	new.run.rwt = paste(n.run, n.rwt, vers[1], sep=".")
	## Destination directories
	#d.cwd = getwd(); on.exit(setwd(d.cwd))
	d.run = file.path(d.base,paste0("Run",n.run))
	d.mpd = file.path(d.run,paste("MPD",new.run.rwt,sep="."))
#browser();return()

	## Previous run from which to poach ssfiles
	if (run.rwt.ver[4]==0 && run.rwt.ver[5]==0) {
		## Preliminary run so get files from d.base
		p.run = p.mpd = d.base
	} else {
		p.run = file.path(d.base,paste0("Run",l.run))
		p.mpd = file.path(p.run,paste("MPD",pre.run.rwt,sep="."))
	}
#browser();return()
	if (!all(file.exists(c(p.run,p.mpd))))
		stop("Previous Run directory and/or MPD directory does not exist")
	if (!file.exists(d.run)) dir.create(d.run)
	if (!file.exists(d.mpd)) dir.create(d.mpd)
	p.ss  = setdiff(list.files(p.mpd, pattern="\\.ss$"),"runnumber.ss")
	if (length(p.ss)!=4) {
		.flush.cat("Need 4 distinct ss files----------","\n"); browser(); return() }
	for (fnam in c("starter","forecast")) {
		pnam  = paste0(fnam, ".ss")
		pee   = file.copy(from=file.path(p.mpd,pnam), to=d.mpd, copy.date=T, overwrite=T)
	}
	for (fnam in c("data","control")) {
		from = file.path(p.mpd,paste(fnam,l.run,l.rwt,"ss",sep="."))
		to   = file.path(d.mpd,paste(fnam,n.run,n.rwt,"ss",sep="."))
		if (file.exists(from)) {
			poo   = file.copy(from=from, to=to, copy.date=T, overwrite=T)
		} else {
			.flush.cat("sumtingwong", "\n"); browser(); return()
		}
	}
	d.ss    = file.path(d.mpd, list.files(d.mpd, pattern="\\.ss$"))
	sfile   = d.ss[grep("/starter\\.",d.ss)]
	ffile   = d.ss[grep("/forecast\\.",d.ss)]
	dfile   = d.ss[grep("/data\\.",d.ss)]
	cfile   = d.ss[grep("/control\\.",d.ss)]
#browser();return()
	## Modify starter file to reflect new run
	#starter = readLines(file.path(d.mpd,"starter.ss"))
	#dline   = grep(paste("data", pre.run.rwt,sep="."), starter)
	#cline   = grep(paste("control", pre.run.rwt, sep="."), starter)
	#starter[dline] = sub(paste0("data\\.",l.run,"\\.",l.rwt),paste0("data.",new.run.rwt),starter[dline])
	#starter[cline] = sub(paste0("control\\.",l.run,"\\.",l.rwt),paste0("control.",new.run.rwt),starter[cline])
	#writeLines(starter, con=file.path(d.mpd,"starter.ss"))
#browser();return()
	starter = readLines(sfile)
	dline   = grep("(\\s+)?data\\.", starter)
	cline   = grep("(\\s+)?control\\.", starter)
	starter[dline] = sub(grep("data",p.ss,value=TRUE), basename(dfile), starter[dline])
	starter[cline] = sub(grep("control",p.ss,value=TRUE), basename(cfile), starter[cline])
	writeLines(starter, con=sfile)

	## Modify data file to add cvpro using the Francis method
	#data = readLines(file.path(d.mpd,paste0("data.",new.run.rwt,".ss")))
	data = readLines(dfile)
	if (!is.null(cvpro)){
		data = gsub("^\\s+|\\s+$", "", data)  ## get rid of leading and trailing spaces
		vline = intersect(grep("_index$",data), grep("^#",data,invert=TRUE))
		if (length(vline)==0) {
			.flush.cat("Cannot identify index data----------","\n"); browser(); return() }
		vbits = strsplit(data[vline], split="\\s+")
		if (length(.su(sapply(vbits,function(x){x[3]}))) != length(cvpro)){
			.flush.cat("User inputs for cvpro do not match indices in data file\n")#; browser();return()
		}
		vdump = sapply(1:length(vbits),function(i){vbits[[i]][5] <<- sqrt(as.numeric(vbits[[i]][5])^2 + cvpro[vbits[[i]][3]]^2) })
#browser();return()
		vnew  = sapply(vbits,function(x){ paste0(x,collapse=" ")})
		data[vline] = vnew
		#writeLines(data, con=file.path(d.mpd,paste0("data.",new.run.rwt,".ss")))
		writeLines(data, con=dfile)
	}
	## Modify control file to add Francis reweights
	#control = readLines(file.path(d.mpd,paste0("control.",new.run.rwt,".ss")))
	control = readLines(cfile)
	if (!is.null(w)){
		control = gsub("^\\s+|\\s+$", "", control)  ## get rid of leading and trailing spaces
		fline = intersect(grep("vadj_af",control), grep("^#",control,invert=TRUE))
		if (length(fline)==0) {
			.flush.cat("Cannot identify variance adjustment factors for AFs----------","\n"); browser(); return() }
#browser();return()
		if (length(fline)!=length(w)){
			.flush.cat("Francis reweights do not match control lines\n"); browser();return()
		}
		fbits = strsplit(control[fline], split="\\s+")
		wdump = sapply(1:length(fbits),function(i){fbits[[i]][3] <<- w[fbits[[i]][2]] }) #w[i]})
		fnew  = sapply(fbits,function(x){ paste0(x,collapse=" ")})
		control[fline] = fnew
		#writeLines(control, con=file.path(d.mpd,paste0("control.",new.run.rwt,".ss")))
		writeLines(control, con=cfile)
	}
	if (modify) {
		.flush.cat(paste0("***Edit the ss files in\n\t",d.mpd,"\nbefore proceeding with an MPD fit in SS."),"\n\n")
	} else {
		p.all  = list.files(p.mpd,full.names=F)
		p.copy = grep("\\.png$|\\.ss$|\\.bak$|french",p.all,value=T,invert=T)
		n.copy = sub(pre.run.rwt,new.run.rwt,p.copy)  ## change names of files sporting previous run.rwt.ver
		file.copy(from=file.path(p.mpd,p.copy), to=file.path(d.mpd,n.copy), copy.date=T, overwrite=T)
#browser();return()
		##file.copy sumting
	}
#	if (!is.null(w) || !is.null(cvpro)){
#	}
	setwd(d.mpd)
	return(invisible(list(starter=starter, data=data, control=control)))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~prepMPD


## quickDT------------------------------2024-08-26
##  Produce quick decision tables for Bt/LRP, Bt/USR
##  Designed for POP 2023; cannot handle composite base case
##  Added code to deal with inconsistencies between 
##   multi- and single-area model output.
##  Function is largel deprecated; use 'agileDT'
## ---------------------------------------------RH
quickDT <- function(compo, vals=c("BtLRP","BtUSR","BtBmsy","utumsy"),
   areas=c("5ABC","3CD","5DE"), cp=paste0("CC.", pad0(1:3,2)) )
{
	if (all(is.na(compo$xavgCP))) { ## not a multi-area model
		xavgCP = compo$avgCP
		areas  = "BC"
		xavgPJ = compo$avgPJ
#browser();return()
	} else {  ## multi-area model
		xavgCP = compo$xavgCP
		xavgPJ = compo$xavgPJ
	}
	cc = as.vector(t(xavgCP[1,1:length(areas),(1:length(cp))+1]))
	for (i in 1:length(vals)) {
		ii  = vals[i]
		num = substring(ii,1,2)
		den = substring(ii,3)
		if (all(areas=="BC")) {  ## deal with inconsistencies in naming arrays
			pj  = xavgPJ[,,ii,cp,drop=F]
			dimnames(pj)[3] = areas
			names(dimnames(pj))[3]="area"
			names(dimnames(pj))[2]="year"
#browser();return()
		} else {
			pj  = xavgPJ[,,ii,,cp]
		}
#browser();return()
		prb = apply(pj, c("proj","year","area"), function(x,num){
			xx=x[!is.na(x)]; xxx=if(num=="Bt") xx>1 else xx<1; sum(xxx)/length(xxx)
		}, num=num) ## over dims 4=catch policy, 2=year, 3=area
		dns = dimnames(prb); dns[[2]] = c("catch",dns[[2]])
		prb2 = array(NA, dim=sapply(dns,length), dimnames=dns)
		prb2[,-1,] = prb
		prb2[,1,]  = cc  ## insert catches for each policy by area
		num2 = paste0("$",sub("([Bu])","\\1_",num),"$")
		den2 = paste0("$",sub("msy","_\\\\text{MSY}",sub("USR","0.8Bmsy",sub("LRP","0.4Bmsy",den))),"$")
		tabcap = paste0("Probability ", num2, ifelse(num=="Bt"," > "," < "), den2)
		tex = texArray(prb2, sigdig=4, zero="0", use.row.names=F, name.row.names="", add.header.column=T, new.header=areas, outnam=ii, table.caption=tabcap)
		cmd1 = paste0("pdflatex.exe -quiet ", ii, "+.tex")
		cmd2 = paste0("latexmk.exe  -c ", ii, "+.tex")
#browser();return()
		system(command=cmd1, intern=TRUE, wait=TRUE)
		system(command=cmd2, intern=TRUE, wait=TRUE)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~quickDT


## repeatMPD----------------------------2024-04-03
## Repeat MPDs for axes of uncertainty to visualise likelihood density.
## Run from anywhere but specify the MPD directory as the starting point.
## Use 'plotSS.profile.r' to gather and plot the various MPDs.
## ---------------------------------------------RH
repeatMPD <- function(
   Pfix=list(FFF=seq(0.05,0.06,0.01), MMM=seq(0.05,0.06,0.01), RRR=NULL, SSS=NULL),
   A=60, dir.mpd=getwd(), dir.par, prefix="control.MMM.", clean=FALSE, 
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
				dc = c(dc, "")  ## decimal part as a character
				dp = c(dp, 0)   ## number of decimal places excluding zeroes
			}
#browser();return()
		}
		## Automatically determine minimum number of decimal places to yield unique values
		if (max(dp)<=0) {
			udp =0
		} else {
			udp = max(dp) + 1
			for (i in max(dp):1){
				#print(length(unique(round(x,i)))==length(x))
				if (length(unique(round(x,i))) == length(x) && !any(round(x,i)==0))
					udp = udp - 1
			}
		}
#browser() ;return()
		return(list(dc=dc, dp=dp, mindp=udp))
	}
	cleanup <- function(){
		junkpat = c("^admodel","^ss3\\.","\\.log$","\\.sso$","\\.ss_new$","\\.png$","\\.bak$","\\.old$","runnumber\\.ss")
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

	## SS3 changed its output files 'ss.xxx' to 'ss3.xxx'; also 'data.ss_new' now appears as 'data_echo.ss_new'  (RH 240327)
	#keep  = c("admodel.hes",paste0(c("Report","CompReport","covar","warning"),".sso"), paste("ss",c("cor","eva","par","rep","std","mpd"),sep="."), "data.ss_new") ## MPD files to save
	keep  = c("admodel.hes",paste0(c("Report","CompReport","covar","warning"),".sso"), paste("ss3",c("cor","eva","par","rep","std","mpd"),sep="."), "data_echo.ss_new") ## MPD files to save
	Nprof = length(Y[[1]]) * length(A)
	Iprof = 0

	if (!file.exists(dir.mpd))
		stop ("Source MPD directory does not exist")
	dir.nll = sub("MPD","NLL",dir.mpd) ## (RH 240327)

	if (missing(dir.par)){
		fornow  = lapply(Y,function(x){ paste0("(", paste0(range(x),collapse="-"),")")})
		dir.sub =  paste0(c("Prof",paste0(substring(names(fornow),1,1),fornow)),collapse=".")
		#dir.sub =  paste0("pro.M(",paste0(range(M),collapse="-"),").A(",paste0(A,collapse=","),")")
		dir.par = file.path(dir.nll,dir.sub)
	} else {
		if (basename(dir.par)==dir.par || !grepl("^\\./|^[[:alpha:]]:/",dir.par))  ## accept proper relative paths or those starting from a root
			dir.par = file.path(dir.nll,sub("^/","",dir.par))
	}
	if (!file.exists(dir.nll)) dir.create(dir.nll)
	if (!file.exists(dir.par)) dir.create(dir.par)

	## Copy input files from MPD directory to NLL directory
	dfile = paste0(sub("\\.v[0-9]$","",sub("MPD|NLL|Like","data",basename(dir.mpd))),".ss")
	cfile = paste0(sub("\\.v[0-9]$","",sub("MPD|NLL|Like","control",basename(dir.mpd))),".ss")
	sfile = "starter.ss"
	ffile = "forecast.ss"
	if (!all(file.exists(file.path(dir.mpd,c(sfile,dfile,cfile,ffile)))))
		stop ("Cannot locate all four input files")
	file.copy(from=file.path(dir.mpd, sfile), to=file.path(dir.nll,sub("starter","starter.base", sfile)), copy.date=TRUE)
	for (ifile in c(dfile,cfile,ffile))
		file.copy(from=file.path(dir.mpd, ifile), to=file.path(dir.nll, ifile), copy.date=TRUE)
#browser() ;return()

	for (a in A){
		## Hangover from using multiple plus-class ages A
		##  Automation won't work unless further coding is added to deal with multiple MPDs 
		## Look for ss files
#		afile = file.path(dir.mpd,paste0(prefix,"A",pad0(a,ncA),".ss"))
#		sfile = file.path(dir.mpd,c("starter.base.ss", "forecast.ss", paste0(sub("\\.v[0-9]$","",sub("MPD|NLL|Like","data",basename(dir.mpd))),".ss")))
#		for (bfile in c(afile,sfile)) {
#			if(!file.exists(bfile)) {
#				if (file.exists(sub("/Like\\.","/MPD.",bfile)))
#					file.copy(sub("/Like\\.","/MPD.",bfile), dir.nll, copy.date=T)
#				else
#					stop(paste0("Cannot find file\n\t'", bfile, "'"))
#			}
#		}

		## Create the templates for the LL profile analyses (RH 240327)
		afile = file.path(dir.nll,paste0(prefix,"A",pad0(a,ncA),".ss"))
		if (!file.exists(afile)) {
			aline = readLines(file.path(dir.nll,cfile))
			if (grepl("RRR",prefix)) {
				zrrr  = grep("ln\\(r0", aline)[1]
				rrr   = strsplit(aline[zrrr],"\\s+")[[1]]
				rrr[3] = "RRR"; rrr[6]="0"; rrr[7]="-1"
				aline[zrrr] = paste0(rrr,collapse=" ")
				## If fixing R0, at least one parameter must be estimated in phase 1
				zfm  = grep("mgp[FM]_gp1_M$", aline)
				fm   = strsplit(aline[zfm],"\\s+")
				fm   = lapply(fm,function(x) {x[7]=1; return(x)})
				aline[zfm[1]] = paste0(fm[[1]],collapse=" ")
				aline[zfm[2]] = paste0(fm[[2]],collapse=" ")
#browser();return()
			}
			if (grepl("FFF",prefix)) {
				zfff  = grep("mgpF_gp1_M$", aline)[1]
				fff   = strsplit(aline[zfff],"\\s+")[[1]]
				fff[3] = "FFF"; fff[6]="0"; fff[7]="-4"
				aline[zfff] = paste0(fff,collapse=" ")
			}
			if (grepl("MMM",prefix)) {
				zmmm  = grep("mgpM_gp1_M$", aline)[1]
				mmm   = strsplit(aline[zmmm],"\\s+")[[1]]
				mmm[3] = "MMM"; mmm[6]="0"; mmm[7]="-4"
				aline[zmmm] = paste0(mmm,collapse=" ")
			}
			if (grepl("SSS",prefix)) {
				zsss  = grep("spr_sigma_r$", aline)[1]
				sss   = strsplit(aline[zsss],"\\s+")[[1]]
				sss[3] = "SSS"; sss[6]="0"; sss[7]="-1"; sss[1]="0"; sss[2]="2"
				aline[zsss] = paste0(sss,collapse=" ")
#browser() ;return()
			}
			writeLines(text=aline, con=afile)
		}
		aline = readLines (afile)

		for (i in 1:Nprof){
			Iprof = Iprof + 1
			yval = sapply(Y,function(x){x[i]})
			ychr = sapply(Ydec,function(x){x[i]})
			ychr = sapply(1:length(ychr),function(f){
				if(ipY[f]==0) ochr = show0(x=round(ychr[f],dpY[f]), n=dpY[f]) 
				else          ochr = pad0(x=ychr[f], n=ipY[f]+dpY[f], f=dpY[f])
#browser() ;return()
			return(ochr)
			})

			ipref = paste0("control.",substring(gsub("control|\\.","",prefix),1,1), paste0(sub("^0\\.","",ychr),collapse=""),".A",pad0(a,ncA),collapse="")
			.flush.cat(paste0("Processing run '", ipref, "' ..."), "\n")
			ifile = paste0(ipref,".ss")
			## for some reason, certain combos of M & A do not converge unless M is nudged
			if (strSpp=="WWR" && M==0.03 && A==45) smudge = 0.00001 ## necessary for WWR
			else smudge = 0 
			iline = aline
			## Replace dummy placeholder (RRR or FFF or MMM) with fixed value
			for (l in 1:length(yval))
				iline = gsub(Yvar[l], yval[l]+smudge, iline)
#browser() ;return()
			writeLines(iline, con=file.path(dir.nll,ifile))

			## Modify starter file to reflect new run
			starter = readLines(file.path(dir.nll,"starter.base.ss"))
			#dline   = grep(paste("data",pre.run.rwt,sep="."),starter)
			cline   = grep("^control",starter)
			#starter[dline] = sub(paste0("data\\.",l.run,"\\.",l.rwt),paste0("data.",new.run.rwt),starter[dline])
			starter[cline] = sub("^control.+\\.ss",paste0(ipref,".ss"),starter[cline])
			writeLines(starter, con=file.path(dir.nll,"starter.ss"))
#browser();return()

			#expr = paste("mess = shell(cmd=\"awatea -ind ",mfile,argsMPD,"\", wait=TRUE, intern=TRUE)",sep="")
			setwd(dir.nll)
			if (clean) cleanup()
			expr = paste("mpd = shell(cmd=\"ss\", wait=TRUE, intern=TRUE)",sep="")
			.flush.cat("   ", expr, "\n")
			eval(parse(text=expr))
			writeLines(mpd ,con=file.path(dir.nll,"ss.mpd"))

			for (jfile in keep){
				if (file.exists(file.path(dir.nll,jfile))){
					if (jfile=="ss3.par"){
						## r4ss' SS_profile seems to rename ss.par; repeat here but still copy ss.par below
						kfile = paste0(jfile,"_",Iprof,".sso")
						file.copy(from=file.path(dir.nll,jfile), to=file.path(dir.par,kfile), overwrite=TRUE, copy.date=TRUE)
					}
					#else {
					jbits = strsplit(jfile,"\\.")[[1]]
					#jbits = c(jbits,"sumtingwong")  ## just for testing
					jlen  = length(jbits)-1
					jbits[jlen] = paste0(jbits[jlen],Iprof)
					kfile = paste0(jbits,collapse=".")
					#}
					file.copy(from=file.path(dir.nll,jfile), to=file.path(dir.par,kfile), overwrite=TRUE, copy.date=TRUE)
				}
			}
#browser();return()
			if (clean) cleanup()
			rubbish = gc(verbose=FALSE)
		} ## end i loop (number of profiles)
	} ## end a loop (maximum age for plus class)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~repeatMPD


## runSweave----------------------------2025-08-22
##  Run Sweave code to build pdfs for MPD and 
##  MCMC runs (not Results appendix)
## ---------------------------------------------RH
runSweave <- function (d.model=getwd(), d.sweave, type="MPD", figs.only=FALSE, debug=FALSE)
{
	if (!grepl(ver,getwd())) {.flush.cat("Appears to be version msismatch: ",ver," not in ",basename(getwd()),"\n",sep=""); browser();return()}

	if (missing(d.sweave))
		d.sweave = "C:/Users/haighr/Files/Projects/R/Develop/PBSsynth/Authors/Rcode/develop"
	f.sweave = paste0(basename(d.model), ifelse(figs.only,".figs",""), ifelse(debug,".debug",""), ".Rnw")
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


## tabDQs-------------------------------2025-09-25
##  Tables of Derived Quantities
## ---------------------------------------------RH
tabDQs <- function(xavgRP, xavgTS, csv=TRUE)
{
	areaRPL = areaRPL.pjs = areaRPQ = areaRPQ.pjs = list()
	for (a in dimnames(xavgRP)$area) {
		areaRP    = xavgRP[,,a]
		#areaTS    = xavgTS[,as.character(c(startYear,prevYear,currYear)),c("B","u"),a]
		areaTS    = xavgTS[,as.character(startYear:currYear),c("B","u","fR"),a]
		areaRPtab = data.frame(
			B0        = areaRP[,"B0"],
			Bcurr     = areaTS[,as.character(currYear),"B"],
			BcurrB0   = areaTS[,as.character(currYear),"B"] / areaTS[,as.character(startYear),"B"],
			Uprev     = areaTS[,as.character(prevYear),"u"],
			Umax      = apply(areaTS[,,"u"],1,max,na.rm=T), ## for each mcmc sample across the time series,
			MSY       = areaRP[,"MSY"],
			Bmsy      = areaRP[,"Bmsy"],
			LRP       = 0.4 * areaRP[,"Bmsy"],
			USR       = 0.8 * areaRP[,"Bmsy"],
			BcurrBmsy = areaTS[,as.character(currYear),"B"] / areaRP[,"Bmsy"],
			BmsyB0    = areaRP[,"Bmsy"] / areaRP[,"B0"],
			Umsy      = areaRP[,"umsy"],
			UprevUmsy = areaTS[,as.character(prevYear),"u"] / areaRP[,"umsy"],
			pR        = areaTS[,as.character(prevYear),"fR"],
			pVB       = areaRP[,"pVB"]
		)
#browser();return()
		areaRPtab.pjs = data.frame(
			B0        = areaRP[,"B0"],
			Bcurr     = areaTS[,as.character(currYear),"B"],
			BcurrB0   = areaTS[,as.character(currYear),"B"] / areaTS[,as.character(startYear),"B"],
			Fprev     = -log(1-areaTS[,as.character(prevYear),"u"]),
			Uprev     = areaTS[,as.character(prevYear),"u"],
			MSY       = areaRP[,"MSY"],
			Bmsy      = areaRP[,"Bmsy"],
			BcurrBmsy = areaTS[,as.character(currYear),"B"] / areaRP[,"Bmsy"],
			BmsyB0    = areaRP[,"Bmsy"] / areaRP[,"B0"],
			Fmsy      = areaRP[,"Fmsy"],
			Umsy      = areaRP[,"umsy"],
			UprevUmsy = areaTS[,as.character(prevYear),"u"] / areaRP[,"umsy"]
		)
		areaRPL[[a]]     = areaRPtab
		areaRPL.pjs[[a]] = areaRPtab.pjs
		areaRPQ[[a]]     = apply(areaRPtab, 2, function(x){quantile(x, probs=quants5, na.rm=T)})
		areaRPQ.pjs[[a]] = apply(areaRPtab.pjs, 2, function(x){quantile(x, probs=quants5, na.rm=T)})
		write.csv(areaRPQ[[a]], paste0("subarea", a, ifelse(exists("istock"),paste0(".",istock),""), ".csv"))
	}
	return(list(areaRPL=areaRPL, areaRPL.pjs=areaRPL.pjs, areaRPQ=areaRPQ, areaRPQ.pjs=areaRPQ.pjs))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tabDQs


## weightAF-----------------------------2022-02-24
##  Weight age frequencies using harmonic mean ratio method
## ---------------------------------------------RH
weightAF <- function(replist, fleets, abase="agedbase", afld="Nsamp_adj",
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
	wsub = 1/((Rmax/10^floor(log10(Rmax)))*M^0.7)   ## WTF? some sort of discount factor?
	wadj = wmean / wsub

	w.mcallister = list(agedat=fltbase, hmean=hmean, amean=amean, wmean=wmean, Rmax=Rmax, M=M, wsub=wsub, wadj=wadj)
#browser();return()
	ttput(w.mcallister)
	return(w.mcallister)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~weightAF

