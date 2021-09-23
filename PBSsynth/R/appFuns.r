## ===============================================
## SS RESULTS APPENDIX FUNCTIONS
## --------------------------===
## plotSS.compo..........Make composite figures
## plotSS.senso..........Make sensitivity figures
## tabSS.compo...........Make base case tables
## tabSS.decision........Make decision tables
## tabSS.senso...........Make sensitivity tables
## ===============================================

## plotSS.compo-------------------------2021-07-28
## Make Composite Figures
## ---------------------------------------------RH
plotSS.compo = function(compo, spp.code="YMR", istock="YMR", 
   subset=NULL, redo.figs=FALSE, redo.panels=FALSE, 
   ptypes, lang, pngres=400, PIN=c(9,9))
{
	on.exit(gc(verbose=FALSE))
	unpackList(stock[[istock]][["Controls"]])
	unpackList(compo)
	unpackList(tcall(data.compo.figs))  ## created in Rnw
	createFdir(lang)

	modYrsChar  = as.character(modYrs)
	## Diagnostics for select parameters
	P.names     = colnames(avgPA)
	P.runs      = as.numeric(sapply(strsplit(rownames(avgPA),"\\."),function(x){x[1]}))
	P.run.ord   = unique(P.runs)
	P.run.nmc   = table(P.runs)[as.character(P.run.ord)]
	P.run.num   = rep(1:length(P.run.nmc), P.run.nmc)
	use.run.rwt = is.element(P.runs, P.run.ord) ## default use all runs but a subset might be used for management advice

#redoFigs=FALSE
#if (redoFigs) {

sumting=F
if (sumting) {
	## Spawning biomass (envelope)
	.flush.cat("SB envelope","\n")
	catpol = avgPJ[,,"Bt",c("CC.01","CC.05","CC.08")]  ## YMR RPR wanted to see projections using catches of (0,1250,2500)
	plotSS.pmcmc(Bt.mcmc, yrs=modYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=expression(paste("Spawning Biomass ", group("(",italic(B)[italic(t)],")"))), outnam=paste0(prefix,"compo.Bt"), xyType="envelope", ptypes=ptypes, pyrs=proYrs, catpol=catpol)

	## Spawning biomass relative to that at MSY (envelope)
	.flush.cat("BtBmsy envelope","\n")
	catpol = avgPJ[,,"BtBmsy",c("CC.01","CC.05","CC.08")]  ## YMR RPR wanted to see projections using catches of (0,1250,2500)
	plotSS.pmcmc(BtBmsy.mcmc, yrs=modYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=expression(italic(B)[italic(t)] / italic(B)[MSY]), outnam=paste0(prefix,"compo.BtBmsy"), xyType="envelope", ptypes=ptypes, pyrs=proYrs, LRP=0.4, USR=0.8, catpol=catpol)

	## Spawning biomass relative to that at B0 (envelope)
	.flush.cat("BtB0 envelope","\n")
	catpol = avgPJ[,,"BtB0",c("CC.01","CC.05","CC.08")]  ## YMR RPR wanted to see projections using catches of (0,1250,2500)
	plotSS.pmcmc(BtB0.mcmc, yrs=modYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=expression(paste("Depletion ", group("(",italic(B)[italic(t)] / italic(B)[0],")"))), outnam=paste0(prefix,"compo.BtB0"), xyType="envelope", ptypes=ptypes, pyrs=proYrs, LRP=0.2, USR=0.4, catpol=catpol)
} ## end sumting

	## Exploitation rate (quantile boxes)
	.flush.cat("ut boxes","\n")
	plotSS.pmcmc(ut.mcmc, yrs=modYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=expression(paste("Exploitation rate ", group("(",italic(u)[italic(t)],")"))), outnam=paste0(prefix,"compo.ut"), ptypes=ptypes, USR=NULL) #mean(apply(ut.mcmc[,as.character(modYrs)],1,median)))

	## Recruitment (quantile boxes)
	.flush.cat("Rt boxes","\n")
	plotSS.pmcmc(Rt.mcmc, yrs=modYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=expression(paste("Recruitment ", group("(",italic(R)[italic(t)],")"))), outnam=paste0(prefix,"compo.Rt"), ptypes=ptypes, pyrs=proYrs, USR=NULL) ##mean(apply(Rt.mcmc[,as.character(modYrs)],1,median)))
browser();return()

	## Exploitation rate relative to that at MSY (quantile boxes)
	.flush.cat("utumsy boxes","\n")
	plotSS.pmcmc(utumsy.mcmc, yrs=modYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=expression(paste("Exploitation rate relative to ", italic(u)[MSY], group("(",italic(u)[italic(t)] / italic(u)[MSY],")"))), outnam=paste0(prefix,"compo.utumsy"), ptypes=ptypes, USR=1)

	## Recruitment deviations (quantile boxes)
	.flush.cat("Rtdev boxes","\n")
	plotSS.pmcmc(Rtdev.mcmc, yrs=modYrs, lang=lang, cex.axis=1.2, cex.lab=1.5, yLab=expression(paste("Recruitment deviations", group("(",delta~italic(R)[italic(t)],")"))), outnam=paste0(prefix,"compo.Rtdev"), ptypes=ptypes, USR=0, yLim=c(-2.75,3.5))
#browser();return()

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
		fout = fout.e = paste0(prefix,"compo.", iii, ".traces")
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=7, horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
				#colnames(P.ii) = paste(iii,paste0(switch(l, 'e'="R",'f'="E"),names(P.i)),sep="_")
				colnames(P.ii) = paste(iii,base.lab,sep=" ")
				panelTraces(mcmc=P.ii, mpd=ampdPA[,i], xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, same.limits=ifelse(i%in%c(1),F,T),lang=l)
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
#browser();return()

		## Split chains
		fout = fout.e = paste0(prefix,"compo.", ii, ".chains")
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=7, horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))  ## mar and oma ignored, fixed in call to `mochaLatte'
				panelChains(mcmc=P.ii, axes=TRUE, pdisc=0, between=list(x=0, y=0), col.trace=c("red","blue","black"), xlab="Parameter Value", ylab="Cumulative Frequency", cex.axis=1.2, cex.lab=1.4, yaxt="n", lang=l)
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
#browser();return()

		## ACF plots
		fout = fout.e = paste0(prefix,"compo.", ii, ".acfs")
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), width=8, height=8, horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
				plotACFs(P.ii, lag.max=60, lang=l)
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}

	## BtBmsy UtUmsy snail trail plot
	## ------------------------------
	z = use.run.rwt
	BoverBmsy = avgTS[z,,"BtBmsy"]
	UoverUmsy = avgTS[z,,"utumsy"]
	outnam    = paste0(prefix,"compo.snail")
	plotSnail(BoverBmsy, UoverUmsy, yrs=modYrs, p=tcall(quants3)[c(1,3)], xLim=NULL, yLim=NULL, ngear=length(gseries), assYrs=assYrs, outs=F, Cnames="Trawl+", ptypes=ptypes, outnam=outnam, lang=lang)

#	fout = fout.e = paste0(istock,".compo.snail")
#	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
#		changeLangOpts(L=l)
#		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
#		for (p in ptypes) {
#			if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=5, horizontal=FALSE, paper="special")
#			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2]*0.75)
#			par(mfrow=c(1,1), mar=c(3,3.75,0.5,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
#			## see 'assYrs' in run-masterMCMC.Snw
#			plotSnail(BoverBmsy, UoverUmsy, p=tcall(quants3)[c(1,3)], xLim=NULL, yLim=NULL, ngear=ngear, currYear=currYear, assYrs=assYrs, outs=F, Cnames="Trawl+", lang=l)
#			if (p %in% c("eps","png")) dev.off()
#		} ## end p (ptypes) loop
#	}; eop()
#browser();return()
#}

	## Prepare base composite and components for function 'compBmsy'
	## ------------------------------------------------------------
	Bbase = list()
	L1    = toupper(ifelse(spp.code=="REBS", istock, spp.code))
	Bbase[[L1]] = list()
	ibase  = toupper(istock)

	z  = use.run.rwt
	zz = sapply(split(z,P.runs),unique)[as.character(P.run.ord)]
	Ubase = sum(zz)
	runlab = base.lab[zz]

	Btemp = data.frame(run=P.runs, BcurrBmsy=avgTS[,as.character(currYear),"BtBmsy"])[z,]
	Bcurr.Bmsy = Btemp[,"BcurrBmsy"]
	Bbase[[L1]][[ibase]] = Bcurr.Bmsy ## composite from multiple bases
	for (i in 1:length(P.run.ord[zz])) {
		ii = P.run.ord[i]
		if (!is.element(ii, Btemp$run)) next
		iBtemp = Btemp[is.element(Btemp$run,ii),]
		iii = runlab[i]
		Bbase[[L1]][[paste0("R",ii)]] = iBtemp$BcurrBmsy
	}

	L1nam = ifelse(L1=="BSR","REBS North", ifelse(L1=="RER","REBS South", L1))
	#Mnams = c(paste0(L1nam," Composite"), gsub("_| AE=3","",base.runs.lab))
	Mnams = c(paste0(L1nam," Composite"), runlab)
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
		bord=c(1:N); medcol=c(ifelse(istock==names(stock)[1],"blue","red"),rep("grey30",Ubase)); 
		boxfill=c(ifelse(istock==names(stock)[1],"aliceblue","mistyrose"), rep("grey95",Ubase))
	}
	boxlim = c(0, max(sapply(Bbase[[L1]],quantile,tcall(quants5)[5])) )
	boxlim[2] = boxlim[2] + diff(boxlim) * 0.04
#ptypes="png"; scenario=0; so("compBmsy.r")

	if ("win" %in% ptypes) resetGraph()
	par(mfrow=c(1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
	mess =  paste0("list(",paste0(paste0(ptypes,"=TRUE"),collapse=","),")")

	out = compBmsy(Bspp=Bbase, spp=L1, boxwidth=0.5, medcol=medcol, boxfill=boxfill, boxlim=boxlim, whisklwd=2, staplelwd=2, Mnams=Mnams[bord], width=9, height=6, figgy=eval(parse(text=mess)), pngres=pngres, spplabs=F, t.yr=currYear, left.space=9, top.space=1, fout=paste0(prefix,"compo", ifelse(Ubase>1,".stock.","."), "status"), calcRat=F, lang=lang)
#}
#browser();return()

	## Make quantile plots of component-run parameters and quantities
	## --------------------------------------------------------------
	if (redo.panels) {
		so("mochaLatte.r","awatea")
		so("panelBoxes.r","awatea")

		P.collect = c(1:11) ## YMR: collect all parameters
		P.pars    = data.frame( avgPA[,P.collect] )
		colnames(P.pars) = colnames(avgPA)[P.collect]

		Q.pars   = data.frame(
			Bcurr      = avgRP[,"Bcurr"],
			B0         = avgRP[,"B0"],
			Bcurr.B0   = avgRP[,"Bcurr"]/avgRP[,"B0"],
			MSY        = avgRP[,"MSY"],
			Bmsy       = avgRP[,"Bmsy"],
			Bmsy.B0    = avgRP[,"Bmsy"]/avgRP[,"B0"],
			Ucurr      = avgRP[,"ucurr"],
			Umsy       = avgRP[,"umsy"],
			Umax       = apply(avgTS[,,"ut"],1,max) ## for each mcmc sample across the time series
			#Ucurr.Umsy = avgRP[,"ucurr"]/avgRP[,"umsy"]
		)
		nchains = length(P.run.ord)

		names(Q.pars) = sub("\\.", " / ", gsub("U","u", 
			gsub("Ucurr",paste0("U",currYear), 
			gsub("Bcurr",paste0("B",currYear), names(Q.pars)))))

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

		fout = fout.e = paste0(prefix,"compo.pars.qbox")
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=8, height=9)
				pflds  = setdiff(colnames(P.pars), "N")
				#pflds  = c("q_6",pflds)
				panelBoxes(P.pars[,pflds], nchains=nchains, xlab=linguaFranca("Base Runs",l), ylab=linguaFranca("Parameter estimates",l), cex.strip=1.2, cex.axis=1.1, cex.lab=1.5, outline=FALSE, xlim=xlim, boxfill=boxfill, xfac=paste0("B",names(run.num[use.num])), mar=c(0,3.8,0,0), oma=c(4,2,0.5,1))  ## RH 200508 (for subsets of B)
				if (p %in% c("png","eps")) dev.off()
			}
		}; eop()
#browser();return()

		fout = fout.e = paste0(prefix,"compo.rfpt.qbox")
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=8, height=8)
				panelBoxes(Q.pars, nchains=nchains, xlab=linguaFranca("Base Runs",l), ylab=linguaFranca("Derived Quantities",l), cex.strip=1.2, cex.axis=1.1, cex.lab=1.5, outline=FALSE, xlim=xlim, boxfill=boxfill, xfac=paste0("B",names(run.num[use.num])), mar=c(0,3.8,0,0), oma=c(4,2,0.5,1))  ## RH 200508 (for subsets of B)
				if (p %in% c("png","eps")) dev.off()
			}
		}; eop()
	}
#browser();return()

	## Prepare stock base composites for function 'compBmsy'
	## If comparing two stocks (e.g., RSR,POP)
	if (length(stock)>2) { ## disable for now (fix later)
		if (istock==names(stock)[1]){
			Bstock = list()
			Bstock[[toupper(spp.code)]] = list()
		}
		ibase  = istock
		Bstock[[1]][[ibase]] = Bcurr.Bmsy ## composite from multiple bases
		tput(Bstock)
		## May have more than on stock so need to save Bstock somewhere for retrieval later
		if (istock==rev(names(stock))[1]) {
			tget(Bstock)
			Mnams = paste0(spp.code,"\n",sapply(stock,function(x){x$area.name}))
			N = length(Mnams)
			bord=c(1:N); medcol=c("blue","red")[1:N]; boxfill=c("aliceblue","mistyrose")[1:N] ## will need more if Nstock>2
			if ("win" %in% ptypes) resetGraph()
			mess =  paste0("list(",paste0(paste0(ptypes,"=TRUE"),collapse=","),")")
			out = compBmsy(Bspp=Bstock, spp=c(spp.code), boxwidth=0.5, medcol=medcol, boxfill=boxfill, whisklwd=2, staplelwd=2, Mnams=Mnams[bord], width=9, height=4, figgy=eval(parse(text=mess)), pngre=pngres, spplabs=F, t.yr=currYear, left.space=6, top.space=1, fout=paste0(spp.code,".base", ifelse(Ubase>1,".composite.","."), "status"), calcRat=F, lang=lang)
		}
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.compo


## plotSS.senso-------------------------2021-09-16
## Make Sensitivity Figures
## ---------------------------------------------RH
plotSS.senso = function(senso, ptypes, spp.code="YMR", istock="YMR",
	subset=NULL, redo.figs=FALSE, redo.panels=FALSE, lang, useRlow=FALSE,
	png=F, pngres=400, PIN=c(9,9))
{
	on.exit(gc(verbose=FALSE))
	unpackList(stock[[istock]][["Controls"]])
	unpackList(senso)
	#unpackList(tcall(data.compo.figs))  ## created in Rnw
	createFdir(lang)

	col.senso = c("green4","blue","red","purple","orange","skyblue","gold3","darkolivegreen","lightseagreen","hotpink","brown","cyan","darkorchid1","chartreuse3")
	lty.senso = c("22", "44", "66", "13", "73", "1551", "1343", "2262", "3573", "654321", "12345678", "88", "17", "5937")
	#modYrs      = startYear:currYear
	modYrsChar  = as.character(modYrs)
	## Diagnostics for select parameters
	P.names     = colnames(senPA) ## paremeter names
	S.mpd       = smpdPA; storage.mode(S.mpd)="double"     ## For some reason, this matrix is stored as a list (maybe due to NA values?)
	S.runs      = as.numeric(sapply(strsplit(rownames(senPA),"\\."),function(x){x[1]})) ## vector of Runs
	S.run.rwt   = sapply(strsplit(rownames(senPA),"\\."),function(x){paste0(x[1],".",x[2])}) ## vector of Runs and Rwts
	S.run.ord   = unique(S.runs)                           ## unique Run numbers (1st is central run)
	S.run.nmc   = table(S.runs)[as.character(S.run.ord)]   ## vector of sensitivity runs ordered
	S.run.num   = rep(1:length(S.run.nmc), S.run.nmc) - 1  ## 1st run is Central Run
	S.num       = unique(S.run.num)
	## Look into subsetting later
	use.run.rwt = unique(S.run.rwt)
	#is.element(P.runs, P.run.ord) ## default use all runs but a subset might be used for management advice

	## Create Sensitivity labels
	S.prefix    = paste0("S",pad0(S.num,2)," (R", pad0(S.run.ord,2) ,") ")
	iCR         = grep(exp.run.rwt,unique(B.index))
	S.prefix[1] = paste0("B",iCR, " (R", strsplit(exp.run.rwt,"\\.")[[1]][1], ") ")
	S.labels    = paste0(S.prefix, gsub("\\_"," ", c("Central Run",sen.lab)))
	S.labels    = sub("\\\\pc", "%", S.labels)  ## just in case
#browser();return()

	## Calculate quantiles for later (eventually make this function global)
	calcQs = function (dat, ivec, ovec) {
		F.data = as.data.frame(dat)
		F.list = split(F.data, ivec)
		rr.ord = match(ovec, substring(names(F.list),1,2))  ## oder not so important for sensitivities
#browser();return()
		F.list = F.list[rr.ord]
		F.qnts = lapply(F.list,function(x){
			z = apply(x,2,function(xx){!any(is.na(xx))})
			out = apply(x[z],2,quantile,quants5)
			return(out)
		}) ## lapply does not sort like split does
		return(F.qnts)
	}
	P.qnts = calcQs(senPA, ivec=S.run.rwt, ovec=S.run.ord)
	D.qnts = calcQs(senTS[,,"BtB0"], ivec=S.run.rwt, ovec=S.run.ord)
	U.qnts = calcQs(senTS[,,"ut"], ivec=S.run.rwt, ovec=S.run.ord)
	R.qnts = calcQs(senTS[,,"Rt"], ivec=S.run.rwt, ovec=S.run.ord)
	RD.qnts = calcQs(senTS[,,"Rtdev"], ivec=S.run.rwt, ovec=S.run.ord)
#browser();return()

	L1 = if (spp.code %in% c("REBS")) istock else spp.code  ## RH 200416

sumting=F
#if (sumting) {

	## Diagnostics for select parameters
	if (redo.figs) {
		for (i in c(1)){
			ii   = P.names[i]
			iii  = gsub("[_|[:space:]]","",ii)
			P.i  = split(senPA[,i], S.runs)[as.character(S.run.ord)]
			## splits by 1:length(run.num) in 'gather.compo.case.r' so retains correct order
			P.ii = data.frame(P.i)
			colnames.e = paste0(ii, ": ", S.labels)
			colnames.f = linguaFranca(colnames.e, "f")
			fout = fout.e = paste0(prefix,"senso.", ii, ".traces")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				colnames(P.ii) = switch(l, 'e'=colnames.e, 'f'=colnames.f)
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=7, horizontal=FALSE,  paper="special")
					else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
					##par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0)) ## ignored
					panelTraces(mcmc=P.ii, mpd=S.mpd[,i], xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, cex.strip=1.0, same.limits=ifelse(ii%in%c("LN(R0)"),F,T), lang=l, mar=c(0,2.5,0,0), oma=c(3.2,2,0.5,0.5))
					if (p %in% c("eps","png")) dev.off()
				} ## end p (ptypes) loop
			}; eop()
#browser();return()

			fout = fout.e = paste0(prefix,"senso.", ii, ".chains")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				colnames(P.ii) = switch(l, 'e'=colnames.e, 'f'=colnames.f)
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(fout,".eps"), width=6.25, height=7, horizontal=FALSE,  paper="special")
					else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
					##par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))  ## mar and oma ignored, fixed in call to `mochaLatte'
					panelChains(mcmc=P.ii, axes=TRUE, pdisc=0, between=list(x=0, y=0), col.trace=c("red","blue","black"), xlab="Parameter Value", ylab="Cumulative Frequency", cex.axis=1.2, cex.lab=1.4, cex.strip=ccex.strip, lang=l, mar=c(2,0,0,0), oma=c(3,4,0.5,0.5)) #, yaxt="n"
					if (p %in% c("eps","png")) dev.off()
				} ## end p (ptypes) loop
			}; eop()
#browser();return()

			fout = fout.e = paste0(prefix,"senso.", ii, ".acfs")
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				colnames(P.ii) = switch(l, 'e'=colnames.e, 'f'=colnames.f)
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(fout,".eps"), width=8, height=8, horizontal=FALSE,  paper="special")
					else if (p=="png") png(paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
					plotACFs(P.ii, lag.max=60, lang=l)
					if (p %in% c("eps","png")) dev.off()
				} ## end p (ptypes) loop
			}; eop()
		}
	}
#browser();return()
#} ## sumting

	## Make quantile plots of component-run parameters and quantities
	## Note: When outliers are plotted, these quantile plots take forever when sent to a png file.
	if (redo.panels) {
		#so("mochaLatte.r","awatea")
		#so("panelBoxes.r","awatea")

		#P.pool = P.mcmc
		#P.sens = S.mcmc ## created in 'gather.sens.runs.'
		### BOR central is 2nd of 3 but need to reduce  by 1 to get right numbe rof rows (1001:2000)
		##i.cent = ((grep(exp.run.rwt,stock[[istock]]$Base$run.rwts)-1)*nmcmc)+(1:nmcmc)  ## a little bit arbitrary
		#i.cent = (1:length(B.index))[is.element(B.index,exp.run.rwt)]  ## RH 200424: Should be more stable when central run is not necessarily 'central'
		#P.cent = P.pool[i.cent,]  ## grab the central run (Run05 for WWR 2019, Run02 for BOR 2019)
		#pflds  = intersect(names(P.cent),names(P.sens))  ## get fields in common
		#P.cent.sens = rbind(P.cent[,pflds],P.sens[,pflds])

#browser();return()
		#qflds  = c("Bcurr", "B0", "Bcurr.B0", "MSY", "Bmsy", "Bmsy.B0", "Umsy")
		#Q.sens = iQ.sens[[1]][,qflds] 
		#Q.cent = data.frame(Bcurr, B0, Bcurr.B0, MSY, Bmsy, Bmsy.B0, Umsy)
		#for (g in 1:Ngear) {
		#	Q.sens.add = iQ.sens[[g]][,c("Ucurr","Umax")]
		#	colnames(Q.sens.add) = paste0(c("Ucurr_", "Umax_"),gseries[g])
		#	Q.sens = cbind(Q.sens, Q.sens.add)
		#	Q.cent.add = data.frame(Ucurr[[g]], Umax[[g]])
		#	colnames(Q.cent.add) = paste0(c("Ucurr_", "Umax_"),gseries[g])
		#	Q.cent = cbind(Q.cent, Q.cent.add)
		#}
		#Q.cent = Q.cent[i.cent,]
		#Q.cent.sens = rbind(Q.cent,Q.sens)
		#nchains = nrow(Q.cent.sens)/nmcmc

		P.collect   = c(1:6) ## YMR: collect  all parameters except M
		P.cent.sens = data.frame( senPA[,P.collect] )
		colnames(P.cent.sens) = colnames(senPA)[P.collect]

		Q.cent.sens   = data.frame(
			Bcurr      = senRP[,"Bcurr"],
			B0         = senRP[,"B0"],
			Bcurr.B0   = senRP[,"Bcurr"]/senRP[,"B0"],
			MSY        = senRP[,"MSY"],
			Bmsy       = senRP[,"Bmsy"],
			Bmsy.B0    = senRP[,"Bmsy"]/senRP[,"B0"],
			Ucurr      = senRP[,"ucurr"],
			Umsy       = senRP[,"umsy"],
			Umax       = apply(senTS[,,"ut"],1,max) ## for each mcmc sample across the time series
			#Ucurr.Umsy = senRP[,"ucurr"]/senRP[,"umsy"]
		)
		nchains = length(S.run.ord)

		names(Q.cent.sens) = sub("\\.", " / ", gsub("U","u", 
			gsub("Ucurr",paste0("U",currYear), 
			gsub("Bcurr",paste0("B",currYear), names(Q.cent.sens)))))

		if (spp.code %in% c("BOR","WWR","POP")) {
			sflds = c("R_0","h","q_1","mu_1","q_7","mu_7")
		} else if (spp.code %in% c("REBS")) {
			if (istock=="BSR")
				sflds = c("R_0","mu_3","q_1","mu_1","q_2","mu_2")
			if (istock=="RER")
				sflds = c("R_0","q_1","q_2","q_3","mu_4","mu_5")
		} else 
			P.pars  = P.cent.sens

		## Sometimes want to exclude sensitivities  ## RH 200416
		verboten = NULL
		if (spp.code=="BOR") {
			verboten = "no_CPUE"
			zap.runs = match(verboten,sen.lab)
			#zap.runs = grep("9",c(sen.run.nums))
		}
		if (spp.code=="YMR") {
			verboten = c("estimate_M","start_Rdevs_in_1970")
			zap.runs = match(verboten,sen.lab) ## + 1  ## to account for central run
		}
		if (!is.null(verboten)) {
#browser();return()
			## Not really robust to mcsub choices
			nmcmc    = sen.mcsubs[zap.runs]
			bad.rows = lapply(1:length(zap.runs),function(i){ (nmcmc[[i]][1]:nmcmc[[i]][2]) + (zap.runs[i] * nmcmc[[i]][2]) })
			zap.rows = unlist(bad.rows)
			#P.pars[zap.rows,] = NA
		}
		fout = fout.e = paste0(prefix,"senso.pars.qbox")
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="png")
					png(filename=paste0(fout,".png"), units="in", res=pngres, width=8, height=8)
				#boxfill = c("gainsboro","green4","blue","red","purple","orange","skyblue","gold3","salmon3","red2","hotpink","darkred","dodgerblue") ## colours used in the trajectory plots (REBS)
				boxfill = c("gainsboro",col.senso) ## colours used in the trajectory plots (YMR)
				panelBoxes(P.pars, xlim=c(0.25,nchains+0.75), xlab=linguaFranca("Sensitivity Runs",l), ylab=linguaFranca("Parameter estimates",l), nchains=nchains, xfac=c("CR", paste0("S",pad0(1:(nchains-1),2))), boxfill=boxfill, cex.strip=1.2, cex.axis=1.2, cex.lab=1.5, outline=FALSE)
				if (p %in% c("png","eps")) dev.off()
			}
		}; eop()
#browser();return()

		Q.pars = Q.cent.sens
		if (!is.null(verboten)) {
			Q.pars[zap.rows,] = NA
			fn.ylim = function(x){extendrange(sapply(split(x,names(x)), quantile,quants5[c(1,5)], na.rm=TRUE))}
		} else {
			## Sometimes the 95pc limits are outrageously high so use this function:
			fn.ylim = function(x,i=rep(1:12,each=2000)){ ## hardwire index i for now
				xr = range(x)
#browser();return()
				xx = split(x,i);
				if (any(sapply(xx, function(xxx){quantile(xxx,0.99) > 10*quantile(xxx,0.75)}))){
					xr[2] = quantile(x,0.95)
				}
				return(xr)
			}
		}
#browser();return()
		
		fout = fout.e = paste0(prefix,"senso.rfpt.qbox")
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			LQ.pars = Q.pars
			names(LQ.pars) = linguaFranca(names(LQ.pars),l)
			for (p in ptypes) {
				if (p=="png") png(filename=paste0(fout,".png"), units="in", res=pngres, width=8, height=8)
				#boxfill = c("gainsboro","green4","blue","red","purple","orange","skyblue","gold3","salmon3","red2","hotpink","darkred","dodgerblue") ## colours used in the trajectory plots (REBS)
				boxfill = c("gainsboro",col.senso) ## colours used in the trajectory plots (YMR)
				panelBoxes(LQ.pars, xlim=c(0.25,nchains+0.75), xlab=linguaFranca("Sensitivity Runs",l), ylab=linguaFranca("Derived Quantities",l), nchains=nchains, xfac=c("CR", paste0("S",pad0(1:(nchains-1),2))), boxfill=boxfill, cex.strip=1.2, cex.axis=1.1, cex.lab=1.5, outline=FALSE, fn.ylim=fn.ylim)
				if (p %in% c("png","eps")) dev.off()
			}
		}; eop()
	}
#browser();return()

	## Make plots of median trajectories
	## ---------------------------------
	ii   = as.character(startYear:currYear)
	bb   = list('BtB0'=D.qnts, 'U'=U.qnts, 'R'=R.qnts, 'RD'=RD.qnts)
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
#browser();return()
	bdat$label = S.labels

	tlty = c("solid", lty.senso)
	tcol = c("black",col.senso)

	#so("plotTraj.r","awatea")
	if (redo.figs) {
		traj.names = c("BtB0","U","R","RD")
		Ntraj = length(traj.names)
		for (k in 1:Ntraj){
			kk = traj.names[k]
			logR = FALSE

			traj.meds = bdat[[kk]]
			if (is.null(traj.meds)) next
			fout = fout.e = paste0(prefix,"senso.traj.",kk)
			traj = traj.meds[intersect(rownames(traj.meds),ii),]
			Nruns = ncol(traj)
			Nyrs  = nrow(traj)
			x     = as.numeric(rownames(traj))
			xlim  = range(x)
			ylim = range(traj,na.rm=TRUE)
			y0=ifelse(kk=="RD",F,T);   if (y0) ylim[1] = 0; if (kk=="RD") ylim[1]=-3; if (kk=="BtB0") ylim[1]=-0.2
			logR=F; if (kk=="R" && logR) ylim[1]=1
			## tcol and legtxt need to be reset inside language loop because they are altered when bdat is supplied
			tcol    = rep(tcol,Nruns)[1:Nruns]
			tlty    = rep(tlty,Nruns)[1:Nruns]
			legtxt  = bdat$label
			for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
				changeLangOpts(L=l)
				fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				for (p in ptypes) {
					if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
					else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
					par(mfrow=c(1,1), mar=c(3,3.5,0.5,0), oma=c(0,0,0,1), mgp=c(2,0.5,0))
					plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="",log=ifelse(kk=="R"&&logR,"y",""))
					if (kk %in% c("BtB0"))
						abline(h=c(0.2,0.4,1), lty=5, col="grey20") #col=c("salmon","darkorchid","navy"))
					if (kk %in% c("RD"))
						abline(h=c(0), lty=5, col="grey20") #col=c("salmon","darkorchid","navy"))
					for (j in Nruns:1) {
#.flush.cat(j," ",jj," ", tlty[j], "\n")
						jj = colnames(traj)[j]
						y  = traj[,jj]
						lines(x,y, col=tcol[j], lwd=ifelse(j==1,3,2), lty=tlty[j])
					}
					#lines(x=as.numeric(names(bline)), y=if (jj=="R" && logR) log10(bline) else bline, col="black", lty=1, lwd=3)
					mtext(linguaFranca("Year",l), side=1, line=1.75, cex=ifelse(Nruns==1,1.5,1))
					ylab = switch(kk, 'BtB0'="Spawning Biomass Depletion", 'B'="Spawning Biomass", 'VB'="Vulnerable Biomass", 'R'="Recruitment", 'RD'="Recruitment Deviation", 'U'="Exploitation Rate", "Unknown")
					#mtext(linguaFranca(ylabs[[jj]],l), side=2, line=1.8, cex=ifelse(Ntraj==1,1.5,1.2))
					mtext(linguaFranca(ylab,l), side=2, line=1.8, cex=ifelse(Ntraj==1,1.5,1.2))
					#if (k==1){
						legtxt = gsub("_"," ",legtxt)
						addLegend(ifelse(kk%in%c("R","U"),0.025,0.025), ifelse(kk%in%c("BtB0","RD"),0.01,0.975), col=tcol, seg.len=5, legend=linguaFranca(legtxt,l), bty="o", box.col="grey", bg=ifelse(kk%in%c("BtB0"),"white","transparent"), xjust=ifelse(kk%in%c("R","U"),0,0), yjust=ifelse(kk%in%c("BtB0","RD"),0,1), lwd=2, lty=tlty)
					#}
					#if (Ntraj==1) {
						axis(1,at=intersect(seq(1900,2500,5),x),labels=FALSE,tcl=-0.2)
						axis(1,at=intersect(seq(1900,2500,10),x),labels=FALSE)
					#}
					box()
					if (p %in% c("png","eps")) dev.off()
				}
			}; eop()
		}
	}
#browser();return()

	## Prepare sensitivity runs for function 'compBmsy'
	## ------------------------------------------------
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

		Stemp = data.frame(run=S.runs[zzz], BcurrBmsy=senTS[,as.character(currYear),"BtBmsy"])[zzz,]
		Bcurr.Bmsy = Stemp[is.element(Stemp$run, central.run),"BcurrBmsy"]
		medCR = median(Bcurr.Bmsy)
		Bsens[[L1]][[runlab[1]]] = Bcurr.Bmsy ## central run from composite for snesitivity comparison
		for (i in 2:length(S.run.ord[zz])) {
			ii = S.run.ord[i]
			if (!is.element(ii, Stemp$run)) next
			iStemp = Stemp[is.element(Stemp$run,ii),]
			iii = runlab[i]
			Bsens[[L1]][[iii]] = iStemp$BcurrBmsy
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

#ptypes="png"; scenario=0; so("compBmsy.r")
		if ("win" %in% ptypes) resetGraph()
		mess =  paste0("list(",paste0(paste0(ptypes,"=TRUE"),collapse=","),")")

		Bspp = Bsens
		#Bspp = Bsens["BtBmsy"]; names(Bspp)=L1
		#if (spp.code=="BOR")
		#	Bspp[[L1]] = Bspp[[L1]][c(0,iisens)+1]  ## BOR 2021: remove no_CPUE
		#if (useRlow)
		#	Bspp[[L1]][[paste0(names(Bspp$BOR)[1],"_low")]] = Bspp[[L1]][[1]][Rlow]
		boxlim    = c(0, max(sapply(Bspp[[L1]],quantile,quants5[5])) )
		boxlim[2] = boxlim[2] + diff(boxlim) * 0.04
		left.space = 
			if (istock %in% c("BCN")) c(7.5,10)
			else if (istock %in% c("YMR")) c(12,15)
			else c(10,10)
		out = compBmsy(Bspp=Bsens, spp=L1, boxwidth=0.5, medcol=medcol, boxfill=boxfill, boxlim=boxlim, whisklwd=2, staplelwd=2, Mnams=Mnams[bord], width=9, height=6, figgy=eval(parse(text=mess)), pngres=pngres, spplabs=F, t.yr=currYear, left.space=left.space, top.space=1.5, fout=paste0(prefix,"senso", ifelse(Ubase>1,".stock.","."), "status",ifelse(useRlow,"+","")), calcRat=F, lang=lang, ratios=c(0.4,0.8,medCR), rlty=c(2,2,3), rlwd=c(2,2,1))
	}
#browser();return()

## Collect BtBmsy for discussions in Sweave ??
	#BtBmsy.med = sapply(Bsens[[1]],function(x){median(x)})
	#BtBmsy.sens.MCMC = list()
	#for (i in 1:(Nsens)) {
	#	iBmcmc = stock[[istock]]$Sens$Bmcmc.sens[[i]]
	#	ii = names(Bmcmc.sens[[i]][[toupper(spp.code)]])
	#	Bsens[[1]][[ii]] = iBmcmc[[toupper(spp.code)]][[ii]]$BoverBmsy[,currYearChar]
	#	BtBmsy.sens.MCMC[[ii]] = iBmcmc[[toupper(spp.code)]][[ii]]$BoverBmsy
	#}
	#BtBmsy.med = sapply(BtBmsy.sens.MCMC,function(x){apply(x,2,median)})

	## Identify lowest median stock status for sensitivities
	#lowSS      = which(BtBmsy.med[-1]==min(BtBmsy.med[-1]))
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.senso


## tabSS.compo--------------------------2021-08-24
## Make Base Case Tables
## ---------------------------------------------RH
tabSS.compo = function(istock="YMR", prefix="ymr.", compo,
  useRlow=FALSE, qRlow=0.25, sigdig=4)
{
	on.exit(gc(verbose=FALSE))
	unpackList(stock[[istock]][["Controls"]])
	unpackList(compo)

	modYrs      = startYear:currYear
	modYrsChar  = as.character(modYrs)
	## Diagnostics for select parameters
	P.mpd       = ampdPA; storage.mode(P.mpd)="double"     ## For some reason, this matrix is stored as a list (maybe due to NA values?)
	P.names     = colnames(avgPA)
	P.runs      = as.numeric(sapply(strsplit(rownames(avgPA),"\\."),function(x){x[1]}))
	P.run.rwt   = sapply(strsplit(rownames(avgPA),"\\."),function(x){paste0(x[1],".",x[2])}) ## vector of Runs and Rwts
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
	tabPmed = t(apply(avgPA,2,quantile,tcall(quants5)))  ## arrays so cannot use sapply
	tabPmed = tabPmed[P.ord,]
	tabPmed = formatCatch(tabPmed,N=sigdig)

	names.pars = rownames(tabPmed)
	fixsub = grep("_",names.pars)
	names.pars[fixsub] = paste0(sub("_","~(",names.pars[fixsub]),")")
	names.pars = sub("varL\\(","\\\\log v(\\\\mathrm{L}",names.pars)
	names.pars[fixsub] = sub("\\(","_{", sub("\\)","}", names.pars[fixsub]))
	names.pars = gsub("mu","\\\\mu",names.pars)
	names.pars = gsub("Delta","\\\\Delta",names.pars)
	names.pars = sub("LN\\(R0)","\\\\log R_{0}",names.pars)
	names.pars[fixsub] = sub("\\(","(\\\\text{", sub("\\)","})", names.pars[fixsub]))
	#names.pars = gsub("log","\\\\mathrm{log}",gsub("v_","v_{",gsub("L","L}",names.pars)))

	rownames(tabPmed) =  paste(rep("$",nrow(tabPmed)),names.pars,rep("$",nrow(tabPmed)),sep="")
	#colnames(tabPmed) =  gsub("\\%","\\\\%",colnames(tabPmed))

	xtab.compo.pars = xtable(tabPmed, align="lrrrrr",
		label   = paste0("tab:",prefix,"base.pars"), digits = if (exists("formatCatch")) NULL else sigdig,
		caption = paste0("Composite base case: the ", texThatVec(tcall(quants5)), " quantiles for pooled model parameters (defined in \\AppEqn) from MCMC estimation of ", "\\numberstringnum{", NrefM, "} component model runs of \\Nmcmc{} samples each.") )
	xtab.compo.pars.out = capture.output(print(xtab.compo.pars,
		caption.placement="top", sanitize.rownames.function=function(x){x}, add.to.row =list(pos=list(-1), command=c("\\\\[-1.0ex]")) ) )
	tput(xtab.compo.pars.out)
#browser();return()
# \\hline\\\\[-2.2ex] CC & 2023 & 2024 & 2025 & 2026 & 2027 & 2028 & 2029 & 2030 & 2031 \\\\[0.2ex]\\hline\\\\[-1.5ex]
	##----------Table 2-----------
	## MCMC Derived Parameters
	##----------------------------
	yrsUse = modYrs[c(1,length(modYrs)+c(-1,0))]
	yrsChr = as.character(yrsUse)

	B.mcmc = data.frame (
		B0         = avgRP[,"B0"],
		Bcurr      = avgRP[,"Bcurr"],
		Bcurr.B0   = avgRP[,"Bcurr"] / avgRP[,"B0"],
		ucurr      = avgRP[,"ucurr"],
		umax       = apply(avgTS[,,"ut"],1,max), ## for each mcmc sample across the time series
		MSY        = avgRP[,"MSY"],
		Bmsy       = avgRP[,"Bmsy"],
		LRP        = avgRP[,"LRP"],
		USR        = avgRP[,"USR"],
		Bcurr.Bmsy = avgRP[,"Bcurr"] / avgRP[,"Bmsy"],
		Bmsy.B0    = avgRP[,"Bmsy"] / avgRP[,"B0"],
		umsy       = avgRP[,"umsy"],
		ucurr.umsy = avgRP[,"ucurr"] / avgRP[,"umsy"]
	)
	tabmsy = t(apply(B.mcmc,2,quantile,tcall(quants5))) 
	tabmsy = formatCatch(tabmsy, N=sigdig)
	#colnames(tabmsy) =  gsub("\\%","\\\\%",colnames(tabmsy))
	names.msy = rownames(tabmsy)
	names.msy =
		gsub("_Trawl", "~(\\\\mathrm{trawl})",
		gsub("_Other", "~(\\\\mathrm{other})",
		gsub("LRP",  paste0("0.4B_{\\\\mathrm{MSY}}"),
		gsub("USR",  paste0("0.8B_{\\\\mathrm{MSY}}"),
		gsub("VB",  "V",
		gsub("[Uu]max", "u_\\\\mathrm{max}",
		gsub("msy", "_\\\\mathrm{MSY}",
		gsub("curr",  paste0("_{",currYear,"}"),
		gsub("[Uu]msy",  "u_\\\\mathrm{MSY}",
		gsub("[Uu]curr", paste0("u_{",prevYear,"}"),
		gsub("B0", "B_{0}",
		gsub("\\.", "/",
		gsub("currB","curr/B",
		gsub("curru","curr/u",
		gsub("msyB","msy/B",
		names.msy)))))))))))))))
	rownames(tabmsy) =  paste(rep("$",nrow(tabmsy)),names.msy,rep("$",nrow(tabmsy)),sep="")
	
	cap.msy = paste0(
		"Composite base case: the ", texThatVec(tcall(quants5)), " quantiles of MCMC-derived quantities from \\Nbase", 
		" samples pooled from ", NrefM, " component runs. Definitions are: ",
		"$B_0$ -- unfished equilibrium spawning biomass (mature females), ",
		#"$V_0$ -- unfished equilibrium vulnerable biomass (males and females), ",
		"$B_{", currYear, "}$ -- spawning biomass at the end of ", currYear, ", ",
		#"$V_{", currYear, "}$ -- vulnerable biomass in the middle of ", currYear, ", ",
		"$u_{", prevYear, "}$ -- exploitation rate (ratio of total catch to vulnerable biomass) in the middle of ", prevYear, ", ",
		"$u_\\mathrm{max}$ -- maximum exploitation rate (calculated for each sample as the maximum exploitation rate from ",
		modYrs[1], "-", currYear, "), ",
		"$B_\\mathrm{MSY}$ -- equilibrium spawning biomass at MSY (maximum sustainable yield), ",
		"$u_\\mathrm{MSY}$ -- equilibrium exploitation rate at MSY, ",
		#"$V_\\mathrm{MSY}$ -- equilibrium vulnerable biomass at MSY. ",
		"All biomass values (and MSY) are in tonnes. ", refCC.sentence
	)
	
	xtab.compo.rfpt = xtable(tabmsy, align="lrrrrr",
		label=paste0("tab:",prefix,"base.rfpt"), digits=if (exists("formatCatch")) NULL else sigdig, caption=cap.msy )
	xtab.compo.rfpt.out = capture.output( print(xtab.compo.rfpt,
		caption.placement="top", sanitize.rownames.function=function(x){x}, add.to.row=list(pos=list(-1,3,11), command=c("\\\\[-1.0ex]", "\\hdashline \\\\[-1.75ex]", "\\hdashline \\\\[-1.75ex]")), hline.after=c(-1,0,5,nrow(xtab.compo.rfpt)) ) )
	tput(xtab.compo.rfpt.out)

#browser();return()

	## Component run likelihoods
	## -------------------------
	LL.compo =data.frame(
		Run = c(77L, 71L, 75L, 72L, 76L),
		M = c(0.04, 0.045, 0.05, 0.055, 0.06),
		CPUE = c(-18.44, -18.29, -18.06, -17.77, -17.41),
		QCS = c(1.275, 1.060, 0.870, 0.703, 0.555),
		WCVI = c(7.863, 7.898, 7.920, 7.934, 7.939),
		WCHG = c(20.36, 20.00, 19.68, 19.40, 19.14),
		GIG = c(14.33, 14.41, 14.48, 14.54, 14.57),
		Index = c(25.39, 25.08, 24.89, 24.80, 24.80),
		AF = c(453.6, 456.1, 456.4, 457.2, 457.5),
		Recruit = c(47.46, 43.51, 41.93, 40.80, 39.93),
		Total = c(638.5, 636.6, 635.1, 634.7, 634.0)
	)
#browser();return()
	xtab.cruns.ll = xtable(formatCatch(LL.compo), align=paste0("l", paste0(rep("r",dim(LL.compo)[2]),collapse="")),
		label   = paste0("tab:",prefix,"log.likes"), digits=NULL, 
		caption = "Log likelihood (LL) values reported by component base runs for survey indices, age composition (AF), recruitment, and total (not all LL components reported here)")
	xtab.cruns.ll.out = capture.output(print(xtab.cruns.ll,  include.rownames=FALSE,
		caption.placement="top", sanitize.rownames.function=function(x){x}, add.to.row =list(pos = list(-1), command = c("\\\\[-1.0ex]")) ) )
	tput(xtab.cruns.ll.out)

	## PJS wants to see parmeter comparisons among component runs (of course he does)
	## ----------------------------------------------------------
	## Calculate quantiles for later (eventually make this function global)
	calcQs = function (dat, ivec, ovec, qval=quants3) {
		F.data = as.data.frame(dat)
		F.list = split(F.data, ivec)
		rr.ord = match(ovec, substring(names(F.list),1,2))  ## oder not so important for sensitivities
		F.list = F.list[rr.ord]
		F.qnts = lapply(F.list,function(x){
			z = apply(x,2,function(xx){!any(is.na(xx))})
			out = apply(x[z],2,quantile,qval)
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
		t(t(apply(mat, 1, function(x) {
			paste0(c("|", c(x[1], "| ", x[2], " (", x[3],",",x[4],")")), collapse="")
			#paste0(c(x[1], " , ", x[2], " (", x[3],",",x[4],")"), collapse="")
		})))
	})
	P.tab = array(NA,dim=rev(dim(ampdPA)), dimnames=list(names.pars,names(P.qnts)))
	for (i in names(P.tabs)){
		P.tab[,i] = P.tabs[[i]]
	}
	rownames(P.tab) =  paste(rep("$",nrow(P.tab)),  gsub("\\s*~\\(.*?\\)$","",rownames(P.tab)), rep("$",nrow(P.tab)),sep="")
	colnames(P.tab) = B.labs

	xtab.cruns.pars = xtable(P.tab, align=paste0("l", paste0(rep("c",dim(P.tab)[2]),collapse="")),
		label   = paste0("tab:",prefix,"runs.pars"), digits = if (exists("formatCatch")) NULL else sigdig,
		caption = paste0("Component base case runs: model parameter MPDs (delimited by `|') and MCMC medians (with 0.05 and 0.95 quantile limits) for each of the ",
		"\\numberstringnum{", NrefM, "} component model runs of \\Nmcmc{} samples each.") )
	xtab.cruns.pars.out = capture.output(print(xtab.cruns.pars,
		caption.placement="top", sanitize.rownames.function=function(x){x}, add.to.row =list(pos = list(-1), command = c("\\\\[-1.0ex]")) ) )
	tput(xtab.cruns.pars.out)

	## PJS ALSO wants to see comparisons among derived quantities for component runs (of course he does)
	## ----------------------------------------------------------
	Q.qnts = calcQs(B.mcmc, ivec=P.run.rwt, ovec=P.run.ord)
	Q.med  = lapply(Q.qnts,function(x){apply(x,2,function(xx){c(xx[2],xx[1],xx[3])})})
#browser();return()

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
	rownames(Q.tab) = names.msy
	rownames(Q.tab) = paste(rep("$",nrow(Q.tab)),  gsub("\\s*~\\(.*?\\)$","",rownames(Q.tab)), rep("$",nrow(Q.tab)),sep="")
	colnames(Q.tab) = B.labs

	xtab.cruns.rfpt = xtable(Q.tab, align=paste0("l", paste0(rep("r",dim(Q.tab)[2]),collapse="")),
		label   = paste0("tab:",prefix,"runs.rfpt"), digits = if (exists("formatCatch")) NULL else sigdig,
		caption = paste0("Component base case runs: MCMC median (with 0.05 and 0.95 quantile limits) for derived model quantities for each of the ",
		"\\numberstringnum{", NrefM, "} component model runs of \\Nmcmc{} samples each.") )
	xtab.cruns.rfpt.out = capture.output(print(xtab.cruns.rfpt,
		caption.placement="top", sanitize.rownames.function=function(x){x}, add.to.row=list(pos=list(-1,3,11), command=c("\\\\[-1.0ex]", "\\hdashline \n", "\\hdashline \n")), hline.after=c(-1,0,5,nrow(xtab.cruns.rfpt)) ) )
	tput(xtab.cruns.rfpt.out)
#browser();return()

	save("xtab.compo.pars.out", "xtab.compo.rfpt.out", "xtab.cruns.pars.out", "xtab.cruns.rfpt.out", "xtab.cruns.ll.out", file=paste0(prefix,"compo.tabs", ifelse(useRlow, paste0(".Rlow(q=",pad0(qRlow*100,2),")"),""), ".rda") )
#browser();return()
return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tabSS.compo


## tabSS.decision-----------------------2021-08-24
## Make Decision Tables (Probabilities)
## ---------------------------------------------RH
tabSS.decision = function(istock="YMR", prefix="ymr.", compo,
  useRlow=FALSE, qRlow=0.25, decdig=2)
{
	on.exit(gc(verbose=FALSE))
	require(xtable)
	unpackList(stock[[istock]][["Controls"]])
	unpackList(compo)

	proYears    = (currYear+1):projYear
	refYears    = as.character(currYear:projYear)  ## include current year in decision tables

	## Diagnostics for select parameters
	P.mpd       = ampdPA; storage.mode(P.mpd)="double"     ## For some reason, this matrix is stored as a list (maybe due to NA values?)
	P.names     = colnames(avgPA)
	P.runs      = as.numeric(sapply(strsplit(rownames(avgPA),"\\."),function(x){x[1]}))
	P.run.rwt   = sapply(strsplit(rownames(avgPA),"\\."),function(x){paste0(x[1],".",x[2])}) ## vector of Runs and Rwts
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
	#lab.recentCatchYears = texThatVec(rev(rev(rownames(catch))[1:num.recentCatchYears]))
	recentCatchMean = mean(recentCatch); num.recentCatchYears=length(recentCatch)
	recentHarvMean  = 0.01 ## dummy for now

	lab.recentCatchYears = texThatVec(names(recentCatch))
	refCC.sentence = paste0(" For reference, the average catch over the last ", num.recentCatchYears, " years (", lab.recentCatchYears, ") was ", round(recentCatchMean, dig=0), "~t. ")
	refHR.sentence = paste0(" For reference, the average harvest rate over the last ", num.recentCatchYears, " years (", lab.recentCatchYears, ") was ", round(recentHarvMean, dig=3), ". ")

	## Determine subset vector for low recruitment events
	R.mcmc = avgTS[,,"Rt"]
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
#browser();return()

#	## Subset vector objects for composite base case (same code as in 'make.compo.figs.r')
#	vec.nams = c("B0", "Bcurr", "Bcurr.B0", "MSY", "Bmsy", "LRP", "USR", "Bcurr.Bmsy", "Bmsy.B0", "VBmsy", "Umsy")
#	for (v in vec.nams)
#		eval(parse(text=paste0(v, " = ", v, "[Rsub]")))
#	vec.list = c("VB0", "VBcurr", "VBcurr.VB0", "Ucurr", "Umax", "VBmsy.VB0", "Ucurr.Umsy")
#	for (vv in vec.list)
#		eval(parse(text = paste0(vv, " = lapply(", vv, ",function(x){ x[Rsub] })")))
#
#	## Subset data.frames for composite base case
#	df.nams = c("B.mcmc", "R.mcmc", "P.mcmc", "BoverBmsy")
#	for (d in df.nams)
#		eval(parse(text=paste0(d, " = ", d, "[Rsub,]")))
#	df.list = c("U.mcmc", "VB.mcmc","UoverUmsy")
#	for (dd in df.list)
#		eval(parse(text = paste0(dd, " = lapply(", dd, ",function(x){ x[Rsub,] })")))
#
#	B.nams = names(B.pols)
#	B.pols = lapply(names(B.pols),function(i) { A=B.pols[[i]]; Asub=lapply(names(A),function(j) { AA=A[[j]]; AA[Rsub,] }); names(Asub)=names(A); return(Asub) })
#	R.pols = lapply(names(R.pols),function(i) { A=R.pols[[i]]; Asub=lapply(names(A),function(j) { AA=A[[j]]; AA[Rsub,] }); names(Asub)=names(A); return(Asub) })
#	U.pols = lapply(names(U.pols),function(i) { A=U.pols[[i]]; Asub=lapply(names(A),function(j) { AA=A[[j]]; AA[Rsub,] }); names(Asub)=names(A); return(Asub) })
#	names(B.pols) = names(R.pols) = names(U.pols) = B.nams
#
#	B.proj = lapply(1:length(policies), function(i){ B.pols[[names(policies)[i]]][[as.character(policies[i])]] })
#	R.proj = lapply(1:length(policies), function(i){ R.pols[[names(policies)[i]]][[as.character(policies[i])]] })
#	U.proj = lapply(1:length(policies), function(i){ U.pols[[names(policies)[i]]][[as.character(policies[i])]] })
#	names(B.proj) = names(R.proj) = names(U.proj) = names(policies)
#browser();return()

	## -----------------------------------------------
	## 'refProbs3GenList' contains probabilities for the main
	##   decision tables (usually only use MSY-based RefPts).
	## Will include Bmsy and B0 RefPts already calculated above,
	##   e.g., refProbs3GenList$'0.4Bmsy' - refProbs$LRP = matrix of 0's
	## -----------------------------------------------
	so("findTarget.r","synth")

	## This has been taken care of above
	#if (useRlow){
	#	Bmsy   = Bmsy[Rlow]
	#	Bcurr  = Bcurr[Rlow]
	#	B0     = B0[Rlow]
	#	B.mcmc = B.mcmc[Rlow,]
	#	Umsy   = Umsy[Rlow]
	#	Ucurr  = lapply(Ucurr,function(x){x[Rlow]})
	#	B.pols = lapply(actpol, function(p){lapply(1:length(B.pols[[p]]),function(i,x){ x[[i]][Rlow,] }, x=B.pols[[p]])})
	#	U.pols = lapply(actpol, function(p){lapply(1:length(U.pols[[p]]),function(i,x){ x[[i]][Rlow,] }, x=U.pols[[p]])})
	#	names(B.pols) = names(U.pols) = actpol
	#	junk = lapply(actpol,function(p){names(B.pols[[p]]) <<- policy[[p]]$projPolicy})
	#	junk = lapply(actpol,function(p){names(U.pols[[p]]) <<- policy[[p]]$projPolicy})
	#}
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
	CP = round(sapply(avgCP[1,1,],mean,na.rm=T)[cp])  ## should already have been rounded in 'gatherMCMC'

	## create table for B < Bmsy (or whatever)
	## ---------------------------------------
	for (k in c("CC","HR")){
		BRPpList[[k]] = list()
		for(i in 1:length(cp)) {
			ii = cp[i]
			kk = ifelse (grepl("CC",ii), "CC", "HR")
			Bmat = avgPJ[,as.character(proYears),"Bt",ii]  ## subset years until/if we can project further
			Btab = data.frame(Bcurr,Bmat)
			colnames(Btab)= currYear:projYear
			BRPp.temp = sapply(Tlst,function(x){findTarget(Btab, ratio=x$ratio, target=x$target, retVal="p.hi", op=">", yrG=Ngen*gen1)}, simplify=FALSE)
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
	for (k in c("CC","HR")){
		URPpList[[k]] = list()
		for(i in 1:length(cp)) {
			ii = cp[i]
			kk = ifelse (grepl("CC",ii), "CC", "HR")
			Umat = avgPJ[,as.character(proYears),"ut",ii]  ## subset years until/if we can project further
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
	#{
	#	## Also calculate number of years to reach reference target points with a specified confidence
	#	Ttab.temp = Utab.temp = list()
	#	for (j in c(0.50,0.65,0.80,0.95)) {
	#		jj = formatC(j,digits=2,format="f")
	#		Ttab.temp[[jj]]  = pmin(sapply(Tlst, function(x){sapply(B.pols[[i]], findTarget, ratio=x$ratio, target=x$target, conf=j,  retVal="N", op=">", yrG=Ngen*gen1)}), Ngen*gen1)
	#		Utab.temp[[jj]]  = pmin(sapply(Ulst, function(x){sapply(U.pols[[i]], findTarget, ratio=x$ratio, target=x$target, conf=j,  retVal="N", op="<", yrG=Ngen*gen1)}), Ngen*gen1)
	#	}
	#	Ttab.conf[[i]] = Ttab.temp
	#	Utab.conf[[i]] = Utab.temp
	#}
	##=====TABLE: DEFAULT DFO MSY REFERENCE POINTS=====
	## Just use and adapt code from runSweaveMCCMC.Rnw
	projYearsNum = length(proYears)
	for.area = paste0(" (",gsub(")","]",gsub("\\(","[",name)),")")
	maxCatSentence = ""  ## only ever used for the ROL assessment (see run-masterMCMC.Snw)
	Nmcmc = sum(Rsub)

	##=====GMU -- Guidance for Setting TAC=====

	adjUyrs = function(tab){
		## adjust the exploitation years to reflect catch before start of current year
		colnames(tab) = as.character(as.numeric(colnames(tab))-1)
		return(tab)
	}

	##--- LRP P(Bt>0.4Bmsy) ----------------
	if (exists("BRPpList") && !is.null(BRPpList$CC$'0.4Bmsy')) {
		xtabLRP = texArray(BRPpList$CC$'0.4Bmsy'[,refYears], table.label=paste0("tab:",prefix,"gmu.LRP.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the limit reference point $0.4 \\Bmsy$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{constant catch} strategies (in tonnes). Values are P$(B_t > 0.4 \\Bmsy)$, i.e.~the probability of the spawning biomass (mature females) at the start of year $t$ being greater than the limit reference point. The probabilities are the proportion (to two decimal places) of the ", Nmcmc, " MCMC samples for which $B_t > 0.4 \\Bmsy$. ", refCC.sentence, maxCatSentence)) )
		gmu.LRP.CCs.out = xtabLRP$tabfile
		tput(gmu.LRP.CCs.out)
	}
#browser();return()
	if (exists("BRPpList") && !is.null(BRPpList$HR$'0.4Bmsy')) {
		xtabLRP2 = texArray(BRPpList$HR$'0.4Bmsy'[,refYears], table.label=paste0("tab:",prefix,"gmu.LRP.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the limit reference point $0.4 \\Bmsy$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{harvest rate} strategies. Values are P$(B_t > 0.4 \\Bmsy)$, i.e.~the probability of the spawning biomass (mature females) at the start of year $t$ being greater than the limit reference point. The probabilities are the proportion (to two decimal places) of the ", Nmcmc, " MCMC samples for which $B_t > 0.4 \\Bmsy$. ", refHR.sentence, maxCatSentence)) )
		gmu.LRP.HRs.out = xtabLRP2$tabfile
		tput(gmu.LRP.HRs.out)
	}
	##--- USR P(Bt>0.8Bmsy) ----------------
	if (exists("BRPpList") && !is.null(BRPpList$CC$'0.8Bmsy')) {
		xtabUSR = texArray(BRPpList$CC$'0.8Bmsy'[,refYears], table.label=paste0("tab:",prefix,"gmu.USR.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the upper stock reference point $0.8 \\Bmsy$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{constant catch} strategies (in tonnes), such that values are P$(B_t > 0.8 \\Bmsy)$. ", refCC.sentence, maxCatSentence)) )
		gmu.USR.CCs.out = xtabUSR$tabfile
		tput(gmu.USR.CCs.out)
	}
	if (exists("BRPpList") && !is.null(BRPpList$HR$'0.8Bmsy')) {
		xtabUSR2 = texArray(BRPpList$HR$'0.8Bmsy'[,refYears], table.label=paste0("tab:",prefix,"gmu.USR.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the upper stock reference point $0.8 \\Bmsy$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.8 \\Bmsy)$. ", refHR.sentence, maxCatSentence)) )
		gmu.USR.HRs.out = xtabUSR2$tabfile
		tput(gmu.USR.HRs.out)
	}
	##--- Bmsy P(Bt>Bmsy) ------------------
	if (exists("BRPpList") && !is.null(BRPpList$CC$'Bmsy')) {
		xtabBmsy = texArray(BRPpList$CC$'Bmsy'[,refYears], table.label=paste0("tab:",prefix,"gmu.Bmsy.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the reference point $\\Bmsy$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{constant catch} strategies (in tonnes), such that values are P$(B_t > \\Bmsy)$. ", refCC.sentence, maxCatSentence)) )
		gmu.Bmsy.CCs.out = xtabBmsy$tabfile
		tput(gmu.Bmsy.CCs.out)
	}
	if (exists("BRPpList") && !is.null(BRPpList$HR$'Bmsy')) {
		xtabBmsy2 = texArray(BRPpList$HR$'Bmsy'[,refYears], table.label=paste0("tab:",prefix,"gmu.Bmsy.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the reference point $\\Bmsy$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > \\Bmsy)$. ", refHR.sentence, maxCatSentence)) )
		gmu.Bmsy.HRs.out = xtabBmsy2$tabfile
		tput(gmu.Bmsy.HRs.out)
	}
	##--- Umsy P(ut>umsy) ------------------
	if (exists("URPpList") && !is.null(URPpList$CC$'umsy')) {
		xtabUmsy = texArray(adjUyrs(URPpList$CC$'umsy'[,refYears]), table.label=paste0("tab:",prefix,"gmu.umsy.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the reference point $\\umsy$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{constant catch} strategies, such that values are P$(u_t < \\umsy)$. ", refCC.sentence, maxCatSentence)) )
		gmu.umsy.CCs.out = xtabUmsy$tabfile
		tput(gmu.umsy.CCs.out)
	}
	if (exists("URPpList") && !is.null(URPpList$HR$'umsy')) {
		xtabUmsy2 = texArray(adjUyrs(URPpList$HR$'umsy'[,refYears]), table.label=paste0("tab:",prefix,"gmu.umsy.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the reference point $\\umsy$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{harvest rate} strategies, such that values are P$(u_t < \\umsy)$. ", refHR.sentence, maxCatSentence)) )
		gmu.umsy.HRs.out = xtabUmsy2$tabfile
		tput(gmu.umsy.HRs.out)
	}
	##--- Bcurr P(Bt>Bcurr) ----------------
	if (exists("BRPpList") && !is.null(BRPpList$CC$'Bcurr')) {
		xtabBcurr = texArray(BRPpList$CC$'Bcurr'[,refYears], table.label=paste0("tab:",prefix,"gmu.Bcurr.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the reference point $B_{\\currYear}$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > B_{\\currYear})$. ", refCC.sentence, maxCatSentence)) )
		gmu.Bcurr.CCs.out = xtabBcurr$tabfile
		tput(gmu.Bcurr.CCs.out)
	}
	if (exists("BRPpList") && !is.null(BRPpList$HR$'Bcurr')) {
		xtabBcurr2 = texArray(BRPpList$HR$'Bcurr'[,refYears], table.label=paste0("tab:",prefix,"gmu.Bcurr.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the reference point $B_{\\currYear}$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > B_{\\currYear})$. ", refHR.sentence, maxCatSentence)) )
		gmu.Bcurr.HRs.out = xtabBcurr2$tabfile
		tput(gmu.Bcurr.HRs.out)
	}
	##--- ucurr P(ut>ucurr) ----------------
	if (exists("URPpList") && !is.null(URPpList$CC$'ucurr')) {
		xtabUcurr = texArray(adjUyrs(URPpList$CC$'ucurr'[,refYears]), table.label=paste0("tab:",prefix,"gmu.ucurr.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the reference point $u_{\\prevYear}$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{constant catch} strategies, such that values are P$(u_t < u_{\\prevYear})$. ", refCC.sentence, maxCatSentence)) )
		gmu.ucurr.CCs.out = xtabUcurr$tabfile
		tput(gmu.ucurr.CCs.out)
	}
	if (exists("URPpList") && !is.null(URPpList$HR$'ucurr')) {
		xtabUcurr2 = texArray(adjUyrs(URPpList$HR$'ucurr'[,refYears]), table.label=paste0("tab:",prefix,"gmu.ucurr.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for the reference point $u_{\\prevYear}$ featuring current- and ",projYearsNum,"-year projections for a range of \\itbf{harvest rate} strategies, such that values are P$(u_t < u_{\\prevYear})$. ", refHR.sentence, maxCatSentence)) )
		gmu.ucurr.HRs.out = xtabUcurr2$tabfile
		tput(gmu.ucurr.HRs.out)
	}
	##--- Alt.LRP P(Bt>0.2B0) --------------
	if(exists("BRPpList") && !is.null(BRPpList$CC$'0.2B0')) {
		xtab20B0 = texArray(BRPpList$CC$'0.2B0'[,refYears], table.label=paste0("tab:",prefix,"gmu.20B0.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for an alternative reference point $0.2 B_0$ featuring current- and ",projYearsNum," year projections for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.2 B_0)$. ", refCC.sentence, maxCatSentence)) )
		gmu.20B0.CCs.out = xtab20B0$tabfile
		tput(gmu.20B0.CCs.out)
	}
	if(exists("BRPpList") && !is.null(BRPpList$HR$'0.2B0')) {
		xtab20B02 = texArray(BRPpList$HR$'0.2B0'[,refYears], table.label=paste0("tab:",prefix,"gmu.20B0.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for an alternative reference point $0.2 B_0$ featuring current- and ",projYearsNum," year projections for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.2 B_0)$. ", refHR.sentence, maxCatSentence)) )
		gmu.20B0.HRs.out = xtab20B02$tabfile
		tput(gmu.20B0.HRs.out)
	}
	##--- Alt.USR P(Bt>0.4B0) --------------
	if(exists("BRPpList") && !is.null(BRPpList$CC$'0.4B0')) {
		xtab40B0 = texArray(BRPpList$CC$'0.4B0'[,refYears], table.label=paste0("tab:",prefix,"gmu.40B0.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for an alternative reference point $0.4 B_0$ featuring current- and ",projYearsNum," year projections for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.4 B_0)$. ", refCC.sentence, maxCatSentence)) )
		gmu.40B0.CCs.out = xtab40B0$tabfile
		tput(gmu.40B0.CCs.out)
	}
	if(exists("BRPpList") && !is.null(BRPpList$HR$'0.4B0')) {
		xtab40B02 = texArray(BRPpList$HR$'0.4B0'[,refYears], table.label=paste0("tab:",prefix,"gmu.40B0.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0,
			table.caption=eval(paste0(name, ": decision table for an alternative reference point $0.4 B_0$ featuring current- and ",projYearsNum," year projections for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.4 B_0)$. ", refHR.sentence, maxCatSentence)) )
		gmu.40B0.HRs.out = xtab40B02$tabfile
		tput(gmu.40B0.HRs.out)
	}
#	##=====COSEWIC -- Reference Criteria=====
	refYears     = as.character(currYear:pgenYear)
	pgenYearsNum = length(refYears)
	allYrs       = as.character(currYear + seq(0, pgenYearsNum-1, 1))
	Npyrs        = length(currYear:projYear)-1
	shortYrs     = as.character(currYear + seq(0,Npyrs,1))
	shortYrs     = intersect(shortYrs,allYrs)

	##--- A2 criterion 0.5B0 COSEWIC short --------
	if (exists("BRPpList") && !is.null(BRPpList$CC$'0.5B0')) {
		mess = paste0(name, ": decision table for COSEWIC reference criterion A2 `Endangered' featuring current-", Npyrs, "-year projections and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.5 B_0)$.", refCC.sentence)
		cosewic.50B0.CCs.out = texArray(BRPpList$CC$'0.5B0'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50B0.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
		tput(cosewic.50B0.CCs.out)
	}
	##--- A2 criterion 0.7B0 COSEWIC short --------
	if (exists("BRPpList") && !is.null(BRPpList$CC$'0.7B0')) {
		mess = paste0(name, ": decision table for COSEWIC reference criterion A2 `Threatened' featuring current-", Npyrs, "-year projections and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.7 B_0)$.", refCC.sentence)
		cosewic.70B0.CCs.out = texArray(BRPpList$CC$'0.7B0'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.70B0.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
		tput(cosewic.70B0.CCs.out)
	}
	##--- A2 decline <=30pct COSEWIC short --------
	if (exists("BRPpList") && !is.null(BRPpList$CC$'0.3Gen')) {
		mess = paste0(name, ": probability of satisfying the A2 criterion of $\\leq 30 \\%$ decline from ", Ngen, " generations (", Ngen*gen1, " years) earlier featuring current- and ", Npyrs, "-year projections and for a range of \\itbf{constant catch} strategies. ", refCC.sentence, maxCatSentence)
		cosewic.30Gen.CCs.out = texArray(BRPpList$CC$'0.3Gen'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.30Gen.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
		tput(cosewic.30Gen.CCs.out)
	}
	##--- A2 decline <=50pct COSEWIC short --------
	if (exists("BRPpList") && !is.null(BRPpList$CC$'0.5Gen')) {
		mess = paste0(name, ": probability of satisfying the A2 criterion of $\\leq 50 \\%$ decline from ", Ngen, " generations (", Ngen*gen1, " years) earlier featuring current- and ", Npyrs, "-year projections and for a range of \\itbf{constant catch} strategies. ", refCC.sentence, maxCatSentence)
		cosewic.50Gen.CCs.out = texArray(BRPpList$CC$'0.5Gen'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50Gen.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
		tput(cosewic.50Gen.CCs.out)
	}
#browser();return()
	save(
	"gmu.LRP.CCs.out", "gmu.USR.CCs.out", "gmu.Bmsy.CCs.out", "gmu.umsy.CCs.out",
	"gmu.Bcurr.CCs.out", "gmu.ucurr.CCs.out", "gmu.20B0.CCs.out", "gmu.40B0.CCs.out",
	"cosewic.50B0.CCs.out", "cosewic.70B0.CCs.out", "cosewic.30Gen.CCs.out", "cosewic.50Gen.CCs.out",
	file=paste0(prefix,"decision.tables", ifelse(useRlow, paste0(".Rlow(q=",pad0(qRlow*100,2),")"),""), ".rda")
	)


#	##=====GMU -- Guidance for Rebuilding=====
#
#	## Note: will need to modify labels when more than one stock but leave as is for BOR 2019.
#
#	refYears = as.character(currYear:pgenYear)
#	pgenYearsNum = length(refYears)
#	if (pgenYearsNum > 11) {
#		allYrs = as.character(currYear + seq(0, pgenYearsNum-1, 1))
#		Npyrs  = length(currYear:projYear)-1
#
#		shortYrs = as.character(currYear + seq(0,Npyrs,1))
#		shortYrs  = intersect(shortYrs,allYrs)
#	
#		## Assumes 3 generations:
#		## longYrs  = as.character(currYear + c(seq(0,round(2*gen1),round(gen1/4)),seq(2*gen1,3*gen1,gen1/2)[-1]))
#		## longYrs  = intersect(longYrs,allYrs)
#
#		## Usual shuffling shit:
#		yrs5  = intersect(seq(2000,3000,5),currYear:pgenYear)
#		yrs10 = rev(seq(rev(yrs5)[1],currYear,-10))
#		yrs5  = rev(seq(rev(yrs5)[1],currYear,-5))
#		longYrs = as.character(c(currYear, yrs5[1:5],setdiff(yrs10,yrs5[1:5])))
#
##browser();return()
#
#		##--- LRP Rebuild P(Bt>0.4Bmsy) ----------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.4Bmsy')) {
#			mess = paste0(name, ": decision table for the limit reference point $0.4 \\Bmsy$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.4 \\Bmsy)$.", refCC.sentence)
#			gmu.LRP.CCl.out = texArray(BRPpList$CC$'0.4Bmsy'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.LRP.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(gmu.LRP.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.4Bmsy')) {
#			mess = paste0(name, ": decision table for the limit reference point $0.4 \\Bmsy$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.4 \\Bmsy)$.", refHR.sentence)
#			gmu.LRP.HRl.out = texArray(BRPpList$HR$'0.4Bmsy'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.LRP.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(gmu.LRP.HRl.out)
#		}
#		##--- USR Rebuild P(Bt>0.8Bmsy) ----------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.8Bmsy')) {
#			mess = paste0(name, ": decision table for the upper stock reference $0.8 \\Bmsy$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.8 \\Bmsy)$.", refCC.sentence)
#			gmu.USR.CCl.out = texArray(BRPpList$CC$'0.8Bmsy'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.USR.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(gmu.USR.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.8Bmsy')) {
#			mess = paste0(name, ": decision table for the upper stock reference $0.8 \\Bmsy$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.8 \\Bmsy)$.", refHR.sentence)
#			gmu.USR.HRl.out = texArray(BRPpList$HR$'0.8Bmsy'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.USR.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(gmu.USR.HRl.out)
#		}
#		##--- Bmsy Rebuild P(Bt>Bmsy) ----------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'Bmsy')) {
#			mess = paste0(name, ": decision table for biomass at maximum sustainable yield $\\Bmsy$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > \\Bmsy)$.", refCC.sentence)
#			gmu.Bmsy.CCl.out = texArray(BRPpList$CC$'Bmsy'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.Bmsy.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(gmu.Bmsy.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'Bmsy')) {
#			mess = paste0(name, ": decision table for biomass at maximum sustainable yield $\\Bmsy$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > \\Bmsy)$.", refHR.sentence)
#			gmu.Bmsy.HRl.out = texArray(BRPpList$HR$'Bmsy'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.Bmsy.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(gmu.Bmsy.HRl.out)
#		}
#		##--- umsy Rebuild P(ut>umsy) ----------
#		if (exists("URPpList") && !is.null(URPpList$CC$'umsy')) {
#			mess = paste0(name, ": decision table for harvest rate at maximum sustainable yield $\\umsy$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(u_t < \\umsy)$.", refCC.sentence)
#			gmu.umsy.CCl.out = texArray(URPpList$CC$'umsy'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.umsy.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(gmu.umsy.CCl.out)
#		}
#		if (exists("URPpList") && !is.null(URPpList$HR$'umsy')) {
#			mess = paste0(name, ": decision table for harvest rate at maximum sustainable yield $\\umsy$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(u_t < \\umsy)$.", refHR.sentence)
#			gmu.umsy.HRl.out = texArray(URPpList$HR$'umsy'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.umsy.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(gmu.umsy.HRl.out)
#		}
#		##--- Bcurr Rebuild P(Bt>Bcurr) ----------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'Bcurr')) {
#			mess = paste0(name, ": decision table for comparing projected biomass to current biomass $B_{", currYear, "}$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > B_{", currYear, "})$.", refCC.sentence)
#			gmu.Bcurr.CCl.out = texArray(BRPpList$CC$'Bcurr'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.Bcurr.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(gmu.Bcurr.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'Bcurr')) {
#			mess = paste0(name, ": decision table for comparing projected biomass to current biomass $B_{", currYear, "}$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > B_{", currYear, "})$.", refHR.sentence)
#			gmu.Bcurr.HRl.out = texArray(BRPpList$HR$'Bcurr'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.Bcurr.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(gmu.Bcurr.HRl.out)
#		}
#		##--- ucurr Rebuild P(ut>ucurr) ----------
#		if (exists("URPpList") && !is.null(URPpList$CC$'ucurr')) {
#			mess = paste0(name, ": decision table for comparing projected harvest rate to current harvest rate $u_{", currYear-1, "}$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(u_t < u_{", currYear-1, "})$.", refCC.sentence)
#			gmu.ucurr.CCl.out = texArray(URPpList$CC$'ucurr'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.ucurr.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(gmu.ucurr.CCl.out)
#		}
#		if (exists("URPpList") && !is.null(URPpList$HR$'ucurr')) {
#			mess = paste0(name, ": decision table for comparing projected harvest rate to current harvest rate $u_{", currYear-1, "}$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(u_t < u_{", currYear-1, "})$.", refHR.sentence)
#			gmu.ucurr.HRl.out = texArray(URPpList$HR$'ucurr'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.ucurr.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(gmu.ucurr.HRl.out)
#		}
#		##--- 0.2B0 Rebuild P(Bt>0.2B0) ----------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.2B0')) {
#			mess = paste0(name, ": decision table for alternative limit reference point $0.2B_0$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.2B_0)$.", refCC.sentence)
#			gmu.20B0.CCl.out = texArray(BRPpList$CC$'0.2B0'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.20B0.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(gmu.20B0.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.2B0')) {
#			mess = paste0(name, ": decision table for alternative limit reference point $0.2B_0$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.2B_0)$.", refHR.sentence)
#			gmu.20B0.HRl.out = texArray(BRPpList$HR$'0.2B0'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.20B0.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(gmu.20B0.HRl.out)
#		}
#		##--- 0.4B0 Rebuild P(Bt>0.4B0) ----------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.4B0')) {
#			mess = paste0(name, ": decision table for alternative upper stock reference $0.4B_0$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.4B_0)$.", refCC.sentence)
#			gmu.40B0.CCl.out = texArray(BRPpList$CC$'0.4B0'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.40B0.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(gmu.40B0.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.4B0')) {
#			mess = paste0(name, ": decision table for alternative upper stock reference $0.4B_0$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.4B_0)$.", refHR.sentence)
#			gmu.40B0.HRl.out = texArray(BRPpList$HR$'0.4B0'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"gmu.40B0.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(gmu.40B0.HRl.out)
#		}
#	}
#	##=====COSEWIC -- Reference Criteria=====
#	if (pgenYearsNum > 11) {
#		##--- A2 decline <=50pct COSEWIC short --------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.5Gen')) {
#			mess = paste0(name, ": decision table for probabilities of satisfying the A2  criterion of $\\leq 50 \\%$ decline over ", Ngen, " generations (", Ngen*gen1, " years) for ", Npyrs, "-year projections and for a range of \\itbf{constant catch} strategies. ", refCC.sentence, maxCatSentence)
#			cosewic.50Gen.CCs.out = texArray(BRPpList$CC$'0.5Gen'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50Gen.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(cosewic.50Gen.CCs.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.5Gen')) {
#			mess = paste0(name, ": decision table for probabilities of satisfying the A2 criterion of $\\leq 50 \\%$ decline over ", Ngen, " generations (", Ngen*gen1, " years) for ", Npyrs, "-year projections and for a range of \\itbf{harvest rate} strategies. ", refHR.sentence, maxCatSentence)
#			cosewic.50Gen.HRs.out = texArray(BRPpList$HR$'0.5Gen'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50Gen.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(cosewic.50Gen.HRs.out)
#		}
#		##--- A2 decline <=30pct COSEWIC short --------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.3Gen')) {
#			mess = paste0(name, ": decision table for probabilities of satisfying the A2 criterion of $\\leq 30 \\%$ decline over ", Ngen, " generations (", Ngen*gen1, " years) for ", Npyrs, "-year projections and for a range of \\itbf{constant catch} strategies. ", refCC.sentence, maxCatSentence)
#			cosewic.30Gen.CCs.out = texArray(BRPpList$CC$'0.3Gen'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.30Gen.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(cosewic.30Gen.CCs.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.3Gen')) {
#			mess = paste0(name, ": decision table for probabilities of satisfying the A2 criterion of $\\leq 30 \\%$ decline over ", Ngen, " generations (", Ngen*gen1, " years) for ", Npyrs, "-year projections and for a range of \\itbf{harvest rate} strategies. ", refHR.sentence, maxCatSentence)
#			cosewic.30Gen.HRs.out = texArray(BRPpList$HR$'0.3Gen'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.30Gen.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(cosewic.30Gen.HRs.out)
#		}
#		##--- A2 criterion 0.5B0 COSEWIC short --------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.5B0')) {
#			mess = paste0(name, ": decision table for reference criterion $0.5 B_0$ for ", Npyrs, "-year projections and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.5 B_0)$.", refCC.sentence)
#			cosewic.50B0.CCs.out = texArray(BRPpList$CC$'0.5B0'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50B0.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(cosewic.50B0.CCs.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.5B0')) {
#			mess = paste0(name, ": decision table for reference criterion $0.5 B_0$ for ", Npyrs, "-year projections and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.5 B_0)$.", refHR.sentence)
#			cosewic.50B0.HRs.out = texArray(BRPpList$HR$'0.5B0'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50B0.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(cosewic.50B0.HRs.out)
#		}
#		##--- A2 criterion 0.7B0 COSEWIC short --------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.7B0')) {
#			mess = paste0(name, ": decision table for reference criterion $0.7 B_0$ for ", Npyrs, "-year projections and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.7 B_0)$.", refCC.sentence)
#			cosewic.70B0.CCs.out = texArray(BRPpList$CC$'0.7B0'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.70B0.CCs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(cosewic.70B0.CCs.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.7B0')) {
#			mess = paste0(name, ": decision table for reference criterion $0.7 B_0$ for ", Npyrs, "-year projections and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.7 B_0)$.", refHR.sentence)
#			cosewic.70B0.HRs.out = texArray(BRPpList$HR$'0.7B0'[,shortYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.70B0.HRs"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(cosewic.70B0.HRs.out)
#		}
#		##--- A2 decline <=50pct COSEWIC long ---------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.5Gen')) {
#			mess = paste0(name, ": decision table for probabilities of satisfying the A2  criterion of $\\leq 50 \\%$ decline over ", Ngen, " generations (", Ngen*gen1, " years) for selected projection years and for a range of \\itbf{constant catch} strategies. ", refCC.sentence, maxCatSentence)
#			cosewic.50Gen.CCl.out = texArray(BRPpList$CC$'0.5Gen'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50Gen.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(cosewic.50Gen.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.5Gen')) {
#			mess = paste0(name, ": decision table for probabilities of satisfying the A2 criterion of $\\leq 50 \\%$ decline over", Ngen, " generations (", Ngen*gen1, " years) for selected projection years and for a range of \\itbf{harvest rate} strategies. ", refHR.sentence, maxCatSentence)
#			cosewic.50Gen.HRl.out = texArray(BRPpList$HR$'0.5Gen'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50Gen.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(cosewic.50Gen.HRl.out)
#		}
#		##--- A2 decline <=30pct COSEWIC long ---------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.3Gen')) {
#			mess = paste0(name, ": decision table for probabilities of satisfying the A2 criterion of $\\leq 30 \\%$ decline over ", Ngen, " generations (", Ngen*gen1, " years) for selected projection years and for a range of \\itbf{constant catch} strategies. ", refCC.sentence, maxCatSentence)
#			cosewic.30Gen.CCl.out = texArray(BRPpList$CC$'0.3Gen'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.30Gen.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(cosewic.30Gen.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.3Gen')) {
#			mess = paste0(name, ": decision table for probabilities of satisfying the A2 criterion of $\\leq 30 \\%$ decline over ", Ngen, " generations (", Ngen*gen1, " years) for selected projection years and for a range of \\itbf{harvest rate} strategies. ", refHR.sentence, maxCatSentence)
#			cosewic.30Gen.HRl.out = texArray(BRPpList$HR$'0.3Gen'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.30Gen.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(cosewic.30Gen.HRl.out)
#		}
#		##--- A2 criterion 0.5B0 COSEWIC long ---------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.5B0')) {
#			mess = paste0(name, ": decision table for reference criterion $0.5 B_0$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that valuTabes are P$(B_t > 0.5 B_0)$.", refCC.sentence)
#			cosewic.50B0.CCl.out = texArray(BRPpList$CC$'0.5B0'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50B0.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(cosewic.50B0.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.5B0')) {
#			mess = paste0(name, ": decision table for reference criterion $0.5 B_0$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.5 B_0)$.", refHR.sentence)
#			cosewic.50B0.HRl.out = texArray(BRPpList$HR$'0.5B0'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.50B0.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(cosewic.50B0.HRl.out)
#		}
#		##--- A2 criterion 0.7B0 COSEWIC long ---------
#		if (exists("BRPpList") && !is.null(BRPpList$CC$'0.7B0')) {
#			mess = paste0(name, ": decision table for reference criterion $0.7 B_0$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{constant catch} strategies, such that values are P$(B_t > 0.7 B_0)$.", refCC.sentence)
#			cosewic.70B0.CCl.out = texArray(BRPpList$CC$'0.7B0'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.70B0.CCl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="CC", tablewidth=6.0)$tabfile
#			tput(cosewic.70B0.CCl.out)
#		}
#		if (exists("BRPpList") && !is.null(BRPpList$HR$'0.7B0')) {
#			mess = paste0(name, ": decision table for reference criterion $0.7 B_0$ for selected projection years over ", Ngen, " generations (", Ngen*gen1, " years) and for a range of \\itbf{harvest rate} strategies, such that values are P$(B_t > 0.7 B_0)$.", refHR.sentence)
#			cosewic.70B0.HRl.out = texArray(BRPpList$HR$'0.7B0'[,longYrs], table.caption=eval(mess), table.label=paste0("tab:",prefix,"cosewic.70B0.HRl"), sigdig=decdig, use.round=T, zero="0", use.row.names=T, name.row.names="HR", tablewidth=6.0)$tabfile
#			tput(cosewic.70B0.HRl.out)
#		}
#	}
#	##=====Time to Targets=====
#	if (pgenYearsNum > 11) {
#		sanicol = function(x){
#			sani = 
#				gsub("Gen", "$B_{t\\\\text{-}\\\\mathrm{G}}$",  ## just call it 'G' and specifiy how many G in table caption
#				gsub("Bcurr",paste0("$B_{", currYear, "}$"),
#				gsub("B0","$B_0$",
#				gsub("Bmsy","$B_\\\\mathrm{MSY}$",
#				gsub("0.8Bmsy","USR",
#				gsub("0.4Bmsy","LRP",
#			x))))))
#			return(sani)
#		}
#		#print(sanicol(colnames(Ttab.conf$'0.50')))
#	
#		##--- Time to achieve target with 50% confidence
#		if (exists("Ttab.conf") && !is.null(Ttab.conf$CC$'0.50')) {
#			mess = paste0(name, ": estimated time (years) for projected biomass $B_t$ to exceed reference points and criteria with a probability of 50\\%, for for a range of \\itbf{constant catch} strategies. An estimated time of 0 means that the condition is satisfied and remains so over the ", Ngen*gen1, "-year projection; an estimated time of ", Ngen*gen1, " means that the condition never becomes satisfied over the ", Ngen, "-generation projection. A further condition is that the probability of satisfying the condition must increase for two consecutive years. Columns respectively correspond to the provisional DFO reference points: LRP~= $0.4\\Bmsy$, USR~= $0.8\\Bmsy$; alternative reference points: $\\Bmsy$, $B_{", currYear, "}$, $0.2B_0$, $0.4B_0$; and COSEWIC reference criteria: 0.5$B_{t\\text{-}\\mathrm{G}}$~= $\\leq 50\\%$ decline over ", Ngen, " generations (G), 0.7$B_{t\\text{-}\\mathrm{G}}$~= $\\leq 30\\%$ decline over ", Ngen, "G, $0.5B_0$, $0.7B_0$.")
#			Ttime50.CC.out = capture.output(print(xtable(Ttab.conf$CC$'0.50', caption=eval(mess), label=paste0("tab:",prefix,"Ttime50.CC"), digits=0), caption.placement="top", table.placement="!ht", sanitize.colnames.function=sanicol) )
#			tput(Ttime50.CC.out)
#		}
#		if (exists("Ttab.conf") && !is.null(Ttab.conf$HR$'0.50')) {
#			mess = paste0(name, ": estimated time (years) for projected biomass $B_t$ to exceed reference points and criteria with a probability of 50\\%, for for a range of \\itbf{harvest rate} strategies. See caption in Table~\\ref{tab:", prefix, "Ttime50.CC} for further details")
#			Ttime50.HR.out = capture.output(print(xtable(Ttab.conf$HR$'0.50', caption=eval(mess), label=paste0("tab:",prefix,"Ttime50.HR"), digits=0), caption.placement="top", table.placement="!ht", sanitize.colnames.function=sanicol) )
#			tput(Ttime50.HR.out)
#		}
#		##--- Time to achieve target with 65% confidence
#		if (exists("Ttab.conf") && !is.null(Ttab.conf$CC$'0.65')) {
#			mess = paste0(name, ": estimated time (years) for projected biomass $B_t$ to exceed reference points and criteria with a probability of 65\\%, for for a range of \\itbf{constant catch} strategies. See caption in Table~\\ref{tab:", prefix, "Ttime50.CC} for further details.")
#			Ttime65.CC.out = capture.output(print(xtable(Ttab.conf$CC$'0.65', caption=eval(mess), label=paste0("tab:",prefix,"Ttime65.CC"), digits=0), caption.placement="top", table.placement="!ht", sanitize.colnames.function=sanicol) )
#			tput(Ttime65.CC.out)
#		}
#		if (exists("Ttab.conf") && !is.null(Ttab.conf$HR$'0.65')) {
#			mess = paste0(name, ": estimated time (years) for projected biomass $B_t$ to exceed reference points and criteria with a probability of 65\\%, for for a range of \\itbf{harvest rate} strategies. See caption in Table~\\ref{tab:", prefix, "Ttime50.CC} for further details.")
#			Ttime65.HR.out = capture.output(print(xtable(Ttab.conf$HR$'0.65', caption=eval(mess), label=paste0("tab:",prefix,"Ttime65.HR"), digits=0), caption.placement="top", table.placement="!ht", sanitize.colnames.function=sanicol) )
#			tput(Ttime65.HR.out)
#		}
#		##--- Time to achieve target with 80% confidence
#		if (exists("Ttab.conf") && !is.null(Ttab.conf$CC$'0.80')) {
#			mess = paste0(name, ": estimated time (years) for projected biomass $B_t$ to exceed reference points and criteria with a probability of 80\\%, for for a range of \\itbf{constant catch} strategies. See caption in Table~\\ref{tab:", prefix, "Ttime50.CC} for further details.")
#			Ttime80.CC.out = capture.output(print(xtable(Ttab.conf$CC$'0.80', caption=eval(mess), label=paste0("tab:",prefix,"Ttime80.CC"), digits=0), caption.placement="top", table.placement="!ht", sanitize.colnames.function=sanicol) )
#			tput(Ttime80.CC.out)
#		}
#		if (exists("Ttab.conf") && !is.null(Ttab.conf$HR$'0.80')) {
#			mess = paste0(name, ": estimated time (years) for projected biomass $B_t$ to exceed reference points and criteria with a probability of 80\\%, for for a range of \\itbf{harvest rate} strategies. See caption in Table~\\ref{tab:", prefix, "Ttime50.CC} for further details.")
#			Ttime80.HR.out = capture.output(print(xtable(Ttab.conf$HR$'0.80', caption=eval(mess), label=paste0("tab:",prefix,"Ttime80.HR"), digits=0), caption.placement="top", table.placement="!ht", sanitize.colnames.function=sanicol) )
#			tput(Ttime80.HR.out)
#		}
#		##--- Time to achieve target with 95% confidence
#		if (exists("Ttab.conf") && !is.null(Ttab.conf$CC$'0.95')) {
#			mess = paste0(name, ": estimated time (years) for projected biomass $B_t$ to exceed reference points and criteria with a probability of 95\\%, for for a range of \\itbf{constant catch} strategies. See caption in Table~\\ref{tab:", prefix, "Ttime50.CC} for further details.")
#			Ttime95.CC.out = capture.output(print(xtable(Ttab.conf$CC$'0.95', caption=eval(mess), label=paste0("tab:",prefix,"Ttime95.CC"), digits=0), caption.placement="top", table.placement="!ht", sanitize.colnames.function=sanicol) )
#			tput(Ttime95.CC.out)
#		}
#		if (exists("Ttab.conf") && !is.null(Ttab.conf$HR$'0.95')) {
#			mess = paste0(name, ": estimated time (years) for projected biomass $B_t$ to exceed reference points and criteria with a probability of 95\\%, for for a range of \\itbf{harvest rate} strategies. See caption in Table~\\ref{tab:", prefix, "Ttime50.CC} for further details.")
#			Ttime95.HR.out = capture.output(print(xtable(Ttab.conf$HR$'0.95', caption=eval(mess), label=paste0("tab:",prefix,"Ttime95.HR"), digits=0), caption.placement="top", table.placement="!ht", sanitize.colnames.function=sanicol) )
#			tput(Ttime95.HR.out)
#		}
#	}
#	save(
#	"gmu.LRP.CCs.out", "gmu.USR.CCs.out", "gmu.Bmsy.CCs.out", "gmu.umsy.CCs.out",
#	"gmu.Bcurr.CCs.out", "gmu.ucurr.CCs.out", "gmu.20B0.CCs.out", "gmu.40B0.CCs.out",
#	"gmu.LRP.CCs.out", "gmu.LRP.HRs.out", "gmu.USR.CCs.out", "gmu.USR.HRs.out",
#	"gmu.Bmsy.CCs.out", "gmu.Bmsy.HRs.out", "gmu.umsy.CCs.out", "gmu.umsy.HRs.out",
#	"gmu.Bcurr.CCs.out", "gmu.Bcurr.HRs.out", "gmu.ucurr.CCs.out", "gmu.ucurr.HRs.out",
#	"gmu.20B0.CCs.out", "gmu.20B0.HRs.out", "gmu.40B0.CCs.out", "gmu.40B0.HRs.out",
#	"gmu.LRP.CCl.out", "gmu.LRP.HRl.out", "gmu.USR.CCl.out", "gmu.USR.HRl.out",
#	"gmu.Bmsy.CCl.out", "gmu.Bmsy.HRl.out", "gmu.umsy.CCl.out", "gmu.umsy.HRl.out",
#	"gmu.Bcurr.CCl.out", "gmu.Bcurr.HRl.out", "gmu.ucurr.CCl.out", "gmu.ucurr.HRl.out",
#	"gmu.20B0.CCl.out", "gmu.20B0.HRl.out", "gmu.40B0.CCl.out", "gmu.40B0.HRl.out",
#	"cosewic.50Gen.CCs.out", "cosewic.50Gen.HRs.out", "cosewic.30Gen.CCs.out", "cosewic.30Gen.HRs.out",
#	"cosewic.50B0.CCs.out", "cosewic.50B0.HRs.out", "cosewic.70B0.CCs.out", "cosewic.70B0.HRs.out",
#	"cosewic.50Gen.CCl.out", "cosewic.50Gen.HRl.out", "cosewic.30Gen.CCl.out", "cosewic.30Gen.HRl.out",
#	"cosewic.50B0.CCl.out", "cosewic.50B0.HRl.out", "cosewic.70B0.CCl.out", "cosewic.70B0.HRl.out",
#	"Ttime50.CC.out", "Ttime50.HR.out", "Ttime65.CC.out", "Ttime65.HR.out",
#	"Ttime80.CC.out", "Ttime80.HR.out", "Ttime95.CC.out", "Ttime95.HR.out",
#	file=paste0(prefix,"decision.tables", ifelse(useRlow, paste0(".Rlow(q=",pad0(qRlow*100,2),")"),""), ".rda") )

##	## Rob Tadey asked for exploitation rate and probability of Ut less than Umsy
##	U.tadey  = t(apply(U.mcmc, 2, quantile, quants3))
##	colnames(U.tadey) = paste0("Q.", formatC(quants3*100, digits=0, width=2, format="d", flag="0"))
##	N.tadey  = apply(U.mcmc<Umsy,2,sum)
##	P.tadey  = N.tadey/nrow(U.mcmc)
##	U.tadey  = data.frame(U.tadey, P.MSY=P.tadey)
##	U.tadey2 = splitTab(U.tadey, np=2, row.label="Year")
##	names.tadey = colnames(U.tadey2)
##	#names.tadey = gsub("Q\\.", "Q$_{0.", 
##	#	gsub("P\\.MSY", "P($u_t<u_\\\\mathrm{MSY}$)", names.tadey))
##	zq = grep("^Q",names.tadey)
##	zp = grep("^P",names.tadey)
##	#names.tadey[zq] = paste0(names.tadey[zq],"}$")
##	names.tadey = gsub("^Q\\.(0)?","",
##		gsub("P\\.MSY", "Prob", names.tadey))
##	names.tadey[zq] = paste0(names.tadey[zq],"\\%")
##	names(U.tadey2) = names.tadey
##
##	z = texArray(U.tadey2, zero="0", outnam="tempfile", sigdig=2, use.round=T, struts=T, table.label=paste0("tab:",prefix,"base.u.tadey"), tablewidth=6.5,
##		table.caption = paste0("Quantiles (0.05, 0.5, 0.95) of annual exploitation rate $u_t$ (harvest rate~= catch divided by vulnerable biomass) from ", startYear, " to current model year ", prevYear, " and projected to ", projYear, " assuming a constant catch of ", catpol, "~t. Prob = P($u_t<u_\\mathrm{MSY}$).") )
##	xtab.compo.U.tadey = (z$tabfile)
##	tput(xtab.compo.U.tadey)
#browser();return()

	#out=texArray(U.tadey2,sigdig=2,use.round=T)
	#xtab.compo.U.tadey = xtable(U.tadey2, label=paste0("tab:",prefix,"base.u.tadey"), digits=decdig,
	#	caption=eval(paste0("Annual exploitation rate (expressed as harvest rate: catch over vulnerable biomass) from ", startYear, " to current model year ", prevYear, " and projected to ", projYear, " assuming a constant catch of ", round(recentCatchMean, dig=0), "~t. ")) )
	#tput(xtab.compo.U.tadey)
return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tabSS.decision


## tabSS.senso--------------------------2021-09-17
## Make Sensitivity Tables
## ---------------------------------------------RH
tabSS.senso = function(istock="YMR", prefix="ymr.", senso, sigdig=4)
{
	on.exit(gc(verbose=FALSE))
	unpackList(stock[[istock]][["Controls"]])
	unpackList(senso)

	## Use same preliminary code from 'plotSS.senso'
	## ---------------------------------------------
	modYrs      = startYear:currYear
	modYrsChar  = as.character(modYrs)
	## Diagnostics for select parameters
	P.names     = colnames(senPA) ## paremeter names
	S.mpd       = smpdPA; storage.mode(S.mpd)="double"     ## For some reason, this matrix is stored as a list (maybe due to NA values?)
	S.runs      = as.numeric(sapply(strsplit(rownames(senPA),"\\."),function(x){x[1]})) ## vector of Runs
	S.run.rwt   = sapply(strsplit(rownames(senPA),"\\."),function(x){paste0(x[1],".",x[2])}) ## vector of Runs and Rwts
	S.run.ord   = unique(S.runs)                           ## unique Run numbers (1st is central run)
	S.run.nmc   = table(S.runs)[as.character(S.run.ord)]   ## vector of sensitivity runs ordered
	S.run.num   = rep(1:length(S.run.nmc), S.run.nmc) - 1  ## 1st run is Central Run
	S.num       = unique(S.run.num)
	## Look into subsetting later
	use.run.rwt = unique(S.run.rwt)
	#is.element(P.runs, P.run.ord) ## default use all runs but a subset might be used for management advice
	if (spp.code %in% c("YMR"))
		P.ord=c(1,seq(2,10,2),seq(3,11,2),12,13)
	else
		P.ord = 1:ncol(smpdPA)

	## Create Sensitivity labels
	S.prefix    = paste0("S",pad0(S.num,2)," (R", pad0(S.run.ord,2) ,") ")
	iCR         = grep(exp.run.rwt,unique(B.index))
	S.prefix[1] = paste0("B",iCR, " (R", strsplit(exp.run.rwt,"\\.")[[1]][1], ") ")
	S.labels    = paste0(S.prefix, gsub("\\_"," ", c("Central Run",sen.lab)))
	S.labels    = sub("\\\\pc", "%", S.labels)  ## just in case

	## Calculate quantiles for later (eventually make this function global)
	calcQs = function (dat, ivec, ovec) {
		F.data = as.data.frame(dat)
		F.list = split(F.data, ivec)
		rr.ord = match(ovec, substring(names(F.list),1,2))  ## oder not so important for sensitivities
		F.list = F.list[rr.ord]
		F.qnts = lapply(F.list,function(x){
			z = apply(x,2,function(xx){!any(is.na(xx))})
			out = apply(x[z],2,quantile,quants5)
			return(out)
		}) ## lapply does not sort like split does
		return(F.qnts)
	}
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
	tab.sens.pars = formatCatch(tabPmed,N=sigdig)
	colnames(tab.sens.pars) = gsub(" +","",S.prefix)

	names.pars = rownames(tab.sens.pars)
	names.pars = sub("Male","{2}",sub("Female","{1}",names.pars))
	fixsub = grep("[^M]_",names.pars)
	names.pars[fixsub] = paste0(sub("_","~(",names.pars[fixsub]),")")
	names.pars = sub("varL\\(","\\\\log v(\\\\mathrm{L}",names.pars)
	names.pars[fixsub] = sub("\\(","_{", sub("\\)","}", names.pars[fixsub]))
	names.pars = gsub("mu","\\\\mu",names.pars)
	names.pars = gsub("Delta","\\\\Delta",names.pars)
	names.pars = sub("LN\\(R0)","\\\\log R_{0}",names.pars)
	names.pars = sub("~\\([[:alpha:]]+(\\+)?\\)$", "", names.pars)  ## keep annotation short 
	#names.pars[fixsub] = sub("\\(","(\\\\mathrm{", sub("\\)","})", names.pars[fixsub]))
	rownames(tab.sens.pars) =  paste(rep("$",nrow(tab.sens.pars)),names.pars,rep("$",nrow(tab.sens.pars)),sep="")

	sen.leg = paste0("Sensitivity runs: ", paste0("S", formatC(S.num[-1], width=0, format="d", flag="0"),"~= ", gsub("\\%","\\\\%",gsub("_"," ",sen.lab)), collapse=", "))
	if (spp.code=="REBS") {  ## need to add Other fishery because it has ages but no CPUE index
		nseries = length(iseries)
		iseries[nseries] = paste0(iseries[nseries],"/Trawl fishery")
		iseries = c(iseries, "Other fishery")
	}
	cap.par = paste0(
		name, ": median values of MCMC samples for the primary estimated parameters, ",
		"comparing the central run to ", Nsens, " sensitivity runs (\\Nmcmc{} samples each). C~=Central, R~= Run, S~= Sensitivity. ",
		#"D~=Diagnostic sensitivity ($M$=0.11). ",
		"Numeric subscripts other than those for $R_0$ and $M$ indicate the following gear types $g$: ",
		texThatVec(paste0(1:length(iseries),"~= ",iseries),simplify=F), ". ", sen.leg
	)
	xtab.sens.pars = xtable(tab.sens.pars, align=paste0("c",paste0(rep("r",Nsens+1),collapse="")),
		label   = paste0("tab:",prefix,"sens.pars"), digits = if (exists("formatCatch")) NULL else sigdig,
		caption = cap.par )
	xtab.sens.pars.out = capture.output( print(xtab.sens.pars, caption.placement="top",
		sanitize.rownames.function=function(x){x},
		add.to.row =list(pos = list(-1), command = c("\\\\[-1.0ex]")),
		size = "\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\footnotesize"
	) )
	tput(xtab.sens.pars.out)
#browser();return()

	## MCMC MSY-based quantities
	## -------------------------
	ls2df = function(x) {
		nn = as.character(substitute(x))
		xx = as.data.frame(x)
		colnames(xx) = paste(nn,colnames(xx),sep="_")
		return(xx)
	}
	#D.qnts = calcQs(senTS[,,"BtB0"], ivec=S.run.rwt, ovec=S.run.ord)
	#U.qnts = calcQs(senTS[,,"ut"], ivec=S.run.rwt, ovec=S.run.ord)
	#R.qnts = calcQs(senTS[,,"Rt"], ivec=S.run.rwt, ovec=S.run.ord)
	#RD.qnts = calcQs(senTS[,,"Rtdev"], ivec=S.run.rwt, ovec=S.run.ord)

	Q.mcmc = data.frame (
		B0         = senRP[,"B0"],
		Bcurr      = senRP[,"Bcurr"],
		Bcurr.B0   = senRP[,"Bcurr"] / senRP[,"B0"],
		ucurr      = senRP[,"ucurr"],
		umax       = apply(senTS[,,"ut"],1,max), ## for each mcmc sample across the time series
		MSY        = senRP[,"MSY"],
		Bmsy       = senRP[,"Bmsy"],
		LRP        = senRP[,"LRP"],
		USR        = senRP[,"USR"],
		Bcurr.Bmsy = senRP[,"Bcurr"] / senRP[,"Bmsy"],
		Bmsy.B0    = senRP[,"Bmsy"] / senRP[,"B0"],
		umsy       = senRP[,"umsy"],
		ucurr.umsy = senRP[,"ucurr"] / senRP[,"umsy"]
	)
	Q.mcmc.sens = split(Q.mcmc, S.run.rwt)
	tabQmed = sapply(Q.mcmc.sens, function(Q){
		sapply(Q, median)
	})
	tab.sens.rfpt = formatCatch(tabQmed,N=sigdig-1)  ## use 3 instead of 4
	colnames(tab.sens.rfpt) = gsub(" +","",S.prefix)

	names.rfpt = rownames(tab.sens.rfpt)
	## Use routine from 'make.base.tabs.r':
	names.rfpt =
		gsub("LRP",  paste0("0.4B_{",currYear,"}"),
		gsub("USR",  paste0("0.8B_{",currYear,"}"),
		gsub("umax", "u_\\\\mathrm{max}",
		gsub("msy", "_\\\\mathrm{MSY}",
		gsub("curr",  paste0("_{",currYear,"}"),
		gsub("umsy",  "u_\\\\mathrm{MSY}",
		gsub("ucurr", paste0("u_{",currYear,"}"),
		gsub("0", "_{0}",
		gsub("\\.", "/",
		names.rfpt)))))))))
	rownames(tab.sens.rfpt) =  paste(rep("$",nrow(tab.sens.rfpt)),names.rfpt,rep("$",nrow(tab.sens.rfpt)),sep="")

	cap.rfpt = paste0(
		name, ": medians of MCMC-derived quantities from the central run and ", Nsens,
		" sensitivity runs (\\Nmcmc{} samples each) from their respective MCMC posteriors. Definitions are: ",
		"$B_0$ -- unfished equilibrium spawning biomass (mature females), ",
		#"$V_0$ -- unfished equilibrium vulnerable biomass (males and females), ",
		"$B_{", currYear, "}$ -- spawning biomass at the end of ",currYear, ", ",
		#"$V_{", currYear, "}$ -- vulnerable biomass in the middle of ", currYear, ", ",
		"$u_{", currYear, "}$ -- exploitation rate (ratio of total catch to vulnerable biomass) in the middle of ", currYear, ", ",
		"$u_\\mathrm{max}$ -- maximum exploitation rate (calculated for each sample as the maximum exploitation rate from ",
		startYear , " - ", currYear, "), ",
		"MSY -- maximum sustainable yield at equilibrium, ",
		"$B_\\mathrm{MSY}$ -- equilibrium spawning biomass at MSY, ",
		"$u_\\mathrm{MSY}$ -- equilibrium exploitation rate at MSY, ",
		#"$V_\\mathrm{MSY}$ -- equilibrium vulnerable biomass at MSY. ",
		"All biomass values (and MSY) are in tonnes. ", sen.leg
	)
#browser();return()

	xtab.sens.rfpt = xtable(tab.sens.rfpt, align=paste0("l",paste0(rep("r",Nsens+1),collapse="")),
		label   = paste0("tab:",prefix,"sens.rfpt"), digits = if (exists("formatCatch")) NULL else sigdig,
		caption = cap.rfpt )
	xtab.sens.rfpt.out = capture.output( print(xtab.sens.rfpt, caption.placement="top",
		sanitize.rownames.function=function(x){x},
		add.to.row=list(pos=list(-1,3,11), command=c("\\\\[-1.0ex]", "\\hdashline \\\\[-1.75ex]", "\\hdashline \\\\[-1.75ex]")),
		hline.after =  c(-1,0,5,nrow(xtab.sens.rfpt)),
		size = "\\usefont{\\encodingdefault}{\\familydefault}{\\seriesdefault}{\\shapedefault}\\footnotesize"
	) )
	tput(xtab.sens.rfpt.out)
#browser();return()

	## Sensitivity run likelihoods
	## ---------------------------
	LL.senso =data.frame(
		Sen.Run =c("S00 (R75)", "S01 (R78)", "S02 (R79)", "S03 (R80)", "S04 (R81)", "S05 (R82)", "S06 (R83)", "S07 (R84)", "S08 (R85)", "S09 (R86)", "S10 (R87)", "S11 (R88)", "S12 (R91)", "S13 (R92)", "S14 (R93)"),
		Label   =c("central run", "add 1997 WCHG index", "estimate M", "drop CPUE", "Tweedie CPUE", "sigmaR=0.6", "sigmaR=1.2", "reduce catch 33%", "increase catch 50%", "upweight QC AF", "start Rdevs in 1970", "no ageing error", "steepness h=0.5", "double 2021 catch", "AE from age readers"),
		CPUE    =c(-18.06, -18.13, -17.02, 0, -22.47, -16.83, -18.5, -18.4, -17.4, -18.53, -8.074, -17.87, -18.02, -18.06, -17.87),
		QCS     =c(0.8701, 0.9094, 0.4254, 0.6586, 0.9573, 0.6097, 0.9335, 1.118, 0.5719, 0.6807, -0.8074, 0.9635, 0.927, 0.8724, 0.7353),
		WCVI    =c(7.92, 7.918, 7.937, 7.939, 7.952, 7.749, 8.026, 7.893, 7.943, 7.842, 8.249, 7.819, 7.879, 7.921, 7.952),
		WCHG    =c(19.68, 19.13, 18.93, 19.25, 19.73, 19.43, 19.68, 20.08, 19.2, 19.99, 16.41, 19.9, 19.79, 19.68, 19.44),
		GIG     =c(14.48, 14.5, 14.58, 14.32, 14.52, 14.2, 14.58, 13.67, 15.45, 14.42, 12.13, 14.5, 14.48, 14.48, 14.41),
		Index   =c(24.89, 24.33, 24.86, 42.17, 20.69, 25.16, 24.72, 24.36, 25.77, 24.41, 27.9, 25.31, 25.05, 24.89, 24.67),
		AF      =c(456.4, 456.5, 457.2, 705.9, 715.1, 482.6, 442.4, 451.9, 459.8, 626.8, 540.1, 387.2, 457.1, 458.8, 443.3),
		Recruit =c(41.93, 41.96, 39.27, 51.79, 52.24, 55.48, 39.61, 42.42, 42.87, 51.04, 21.02, 35.09, 42.31, 42.02, 40.69),
		Total   =c(635.1, 634.7, 635.4, 912.2, 900.5, 675, 618.6, 630.5, 640.3, 817.6, 700.7, 560.2, 636.4, 637.6, 621.1)
	)
#browser();return()
	xtab.sruns.ll = xtable(formatCatch(LL.senso, X=1:2), align=paste0("l", paste0(rep("r",dim(LL.senso)[2]),collapse="")),
		label   = paste0("tab:",prefix,"log.likes"), digits=NULL, 
		caption = "Log likelihood (LL) values reported by central and sensitivity runs for survey indices, age composition (AF), recruitment, and total (not all LL components reported here)")
	xtab.sruns.ll.out = capture.output(print(xtab.sruns.ll,  include.rownames=FALSE,
		caption.placement="top", sanitize.rownames.function=function(x){x}, add.to.row =list(pos = list(-1,1), command = c("\\\\[-0.5ex]", "\\hdashline \\\\[-1.75ex]")) ) )
	tput(xtab.sruns.ll.out)

	save("xtab.sens.pars.out", "xtab.sens.rfpt.out", "xtab.sruns.ll.out", file=paste0(prefix,"senso.tabs.rda"))

	collect=c("tab.sens.pars","tab.sens.rfpt")
	for (i in collect) 
		eval(parse(text=paste0("stock[[istock]][[\"Sens\"]][[\"",i,"\"]] = ",i)))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tabSS.senso
