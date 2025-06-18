## Set the following to preview figures and tables before running the main Sweave routine and after running 'set.controls.r'
## Mirrors code in POP.Rnw so make sure changes in either are updated in the other.

strSpp="418"; assyr=2024; species.code="YTR"; species.name="Yellowtail Rockfish"  ## also put this line in 'initialise.r'

so("convPN.r","synth")
so("set.controls.r","synth");
istock  = "YTR"
unpackList(stock[[istock]][["Controls"]])
appFdir = "C:/Users/haighr/Files/GFish/PSARC24/YTR/Docs/RD/AppF_Results/"  ## where compo, senso, etc. reside

## Cover inconsistent notation
currYr = currYear

ngear   = 1  ## number of commercial fisheries
nfleet  = 7  ## number of fleets (fisheries + surveys)
ncpue   = 0
nsurv   = 6
cvp1    = 0
modYrs  = startYear:currYear
proYrs  = (currYear+1):projYear
assYrs  = c(1996,1997,2015)  ## Modelled current year

area.names  = c("BC")
narea       = length(area.names)  ## number of areas

run.rwt = paste0(pad0(run.num,2),".",pad0(rwt.num,2))
run.rwt.ver = paste0(run.rwt,".v",ver.num)
use.run.rwt = B.index = run.rwt[use.num]
#Nmodel = length(run.rwt.num)

recentCatchMean  = round(sapply(recentCatch,mean))
coastal = is.element (names(recentCatchMean), "CST")
refCC.sentence   = paste0("For reference, the average catch over the last 5 years (", texThatVec(names(recentCatch[[1]])), ") was ",
	paste0(unlist(formatCatch(recentCatchMean[!coastal])),"~t in ", names(recentCatchMean[!coastal]),  collapse=", "), ".")
if (any(coastal))
	refCC.sentence = sub("\\.$", paste0(", and ", unlist(formatCatch(recentCatchMean[coastal])),"~t along the BC coast."), refCC.sentence)
central.run      = as.numeric(strsplit(exp.run.rwt,"\\.")[[1]][1])
central.rwt      = as.numeric(strsplit(exp.run.rwt,"\\.")[[1]][2])
central.mpd.dir  = file.path(base.dir, paste0("Run",pad0(central.run,2)), paste0("MPD.",sub("[a-z]$","",exp.run.rwt)))
central.mcmc.dir = file.path(base.dir, paste0("Run",pad0(central.run,2)), paste0("MCMC.",exp.run.rwt)) #,".nuts4K"))

if (!exists("compo",envir=.GlobalEnv)) {
	#load(paste0(appFdir,"compo.240704.rda"), envir=.GlobalEnv)  ## pre RPR
	load(paste0(appFdir,"compo.241015.rda"), envir=.GlobalEnv)  ## post RPR
}
unpackList(compo)

P.rc    = PBSmodelling:::.findSquare(ncol(ampdPA))  ## Parameter panel setup
R.rc    = PBSmodelling:::.findSquare(nrow(ampdPA))  ## Component run setup
N.mcmc  = nrow(avgTS)

## Combine TS and PJ's AC (average catch -- default projections)
Bt.mcmc     = cbind(avgTS[,,"Bt"],     avgPJ[,,"Bt","AC.00"])
BtBmsy.mcmc = cbind(avgTS[,,"BtBmsy"], avgPJ[,,"BtBmsy","AC.00"])
BtB0.mcmc   = cbind(avgTS[,,"BtB0"],   avgPJ[,,"BtB0","AC.00"])
ut.mcmc     = cbind(avgTS[,,"ut"],     avgPJ[,,"ut","AC.00"])
utumsy.mcmc = cbind(avgTS[,,"utumsy"], avgPJ[,,"utumsy","AC.00"])
Rt.mcmc     = cbind(avgTS[,,"Rt"],     avgPJ[,,"Rt","AC.00"])
Rtdev.mcmc  = cbind(avgTS[,,"Rtdev"],  avgPJ[,,"Rtdev","AC.00"])

if (exists("xavgTS") && exists("xavgPJ")) {  ## collect area values if they exist (multi-area model)
	areaTS = areaPJ = list()
	## Add the whole area to the area list to simplify plotting later on
	for (i in c("Bt","BtBmsy","BtB0","ut","utumsy","Rt","Rtdev")) {
		areaTS[[area.name]][[paste0(i,".mcmc")]] = get(paste0(i,".mcmc"))
		areaPJ[[area.name]][[i]] = avgPJ[,,i,grep("^AC",dimnames(avgPJ)$proj,invert=T)]  ## grab decision-table catch policy scenarios
	}
	areas  = dimnames(xavgTS)$area
	for (a in areas) {
		areaTS[[a]] = list()
		areaTS[[a]][["Bt.mcmc"]]     = cbind(xavgTS[,,"B",a],      xavgPJ[,,"B",a,"AC.00"])
		areaTS[[a]][["BtBmsy.mcmc"]] = cbind(xavgTS[,,"BtBmsy",a], xavgPJ[,,"BtBmsy",a,"AC.00"])
		areaTS[[a]][["BtB0.mcmc"]]   = cbind(xavgTS[,,"D",a],      xavgPJ[,,"D",a,"AC.00"])
		areaTS[[a]][["ut.mcmc"]]     = cbind(xavgTS[,,"u",a],      xavgPJ[,,"u",a,"AC.00"])
		areaTS[[a]][["utumsy.mcmc"]] = cbind(xavgTS[,,"utumsy",a], xavgPJ[,,"utumsy",a,"AC.00"])
		areaTS[[a]][["Rt.mcmc"]]     = cbind(xavgTS[,,"R",a],      xavgPJ[,,"R",a,"AC.00"])
		#areaTS[[a]][["Rtdev.mcmc"]]  = cbind(xavgTS[,,"Rdev",a],   xavgPJ[,,"Rdev",a,"AC.00"])  ## appears to be no area-specific recdev in Report files
		areaTS[[a]][["Rtdev.mcmc"]]  = NA
		for (i in c("Bt","BtBmsy","BtB0","ut","utumsy","Rt")) {
			ii = switch(i, 'Bt'="B", 'BtBmsy'="BtBmsy", 'BtB0'="D", 'ut'="u", 'utumsy'="utumsy", 'Rt'="R")
			areaPJ[[a]][[i]] = xavgPJ[,,ii,a,grep("^AC",dimnames(avgPJ)$proj,invert=T)]  ## grab decision-table catch policy scenarios
		}
		areaPJ[[a]][["Rtdev"]] = NA
	}
}
packList(c("N.mcmc", "P.rc", "R.rc", "Bt.mcmc", "BtBmsy.mcmc", "BtB0.mcmc", "ut.mcmc", "utumsy.mcmc", "Rt.mcmc", "Rtdev.mcmc", "areaTS", "areaPJ"), target="data.compo.figs")

## Tab composition
RUN.NUM     = texThatVec(run.num, simplify=FALSE)
NrefM       = sum(!is.na(stock[[istock]]$Control$axisU))  ## number or reference models
fornow      = base.lab
if (length(fornow)!=NrefM)
	stop ("Sum Ting Wong: Number of base runs do not match labels...")
NaxU = getActDim(stock[[istock]]$Control$axisU)             ## Number of axes of uncertainty
SaxU = paste0(paste0("\\\\item ",fornow),collapse="\\\\\\\\")  ## Construct latex item list of base runs labels

if (!exists("senso",envir=.GlobalEnv)) {
	#load(paste0(appFdir,"senso.240705.rda"), envir=.GlobalEnv)  ## pre RPR
	load(paste0(appFdir,"senso.241015.rda"), envir=.GlobalEnv)  ## post RPR
}
#if (!exists("penso",envir=.GlobalEnv))
#	load(paste0(appFdir,"penso.231010.rda"), envir=.GlobalEnv)
#if (!exists("renso",envir=.GlobalEnv))
#	load(paste0(appFdir,"renso.231012.rda"), envir=.GlobalEnv)
#unpackList(senso)

