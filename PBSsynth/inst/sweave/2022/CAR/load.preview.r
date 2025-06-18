## Set the following to preview figures and tables before running the main Sweave routine and after running 'set.controls.r'
d.code = "C:/Users/haighr/Files/Projects/R/Develop/PBSsynth/Authors/Rcode/develop/CAR2022"
source(file.path(d.code,"set.controls.r")) 
istock = "CAR"
unpackList(stock[[istock]][["Controls"]]) 


ngear   = 2  ## number of commercial fisheries
nfleet  = 8  ## number of fleets (fisheries + surveys)
ncpue   = 1
nsurv   = 6
cvp1    = 0.178
modYrs  = startYear:currYear
proYrs  = (currYear+1):projYear
assYrs  = c(1999,2005,2007,2009)

run.rwt = paste0(pad0(run.num,2),".",pad0(rwt.num,2))
use.run.rwt = B.index = run.rwt[use.num]

refCC.sentence   = "For reference, the average catch over the last 5 years (2017-2021) was 775~t by Trawl and 13.5~t by Other."
central.run      = as.numeric(strsplit(exp.run.rwt,"\\.")[[1]][1])
central.rwt      = as.numeric(strsplit(exp.run.rwt,"\\.")[[1]][2])
central.mpd.dir  = file.path(base.dir, paste0("Run",central.run), paste0("MPD.",exp.run.rwt))
central.mcmc.dir = file.path(base.dir, paste0("Run",central.run), paste0("MCMC.",exp.run.rwt,".nuts4K"))

load("compo.220715.rda")
unpackList(compo)

P.rc    = .findSquare(ncol(ampdPA))  ## Parameter panel setup
R.rc    = .findSquare(nrow(ampdPA))  ## Component run setup
N.mcmc  = nrow(avgTS)

## Combine TS and PJ's AC
Bt.mcmc     = cbind(avgTS[,,"Bt"], avgPJ[,,"Bt","AC.00"])
BtBmsy.mcmc = cbind(avgTS[,,"BtBmsy"], avgPJ[,,"BtBmsy","AC.00"])
BtB0.mcmc   = cbind(avgTS[,,"BtB0"], avgPJ[,,"BtB0","AC.00"])
ut.mcmc     = cbind(avgTS[,,"ut"], avgPJ[,,"ut","AC.00"])
utumsy.mcmc = cbind(avgTS[,,"utumsy"], avgPJ[,,"utumsy","AC.00"])
Rt.mcmc     = cbind(avgTS[,,"Rt"], avgPJ[,,"Rt","AC.00"])
Rtdev.mcmc  = cbind(avgTS[,,"Rtdev"], avgPJ[,,"Rtdev","AC.00"])

packList(c("N.mcmc", "P.rc", "R.rc", "Bt.mcmc", "BtBmsy.mcmc", "BtB0.mcmc", "ut.mcmc", "utumsy.mcmc", "Rt.mcmc", "Rtdev.mcmc"), target="data.compo.figs")

## Tab composition
RUN.NUM     = texThatVec(run.num, simplify=FALSE)
NrefM       = sum(!is.na(stock[[istock]]$Control$axisU))  ## number or reference models
fornow      = base.lab
if (length(fornow)!=NrefM)
	stop ("Sum Ting Wong: Number of base runs do not match labels...")
NaxU = getActDim(stock[[istock]]$Control$axisU)             ## Number of axes of uncertainty
SaxU = paste0(paste0("\\\\item ",fornow),collapse="\\\\\\\\")  ## Construct latex item list of base runs labels


#load("senso.220624.rda")  ## test sensitivities -- S01(R25), S02(R26), S10(R34), S11(R35)
#load("senso.220627.rda")  ## sensitivities -- S01-S13 (R25-R37)
#load("senso.220630.rda")  ## sensitivities -- S01-S14 (R25-R37, R39)
#load("senso.220715.rda")  ## sensitivities -- S01-S13 (R25-R37)
#load("penso.220720.rda")  ## PDO sensitivities -- (R38 v.6,7,8)

load("senso.220912.rda")  ## After RPR, S14 (Francis reweight) added
#unpackList(senso)

