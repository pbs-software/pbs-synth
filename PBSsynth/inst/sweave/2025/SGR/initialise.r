## Set the following once for each new assessment (RH 230919)
if (!exists("run")) stop("Must specify a run before run-specific settings can be assigned")

strSpp="405"; assyr=2025; species.code="SGR"; species.name="Silvergray Rockfish"  ## also put this line in 'load.preview.r'
#strSpp="467"; assyr=2025; species.code="LIN"; species.name="Outside Lingcod" 
#strSpp="418"; assyr=2024; species.code="YTR"; species.name="Yellowtail Rockfish"

if (strSpp %in% c("405","SGR") && assyr==2025)
{
	## All fleets used in the stock assessment (including sensitivites)
	fleets.all = c("BC Trawl Fishery", "QCS Synoptic", "WCHG Synoptic", "WCVI Synoptic", "HS Synoptic", "GIG Historical", "NMFS Triennial", "HBLL North", "HBLL South")

	## Default values -----------------------------
	fleets.use  = c(1:7)
	fleets.lab  = fleets.all[fleets.use]
	nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
	fleets.idx  = c(2:7)              ## fleets (idx, sel, af, gear) relative to subset
	fleets.sel  = c(1:7)
	fleets.af   = c(1:4)              ## link HS and GIG to QCS; link NMFS to WCVI
	fleets.gear = c(1)                ## commercial gear types
	gear.names  = fleets.lab[fleets.gear]
	ngear       =  length(gear.names) ## number of commercial fisheries
	area.names  = c("BC")
	narea       = length(area.names)  ## number of areas
	col.idx     = rep("purple",6)
	bg.idx      = rep("thistle",6)
	maxage.sel  = 30
	assYrs      = c(1999, 2001, 2014) ## technically, modelled current year (not assessment year)
	## Redefine functions that have disappeared into namespaces
	#.findSquare = PBSmodelling:::.findSquare

	## Started by trying to estimate HS selectivity
	if (run==1) {
		fleets.af   = c(1:5)  ## estimate HS selectivity
	}
	## Add in HBLL surveys
	if (run==3) {
		fleets.use  = c(1:9)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(2:9)              ## fleets (idx, sel, af, gear) relative to subset
		fleets.sel  = c(1:9)
		fleets.af   = c(1:4,8,9)          ## link HS and GIG to QCS; link NMFS to WCVI
		col.idx     = rep("purple",8)
		bg.idx      = rep("thistle",8)
	}
	## Use coastwide CPUE for fleet 1
	if (run %in% c(6,7)) {
		fleets.idx  = c(1:7)              ## fleets (idx, sel, af, gear) relative to subset
		col.idx     = rep("purple",7)
		bg.idx      = rep("thistle",7)
	}
	if (run %in% c(8:100)) {
		fleets.all = c("5ABC Trawl Fishery", "5DE Trawl Fishery", "3CD Trawl Fishery", "QCS Synoptic", "WCHG Synoptic", "WCVI Synoptic", "HS Synoptic", "GIG Historical", "NMFS Triennial")
		fleets.use  = c(1:9)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(1:9)              ## fleets (idx, sel, af, gear) relative to subset
		fleets.sel  = c(1:9)
		fleets.af   = c(1:6)              ## link HS and GIG to QCS; link NMFS to WCVI
		fleets.gear = c(1,2,3)            ## commercial gear types
		gear.names  = fleets.lab[fleets.gear]
		ngear       = length(gear.names) ## number of commercial fisheries
		col.idx     = c("blue","red","green4","blue","red","green4","blue","blue","green4")
		bg.idx      = c("cyan","pink","green","cyan","pink","green","cyan","cyan","green")
	}
	if (run %in% c(8)) { ## added three fishery CPUE series but kept region coastwide (one area)
		col.idx     = rep("purple",9)
		bg.idx      = rep("thistle",9)
	}
	if (run %in% c(9:100)) {
		area.names  = c("5ABC", "5DE", "3CD")
		narea       = length(area.names)  ## number of areas
	}
	if (run %in% c(12:14)) { ## remove CPUE series for 5ABC, 5DE, and 3CD;
		fleets.idx  = c(4:9)
		col.idx     = c("blue","red","green4","blue","blue","green4")
		bg.idx      = c("cyan","pink","green","cyan","cyan","green")
	}
	if (run %in% c(15:16,18)) {
		col.idx     = c("blue","red","green4","blue","red","green4","red","blue","green4")
		bg.idx      = c("cyan","pink","green","cyan","pink","green","pink","cyan","green")
	}

	## Sensitivity runs ---------------------------
	## Something
	if (run==101) {
		fleets.use  = c(1:7,9:10)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(2:9)   ## fleets (idx, sel, af, gear) relative to subset
		fleets.sel  = c(1:9)
	}
	## dome-shaped selectivity
	if (run==102) {
		maxage.sel = 35
	}
	## split trawl fleet into BT and MW trawl fleets
	if (run==103) {
		fleets.use  = c(11,12,2:7)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(3:8)     ## fleets (idx, sel, af, gear) relative to subset
		fleets.sel  = c(1:8)
		fleets.af   = c(1,2,3,4,8) ## fix HS selectivity
		fleets.gear = c(1:2)       ## commercial gear types
		gear.names  = fleets.lab[fleets.gear]
		ngear       =  length(gear.names)  ## number of commercial fisheries
	}
} ## end SGR in 2025

if (strSpp %in% c("467","LIN") && assyr==2025)
{
	## All fleets used in the stock assessment (including sensitivites)
	fleets.all = c("Bottom Trawl", "Hook and Line", "Recreational", "Synoptic Survey", "HBLL Survey", "IPHC Survey", "NMFS Triennial", "BT Discard")

	## Default values -----------------------------
	fleets.use  = c(1:8)
	fleets.lab  = fleets.all[fleets.use]
	nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
	fleets.idx  = c(4:7)              ## fleets (idx, sel, af, gear) relative to subset
	fleets.sel  = c(1:8)
	fleets.af   = c(1,2,4,5,6,8)              ## 
	fleets.gear = c(1:3,8)                ## commercial gear types
	gear.names  = fleets.lab[fleets.gear]
	ngear       =  length(gear.names) ## number of commercial fisheries
	area.names  = c("BC")
	narea       = length(area.names)  ## number of areas
	maxage.sel  = 20
	assYrs      = c(2012) ## technically, modelled current year (not assessment year)
} ## end LIN in 2025

if (strSpp %in% c("418","YTR") && assyr==2024)
{
	## All fleets used in the stock assessment (including sensitivites)
	fleets.all = c("BC Trawl Fishery", "QCS Synoptic", "WCVI Synoptic", "WCHG Synoptic", "HS Synoptic", "GIG Historical", "NMFS Triennial", "WCVI Historical", "HBLL North", "HBLL South", "BC BT Fishery", "BC MW Fishery")

	## Default values -----------------------------
	fleets.use  = c(1:7)
	fleets.lab  = fleets.all[fleets.use]
	nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
	fleets.idx  = c(2:7)     ## fleets (idx, sel, af, gear) relative to subset
	fleets.sel  = c(1:7)
	fleets.af   = c(1,2,3,7) ## fix HS selectivity
	fleets.gear = c(1)       ## commercial gear types
	gear.names  = fleets.lab[fleets.gear]
	ngear       =  length(gear.names)  ## number of commercial fisheries
	area.names  = c("BC")
	narea       = length(area.names)  ## number of areas
	maxage.sel  = 25
	assYrs      = c(1996,1997,2015)  ## Modelled current year
	## Redefine functions that have disappeared into namespaces
	.findSquare = PBSmodelling:::.findSquare

	## Started by trying to estimate HS selectivity
	if (run==1) {
		fleets.af   = c(1,2,3,5,7)  ## estimate HS selectivity
	}
	## Sensitivity runs
	## Add HBLL abundance indices
	if (run==4) {
		fleets.use  = c(1:7,9:10)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(2:9)   ## fleets (idx, sel, af, gear) relative to subset
		fleets.sel  = c(1:9)
	}
	## dome-shaped selectivity
	if (run==11) {
		maxage.sel = 35
	}
	## split trawl fleet into BT and MW trawl fleets
	if (run==17) {
		fleets.use  = c(11,12,2:7)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(3:8)     ## fleets (idx, sel, af, gear) relative to subset
		fleets.sel  = c(1:8)
		fleets.af   = c(1,2,3,4,8) ## fix HS selectivity
		fleets.gear = c(1:2)       ## commercial gear types
		gear.names  = fleets.lab[fleets.gear]
		ngear       =  length(gear.names)  ## number of commercial fisheries
	}
} ## end YTR in 2024

if (strSpp %in% c("607","PEL") && assyr==2023)
{
	## Only one fleet specified in Mackenzies's setup
	fleets.all = c("Trawl Fishery", "QCS Synoptic", "HS Synoptic", "WCVI Synoptic")

	## Default values -----------------------------
	fleets.use  = c(1:4)
	fleets.lab  = fleets.all[fleets.use]
	nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
	fleets.idx  = c(2:4)
	fleets.sel  = c(1:4)
	fleets.af   = c(1:4)
	fleets.gear = c(1)
	gear.names  = fleets.lab[fleets.gear]
	ngear       =  length(gear.names)  ## number of commercial fisheries
	area.names  = c("BC")
	narea       = length(area.names)  ## number of areas
	maxage.sel  = 25
} ## end PEL in 2023
if (strSpp %in% c("396","POP") && assyr==2023)
{
	## All fleets used in the stock assessment (including sensitivites)
	fleets.all = c("5ABC Trawl Fishery", "3CD Trawl Fishery", "5DE Trawl Fishery", "QCS Synoptic", "WCVI Synoptic", "WCHG Synoptic", "GIG Historical", "NMFS Triennial", "WCVI Historical", "3CD Midwater Fishery", "5ABC Midwater Fishery", "HS Synoptic")

	## Default values -----------------------------
	fleets.use  = c(1:9)
	fleets.lab  = fleets.all[fleets.use]
	nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
	fleets.idx  = c(4:9)
	fleets.sel  = c(1:9)
	fleets.af   = c(1:8)
	fleets.gear = c(1:3)
	gear.names  = fleets.lab[fleets.gear]
	ngear       =  length(gear.names)  ## number of commercial fisheries
	area.names  = c("5ABC","3CD","5DE")
	narea       = length(area.names)  ## number of areas
	maxage.sel  = 25
	assYrs      = c(2001,2010,2012,2017)


	## 5ABC Area ----------------------------------
	if (run %in% c(3,18,24)) {
		fleets.use  = c(1,4,7)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(2:3)  ## relative to subset
		fleets.sel  = c(1:3)
		fleets.af   = c(1:3)
		fleets.gear = c(1)
		gear.names  = fleets.lab[fleets.gear]
		ngear       =  length(gear.names)  ## number of commercial fisheries
		area.names  = c("5ABC")
		narea       = length(area.names)  ## number of areas
		assYrs      = c(2001,2010,2017)
	}
	## 3CD Area -----------------------------------
	if (run %in% c(4,19,25)) {
		fleets.use  = c(2,5,8,9)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(2:4)  ## relative to subset
		fleets.sel  = c(1:4)
		fleets.af   = c(1:3)
		fleets.gear = c(1)
		gear.names  = fleets.lab[fleets.gear]
		ngear       =  length(gear.names)  ## number of commercial fisheries
		area.names  = c("3CD")
		narea       = length(area.names)  ## number of areas
		assYrs      = c(2012)
	}
	## 5DE Area -----------------------------------
	if (run %in% c(5,20,26)) {
		fleets.use  = c(3,6)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(2)  ## relative to subset
		fleets.sel  = c(1:2)
		fleets.af   = c(1:2)
		fleets.gear = c(1)
		gear.names  = fleets.lab[fleets.gear]
		ngear       =  length(gear.names)  ## number of commercial fisheries
		area.names  = c("5DE")
		narea       = length(area.names)  ## number of areas
		assYrs      = c(2012)
	}
	## Add 3CD and 5ABC midwater fleets
	if (run %in% c(22)) {
		fleets.use  = c(1:11)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(4:9)  ## relative to subset
		fleets.sel  = c(1:11)
		fleets.af   = c(1:8,10)
		fleets.gear = c(1:3,10:11)
		gear.names  = fleets.lab[fleets.gear]
		ngear       =  length(gear.names)  ## number of commercial fisheries
		area.names  = c("5ABC","3CD","5DE")
		narea       = length(area.names)  ## number of areas
	}
	## Add HS synoptic survey
	if (run %in% c(36)) {
		fleets.use  = c(1:9,12)
		fleets.lab  = fleets.all[fleets.use]
		nfleet      = length(fleets.use)  ## number of fleets (fisheries + surveys)
		fleets.idx  = c(4:10)  ## relative to subset
		fleets.sel  = c(1:10)
		fleets.af   = c(1:8)
		fleets.gear = c(1:3)
		gear.names  = fleets.lab[fleets.gear]
		ngear       =  length(gear.names)  ## number of commercial fisheries
		area.names  = c("5ABC","3CD","5DE")
		narea       = length(area.names)  ## number of areas
	}
} ## end POP in 2023
if (strSpp %in% c("437","CAR") && assyr==2022)
{
	## Default values
	fleets.lab = c("Trawl Fishery","Other Fishery","QCS Synoptic","WCVI Synoptic","NMFS Triennial","HS Synoptic","WCHG Synoptic","GIG Historical")
	fleets.idx = fleets.sel = c(1,3:8)
	fleets.af  = c(1,3:5)
	## R34 (S10) -- Add HS and WCHG AF data
	if (run %in% c(34)) {
		if (type=="MPD" && rwt==0)
			fleets.af  = c(1,3:7)
		else
			fleets.af  = c(1,3:7)
	}
	## R35 (S11) -- Use HBLL North and South survey series
	if (run %in% c(35,47)) {
		fleets.lab = c("Trawl Fishery","Other Fishery","QCS Synoptic","WCVI Synoptic","NMFS Triennial","HS Synoptic","WCHG Synoptic","GIG Historical","HBLL North","HBLL South")
		fleets.idx = fleets.sel = c(1,3:10)
		fleets.af  = c(1,3:5,9,10)
	}
	 ## R37 (S13) -- Remove commercial CPUE times series
	if (run %in% c(37)) {
		fleets.idx = c(3:8)
	}
	## R42 (PDO) -- Use PDO index as an abundance index
	if (run %in% c(42)) {
		fleets.lab = c("Trawl Fishery","Other Fishery","QCS Synoptic","WCVI Synoptic","NMFS Triennial","HS Synoptic","WCHG Synoptic","GIG Historical","PDO Winter")
		fleets.idx = fleets.sel = c(1,3:9)
		fleets.af  = c(1,3:5)
	}
	## R46 (PDO) -- Remove NMFS Triennial and GIG Historical
	if (run %in% c(46)) {
		fleets.lab = c("Trawl Fishery","Other Fishery","QCS Synoptic","WCVI Synoptic","HS Synoptic","WCHG Synoptic")
		fleets.idx = fleets.sel = c(1,3:6)
		fleets.af  = c(1,3:4)
	}
	maxage = 40
}
if (strSpp %in% c("440","YMR") && assyr==2011)
{
	fleets.af  = c(1:3); fleets.idx = 2:6
	fleets.lab = c("Trawl Fishery", "GIG Historical", "QCS Synoptic", "QCS Shrimp", "WCHG Synoptic", "WCVI Synoptic")
}
if (strSpp %in% c("440","YMR") && assyr==2021)
{
	fleets.af  = c(1:5); fleets.idx = ifelse(is.element(run,c(50:53,55:79,81:100))|run=="75a",1,2):5
	fleets.lab = c("Trawl+ Fishery","QCS Synoptic","WCVI Synoptic","WCHG Synoptic","GIG Historical")
}
if (strSpp %in% c("REBS","REBSN"))
{
	fleets.af  = c(1,3); fleets.idx = 1:3
	fleets.lab = c("Trawl Fishery","Other Fishery","WCHG Synoptic")
}
#else stop("Species specified has no fleet indices specified")

