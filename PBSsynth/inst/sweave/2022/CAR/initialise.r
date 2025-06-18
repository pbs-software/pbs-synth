## Set the following once for each new assessment
strSpp="437"; assyr=2022

## Specify fleets for now (mainly for runs when Sweave is not run) -- get a better solution later
if (strSpp %in% c("437","CAR") && assyr==2022) {
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
if (strSpp %in% c("440","YMR") && assyr==2011){
	fleets.af  = c(1:3); fleets.idx = 2:6
	fleets.lab = c("Trawl Fishery", "GIG Historical", "QCS Synoptic", "QCS Shrimp", "WCHG Synoptic", "WCVI Synoptic")
}
if (strSpp %in% c("440","YMR") && assyr==2021){
	fleets.af  = c(1:5); fleets.idx = ifelse(is.element(run,c(50:53,55:79,81:100))|run=="75a",1,2):5
	fleets.lab = c("Trawl+ Fishery","QCS Synoptic","WCVI Synoptic","WCHG Synoptic","GIG Historical")
}
if (strSpp %in% c("REBS","REBSN")){
	fleets.af  = c(1,3); fleets.idx = 1:3
	fleets.lab = c("Trawl Fishery","Other Fishery","WCHG Synoptic")
}
# else stop("Species specified has no fleet indices specified")



