## copySSfiles -------------------------2026-02-18
##  Copy SS3 input files to folders suitable for archiving on GitHub
## ---------------------------------------------RH
copySSfiles <- function (strSpp="405", assyr=2025)
{
	rwts  = c("00","01")  ## will have to be more complicated if different number of rweights by run
	if (strSpp=="405" && assyr==2025) {
		dir.from = "C:/Users/haighr/Files/GFish/PSARC25/SGR/Data/SS3/SGR2025"
		dir.to   = "C:/Users/haighr/Files/Projects/R/Develop/PBSsynth/Authors/input/2025/SGR"
		areas = list()
		areas[["BC_coast"]] = list()
		areas[["BC_coast"]][["nbase"]] = 1
		areas[["BC_coast"]][["runs"]]  = c(29,33,23,31,seq(35,47,2),51)
		areas[["BC_coast"]][["vers"]]  = c(2,1,3,2,rep(1,8))
		areas[["BC_3area"]] = list()
		areas[["BC_3area"]][["nbase"]] = 1
		areas[["BC_3area"]][["runs"]]  = c(28,32,21,30,seq(34,46,2),50,48,49)
		areas[["BC_3area"]][["vers"]]  = c(2,1,3,2,rep(1,10))
	}
	## Generic code
	## ------------
	for (a in 1:length(areas)) {
		aa = names(areas)[a]
		.flush.cat("Copying files for model ", aa,"\n")
		adir.to = file.path(dir.to, aa)
		if (!dir.exists(adir.to))
			dir.create(adir.to)
		nbase = areas[[aa]][["nbase"]]
		runs = areas[[aa]][["runs"]]
		vers = areas[[aa]][["vers"]]
		for (i in 1:length(runs)) {
			ii = pad0(runs[i],2)
			vv = vers[i]
			ipref = if (i %in% 1:nbase) paste0("B",i) else paste0("S", pad0(i-nbase,2))
			.flush.cat("   ", ipref,"\n")
			idir  = paste0(ipref, ".R", ii, "v", vv)
			idir.to = file.path(adir.to, idir)
			if (!dir.exists(idir.to))
				dir.create(idir.to)
			for (j in 1:length(rwts)) {
				jj = rwts[j]
				jdir = paste0(ii, ".", jj)
				jdir.to = file.path(idir.to, jdir)
				if (!dir.exists(jdir.to))
					dir.create(jdir.to)
				jdir.from = paste0(dir.from, "/Run", ii, "/MPD.", jdir, ".v", vv)
				if (!dir.exists(jdir.from))
					stop ("Source directory specified incorrectly")
				ssfiles = paste0(c("starter", "forecast", paste0(c("data.","control."), paste0(c(ii,jj),collapse="."))),".ss")
				file.copy(from=file.path(jdir.from,ssfiles), to=jdir.to, overwrite=TRUE, copy.mode=FALSE, copy.date=TRUE)
			}
		}
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~copySSfiles
copySSfiles()

