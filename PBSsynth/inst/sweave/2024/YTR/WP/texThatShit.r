## Function to deal with lists of g that vary by stock
texThatShit = function(g, gpad=c("trawl survey series","(bottom trawl)"), T, U, Ngear=1)
{
	gN    = length(g)
	gtab  = as.character()
	gnam  = list()
	onebs = "\\\\"
	twobs = "\\\\\\\\"
	for (i in 1:gN) {
		ii = names(g)[i]
		ivec = g[[i]]
		gvec = gsub("QCS", "Queen Charlotte Sound (QCS)",
			gsub("WCVI", "West Coast Vancouver Island (WCVI)",
			gsub("WCHG", "West Coast Haida Gwaii (WCHG)",
			gsub("GIG",  "Goose Island Gully (GIG)",
			gsub("HS",   "Hecate Strait (HS)",
			gsub("CPUE", "catch per unit effort (CPUE)",
			gsub("Trawl Fishery\\+", "All fisheries but dominated by trawl",
			ivec)))))))
		gnam[[ii]] = gvec
		if (gN>1)
			gtab = c(gtab, paste0(" & ", names(g)[i], ":", twobs))
		iN   = length(ivec)
		if (is.null(gpad))
			ipad = rep("", iN)
		else 
			ipad = paste0(c(rep(gpad[1], Ngear), rep(gpad[2],iN-Ngear)))
		gtab = c(gtab, paste0(" & ~~", 1:iN, " -- ", ivec, ipad, twobs))
#browser();return()
	}
	texYrs = function(z) {
		zL  = as.character(substitute(z))
		zN  = length(z)
		ztab = as.character()
		for (i in 1:zN) {
			if (gN>1)
				ztab = c(ztab, paste0(" & ", names(z)[i], ":", twobs))
			ivec = z[[i]]           ## not a vector but a list
			iN   = length(ivec)
			for (j in 1:iN) {
				jvec = ivec[[j]]     ## years for g abundance/age data
				if (all(is.na(jvec))) next
				jN = length(jvec)
				if (jN>1) {
					jdiff = diff(jvec)
					if (all(jdiff==1) && length(jvec)>2)
						jazz = paste0(jvec[1],", ..., ",rev(jvec)[1])
					else {
						sYrs = jvec[c(2,jdiff)>1]
						eYrs = jvec[c(jdiff,2)>1]
						cYrs = as.character(sYrs)
						zYrs = sYrs!=eYrs
						if (any(zYrs))
							cYrs[zYrs] = paste(sYrs[zYrs],eYrs[zYrs],sep=":")
#if(any(sYrs==eYrs)) {browser();return()}
						jazz = paste0(cYrs,collapse=", ")
					}
				} else {
					jazz = paste0(jvec,collapse=", ")
				}
#browser();return()
				ztab = c(ztab, paste0(" & ~~${", onebs, "bf ", zL ,"}_{", j, "}$ = ",onebs,"{", jazz, onebs, "}", twobs))
			}
		}
		return(ztab)
	}
	ttab = texYrs(T)
	utab = texYrs(U)
	return(list(gnam=gnam,gtab=gtab,ttab=ttab,utab=utab))
}

#out=texThatShit(g=glist, T=Tlist, U=Ulist, Ngear=2, gpad=c(" trawl survey series"," index"))
#texShit    = texThatShit(g=glist, T=Tlist, U=Ulist, Ngear=Ngear, gpad=c(" (commercial data)", " trawl survey series"))


