##==========================================================
## PBS Stock Synthesis plotting functions:
## ---------------------------------------
## make.multifig.........Make a multi-panel figure (mod.r4ss)
## plotSS.comps..........Plot age proportions and model fit (mod.r4ss)
## plotSS.dmcmc..........Plot diagnostics (traces, split chains, ACFs) for MCMCs
## plotSS.index..........Plot SS model fit to abundance index series
## plotSS.pairs..........Pairs|density plot comparison among parameters
## plotSS.pars...........Plot parameter fits and priors
## plotSS.pmcmc..........Plot parameter quantile boxplots for MCMCs
## plotSS.profile........Plot SS likelihood profiles
## plotSS.rdevs..........Plot reruitment deviations
## plotSS.rmcmc..........Plot routine output for MCMCs
## plotSS.selex..........Plot selectivity curves with maturity ogive (mod.r4ss)
## plotSS.stdres.........Plot standardised residuals (mod.PBSawatea)
## plotSS.stock.recruit..Plot stock-recruitment function (based on MPDs)
## plotSS.ts.............Plot time series from MPD output from SS (mod.r4ss)
## plt.ageResids.........Plot age residuals by age class (mod.PBSawatea)
## plt.yearResids........Plot age residuals by year (mod.PBSawatea)
## plt.cohortResids......Plot age residuals by cohort (mod.PBSawatea)
## plt.selectivity.......Transferred selectivity code from PBSscape.r (mod.PBSawatea)
##==========================================================


## make.multifig------------------------2020-08-25
##  Make a multi-panel figure
##  Source: R package 'r4ss' v.1.39.1
## ----------------------------------------r4ss|RH
make.multifig = function (ptsx, ptsy, yr, linesx=0, linesy=0, ptsSD=0, 
   sampsize=0, effN=0, showsampsize=TRUE, showeffN=TRUE, 
   sampsize_label="N=", effN_label="effN=", sampsizeround=1, 
   maxrows=6, maxcols=6, rows=1, cols=1, fixdims=TRUE, 
   main="", cex.main=1, xlab="", ylab="", size=1, 
   cexZ1=1.5, bublegend=TRUE, maxsize=NULL, do.sqrt=TRUE, 
   minnbubble=8, allopen=TRUE, xbuffer=c(0.1,0.1), 
   ybuffer=c(0,0.15), yupper=NULL, ymin0=TRUE, xlas=0, ylas=NULL, 
   axis1=NULL, axis2=NULL, axis1labs=NULL, linepos=1, 
   type="o", polygons=TRUE, bars=FALSE, barwidth="default", 
   ptscex=1, ptscol=1, ptscol2=1, 
   colvec=c(rgb(1,0,0,0.7), rgb(0,0,1,0.7), rgb(0.1,0.1,0.1,0.7)),
   linescol=c(rgb(0,0.8,0,0.7), rgb(1,0,0,0.7), rgb(0,0,1,0.7)),
   lty=1, lwd=2, pch=1, nlegends=3, 
   legtext=list("yr", "sampsize", "effN"), legx="default", 
   legy="default", legadjx="default", legadjy="default", 
   legsize=c(1.2, 1), legfont=c(2, 1), venusmars=TRUE, 
   sampsizeline=FALSE, effNline=FALSE, sampsizemean=NULL, 
   effNmean=NULL, ipage=0, scalebins=FALSE, sexvec=NULL, 
   multifig_colpolygon=c("grey60","grey80","grey70"), 
   multifig_oma=NULL, ...) 
{
#polygons=F; bars=T
	twosex <- TRUE
	if (is.null(sexvec)) {
		twosex <- FALSE
	}
	if (length(unique(sexvec)) == 1) {
		twosex <- FALSE
	}
	male_mult <- 1
	if (twosex) {
		male_mult <- -1
	}
	yrvec <- sort(unique(yr))
	npanels <- length(yrvec)
	nvals <- length(yr)
	nrows <- min(ceiling(sqrt(npanels)), maxrows)
	ncols <- min(ceiling(npanels/nrows), maxcols)
	rc=.findSquare(npanels); nrows=rc[1]; ncols=rc[2]  ## (RH 200825)
	if (fixdims) {
		nrows <- maxrows
		ncols <- maxcols
	}
	allbins.obs <- sort(unique(ptsx))
	if (scalebins) {
		if (diff(range(allbins.obs)) > 0 && length(allbins.obs) > 
			2 && length(unique(diff(allbins.obs)))) {
			diffs <- diff(allbins.obs)
			diffs <- c(diffs, diffs[length(diffs)])
			bin.width.table <- data.frame(bin=allbins.obs, width=diffs)
		}
		else {
			scalebins <- FALSE
			warning("Setting scalebins=FALSE. Bins are equal length or too few.")
		}
	}
	npages <- ceiling(npanels/nrows/ncols)
	doSD <- length(ptsSD) == length(ptsx) & max(ptsSD) > 0
	if (doSD) {
		polygons <- FALSE
	}
	if (length(linesx) == 1 | length(linesy) == 1) {
		linepos <- 0
		linesx <- ptsx
		linesy <- ptsy
	}
	anyscaled <- FALSE
	if (bars & barwidth == "default") 
		barwidth <- 400/max(table(yr) + 2)/ncols
	if (length(size) == 1) {
		size <- rep(size, length(yr))
	}
	bub <- diff(range(size, na.rm=TRUE)) != 0
	xrange <- range(c(ptsx, linesx, ptsx, linesx))
	if (ymin0) {
		yrange <- c(0, max(ptsy, linesy))
	}
	else {
		yrange <- range(c(ptsy, linesy, ptsy, linesy))
	}
	yrange <- c(min(yrange[1], yupper), min(yrange[2], yupper))
	xrange_big <- xrange + c(-1, 1) * xbuffer * diff(xrange)
	yrange_big <- yrange + c(-1, 1) * ybuffer * diff(yrange)
	if (twosex & !bub) {
		yrange_big <- range(-yrange, yrange) + c(-1, 1) * ybuffer * diff(yrange)
	}
	yaxs_lab <- pretty(yrange)
	maxchar_yaxs <- max(nchar(yaxs_lab))
	if (is.null(ylas)) {
		if (maxchar_yaxs < 6) {
			ylas <- 1
		}
		else {
			ylas <- 0
		}
	}
	if (is.null(axis1)) {
		axis1 <- pretty(xrange)
	}
	if (is.null(axis1labs)) {
		axis1labs <- axis1
	}
	if (is.null(axis2)) {
		axis2 <- pretty(yrange)
	}
	if (length(sampsize) == 1) {
		sampsize <- 0
	}
	if (length(effN) == 1) {
		effN <- 0
	}
	par_old <- par()
	if (is.null(multifig_oma)) {
		if (main == "") {
			multifig_oma <- c(4, 5, 1, 1) + 0.1  ## (RH 200825)
		}
		else {
			multifig_oma <- c(5, 5, 5, 1) + 0.1
		}
	}
	#par(mfcol=c(nrows, ncols), mar=rep(0, 4), oma=multifig_oma, ...)
	par(mfrow=c(nrows, ncols), mar=rep(0, 4), oma=multifig_oma, ...)  ## (RH 200825)
	panelrange <- 1:npanels
	if (npages > 1 & ipage != 0) {
		panelrange <- intersect(panelrange, 1:(nrows * ncols) + nrows * ncols * (ipage - 1))
	}
	for (ipanel in panelrange) {
		yr_i <- yrvec[ipanel]
		sexvec_i <- sexvec[yr == yr_i]
		ptsx_i0 <- ptsx[yr == yr_i & sexvec == 0]
		ptsx_i1 <- ptsx[yr == yr_i & sexvec == 1]
		ptsx_i2 <- ptsx[yr == yr_i & sexvec == 2]
		ptsy_i0 <- ptsy[yr == yr_i & sexvec == 0]
		ptsy_i1 <- ptsy[yr == yr_i & sexvec == 1]
		ptsy_i2 <- ptsy[yr == yr_i & sexvec == 2] * male_mult
		if (doSD) {
			ptsSD_i0 <- ptsSD[yr == yr_i & sexvec == 0]
			ptsSD_i1 <- ptsSD[yr == yr_i & sexvec == 1]
			ptsSD_i2 <- ptsSD[yr == yr_i & sexvec == 2]
		}
		linesx_i0 <- linesx[yr == yr_i & sexvec == 0]
		linesx_i1 <- linesx[yr == yr_i & sexvec == 1]
		linesx_i2 <- linesx[yr == yr_i & sexvec == 2]
		linesy_i0 <- linesy[yr == yr_i & sexvec == 0]
		linesy_i1 <- linesy[yr == yr_i & sexvec == 1]
		linesy_i2 <- linesy[yr == yr_i & sexvec == 2] * male_mult
		linesy_i0 <- linesy_i0[order(linesx_i0)]
		linesx_i0 <- sort(linesx_i0)
		linesy_i1 <- linesy_i1[order(linesx_i1)]
		linesx_i1 <- sort(linesx_i1)
		linesy_i2 <- linesy_i2[order(linesx_i2)]
		linesx_i2 <- sort(linesx_i2)
		z_i0 <- size[yr == yr_i & sexvec == 0]
		z_i1 <- size[yr == yr_i & sexvec == 1]
		z_i2 <- size[yr == yr_i & sexvec == 2]
		scaled <- FALSE
		if (scalebins) {
			getwidths <- function(ptsx) {
				if (length(ptsx) > 0) {
					widths <- rep(NA, length(ptsx))
					for (ibin in 1:length(ptsx)) {
						widths[ibin] <- bin.width.table$width[bin.width.table$bin == ptsx[ibin]]
					}
				}
				else {
					widths <- NULL
				}
				return(widths)
			}
			widths_i0 <- getwidths(ptsx_i0)
			widths_i1 <- getwidths(ptsx_i1)
			widths_i2 <- getwidths(ptsx_i2)
			ptsy_i0 <- ptsy_i0/widths_i0
			ptsy_i1 <- ptsy_i1/widths_i1
			ptsy_i2 <- ptsy_i2/widths_i2
			linesy_i0 <- linesy_i0/widths_i0
			linesy_i1 <- linesy_i1/widths_i1
			linesy_i2 <- linesy_i2/widths_i2
			scaled <- TRUE
		}
		if (scaled) {
			anyscaled <- TRUE
			if (ylab == "Proportion") {
				ylab <- "Proportion / bin width"
			}
		}
		plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=xrange_big, ylim=yrange_big, xaxs="i", yaxs=ifelse(bars, "i", "r"), ...)
		abline(h=0, col="grey")
		if (linepos == 2) {
			lines(linesx_i0, linesy_i0, col=linescol[1], lwd=lwd, lty=lty)
			lines(linesx_i1, linesy_i1, col=linescol[2], lwd=lwd, lty=lty)
			lines(linesx_i2, linesy_i2, col=linescol[3], lwd=lwd, lty=lty)
		}
		if (bub) {
			if (length(z_i0) > 0) {
				bubble3(x=ptsx_i0, y=ptsy_i0, z=z_i0, col=rep(colvec[3], length(z_i0)), cexZ1=cexZ1, legend.yadj=1.5, legend=bublegend, legendloc="topright", maxsize=maxsize, minnbubble=minnbubble, allopen=allopen, add=TRUE)
			}
			if (length(z_i1) > 0) {
				bubble3(x=ptsx_i1, y=ptsy_i1, z=z_i1, col=rep(colvec[1], length(z_i1)), cexZ1=cexZ1, legend.yadj=1.5, legend=bublegend, legendloc="topright", maxsize=maxsize, minnbubble=minnbubble, allopen=allopen, add=TRUE)
			}
			if (length(z_i2) > 0) {
				bubble3(x=ptsx_i2, y=abs(ptsy_i2), z=z_i2, col=rep(colvec[2], length(z_i2)), cexZ1=cexZ1, legend.yadj=1.5, legend=bublegend, legendloc="topright", maxsize=maxsize, minnbubble=minnbubble, allopen=allopen, add=TRUE)
			}
			if (linepos == 0) 
				effNline <- 0
			if (effNline > 0 && length(effN) > 0) {
				effN_i1 <- effN[yr == yr_i]
				effN_i1_vec <- unlist(lapply(split(effN_i1, ptsy_i1), unique))
				ptsy_i1_vec <- sort(unique(ptsy_i1))
				lines(effNline * effN_i1_vec, ptsy_i1_vec, col="green3")
				if (!is.null(effNmean)) {
					lines(rep(effNline * effNmean, length(ptsy_i1_vec)), ptsy_i1_vec, col="green3", lty=2)
				}
			}
			if (sampsizeline > 0 && length(sampsize) > 0) {
				sampsize_i1 <- sampsize[yr == yr_i]
				sampsize_i1_vec <- unlist(lapply(split(sampsize_i1, ptsy_i1), unique))
				ptsy_i1_vec <- sort(unique(ptsy_i1))
				lines(sampsizeline * sampsize_i1_vec, ptsy_i1_vec, col=2)
				if (!is.null(sampsizemean)) {
					lines(rep(sampsizeline * sampsizemean, length(ptsy_i1_vec)), ptsy_i1_vec, col=2, lty=3)
				}
			}
		}
		else {
			if (length(ptsx_i0) > 0) {
				if (bars) {
					drawBars(ptsx_i0, ptsy_i0, col="slategray3", fill="slategray3")  ## (RH 200825)
				} else {
					if (polygons)
						polygon(c(ptsx_i0[1], ptsx_i0, tail(ptsx_i0,1)), c(0,ptsy_i0, 0), col=multifig_colpolygon[1])
					points(ptsx_i0, ptsy_i0, type=type, lwd=1,  pch=16, cex=0.7, col=ptscol)
				}
			}
			if (length(ptsx_i1) > 0) {
				if (bars) {
					drawBars(ptsx_i1, ptsy_i1, col=lucent("pink",0.75), fill=lucent("pink",0.75))  ## (RH 200826) -- females
				} else {
					if (polygons)
						polygon(c(ptsx_i1[1], ptsx_i1, tail(ptsx_i1,1)), c(0,ptsy_i1, 0), col=multifig_colpolygon[1])
					points(ptsx_i1, ptsy_i1, type=type, lwd=1,  pch=16, cex=0.7, col=ptscol)
				}
			}
			if (length(ptsx_i2) > 0) {
				if (bars) {
					drawBars(ptsx_i2, ptsy_i2, col=lucent("skyblue",0.75), fill=lucent("skyblue",0.75))  ## (RH 200826) -- males
				} else {
					if (polygons)
						polygon(c(ptsx_i2[1], ptsx_i2, tail(ptsx_i2,1)), c(0,ptsy_i2, 0), col=multifig_colpolygon[1])
					points(ptsx_i2, ptsy_i2, type=type, lwd=1,  pch=16, cex=0.7, col=ptscol)
				}
			}
			if (doSD) {
				old_warn <- options()$warn
				options(warn=-1)
				if (length(ptsx_i0) > 0) {
					arrows(x0=ptsx_i0, y0=qnorm(p=0.05, mean=ptsy_i0, sd=ptsSD_i0), x1=ptsx_i0, y1=qnorm(p=0.95, mean=ptsy_i0, sd=ptsSD_i0), length=0.01, angle=90, code=3, col=ptscol)
				}
				if (length(ptsx_i1) > 0) {
					arrows(x0=ptsx_i1, y0=qnorm(p=0.05, mean=ptsy_i1, sd=ptsSD_i1), x1=ptsx_i1, y1=qnorm(p=0.95, mean=ptsy_i1, sd=ptsSD_i1), length=0.01, angle=90, code=3, col=ptscol)
				}
				if (length(ptsx_i2) > 0) {
					arrows(x0=ptsx_i2, y0=qnorm(p=0.05, mean=ptsy_i2, sd=ptsSD_i2), x1=ptsx_i2, y1=qnorm(p=0.95, mean=ptsy_i2, sd=ptsSD_i2), length=0.01, angle=90, code=3, col=ptscol)
				}
				options(warn=old_warn)
			}
		}
		if (linepos == 1) {
			lines(linesx_i0, linesy_i0, col=linescol[1], lwd=lwd, lty=lty)
			lines(linesx_i1, linesy_i1, col=linescol[2], lwd=lwd, lty=lty)
			lines(linesx_i2, linesy_i2, col=linescol[3], lwd=lwd, lty=lty)
		}
#browser();return()
		usr <- par("usr")
		for (i in 1:nlegends) {
			text_i <- ""
			text_i2 <- ""
			legtext_i <- legtext[[i]]
			if (length(legtext_i) == 1) {
				if (legtext_i == "yr") {
					text_i <- yr_i
				}
				for (sex in sort(unique(sexvec_i))) {
					if (legtext_i == "sampsize" & showsampsize) {
						vals <- unique(sampsize[sexvec == sex & yr == yr_i])
						if (length(vals) > 1) {
							warning("sampsize values are not all equal", "--choosing the first value: ", vals[1], "\n", "  yr=", yr_i, ", and all sampsize values: ", paste(vals, collapse=","), sep="")
							vals <- vals[1]
						}
						text_i <- paste(sampsize_label, round(vals, sampsizeround), sep="")
						if (twosex & sex == 2) {
							text_i2 <- paste(sampsize_label, round(vals, sampsizeround), sep="")
						}
					}
					if (legtext_i == "effN" & showeffN) {
						vals <- unique(effN[sexvec == sex & yr == yr_i])
						if (length(vals) > 1) {
							warning("effN values are not all equal", "--choosing the first value: ", vals[1], "\n", "  yr=", yr_i, ", and all effN values: ", paste(vals, collapse=","), sep="")
							vals <- vals[1]
						}
						text_i <- paste(effN_label, round(vals, sampsizeround), sep="")
						if (twosex & sex == 2) {
							text_i2 <- paste(effN_label, round(vals, sampsizeround), sep="")
						}
					}
				}
			}
			if (length(legtext_i) == nvals) {
				text_i <- legtext_i[yr == yr_i][1]
			}
			if (length(legtext_i) == 1) {
				text_i <- text_i
			}
			if (legx[1] == "default") {
				textx <- ifelse(i == 1, usr[1]+(0.015*diff(usr[1:2])), usr[2]-(0.025*diff(usr[1:2])))  ## (RH 200825)
			}
			else {
				textx <- legx[i]
			}
			if (legy[1] == "default") {
				texty <- usr[4]
				texty2 <- usr[3]
			}
			else {
				texty <- legy[i]
				texty2 <- -legy[i]
			}
			if (legadjx[1] == "default") {
				adjx <- ifelse(i == 1, -0.1, 1)
			}
			else {
				adjx <- legadjx[i]
			}
#if (i==2) {browser();return()}
			if (legadjy[1] == "default") {
				adjy <- ifelse(i < 3, 1.3, 1.3 + 1.3 * (i - 2))
			}
			else {
				adjy <- legadjy[i]
			}
			text(x=textx, y=texty, labels=text_i, adj=c(adjx, 
				adjy), cex=legsize[i], font=legfont[i])
			if (text_i2 != text_i & text_i2 != "") {
				text(x=textx, y=texty2, labels=text_i, adj=c(adjx, -adjy), cex=legsize[i], font=legfont[i])
			}
			if (twosex & !bub & venusmars) {
				pu <- par("usr")
				xval <- pu[2]
				if (length(ptsx_i0) > 0) {
					text(xval, 0.5 * yrange[2], "\\VE+\\MA", vfont=c("serif", "plain"), cex=2, col=linescol[1], pos=2)
				}
				if (length(ptsx_i1) > 0) {
					text(xval, 0.5 * yrange[2], "\\VE", vfont=c("serif", "plain"), cex=2, col=linescol[2], pos=2)
				}
				if (length(ptsx_i2) > 0) {
					text(xval, -0.5 * yrange[2], "\\MA", vfont=c("serif", "plain"), cex=2, col=linescol[3], pos=2)
				}
			}
		}
		mfg <- par("mfg")
		if (mfg[1] == mfg[3] | ipanel == npanels) {
			axis(side=1, at=axis1, labels=axis1labs, las=xlas)
		}
		if (mfg[2] == 1) {
			axis(side=2, at=axis2, las=ylas)
			if (twosex) {
				axis(side=2, at=-axis2[axis2 > 0], labels=format(axis2[axis2 > 0]), las=ylas)
			}
		}
		box()
		if (npanels == 1 | ipanel%%(nrows * ncols) == 1) {
			fixcex <- 1
			if (max(nrows, ncols) == 2) {
				fixcex <- 1/0.83
			}
			if (max(nrows, ncols) > 2) {
				fixcex <- 1/0.66
			}
			if (npanels > 1) {
				title(main=main, line=c(2, 0, 3, 3), outer=TRUE, cex.main=cex.main * fixcex)
				title(xlab=xlab, outer=TRUE, cex.lab=fixcex)
				title(ylab=ylab, line=ifelse(ylas %in% 1:2, max(3, 2 + 0.4 * maxchar_yaxs), 3.5), outer=TRUE, cex.lab=fixcex)
			}
			else {
				title(main=main, xlab=xlab, ylab=ylab, outer=TRUE, cex.main=cex.main)
			}
		}
	}
	par(mfcol=par_old$mfcol, mar=par_old$mar, oma=par_old$oma)
	if (anyscaled) {
		cat("Note: compositions have been rescaled by dividing by binwidth\n")
	}
	return(list(npages=npages, npanels=npanels, ipage=ipage))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~make.multifig


## plotSS.comps-------------------------2020-08-25
##  Plot age proportions and model fit
##  Source: R package 'r4ss' v.1.39.1
## ----------------------------------------r4ss|RH
plotSS.comps = function (replist, subplots=c(1:21, 24), kind="LEN", sizemethod=1, 
   aalyear=-1, aalbin=-1, plot=TRUE, print=FALSE, fleets="all", 
   fleetnames="default", sexes="all", yupper=0.4, datonly=FALSE, 
   samplesizeplots=TRUE, compresidplots=TRUE, bub=FALSE, 
   showyears=TRUE, showsampsize=TRUE, showeffN=TRUE, aggregates_by_mkt=FALSE, 
   sampsizeline=FALSE, effNline=FALSE, minnbubble=3, pntscalar=NULL, 
   scalebubbles=FALSE, cexZ1=1.5, bublegend=TRUE,
   #colvec=c(rgb(1,0,0,0.7), rgb(0,0,1,0.7), rgb(0.1,0.1,0.1,0.7)),
   #linescol=c(rgb(0,0.5,0,0.7), rgb(0.8,0,0,0.7), rgb(0,0,0.8,0.7)),
   colsex = c("orange","limegreen","slategray3"),
   colvec=lucent(colsex,0.75),
   #linescol=rgb(RYB2RGB(1-RGB2RYB(col2rgb(colsex)))),
   linescol=rep(lucent("black",0.75),3),
   xlas=0, ylas=NULL, axis1=NULL, axis2=NULL, 
   axis1labs=NULL, sizebinlabs=NULL, blue=lucent("limegreen",0.75), red=lucent("purple",0.75),
   #blue=rgb(0,0,1,0.7), red=rgb(1,0,0,0.7),
   pwidth=6.5, pheight=5, punits="in", ptsize=10, res=400,
   plotdir="default", cex.main=1, linepos=1, fitbar=FALSE, do.sqrt=TRUE,
   smooth=TRUE, cohortlines=c(), labels=c("Length (cm)", 
   "Age (yr)", "Year", "Observed sample size", "Effective sample size",
   "Proportion", "cm", "Frequency", "Weight", "Length", "(mt)",
   "(numbers x1000)", "Stdev (Age)", "Conditional AAL plot, ", "Size bin"),
   printmkt=TRUE, printsex=TRUE, maxrows=4, maxcols=4, maxrows2=2, maxcols2=4,
   rows=1, cols=1, andre_oma=c(3,0,3,0), andrerows=3, fixdims=TRUE,
   fixdims2=FALSE, maxneff=5000, verbose=TRUE, scalebins=FALSE, 
   addMeans=TRUE, mainTitle=FALSE, outnam, lang="e", ...) 
{
	oldpar = par(no.readonly=TRUE)
	fart = function(opar) { if (any("windows"%in%names(dev.list()))) par(opar); eop() }
	on.exit(fart(oldpar))
	changeLangOpts(L=lang)
	if (missing(outnam))
		outnam = NULL
	ttput(outnam)

	if (!is.element(subplots, c(1:21, 24))) {
		mess = c(" \n","Choose another subplot from:",
			"  1 = Multipanel (by year) fits to age proportions for multiple fleets",
			"  2 = Single panel (by fleet) bubble plots of age proportions by year",
			"  3 = ???",
			"  4 = ???",
			"  5 = ???",
			"  6 = ???",
			"  7 = ???",
			"  8 = ???",
			"  9 = ???",
			" 11 = ???",
			" 12 = ???",
			" 13 = ???",
			" 14 = ???",
			" 15 = ???",
			" 21 = Multipanel (by fleet) fits to combined age proportions across years",
			" 24 = Multipanel (by fleet) bubble plots of age proportions by year"
		)
		stop(paste0(mess, collapse="\n"))
	}
	pngfun <- function(file, caption=NA, lang="e") {
		addsize = ifelse(length(.su(dbase$Yr))>9, 1, 0)
#browser();return()
		createFdir(lang, dir=plotdir)
		changeLangOpts(L=lang)
		fout = switch(lang, 'e' = file.path(plotdir, file), 'f' = file.path(plotdir,"french", file) )
		clearFiles(fout)
		png(filename=fout, width=pwidth+addsize, height=pheight+addsize, units=punits, res=res, pointsize=ptsize)
		plotinfo <- rbind(plotinfo, data.frame(file=file, caption=caption))
		return(plotinfo)
	}
	plotinfo <- NULL
	SS_versionNumeric <- replist$SS_versionNumeric
	lendbase          <- replist$lendbase
	sizedbase         <- replist$sizedbase
	agedbase          <- replist$agedbase
	condbase          <- replist$condbase
	ghostagedbase     <- replist$ghostagedbase
	ghostlendbase     <- replist$ghostlendbase
	ladbase           <- replist$ladbase
	wadbase           <- replist$wadbase
	tagdbase1         <- replist$tagdbase1
	tagdbase2         <- replist$tagdbase2
	nfleets           <- replist$nfleets
	nseasons          <- replist$nseasons
	seasfracs         <- replist$seasfracs
	FleetNames        <- replist$FleetNames
	nsexes            <- replist$nsexes
	accuage           <- replist$accuage
	Age_tuning        <- replist$Age_comp_Eff_N_tuning_check
	titles            <- NULL
	titlemkt          <- ""
	if (plotdir == "default") {
		plotdir <- replist$inputs$dir
	}
	if (fleets[1] == "all") {
		fleets <- 1:nfleets
	}
	else {
		if (length(intersect(fleets, 1:nfleets)) != length(fleets)) {
			stop("Input 'fleets' should be 'all' or a vector of values between 1 and nfleets.")
		}
	}
	if (fleetnames[1] == "default") {
		fleetnames <- FleetNames
	}
	if (sexes[1] == "all") {
		sexes <- 0:nsexes
	}
	if (nsexes == 1) {
		sexes <- 0:nsexes
	}
	if (nsexes == 1 | length(sexes) > 1) {
		titlesex <- ""
		filesex <- ""
	}
	if (nsexes > 1 & length(sexes) == 1) {
		if (sexes == 0) {
			titlesex <- "sexes combined, "
			filesex <- "sex0"
		}
		if (sexes == 1) {
			titlesex <- "female, "
			filesex <- "sex1"
		}
		if (sexes == 2) {
			titlesex <- "male, "
			filesex <- "sex2"
		}
	}
	titlesex <- ifelse(printsex, titlesex, "")
	if (kind == "LEN") {
		dbase_kind <- lendbase
		kindlab=labels[1]
		if (datonly) {
			filenamestart <- "comp_lendat_"
			titledata <- "Length comp data, "
		}
		else {
			filenamestart <- "comp_lenfit_"
			titledata <- "Length comps, "
		}
	}
	if (kind == "GSTLEN") {
		dbase_kind       <- ghostlendbase
		kindlab=labels[1]
		if (datonly) {
			filenamestart <- "comp_gstlendat_"
			titledata     <- "Ghost length comp data, "
		}
		else {
			filenamestart <- "comp_gstlenfit_"
			titledata     <- "Ghost length comps, "
		}
	}
	if (kind == "SIZE") {
		dbase_kind       <- sizedbase[sizedbase$method == sizemethod,]
		if (!is.null(sizebinlabs)) {
			kindlab <- labels[15]
			axis1 <- sort(unique(dbase_kind$Bin))
			if (length(sizebinlabs) == length(axis1)) {
				axis1labs <- sizebinlabs
			}
			else {
				axis1labs <- axis1
				warning("Input 'sizebinlabs' differs in length from the unique Bin\n", "  values associated with sizemethod=", sizemethod, ". Using bin values instead.")
			}
		}
		else {
			sizeunits <- unique(dbase_kind$units)
			if (length(sizeunits) > 1) {
				stop("!error with size units in generalized size comp plots:\n", "   more than one unit value per method.\n")
			}
			if (sizeunits %in% c("in", "cm")) {
				kindlab <- paste(labels[10], " (", sizeunits, ")", sep="")
			}
			if (sizeunits %in% c("lb", "kg")) {
				kindlab <- paste(labels[9], " (", sizeunits, ")", sep="")
			}
		}
		if (datonly) {
			filenamestart <- "comp_sizedat_"
			titledata <- "Size comp data, "
		}
		else {
			filenamestart <- "comp_sizefit_"
			titledata <- "Size comps, "
		}
		if (length(unique(sizedbase$method)) > 1) {
			filenamestart <- paste0(filenamestart, "method", sizemethod, "_")
			titledata <- paste0(titledata, " size method ", sizemethod, ", ")
		}
	}
	if (kind == "AGE") {
		dbase_kind <- agedbase
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_agedat_"
			titledata <- "Age comp data, "
		}
		else {
			filenamestart <- "comp_agefit_"
			titledata <- "Age comps, "
		}
	}
	if (kind == "cond") {
		dbase_kind <- condbase
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_condAALdat_"
			titledata <- "Conditional age-at-length data, "
		}
		else {
			filenamestart <- "comp_condAALfit_"
			titledata <- "Conditional age-at-length, "
		}
	}
	if (kind == "GSTAGE") {
		dbase_kind <- ghostagedbase
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_gstagedat_"
			titledata <- "Ghost age comp data, "
		}
		else {
			filenamestart <- "comp_gstagefit_"
			titledata <- "Ghost age comps, "
		}
	}
	if (kind == "GSTcond") {
		dbase_kind <- ghostagedbase
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_gstCAALdat_"
			titledata <- "Ghost conditional age-at-length data, "
		}
		else {
			filenamestart <- "comp_gstCAALfit_"
			titledata <- "Ghost conditional age-at-length comps, "
		}
	}
	if (kind == "L@A") {
		dbase_kind <- ladbase[ladbase$Nsamp_adj != 0, ]
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_LAAdat_"
			titledata <- "Mean length at age data, "
		}
		else {
			filenamestart <- "comp_LAAfit_"
			titledata <- "Mean length at age fit, "
		}
		dbase_kind$SD <- dbase_kind$Lbin_lo/dbase_kind$Nsamp_adj
	}
	if (kind == "W@A") {
		dbase_kind <- wadbase[wadbase$Nsamp_adj != 0, ]
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_WAAdat_"
			titledata <- "Mean weight at age data, "
		}
		else {
			filenamestart <- "comp_WAAfit_"
			titledata <- "Mean weight at age fit, "
		}
		dbase_kind$SD <- dbase_kind$Lbin_lo/dbase_kind$Nsamp_adj
	}
	if (!(kind %in% c("LEN", "SIZE", "AGE", "cond", "GSTAGE", "GSTLEN", "L@A", "W@A"))) {
		stop("Input 'kind' to plotSS.comps needs to be one of the following:\n  ", "'LEN','SIZE','AGE','cond','GSTAGE','GSTLEN','L@A','W@A'.")
	}
	if (nrow(dbase_kind) > 0) {
		if (aggregates_by_mkt) {
			dbase_kind$Part_group <- dbase_kind$Part
		}
		else {
			dbase_kind$Part_group <- -1
		}
	}
	if (any(dbase_kind$SuprPer == "Sup" & dbase_kind$Used == "skip")) {
		cat("Note: removing super-period composition values labeled 'skip'\n", "   and designating super-period values with a '*'\n")
		dbase_kind <- dbase_kind[dbase_kind$SuprPer == "No" | dbase_kind$Used != "skip", ]
		dbase_kind$YrSeasName <- paste(dbase_kind$YrSeasName, ifelse(dbase_kind$SuprPer == "Sup", "*", ""), sep="")
	}
	ageerr_warning <- TRUE
	dbase_kind <- dbase_kind[dbase_kind$Fleet %in% fleets & dbase_kind$sex %in% sexes, ]
	for (f in fleets) {
		if (length(dbase_kind$Obs[dbase_kind$Fleet == f]) > 0) {
			dbasef <- dbase_kind[dbase_kind$Fleet == f, ]
			if (kind %in% c("cond", "GSTcond") && f %in% Age_tuning$Fleet) {
				HarmEffNage <- NULL
				MeanNage    <- NULL
			}
			else {
				HarmEffNage <- NULL
				MeanNage    <- NULL
			}
			dbase_k <- dbasef
			for (j in unique(dbase_k$Part)) {
				dbase <- dbase_k[dbase_k$Part == j, ]
				max_n_ageerr <- max(apply(table(dbase$Yr.S, dbase$Ageerr) > 0, 1, sum))
				if (max_n_ageerr > 1) {
					if (ageerr_warning) {
						cat("Note: multiple samples with different ageing error types within fleet/year.\n", "   Plots label '2005a3' indicates ageing error type 3 for 2005 sample.\n", "   Bubble plots may be misleading with overlapping bubbles.\n")
						ageerr_warning <- FALSE
					}
					dbase$Yr.S <- dbase$Yr.S + dbase$Ageerr/1000
					dbase$YrSeasName <- paste(dbase$YrSeasName, "a", dbase$Ageerr, sep="")
				}
				if (j == 0) 
					titlemkt <- "whole catch, "
				if (j == 1) 
					titlemkt <- "discard, "
				if (j == 2) 
					titlemkt <- "retained, "
				titlemkt <- ifelse(printmkt, titlemkt, "")
				if (datonly | fitbar) 
					bars <- TRUE
				else bars <- FALSE
				title_sexmkt <- paste(titlesex, titlemkt, sep="")
				filename_fltsexmkt <- paste("flt", f, filesex, "mkt", j, sep="")
				if (1 %in% subplots & kind != "cond") {
					caption <- paste(titledata, title_sexmkt, fleetnames[f], sep="")
					if (mainTitle) {
						ptitle <- caption
					}
					else {
						ptitle <- ""
					}
					titles <- c(ptitle, titles)
					tempfun <- function(ipage, l="e", ...) {
						sexvec <- dbase$sex
#browser();return()
						if (!(kind %in% c("GSTAGE", "GSTLEN", "L@A", "W@A"))) {
							if ("DM_effN" %in% names(dbase) && any(!is.na(dbase$DM_effN))) {
								make.multifig(ptsx=dbase$Bin, ptsy=dbase$Obs, 
									yr=dbase$Yr.S, linesx=dbase$Bin, 
									linesy=dbase$Exp, sampsize=dbase$Nsamp_adj, 
									effN=dbase$DM_effN, showsampsize=showsampsize, 
									showeffN=showeffN, sampsize_label="N input=", 
									effN_label="N samp.=", bars=bars, 
									linepos=(1 - datonly) * linepos, 
									nlegends=3, legtext=list(dbase$YrSeasName, "sampsize", "effN"),
									main=ptitle, cex.main=cex.main, xlab=kindlab, 
									ylab=labels[6], maxrows=maxrows, 
									maxcols=maxcols, rows=rows, cols=cols, 
									fixdims=fixdims, ipage=ipage, scalebins=scalebins, 
									colvec=colvec, linescol=linescol, 
									xlas=xlas, ylas=ylas, axis1=axis1, 
									axis2=axis2, axis1labs=axis1labs, 
									sexvec=sexvec, yupper=yupper, lang=l, ...)
							}
							else {
								if (all(dbase$Nsamp_adj==dbase$Nsamp_in)) {
									sslab = "N samps obs = "
									obsN  = 0
								} else {
									sslab = "N samps wtd = "
									obsN = dbase$Nsamp_in
								}
#browser();return()
								make.multifig(ptsx=dbase$Bin, ptsy=dbase$Obs, 
									yr=dbase$Yr.S, linesx=dbase$Bin, 
									linesy=dbase$Exp, sampsize=dbase$Nsamp_adj, 
									effN=dbase$effN, showsampsize=showsampsize, 
									showeffN=showeffN, sampsize_label=sslab, 
									effN_label="N eff.=", bars=bars, 
									linepos=(1 - datonly) * linepos, 
									nlegends=3, legtext=list(dbase$YrSeasName, "sampsize", "effN"),
									main=ptitle, cex.main=cex.main, xlab=kindlab, 
									ylab=labels[6], maxrows=maxrows, 
									maxcols=maxcols, rows=rows, cols=cols, 
									fixdims=fixdims, ipage=ipage, scalebins=scalebins, 
									colvec=colvec, linescol=linescol, 
									xlas=xlas, ylas=ylas, axis1=axis1, 
									axis2=axis2, axis1labs=axis1labs, 
									sexvec=sexvec, yupper=yupper, lang=l, obsN=obsN, ...)
							}
						}
						if (kind == "GSTAGE") {
							make.multifig(ptsx=dbase$Bin, ptsy=dbase$Obs, 
								yr=dbase$Yr.S, linesx=dbase$Bin, 
								linesy=dbase$Exp, sampsize=dbase$Nsamp_adj, 
								effN=dbase$effN, showsampsize=FALSE, 
								showeffN=FALSE, bars=bars,
								linepos=(1-datonly) * linepos, nlegends=3,
								legtext=list(dbase$YrSeasName, "sampsize", "effN"),
								main=ptitle, cex.main=cex.main, xlab=kindlab, 
								ylab=labels[6], maxrows=maxrows, 
								maxcols=maxcols, rows=rows, cols=cols, 
								fixdims=fixdims, ipage=ipage, scalebins=scalebins, 
								colvec=colvec, linescol=linescol, 
								xlas=xlas, ylas=ylas, axis1=axis1, 
								axis2=axis2, axis1labs=axis1labs, 
								sexvec=sexvec, yupper=yupper, lang=l, ...)
						}
						if (kind == "GSTLEN") {
							make.multifig(ptsx=dbase$Bin, ptsy=dbase$Obs, 
								yr=dbase$Yr.S, linesx=dbase$Bin, 
								linesy=dbase$Exp, sampsize=dbase$Nsamp_adj, 
								effN=dbase$effN, showsampsize=FALSE, 
								showeffN=FALSE, bars=bars,
								linepos=(1-datonly) * linepos, nlegends=3,
								legtext=list(dbase$YrSeasName, "sampsize", "effN"),
								main=ptitle, cex.main=cex.main, xlab=kindlab, 
								ylab=labels[6], maxrows=maxrows, 
								maxcols=maxcols, rows=rows, cols=cols, 
								fixdims=fixdims, ipage=ipage, scalebins=scalebins, 
								colvec=colvec, linescol=linescol, 
								xlas=xlas, ylas=ylas, axis1=axis1, 
								axis2=axis2, axis1labs=axis1labs, 
								sexvec=sexvec, lang=l, ...)
						}
						if (kind %in% c("L@A", "W@A")) {
							make.multifig(ptsx=dbase$Bin, ptsy=dbase$Obs, 
								yr=dbase$Yr.S, linesx=dbase$Bin, 
								linesy=dbase$Exp, ptsSD=dbase$SD, 
								sampsize=dbase$Nsamp_adj, effN=0, 
								showsampsize=FALSE, showeffN=FALSE, 
								nlegends=1, legtext=list(dbase$YrSeasName), 
								bars=bars, linepos=(1-datonly) * linepos,
								main=ptitle, cex.main=cex.main, 
								xlab=kindlab, ylab=ifelse(kind=="W@A", labels[9], labels[1]),
								maxrows=maxrows, maxcols=maxcols, rows=rows, cols=cols, 
								fixdims=fixdims, ipage=ipage, scalebins=scalebins, 
								colvec=colvec, linescol=linescol, 
								xlas=xlas, ylas=ylas, axis1=axis1, 
								axis2=axis2, axis1labs=axis1labs, 
								sexvec=sexvec, lang=l, ...)
						}
					} ## end tempfun
					if (plot) 
						tempfun(ipage=0, l=lang, ...)
					if (print) {
						npages <- if (fixdims) ceiling(length(unique(dbase$Yr.S))/maxrows/maxcols) else 1  ## (RH 200825)
						for (ipage in 1:npages) {
							pagetext <- ""
							caption_count <- ""
							if (npages > 1) {
								pagetext <- paste0("_page", ipage)
								caption_count <- paste0(" (plot ", ipage, " of ", npages, ")")
							}
							caption_extra <- ""
							if (ipage == 1) {
								if ("DM_effN" %in% names(dbase) && any(!is.na(dbase$DM_effN))) {
									ipar <- replist$age_data_info$ParmSelect[f]
									Theta <- as.numeric(replist$Dirichlet_Multinomial_pars$Theta[ipar])
									caption_extra <- paste0(".<br><br>'N input' is the input sample size. ", "'N samp.' is the sample size after adjustment by the ", "Dirichlet-Multinomial <i>&#920</i> parameter based on the ", "formula N samp.=1 / (1+<i>&#920</i>) + N * <i>&#920</i> / (1+<i>&#920</i>). ", "<br><br>For this fleet, <i>&#920</i>=", round(Theta, 3), " and the sample size multiplier is approximately ", "<i>&#920</i> / (1+<i>&#920</i>)=", round(Theta/(1 + Theta), 3), "<br><br>For more info, see<br>", "<blockquote>", "Thorson, J.T., Johnson, K.F., ", "Methot, R.D. and Taylor, I.G. 2017. ", "Model-based estimates of effective sample size ", "in stock assessment models using the ", "Dirichlet-multinomial distribution. ", "<i>Fisheries Research</i>", "192: 84-93. ", "<a href=https://doi.org/10.1016/j.fishres.2016.06.005>", "https://doi.org/10.1016/j.fishres.2016.06.005</a>", "</blockquote>")
								}
								else {
									caption_extra <- paste0(".<br><br>'N samp.' is the input sample size ", "after data-weighting adjustment. ", "N eff. is the calculated effective sample size used ", "in the McAllister-Iannelli tuning method.")
								}
							}
							if (!is.null(ttcall(outnam)))
								file = paste0(sub("\\.png$","",outnam),f,".png")
							else
								file <- paste(filenamestart, filename_fltsexmkt, pagetext, ".png", sep="")
#browser();return()
							plotinfo <- pngfun(file=file, caption=paste0(caption, caption_count, caption_extra), lang=lang)
							tempfun(ipage=ipage, l=lang, ...)
							dev.off(); eop()
						}
					}
				}
				if (datonly) {
					z <- dbase$Obs
					if (scalebubbles) {
						z <- dbase$Nsamp_adj * dbase$Obs
					}
					col <- rep("black", 2)
					titletype <- titledata
					filetype <- "bub"
					allopen <- TRUE
				}
				else {
					z <- dbase$Pearson
					col <- rep(colvec[3], 2)
					titletype <- "Pearson residuals, "
					filetype <- "resids"
					allopen <- FALSE
				}
#browser();return()
				if (2 %in% subplots & bub & kind != "cond") {
					if (length(cohortlines) > 0) {
						growdat <- replist$endgrowth
						growdatF <- growdat[growdat$Sex == 1 & growdat$Morph == min(growdat$Morph[growdat$Sex == 1]), ]
						if (nsexes > 1) {
							growdatM <- growdat[growdat$Sex == 2 & growdat$Morph == min(growdat$Morph[growdat$Sex == 2]), ]
						}
					}
					caption <- paste(titletype, title_sexmkt, fleetnames[f], sep="")
					caption <- paste(caption, " (max=", round(max(z), digits=2), ")", sep="")
					if (mainTitle) {
						ptitle <- caption
					}
					else {
						ptitle <- ""
					}
					titles <- c(ptitle, titles)

					tempfun2 <- function(l="e") {
						xvals <- dbase$Yr.S
						xdiff <- 0.1 * sort(unique(diff(sort(unique(dbase$Yr.S)), na.rm=TRUE)))[1]
						if (is.na(xdiff)) {
							xdiff <- 0.1
						}
						cols <- rep(colvec[3], nrow(dbase))
						if (nsexes > 1) {
							xvals[dbase$sex > 0] <- dbase$Yr.S[dbase$sex > 0] - (dbase$sex[dbase$sex > 0] - 1.5) * xdiff
							if (length(unique(dbase$Yr.S)) == 1) {
								xvals[dbase$sex > 0] <- floor(dbase$Yr.S[dbase$sex > 0]) - (dbase$sex[dbase$sex > 0] - 1.5) * xdiff
							}
							cols[dbase$sex > 0] <- colvec[dbase$sex[dbase$sex > 0]]
						}
						r4ss:::bubble3(x=xvals, y=dbase$Bin, z=z, xlab=linguaFranca(labels[3],l), ylab=linguaFranca(kindlab,l), col=cols, cexZ1=cexZ1, legend=linguaFranca(bublegend,l), las=1, main=linguaFranca(ptitle,l), cex.main=cex.main, maxsize=pntscalar, allopen=allopen, minnbubble=minnbubble)
						if (length(cohortlines) > 0) {
							for (icohort in 1:length(cohortlines)) {
								cat("  Adding line for", cohortlines[icohort], "cohort\n")
								if (kind == "LEN") {
									if (nsexes > 1) {
										lines(growdatF$Age_Mid + cohortlines[icohort], growdatF$Len_Mid, col=colvec[1])
										lines(growdatM$Age_Mid + cohortlines[icohort], growdatM$Len_Mid, col=colvec[2])
									}
									else {
										lines(growdatF$Age_Mid + cohortlines[icohort], growdatF$Len_Mid, col=colvec[3])
									}
								}
								if (kind %in% c("AGE", "GSTAGE")) {
									lines(0.5 + c(cohortlines[icohort], cohortlines[icohort] + accuage), c(0, accuage), col=colvec[3], lty=3)
								}
							}
						}
					} ## end tempfun2 
					if (plot) 
						tempfun2(l=lang)
					if (print) {
						pagetext <- ""
						if (npages > 1) {
							pagetext <- paste("_page", ipage, sep="")
							caption <- paste(caption, " (plot ", ipage, " of ", npages, ")", sep="")
						}
						if (length(grep("Pearson", caption)) > 0) {
							caption <- paste(caption, "<br> \nClosed bubbles are positive residuals", "(observed > expected)", "and open bubbles are negative residuals", "(observed < expected).")
						}
						if (!is.null(ttcall(outnam)))
							file = paste0(sub("\\.png$","",outnam),f,".png")
						else
							file <- paste(filenamestart, filetype, filename_fltsexmkt, pagetext, ".png", sep="")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						tempfun2(l=lang)
						dev.off(); eop()
					}
				}
				if (3 %in% subplots & kind == "cond") {
					caption <- paste(titletype, title_sexmkt, fleetnames[f], sep="")
					caption <- paste(caption, " (max=", round(max(z), digits=2), ")", sep="")
					if (mainTitle) {
						ptitle <- caption
					}
					else {
						ptitle <- ""
					}
					titles <- c(ptitle, titles)
					sampsizeline.old <- sampsizeline
					effNline.old <- effNline
					if (is.logical(sampsizeline) && sampsizeline) {
						sampsizeline <- max(dbase$Bin)/max(dbase$Nsamp_adj, na.rm=TRUE)
						if (!datonly && is.logical(effNline) && effNline) {
							sampsizeline <- effNline <- max(dbase$Bin)/max(dbase$Nsamp_adj, dbase$effN, na.rm=TRUE)
							cat("  Fleet ", f, " ", titlesex, "adj. input & effective N in red & green scaled by ", effNline, "\n", sep="")
						}
						else {
							cat("  Fleet ", f, " ", titlesex, "adj. input N in red scaled by ", sampsizeline, "\n", sep="")
						}
					}
					tempfun3 <- function(ipage, l="e", ...) {
						sexvec <- dbase$sex
						col.index <- sexvec
						col.index[col.index == 0] <- 3
						cols <- colvec[col.index]
						yrvec <- dbase$Yr.S + dbase$sex * 1e-06
						make.multifig(ptsx=dbase$Bin, ptsy=dbase$Lbin_mid, 
							yr=yrvec, size=z, sampsize=dbase$Nsamp_adj, 
							showsampsize=showsampsize, effN=dbase$effN, 
							showeffN=FALSE, cexZ1=cexZ1, bublegend=bublegend, 
							nlegends=1, legtext=list(dbase$YrSeasName), 
							bars=FALSE, linepos=0, main=ptitle, 
							cex.main=cex.main, xlab=labels[2], 
							ylab=labels[1], ymin0=FALSE, maxrows=maxrows2, 
							maxcols=maxcols2, fixdims=fixdims, 
							allopen=allopen, minnbubble=minnbubble, 
							ptscol=cols, ipage=ipage, scalebins=scalebins, 
							sampsizeline=sampsizeline, effNline=effNline, 
							sampsizemean=MeanNage, effNmean=HarmEffNage, 
							colvec=colvec, linescol=linescol, xlas=xlas, 
							ylas=ylas, axis1=axis1, axis2=axis2, 
							axis1labs=axis1labs, sexvec=sexvec, lang=l, ...)
					} ## end tempfun3
					if (plot) 
						tempfun3(ipage=0, l=lang, ...)
					if (print) {
						npages <- if (fixdims) ceiling(length(unique(dbase$Yr.S)) * length(unique(dbase$sex))/maxrows2/maxcols2) else 1  ## (RH 200825)
						for (ipage in 1:npages) {
							pagetext <- ""
							if (npages > 1) {
								pagetext <- paste("_page", ipage, sep="")
								caption <- paste(caption, " (plot ", ipage, " of ", npages, ")", sep="")
							}
							if (!is.null(ttcall(outnam)))
								file = paste0(sub("\\.png$","",outnam), f, letters[ipage], ".png")
							else
								file <- paste(filenamestart, filetype, filename_fltsexmkt, pagetext, ".png", sep="")
							plotinfo <- pngfun(file=file, caption=caption, lang=lang)
							tempfun3(ipage=ipage, l=lang, ...)
							dev.off(); eop()
						}
					}
					sampsizeline <- sampsizeline.old
					effNline <- effNline.old
				}
				if ((4 %in% subplots | 5 %in% subplots) & aalyear[1] > 0 & kind == "cond") {
					for (y in 1:length(aalyear)) {
						aalyr <- aalyear[y]
						if (length(dbase$Obs[dbase$Yr == aalyr]) > 0) {
							ydbase <- dbase[dbase$Yr == aalyr, ]
							sexvec <- ydbase$sex
							if (4 %in% subplots) {
								caption <- paste(aalyr, " age-at-length bin, ", title_sexmkt, fleetnames[f], sep="")
								if (mainTitle) {
									ptitle <- caption
								}
								else {
								ptitle <- ""
								}
								titles <- c(ptitle, titles)
								lenbinlegend <- paste(ydbase$Lbin_lo, labels[7], sep="")
								lenbinlegend[ydbase$Lbin_range > 0] <- paste(ydbase$Lbin_lo, "-", ydbase$Lbin_hi, labels[7], sep="")
								tempfun4 <- function(ipage, l="e", ...) {
									make.multifig(ptsx=ydbase$Bin, ptsy=ydbase$Obs, 
										yr=ydbase$Lbin_lo, linesx=ydbase$Bin, 
										linesy=ydbase$Exp, sampsize=ydbase$Nsamp_adj, 
										effN=ydbase$effN, showsampsize=showsampsize, 
										showeffN=showeffN, nlegends=3, 
										legtext=list(lenbinlegend, "sampsize", "effN"),
										bars=FALSE, linepos=linepos, 
										main=ptitle, cex.main=cex.main, 
										xlab=labels[2], ylab=labels[6], 
										maxrows=maxrows, maxcols=maxcols, 
										rows=rows, cols=cols, fixdims=fixdims, 
										ipage=ipage, scalebins=scalebins, 
										xlas=xlas, ylas=ylas, axis1=axis1, 
										axis2=axis2, axis1labs=axis1labs, 
										sexvec=sexvec, yupper=yupper, lang=l, ...)
								} ## end tempfun4
								if (plot) 
									tempfun4(ipage=0, l=lang, ...)
								if (print) {
									npages <- if (fixdims) ceiling(length(unique(ydbase$Yr.S))/maxrows/maxcols) else 1  ## (RH 200825)
									for (ipage in 1:npages) {
										pagetext <- ""
										if (npages > 1) {
											pagetext <- paste("_page", ipage, sep="")
											caption <- paste(caption, " (plot ", ipage, " of ", npages, ")", sep="")
										}
										if (length(grep("Pearson", caption)) > 0) {
											caption <- paste(caption, "<br> \nClosed bubbles are positive residuals", "(observed > expected)", "and open bubbles are negative residuals", "(observed < expected).")
										}
										if (!is.null(ttcall(outnam)))
											file = paste0(sub("\\.png$","",outnam), f, letters[ipage], ".png")
										else
											file <- paste0(filenamestart, filename_fltsexmkt, "_", aalyr, pagetext, ".png")
										plotinfo <- pngfun(file=file, caption=caption, lang=lang)
										tempfun4(ipage=ipage, l=lang, ...)
										dev.off(); eop()
									}
								}
							}
							if (5 %in% subplots) {
								z <- ydbase$Pearson
								col.index <- sexvec
								col.index[col.index == 0] <- 3
								cols <- colvec[col.index]
								x.vec <- ydbase$Bin + ydbase$sex * 1e-06
								caption <- paste(aalyr, " Pearson residuals for A-L key, ", title_sexmkt, fleetnames[f], sep="")
								caption <- paste(caption, " (max=", round(abs(max(z)), digits=2), ")", sep="")
								if (mainTitle) {
									ptitle <- caption
								}
								else {
									ptitle <- ""
								}
								titles <- c(ptitle, titles)
								tempfun5 <- function(l="e") {
									r4ss:::bubble3(x=x.vec, y=ydbase$Lbin_lo, z=z, xlab=linguaFranca(labels[2],l), ylab=linguaFranca(labels[1],l), col=cols, las=1, main=linguaFranca(ptitle,l), cex.main=cex.main, maxsize=pntscalar, cexZ1=cexZ1, legend=linguaFranca(bublegend,l), allopen=FALSE, minnbubble=minnbubble)
								} ## end tempfun5
								if (plot) 
									tempfun5(l=lang)
								if (print) {
									pagetext <- ""
									if (npages > 1) {
										pagetext <- paste("_page", ipage, sep="")
										caption <- paste(caption, " (plot ", ipage, " of ", npages, ")", sep="")
									}
									if (length(grep("Pearson", caption)) > 0) {
										caption <- paste(caption, "<br> \nClosed bubbles are positive residuals", "(observed > expected)", "and open bubbles are negative residuals", "(observed < expected).")
									}
									if (!is.null(ttcall(outnam)))
										file = paste0(sub("\\.png$","",outnam), f, ".png")
									else
										file <- paste0(filenamestart, "yearresids_", filename_fltsexmkt, "_", aalyr, pagetext, ".png")
									plotinfo <- pngfun(file=file, caption=caption, lang=lang)
									tempfun5(l=lang)
									dev.off(); eop()
								}
							}
						}
					}
				}
				if (6 %in% subplots & aalbin[1] > 0) {
					badbins <- setdiff(aalbin, dbase$Lbin_hi)
					goodbins <- intersect(aalbin, dbase$Lbin_hi)
					if (length(goodbins) > 0) {
						if (length(badbins) > 0) {
							cat("Error! the following inputs for 'aalbin' do not match the Lbin_hi values for the conditional age-at-length data:", badbins, "\n", "	   the following inputs for 'aalbin' are fine:", goodbins, "\n")
						}
						for (ibin in 1:length(goodbins)) {
							ilenbin <- goodbins[ibin]
							abindbase <- dbase[dbase$Lbin_hi == ilenbin, ]
							if (nrow(abindbase) > 0) {
								sexvec <- abindbase$sex
								caption <- paste0("Age-at-length ", ilenbin, labels[7], ", ", title_sexmkt, fleetnames[f])
								if (mainTitle) {
									ptitle <- caption
								}
								else {
									ptitle <- ""
								}
								titles <- c(ptitle, titles)
								tempfun6 <- function(ipage, l="e", ...) {
									make.multifig(ptsx=abindbase$Bin, 
										ptsy=abindbase$Obs, yr=abindbase$Yr.S, 
										linesx=abindbase$Bin, linesy=abindbase$Exp, 
										sampsize=abindbase$Nsamp_adj, effN=abindbase$effN, 
										showsampsize=showsampsize, showeffN=showeffN, 
										nlegends=3, legtext=list(abindbase$YrSeasName, 
										  "sampsize", "effN"), bars=bars, 
										linepos=(1 - datonly) * linepos, 
										main=ptitle, cex.main=cex.main, 
										xlab=kindlab, ylab=labels[6], 
										maxrows=maxrows, maxcols=maxcols, 
										rows=rows, cols=cols, fixdims=fixdims, 
										ipage=ipage, scalebins=scalebins, 
										sexvec=sexvec, lang=l, ...)
								} ## end tempfun6
								if (plot) 
									tempfun6(ipage=0, l=lang, ...)
								if (print) {
									npages <- if (fixdims) ceiling(length(unique(abindbase$Yr.S))/maxrows/maxcols) else 1  ## (RH 200825)
									for (ipage in 1:npages) {
										pagetext <- ""
										if (npages > 1) {
											pagetext <- paste0("_page", ipage)
											caption <- paste0(caption, " (plot ", ipage, " of ", npages, ")")
										}
										if (!is.null(ttcall(outnam)))
											file = paste0(sub("\\.png$","",outnam), f, letters[ipage], ".png")
										else
											file <- paste0(filenamestart, filename_fltsexmkt, "_length", ilenbin, labels[7], pagetext, ".png")
										plotinfo <- pngfun(file=file, caption=caption, lang=lang)
										tempfun6(ipage=ipage, l=lang, ...)
										dev.off(); eop()
									}
								}
							}
						}
					}
				}
				if (7 %in% subplots & samplesizeplots & !datonly & !("DM_effN" %in% names(dbase) && any(!is.na(dbase$DM_effN))) & !(kind %in% c("GSTAGE", "GSTLEN", "L@A", "W@A"))) {
					caption <- paste0("N-EffN comparison, ", titledata, title_sexmkt, fleetnames[f])
					if (mainTitle) {
						ptitle <- caption
					}
					else {
						ptitle <- ""
					}
					titles <- c(ptitle, titles)
					lfitfunc <- function(l="e") {
						if (kind == "cond") {
							dbasegood <- dbase[dbase$Obs >= 1e-04 & dbase$Exp < 0.99 & !is.na(dbase$effN) & dbase$effN < maxneff, ]
						}
						else {
							dbasegood <- dbase
						}
						if (nrow(dbasegood) > 0) {
							dbasegood2 <- dbasegood[, c("YrSeasName", "Nsamp_adj", "effN", "Nsamp_in")]
							dbasegood2 <- unique(dbasegood2)
							#plot(dbasegood2$Nsamp_adj, dbasegood2$effN, xlab=linguaFranca(labels[4],l), main=linguaFranca(ptitle,l), cex.main=cex.main, ylim=c(0, 1.15 * max(dbasegood2$effN)), xlim=c(0, 1.15 * max(dbasegood2$Nsamp_adj)), col=colvec[3], pch=19, ylab=linguaFranca(labels[5],l), xaxs="i", yaxs="i")
							plot(dbasegood2$Nsamp_adj, dbasegood2$effN, xlab=linguaFranca(labels[4],l), main=linguaFranca(ptitle,l), cex.main=cex.main, ylim=c(0, 1.04 * max(dbasegood2$effN)), xlim=c(0, 1.10 * max(dbasegood2$Nsamp_adj)), col="black", bg="yellow", pch=21, ylab=linguaFranca(labels[5],l), xaxs="i", yaxs="i")
							if (showyears) {
								par(xpd=TRUE)
								text(x=dbasegood2$Nsamp_adj, y=dbasegood2$effN, linguaFranca(dbasegood2$YrSeasName,l), adj=c(-0.2, 0.5))
								par(xpd=FALSE)
							}
							#abline(0, 1, col="black", lty=1)
							abline(0, 1, col="darkslategrey", lty=1, lwd=2)
							if (smooth & length(unique(dbasegood2$Nsamp_adj)) > 6 & diff(range(dbasegood2$Nsamp_adj)) > 2) {
								old_warn <- options()$warn
								options(warn=-1)
								psmooth <- loess(dbasegood2$effN ~ dbasegood2$Nsamp_adj, degree=1)
								options(warn=old_warn)
								lines(psmooth$x[order(psmooth$x)], psmooth$fit[order(psmooth$x)], lwd=2, col="red", lty="dashed")
								#lmfit = lm(dbasegood2$effN ~ dbasegood2$Nsamp_adj); abline(lmfit,col="blue") ## (RH 210709)
							}
							if (addMeans) {
								col.hmr = "blue"
								abline(v=mean(dbasegood2$Nsamp_adj), lty="22", col=col.hmr)
								#text(x=mean(dbasegood2$Nsamp_adj), y=0, col=col.hmr, linguaFranca("arithmetic mean",l), srt=90, adj=c(-0.1, -0.3))
								text(x=mean(dbasegood2$Nsamp_adj), y=max(dbasegood2$effN), col=col.hmr, linguaFranca("arithmetic mean",l), srt=90, adj=c(1, -0.3))
								abline(h=1/mean(1/dbasegood2$effN), lty="22", col=col.hmr)
								text(x=0, y=1/mean(1/dbasegood2$effN), col=col.hmr, linguaFranca("harmonic mean",l), adj=c(-0.1, -0.3))
#browser();return()
							}
						}
					} ## end lfitfunc
					if (plot) 
						lfitfunc(l=lang)
					if (print) {
						if (!is.null(ttcall(outnam)))
							file = paste0(sub("\\.png$","",outnam), f, ".png")
						else
							file <- paste(filenamestart, "sampsize_", filename_fltsexmkt, ".png", sep="")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						lfitfunc(l=lang)
						dev.off(); eop()
					}
				}
				if (8 %in% subplots & kind %in% c("LEN", "SIZE", "AGE")) {
					kind2 <- tolower(kind)
					if (plot) {
						tmp <- SSMethod.TA1.8(fit=replist, type=kind2, fleet=f, fleetnames=fleetnames, datonly=datonly, printit=verbose)
					}
					if (print) {
						## Needs frenchification
						file <- paste0(filenamestart, "data_weighting_TA1.8_", fleetnames[f], ".png")
						png(filename=file.path(plotdir, file), width=pwidth, height=pheight, units=punits, res=res, pointsize=ptsize)
						tmp <- SSMethod.TA1.8(fit=replist, type=kind2, fleet=f, fleetnames=fleetnames, datonly=datonly, printit=verbose)
						caption <- paste0("Mean ", gsub("len", "length", tolower(kind)), " for ", fleetnames[f], " with 95% confidence intervals", " based on current samples sizes.")
						if (!is.null(replist$Dirichlet_Multinomial_pars)) {
							caption <- paste("WARNING: this figure is based on multinomial likelihood", "and has not been updated to account for Dirichlet-Multinomial", "likelihood and the sample size adjustment associated with", "the estimated log(<i>&#920</i>) parameters.<br><br>", caption)
						}
						if (!datonly) {
							caption <- paste0(caption, "<br>Francis data weighting method TA1.8:")
							if (!is.null(tmp[1])) {
								vals <- paste0("thinner intervals (with capped ends) show ", "result of further adjusting sample sizes ", "based on suggested multiplier ", "(with 95% interval) for ",  kind2, " data from ", fleetnames[f], ":<br>", round(tmp[1], 4), " (", round(tmp[2], 4), "-", round(tmp[3], 4), ")")
							}
							else {
								vals <- "too few points to calculate adjustments."
							}
							caption <- paste(caption, vals, "<br><br>For more info, see<br>", "<blockquote>Francis, R.I.C.C. (2011).", "Data weighting in statistical fisheries stock assessment", "models. <i>Can. J. Fish. Aquat. Sci.</i>", "68: 1124-1138. ", "<a href=https://doi.org/10.1139/f2011-025>", "https://doi.org/10.1139/f2011-025</a>", "</blockquote>")
						}
						plotinfo <- rbind(plotinfo, data.frame(file=file, caption=caption))
						dev.off(); eop()
					}
				}
				if (9 %in% subplots & kind == "cond" & (f %in% condbase$Fleet)) {
					if (plot) {
						SSMethod.Cond.TA1.8(fit=replist, fleet=f, fleetnames=fleetnames, datonly=datonly)
					}
					if (print) {
						## Needs frenchification
						file <- paste(filenamestart, "data_weighting_TA1.8_condAge", fleetnames[f], ".png", sep="")
						png(filename=file.path(plotdir, file), width=pwidth, height=pheight, units=punits, res=res, pointsize=ptsize)
						tmp <- SSMethod.Cond.TA1.8(fit=replist, fleet=f, fleetnames=fleetnames, datonly=datonly)
						caption <- paste0("Mean age from conditional data", " (aggregated across length bins) for ", fleetnames[f], " with 95% confidence intervals ", " based on current samples sizes.")
						if (!datonly) {
							caption <- paste0(caption, "<br>Francis data weighting method TA1.8:")
							if (!is.null(tmp[1])) {
								vals <- paste0("thinner intervals (with capped ends) show ", "result of further adjusting sample sizes ", "based on suggested multiplier ", "(with 95% interval) for ", "conditional age-at-length data from ", fleetnames[f], ":<br>", round(tmp[1], 4), " (", round(tmp[2], 4), "-", round(tmp[3], 4), ")", sep="")
							}
							else {
								vals <- "too few points to calculate adjustments."
							}
							caption <- paste(caption, vals, "<br><br>For more info, see<br>", "<blockquote>Francis, R.I.C.C. (2011).", "Data weighting in statistical fisheries stock assessment", "models. <i>Can. J. Fish. Aquat. Sci.</i>", "68: 1124-1138.</blockquote>")
						}
						plotinfo <- rbind(plotinfo, data.frame(file=file, caption=caption))
						dev.off(); eop()
					}
				}
				if (10 %in% subplots & kind == "cond" & length(unique(dbase$Bin)) > 1) {
					caption1 <- paste(labels[14], title_sexmkt, fleetnames[f], sep="")
					if (mainTitle) {
						ptitle <- caption1
					}
					else {
						ptitle <- ""
					}
					andrefun <- function(ipage=0, l="e") {
						Lens <- sort(unique(dbase$Lbin_lo))
						Yrs <- sort(unique(dbase$Yr.S))
						ymax <- 1.1 * max(dbase$Bin, na.rm=TRUE)
						xmax <- max(condbase$Lbin_hi, na.rm=TRUE)
						xmin <- min(condbase$Lbin_lo, na.rm=TRUE)
						npanels <- length(Yrs)
						npages <- npanels/andrerows
						panelrange <- 1:npanels
						if (npages > 1 & ipage != 0) 
							panelrange <- intersect(panelrange, 1:andrerows + andrerows * (ipage - 1))
						Yrs2 <- Yrs[panelrange]
						par(mfrow=c(andrerows, 2), mar=c(2, 4, 1, 1), oma=andre_oma)
						for (Yr in Yrs2) {
							y <- dbase[dbase$Yr.S == Yr, ]
							Size <- NULL
							Size2 <- NULL
							Obs <- NULL
							Obs2 <- NULL
							Pred <- NULL
							Pred2 <- NULL
							Upp <- NULL
							Low <- NULL
							Upp2 <- NULL
							Low2 <- NULL
							for (Ilen in Lens) {
								z <- y[y$Lbin_lo == Ilen, ]
								if (length(z[, 1]) > 0) {
									weightsPred <- z$Exp/sum(z$Exp)
									weightsObs <- z$Obs/sum(z$Obs)
									ObsV <- sum(z$Bin * weightsObs)
									ObsV2 <- sum(z$Bin * z$Bin * weightsObs)
									PredV <- sum(z$Bin * weightsPred)
									PredV2 <- sum(z$Bin * z$Bin * weightsPred)
									NN <- z$Nsamp_adj[1]
									if (max(z$Obs, na.rm=TRUE) > 1e-04 & !is.na(NN) && NN > 0) {
										Size <- c(Size, Ilen)
										Obs <- c(Obs, ObsV)
										Pred <- c(Pred, PredV)
										varn <- sqrt(PredV2 - PredV * PredV)/sqrt(NN)
										Pred2 <- c(Pred2, varn)
										varn <- sqrt(max(0, ObsV2 - ObsV * ObsV, na.rm=TRUE))/sqrt(NN)
										Obs2 <- c(Obs2, varn)
										Low <- c(Low, ObsV - 1.64 * varn)
										Upp <- c(Upp, ObsV + 1.64 * varn)
										if (NN > 1) {
											Size2 <- c(Size2, Ilen)
											Low2 <- c(Low2, varn * sqrt((NN - 1)/qchisq(0.95, NN)))
										  Upp2 <- c(Upp2, varn * sqrt((NN - 1)/qchisq(0.05, NN)))
										}
									}
								}
							}
							if (length(Obs) > 0) {
								plot(Size, Obs, type="n", xlab="", ylab=linguaFranca("Age",l), xlim=c(xmin, xmax), ylim=c(0, ymax), yaxs="i")
								label <- ifelse(nseasons == 1, floor(Yr), Yr)
								text(x=par("usr")[1], y=0.9 * ymax, labels=linguaFranca(label,l), adj=c(-0.5, 0), font=2, cex=1.2)
								if (length(Low) > 1) 
									polygon(c(Size, rev(Size)), c(Low, rev(Upp)), col="grey95", border=NA)
								if (!datonly) 
									lines(Size, Pred, col=4, lwd=3)
								points(Size, Obs, pch=16)
								lines(Size, Low, lty=3)
								lines(Size, Upp, lty=3)
								if (par("mfg")[1] == 1) {
									title(main=linguaFranca(ptitle,l), xlab=linguaFranca(labels[1],l), outer=TRUE, line=1)
								}
								box()
								ymax2 <- max(Obs2, Pred2, na.rm=TRUE) * 1.1
								plot(Size, Obs2, type="n", xlab=linguaFranca(labels[1],l), ylab=linguaFranca(labels[13],l), xlim=c(xmin, xmax), ylim=c(0, ymax2), yaxs="i")
								if (length(Low2) > 1) 
									polygon(c(Size2, rev(Size2)), c(Low2, rev(Upp2)), col="grey95", border=NA)
								if (!datonly) 
									lines(Size, Pred2, col=4, lwd=3)
								points(Size, Obs2, pch=16)
								lines(Size2, Low2, lty=3)
								lines(Size2, Upp2, lty=3)
								if (!datonly & par("mfg")[1] == 1) {
									legend("topleft", legend=linguaFranca(c("Observed (with 90% interval)", "Expected"),l), bty="n", col=c(1, 4), pch=c(16, NA), lty=c(NA, 1), lwd=3)
								}
								box()
							}
						}
					} ## end andrefun
					if (plot) 
						andrefun(l=lang)
					if (print) {
						npages <- if (fixdims) ceiling(length(unique(dbase$Yr.S))/andrerows) else 1  ## (RH 200825)
						for (ipage in 1:npages) {
							pagetext <- ""
							caption <- caption1
							if (npages > 1) {
								pagetext <- paste("_page", ipage, sep="")
								caption <- paste(caption, " (plot ", ipage, " of ", npages, ")", sep="")
							}
							if (ipage == 1) {
								caption <- paste(caption, "\nThese plots show mean age and std. dev. in conditional A@L.<br>", "Left plots are mean A@L by size-class (obs. and pred.) ", "with 90% CIs based on adding 1.64 SE of mean to the data.<br>", "Right plots in each pair are SE of mean A@L (obs. and pred.) ", "with 90% CIs based on the chi-square distribution.", sep="")
							}
							if (!is.null(ttcall(outnam)))
								file = paste0(sub("\\.png$","",outnam), f, letters[ipage], ".png")
							else
								file <- paste(filenamestart, "Andre_plots", filename_fltsexmkt, pagetext, ".png", sep="")
							plotinfo <- pngfun(file=file, caption=caption, lang=lang)
							andrefun(ipage=ipage, l=lang)
							dev.off(); eop()
						}
					}
				}
			}
		}
	}
	if (21 %in% subplots & kind != "cond") {
		if (nrow(dbase_kind) > 0) {
			dbase_k <- dbase_kind
			for (j in unique(dbase_k$Part_group)) {
				dbase <- dbase_k[dbase_k$Part_group == j, ]
				if (nrow(dbase) > 0) {
					if (j == -1) 
						titlemkt <- ""
					if (j == 0) 
						titlemkt <- "whole catch, "
					if (j == 1) 
						titlemkt <- "discard, "
					if (j == 2) 
						titlemkt <- "retained, "
					titlemkt <- ifelse(printmkt, titlemkt, "")
					if (datonly | fitbar) {
						bars <- TRUE; polygons=FALSE
					}
					else {
						bars <- FALSE; polygons=TRUE
					}
					title_sexmkt <- paste(titlesex, titlemkt, sep="")
					filename_fltsexmkt <- paste(filesex)
					if (j > -1) {
						filename_fltsexmkt <- paste0(filename_fltsexmkt, "mkt", j)
					}
					caption <- paste(titledata, title_sexmkt, "aggregated across time by fleet", sep="")
					if (mainTitle) {
						ptitle <- caption
					}
					else {
						ptitle <- ""
					}
					titles <- c(ptitle, titles)
					Bins <- sort(unique(dbase$Bin))
					nbins <- length(Bins)
					df <- data.frame(Nsamp_adj=dbase$Nsamp_adj, effN=dbase$effN, obs=dbase$Obs * dbase$Nsamp_adj, exp=dbase$Exp * dbase$Nsamp_adj)
					if ("DM_effN" %in% names(dbase) && any(!is.na(dbase$DM_effN))) {
						df$DM_effN <- dbase$DM_effN
					}
					agg <- aggregate(x=df, by=list(bin=dbase$Bin, f=dbase$Fleet, sex=dbase$sex, mkt=dbase$Part), FUN=sum)
					agg <- agg[agg$f %in% fleets, ]
					agg$obs <- agg$obs/agg$Nsamp_adj
					agg$exp <- agg$exp/agg$Nsamp_adj
					for (f in unique(agg$f)) {
						for (mkt in unique(agg$mkt[agg$f == f])) {
							sub <- agg$f == f & agg$mkt == mkt
							agg$Nsamp_adj[sub] <- max(agg$Nsamp_adj[sub])
							if ("DM_effN" %in% names(agg) && any(!is.na(agg$DM_effN))) {
								agg$DM_effN[sub] <- max(agg$DM_effN[sub], na.rm=TRUE)
							}
							else {
								if (any(!is.na(agg$effN[sub]))) {
									agg$effN[sub] <- max(agg$effN[sub], na.rm=TRUE)
								}
								else {
									agg$effN[sub] <- NA
								}
							}
						}
					}
					namesvec <- fleetnames[agg$f]
					max_n_mkt <- max(apply(table(agg$f, agg$mkt) > 0, 1, sum))
					if (max_n_mkt > 0) {
						mktnames <- c("", "(discards)", "(retained)")
						namesvec <- paste(fleetnames[agg$f], mktnames[agg$mkt + 1])
					}
					if (!(kind %in% c("GSTAGE", "GSTLEN", "L@A", "W@A"))) {
						tempfun7 <- function(ipage, l="e", ...) {
							if ("DM_effN" %in% names(agg) && any(!is.na(agg$DM_effN))) {
								make.multifig(ptsx=agg$bin, ptsy=agg$obs, 
									yr=paste(agg$f, agg$mkt), linesx=agg$bin, 
									linesy=agg$exp, sampsize=agg$Nsamp_adj, 
									effN=agg$DM_effN, showsampsize=showsampsize, 
									showeffN=showeffN, sampsize_label="Sum of N input=", 
									effN_label="Sum of N samp.=", bars=bars, 
									linepos=(1 - datonly) * linepos, 
									nlegends=3, legtext=list(namesvec, "sampsize", "effN"),
									main=ptitle, cex.main=cex.main, xlab=kindlab, 
									ylab=labels[6], maxrows=maxrows, 
									maxcols=maxcols, rows=rows, cols=cols, 
									fixdims=fixdims2, ipage=ipage, 
									scalebins=scalebins, colvec=colvec, 
									linescol=linescol, xlas=xlas, ylas=ylas, 
									axis1=axis1, axis2=axis2, axis1labs=axis1labs, 
									sexvec=agg$sex, yupper=yupper, lang=l, ...)
							}
							else {
								make.multifig(ptsx=agg$bin, ptsy=agg$obs, 
									yr=paste(agg$f, agg$mkt), linesx=agg$bin, 
									linesy=agg$exp, sampsize=agg$Nsamp_adj, 
									effN=agg$effN, showsampsize=showsampsize, 
									showeffN=showeffN, sampsize_label="Sum of N samp.=", 
									effN_label="Sum of N eff.=", bars=bars, 
									linepos=(1 - datonly) * linepos, 
									nlegends=3, legtext=list(namesvec, "sampsize", "effN"),
									main=ptitle, cex.main=cex.main, xlab=kindlab, 
									ylab=labels[6], maxrows=maxrows, 
									maxcols=maxcols, rows=rows, cols=cols, 
									fixdims=fixdims2, ipage=ipage, 
									scalebins=scalebins, colvec=colvec, 
									linescol=linescol, xlas=xlas, ylas=ylas, 
									axis1=axis1, axis2=axis2, axis1labs=axis1labs, 
									sexvec=agg$sex, yupper=yupper, lang=l, ...)
							}
						} ## end tempfun7
						if (plot) 
							tempfun7(ipage=0, l=lang, ...)
						if (print) {
							npages <- if (fixdims) ceiling(length(unique(agg$f))/maxrows/maxcols) else 1  ## (RH 200825)
							for (ipage in 1:npages) {
								if (max_n_mkt > 0) {
									caption <- paste0(caption, ".\n <br> ", "Labels 'retained' and 'discard' indicate", " discarded or retained sampled for each fleet.", " Panels without this designation represent", " the whole catch.\n")
								}
								pagetext <- ""
								if (npages > 1) {
									pagetext <- paste("_page", ipage, sep="")
									caption <- paste(caption, "<br> (plot ", ipage, " of ", npages, ")", sep="")
								}
								if (!is.null(ttcall(outnam)))
									file = paste0(sub("\\.png$","",outnam), f, letters[ipage], ".png")
								else
									file <- paste(filenamestart, filename_fltsexmkt, pagetext, "_aggregated_across_time.png", sep="")
								plotinfo <- pngfun(file=file, caption=caption, lang=lang)
								tempfun7(ipage=ipage, l=lang, ...)
								dev.off(); eop()
							}
						}
					}
					else {
					}
				}
			}
		}
	}
	if (22 %in% subplots & kind != "cond" & nseasons > 1) {
		dbasef <- dbase_kind[dbase_kind$Fleet %in% fleets, ]
		if ("DM_effN" %in% names(dbasef) && any(!is.na(dbasef$DM_effN))) {
			warning("Sample sizes in plots by fleet aggregating across years within each season have not yet been updated to reflect Dirichlet-Multinomial likelihood")
		}
		if (nrow(dbasef) > 0) {
			testor <- length(dbasef$sex[dbasef$sex == 1 & dbasef$Pick_sex == 0]) > 0
			testor[2] <- length(dbasef$sex[dbasef$sex == 1 & dbasef$Pick_sex %in% c(1, 3)]) > 0
			testor[3] <- length(dbasef$sex[dbasef$sex == 2]) > 0
			for (k in (1:3)[testor]) {
				if (k == 1) {
					dbase_k <- dbasef[dbasef$sex == 1 & dbasef$Pick_sex == 0, ]
				}
				if (k == 2) {
					dbase_k <- dbasef[dbasef$sex == 1 & dbasef$Pick_sex %in% c(1, 3), ]
				}
				if (k == 3) {
					dbase_k <- dbasef[dbasef$sex == 2, ]
				}
				sex <- ifelse(k == 3, 2, 1)
				for (j in unique(dbase_k$Part)) {
					dbase <- dbase_k[dbase_k$Part == j, ]
					if (nrow(dbase) > 0) {
						if (k == 1) 
							titlesex <- "sexes combined, "
						if (k == 2) 
							titlesex <- "female, "
						if (k == 3) 
							titlesex <- "male, "
						titlesex <- ifelse(printsex, titlesex, "")
						if (j == 0) 
							titlemkt <- "whole catch, "
						if (j == 1) 
							titlemkt <- "discard, "
						if (j == 2) 
							titlemkt <- "retained, "
						titlemkt <- ifelse(printmkt, titlemkt, "")
						if (datonly | fitbar) {
							bars <- TRUE; polygons=FALSE
						}
						else {
							bars <- FALSE; polygons=TRUE
						}
						title_sexmkt <- paste(titlesex, titlemkt, sep="")
						filename_fltsexmkt <- paste("sex", k, "mkt", j, sep="")
						caption <- paste0(titledata, title_sexmkt, "\naggregated within season by fleet")
						if (mainTitle) {
							ptitle <- caption
						}
						else {
							ptitle <- ""
						}
						titles <- c(ptitle, titles)
						Bins <- sort(unique(dbase$Bin))
						nbins <- length(Bins)
						df <- data.frame(Nsamp_adj=dbase$Nsamp_adj, effN=dbase$effN, obs=dbase$Obs * dbase$Nsamp_adj, exp=dbase$Exp * dbase$Nsamp_adj)
						agg <- aggregate(x=df, by=list(bin=dbase$Bin, f=dbase$Fleet, s=dbase$Seas), FUN=sum)
						agg <- agg[agg$f %in% fleets, ]
						if (any(agg$s <= 0)) {
							cat("super-periods may not work correctly in plots of aggregated comps\n")
							agg <- agg[agg$s > 0, ]
						}
						agg$obs <- agg$obs/agg$Nsamp_adj
						agg$exp <- agg$exp/agg$Nsamp_adj
						for (f in unique(agg$f)) {
							for (s in unique(agg$s[agg$f == f])) {
								infleetseas <- agg$f == f & agg$s == s
								agg$Nsamp_adj[infleetseas] <- max(agg$Nsamp_adj[infleetseas])
								agg$effN[infleetseas] <- max(agg$effN[infleetseas])
							}
						}
						agg$fseas <- agg$f + seasfracs[agg$s]
						namesvec <- paste(fleetnames[agg$f], " s", agg$s, sep="")
						tempfun8 <- function(ipage, l="e", ...) {
							if (!(kind %in% c("GSTAGE", "GSTLEN", "L@A", "W@A"))) {
								make.multifig(ptsx=agg$bin, ptsy=agg$obs, 
									yr=agg$fseas, linesx=agg$bin, linesy=agg$exp, 
									sampsize=agg$Nsamp_adj, effN=agg$effN, 
									showsampsize=showsampsize, showeffN=showeffN, 
									bars=bars, linepos=(1 - datonly) * linepos,
									nlegends=3, legtext=list(namesvec, "sampsize", "effN"),
									main=ptitle, cex.main=cex.main, xlab=kindlab, 
									ylab=labels[6], maxrows=maxrows, 
									maxcols=maxcols, rows=rows, cols=cols, 
									fixdims=fixdims2, ipage=ipage, 
									scalebins=scalebins, colvec=colvec, 
									linescol=linescol, xlas=xlas, ylas=ylas, 
									axis1=axis1, axis2=axis2, axis1labs=axis1labs, 
									sexvec=agg$sex, yupper=yupper, lang=l, ...)
							}
						} ## end tempfun8
						if (plot) 
							tempfun8(ipage=0, l=lang, ...)
						if (print) {
							npages <- if (fixdims) ceiling(length(unique(agg$fseas))/maxrows/maxcols) else 1  ## (RH 200825)
							for (ipage in 1:npages) {
								pagetext <- ""
								if (npages > 1) {
									pagetext <- paste("_page", ipage, sep="")
									caption <- paste(caption, " (plot ", ipage, " of ", npages, ")", sep="")
								}
								if (!is.null(ttcall(outnam)))
									file = paste0(sub("\\.png$","",outnam), f, letters[ipage], ".png")
								else
									file <- paste(filenamestart, filename_fltsexmkt, pagetext, "_aggregated_within_season.png", sep="")
								plotinfo <- pngfun(file=file, caption=caption, lang=lang)
								tempfun8(ipage=ipage, l=lang, ...)
								dev.off(); eop()
							}
						}
					}
				}
			}
		}
	}
	if (23 %in% subplots & kind != "cond" & nseasons > 1) {
		for (f in fleets) {
			dbasef <- dbase_kind[dbase_kind$Fleet == f, ]
			if ("DM_effN" %in% names(dbasef) && any(!is.na(dbasef$DM_effN))) {
				warning("Sample sizes in plots by fleet aggregating across seasons within a year have not yet been updated to reflect Dirichlet-Multinomial likelihood")
			}
			if (nrow(dbasef) > 0) {
				testor <- length(dbasef$sex[dbasef$sex == 1 & dbasef$Pick_sex == 0]) > 0
				testor[2] <- length(dbasef$sex[dbasef$sex == 1 & dbasef$Pick_sex %in% c(1, 3)]) > 0
				testor[3] <- length(dbasef$sex[dbasef$sex == 2]) > 0
				for (k in (1:3)[testor]) {
					if (k == 1) {
						dbase_k <- dbasef[dbasef$sex == 1 & dbasef$Pick_sex == 0, ]
					}
					if (k == 2) {
						dbase_k <- dbasef[dbasef$sex == 1 & dbasef$Pick_sex %in% c(1, 3), ]
					}
					if (k == 3) {
						dbase_k <- dbasef[dbasef$sex == 2, ]
					}
					sex <- ifelse(k == 3, 2, 1)
					for (j in unique(dbase_k$Part)) {
						dbase <- dbase_k[dbase_k$Part == j, ]
						if (nrow(dbase) > 0) {
							if (k == 1) 
								titlesex <- "sexes combined, "
							if (k == 2) 
								titlesex <- "female, "
							if (k == 3) 
								titlesex <- "male, "
							titlesex <- ifelse(printsex, titlesex, "")
							if (j == 0) 
								titlemkt <- "whole catch, "
							if (j == 1) 
								titlemkt <- "discard, "
							if (j == 2) 
								titlemkt <- "retained, "
							titlemkt <- ifelse(printmkt, titlemkt, "")
							if (datonly | fitbar) {
								bars <- TRUE; polgons=FALSE
							}
							else {
								bars <- FALSE; polgons=TRUE
							}
							title_sexmkt <- paste(titlesex, titlemkt, sep="")
							filename_fltsexmkt <- paste("flt", f, "sex", k, "mkt", j, sep="")
							Bins <- sort(unique(dbase$Bin))
							nbins <- length(Bins)
							df <- data.frame(Nsamp_adj=dbase$Nsamp_adj, effN=dbase$effN, obs=dbase$Obs * dbase$Nsamp_adj, exp=dbase$Exp * dbase$Nsamp_adj)
							agg <- aggregate(x=df, by=list(bin=dbase$Bin, f=dbase$Fleet, y=floor(dbase$Yr.S)), FUN=sum)
							agg <- agg[agg$f %in% fleets, ]
							agg$obs <- agg$obs/agg$Nsamp_adj
							agg$exp <- agg$exp/agg$Nsamp_adj
							for (f in unique(agg$f)) {
								for (y in unique(agg$y[agg$f == f])) {
									infleetyr <- agg$f == f & agg$y == y
									agg$Nsamp_adj[infleetyr] <- max(agg$Nsamp_adj[infleetyr])
									agg$effN[infleetyr] <- max(agg$effN[infleetyr])
								}
							}
							agg$fy <- agg$f + agg$y/10000
							caption <- paste(titledata, title_sexmkt, fleetnames[f], "\naggregated across seasons within year", sep="")
							if (mainTitle) {
								ptitle <- caption
							}
							else {
								ptitle <- ""
							}
							tempfun9 <- function(ipage, l="e", ...) {
								if (!(kind %in% c("GSTAGE", "GSTLEN", "L@A", "W@A"))) {
									make.multifig(ptsx=agg$bin, ptsy=agg$obs, 
										yr=agg$fy, linesx=agg$bin, linesy=agg$exp, 
										sampsize=agg$Nsamp_adj, effN=agg$effN, 
										showsampsize=showsampsize, showeffN=showeffN, 
										bars=bars, linepos=(1-datonly) * linepos,
										nlegends=3, legtext=list(agg$y, "sampsize", "effN"),
										main=ptitle, cex.main=cex.main, xlab=kindlab, 
										ylab=labels[6], maxrows=maxrows, 
										maxcols=maxcols, rows=rows, cols=cols, 
										fixdims=fixdims2, ipage=ipage, 
										scalebins=scalebins, colvec=colvec, 
										linescol=linescol, xlas=xlas, 
										ylas=ylas, axis1=axis1, axis2=axis2, 
										axis1labs=axis1labs, sexvec=agg$sex, 
										yupper=yupper, lang=l, ...)
								}
							} ## end tempfun9
							if (plot) 
								tempfun9(ipage=0, l=lang, ...)
							if (print) {
								npages <- if (fixdims) ceiling(length(unique(agg$fy))/maxrows/maxcols) else 1  ## (RH 200825)
								for (ipage in 1:npages) {
									pagetext <- ""
									if (npages > 1) {
										pagetext <- paste("_page", ipage, sep="")
										caption <- paste(caption, " (plot ", ipage, " of ", npages, ")", sep="")
									}
									if (!is.null(ttcall(outnam)))
										file = paste0(sub("\\.png$","",outnam), f, letters[ipage], ".png")
									else
										file <- paste(filenamestart, filename_fltsexmkt, pagetext, "_aggregated_across_seasons_within_year.png", sep="")
									pngfun(file=file, caption=caption, lang=lang)
									tempfun9(ipage=ipage, l=lang, ...)
									dev.off(); eop()
								}
							}
						}
					}
				}
			}
		}
	}
	if (24 %in% subplots & kind %in% c("LEN", "AGE")) {
		for (j in unique(dbase_kind$Part_group)) {
			dbase_parts <- dbase_kind[dbase_kind$Part_group == j, ]
			dbase_parts$FleetPart <- dbase_parts$Fleet + 0.1 * dbase_parts$Part
			panel_table <- data.frame(FleetPart=sort(unique(dbase_parts$FleetPart)))
			panel_table$Fleet <- floor(panel_table$FleetPart)
			panel_table$Part <- round(10 * (panel_table$FleetPart - panel_table$Fleet))
			panel_table$Name <- fleetnames[panel_table$Fleet]
			max_n_mkt <- max(apply(table(panel_table$Fleet, panel_table$Part) > 0, 1, sum))
			if (max_n_mkt > 1) {
				mktnames <- c("", "(discards)", "(retained)")
				panel_table$Name <- paste(panel_table$Name, mktnames[panel_table$Part + 1])
			}
			npanels <- nrow(panel_table)
			panelvec <- 1:npanels
			xlim <- range(dbase_parts$Yr.S)
			xaxislab <- sort(unique(floor(dbase_parts$Yr.S)))
			if (length(cohortlines) > 0) {
				growdat <- replist$endgrowth
				growdatF <- growdat[growdat$Sex == 1 & growdat$Morph == min(growdat$Morph[growdat$Sex == 1]), ]
				if (nsexes > 1) {
					growdatM <- growdat[growdat$Sex == 2 & growdat$Morph == min(growdat$Morph[growdat$Sex == 2]), ]
				}
			}
			if (j == -1)
				titlemkt <- ""
			if (j == 0)
				titlemkt <- "whole catch"
			if (j == 1)
				titlemkt <- "discard"
			if (j == 2)
				titlemkt <- "retained"
			titlemkt <- ifelse(printmkt, titlemkt, "")
			caption_base <- paste0(titletype, titlemkt, ", comparing across fleets")
			caption_base <- gsub(", ,", ", ", caption_base)
			if (mainTitle) {
				ptitle <- caption_base
			}
			else {
				ptitle <- ""
			}
			titles <- c(ptitle, titles)
			filenamemkt <- ifelse(j > -1, paste("mkt", j, sep=""), "")
			multifleet.bubble.fun <- function(ipage=0, l="e") {
				par_old <- par()
				par(mfrow=c(min(npanels, maxrows), 1), mar=c(0.5, 0, 0, 0), oma=c(4, 6, ifelse(mainTitle, 3, 1), 1))
				panelrange <- 1:npanels
				npages <- ceiling(npanels/maxrows)
				if (npages > 1 & ipage != 0)
					panelrange <- intersect(panelrange, 1:maxrows + maxrows * (ipage - 1))
				for (ipanel in panelvec[panelrange]) {
					flt <- panel_table$Fleet[ipanel]
					mkt <- panel_table$Part[ipanel]
					dbase <- dbase_parts[dbase_parts$Fleet == flt & dbase_parts$Part == mkt, ]
					max_n_ageerr <- max(apply(table(dbase$Yr.S, dbase$Ageerr) > 0, 1, sum))
					if (max_n_ageerr > 1) {
						if (ageerr_warning) {
							cat("Note: multiple samples with different ageing error types within fleet/year.\n", "   Plots label '2005a3' indicates ageing error type 3 for 2005 sample.\n", "   Bubble plots may be misleading with overlapping bubbles.\n")
							ageerr_warning <- FALSE
						}
						dbase$Yr.S <- dbase$Yr.S + dbase$Ageerr/(1000 * max_n_ageerr)
						dbase$YrSeasName <- paste(dbase$YrSeasName, "a", dbase$Ageerr, sep="")
					}
					xdiff <- 0.1 * sort(unique(diff(sort(unique(dbase$Yr.S)), na.rm=TRUE)))[1]
					if (is.na(xdiff)) {
						xdiff <- 0.1
					}
					xvals <- dbase$Yr.S
					cols <- rep(colvec[3], nrow(dbase))
					if (nsexes > 1) {
						xvals[dbase$sex > 0] <- dbase$Yr.S[dbase$sex > 0] - (dbase$sex[dbase$sex > 0] - 1.5) * xdiff
						if (length(unique(dbase$Yr.S)) == 1) {
							xvals[dbase$sex > 0] <- floor(dbase$Yr.S[dbase$sex > 0]) - (dbase$sex[dbase$sex > 0] - 1.5) * xdiff
						}
						cols[dbase$sex > 0] <- colvec[dbase$sex[dbase$sex > 0]]
					}
					if (datonly) {
						z <- dbase$Obs
						if (scalebubbles) 
							z <- dbase$Nsamp_adj * dbase$Obs
						titletype <- titledata
						filetype <- "bub"
						allopen <- TRUE
					}
					else {
						z <- dbase$Pearson
						titletype <- "Pearson residuals, "
						filetype <- "resids"
						allopen <- FALSE
					}
					ylim <- range(dbase$Bin)
					ylim[2] <- ylim[2] + 0.2 * diff(ylim)
					r4ss:::bubble3(x=xvals, y=dbase$Bin, z=z, col=cols, cexZ1=cexZ1, legend=linguaFranca(bublegend,l), las=1, main="", cex.main=cex.main, maxsize=pntscalar, allopen=allopen, xlim=xlim, ylim=ylim, axis1=FALSE)
					legend("topleft", title=linguaFranca(panel_table$Name[ipanel],l), legend=NA, bty="n", cex=1.5)
					if (length(cohortlines) > 0) {
						for (icohort in 1:length(cohortlines)) {
							cat("  Adding line for", cohortlines[icohort], "cohort\n")
							if (kind == "LEN") {
								lines(growdatF$Age + cohortlines[icohort], growdatF$Len_Mid, col=colvec[1])
								if (nsexes > 1) {
									lines(growdatM$Age + cohortlines[icohort], growdatM$Len_Mid, col=colvec[2])
								}
							}
							if (kind == "AGE") {
								lines(0.5 + c(cohortlines[icohort], cohortlines[icohort] + accuage), c(0, accuage), col="red")
							}
						}
					}
					if (par()$mfg[1] == par()$mfg[3] | ipanel == tail(panelvec, 1)) {
						axis(1, at=xaxislab)
					}
					else {
						axis(1, at=xaxislab, labels=rep("", length(xaxislab)))
					}
					if (par()$mfg[1] == 1) 
						title(main=ptitle, outer=TRUE, xlab=linguaFranca(labels[3],l), ylab=linguaFranca(kindlab,l))
				}
				par(mfcol=par_old$mfcol, mar=par_old$mar, oma=par_old$oma)
			} ## end multifleet.bubble.fun
			if (length(panelvec) > 0) {
				if (plot) 
					multifleet.bubble.fun(ipage=0, l=lang)
				if (print) {
					npages <- if (fixdims) ceiling(nrow(panel_table)/maxrows) else 1  ## (RH 200825)
					for (ipage in 1:npages) {
						pagetext <- ""
						caption <- caption_base
						if (npages > 1) {
							pagetext <- paste("_page", ipage, sep="")
							caption <- paste0(caption, " (plot ", ipage, " of ", npages, ")")
						}
						if (ipage == 1 & length(grep("Pearson", caption)) > 0) {
							caption <- paste(caption, "<br> \nClosed bubbles are positive residuals", "(observed > expected)", "and open bubbles are negative residuals", "(observed < expected).")
						}
						if (!is.null(ttcall(outnam)))
							file = paste0(sub("\\.png$","",outnam), letters[ipage], ".png")
						else
							file <- paste(filenamestart, filenamemkt, pagetext, "_multi-fleet_comparison.png", sep="")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						multifleet.bubble.fun(ipage=ipage, l=lang)
						dev.off(); eop()
					}
				}
			}
		}
		par(mfcol=c(rows, cols), mar=c(5, 4, 4, 2) + 0.1, oma=rep(0, 4))
	}
	if (!is.null(plotinfo))
		plotinfo$category <- "Comp"
	return(invisible(plotinfo))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.comps


## plotSS.francis-----------------------2021-03-10
##  Plot mean age fits using Francis (2011) methodology
##  Modified function 'r4ss::SSMethod.TA1.8' :
##  Apply Francis composition weighting method TA1.8
## ----------------------------------------r4ss|RH
##  Uses method TA1.8 (described in Appendix A of Francis 2011) to do
##  stage-2 weighting of composition data from a Stock Synthesis model.
##  Outputs a multiplier, \emph{w} (with bootstrap 95\% confidence interval),
##  so that \emph{N2y} = \emph{w} x \emph{N1y}, where \emph{N1y} and
##  \emph{N2y} are the stage-1 and stage-2
##  multinomial sample sizes for the data set in year y.  Optionally
##  makes a plot of observed (with confidence limits, based on \emph{N1y})
##  and expected mean lengths (or ages).
##  \cr\cr
##  CAUTIONARY/EXPLANATORY NOTE.
##  The large number of options available in SS makes it very
##  difficult to be sure that what this function does is
##  appropriate for all combinations of options.  The following
##  notes might help anyone wanting to check or correct the code.
##  \enumerate{
##    \item The code first takes the appropriate database (lendbase, sizedbase,
##          agedbase, or condbase) and removes un-needed rows.
##    \item The remaining rows of the database are grouped into individual
##          comps (indexed by vector indx) and relevant statistics (e.g.,
##          observed and expected mean length or age), and ancillary data,
##          are calculated for each comp (these are stored in pldat - one row
##          per comp).
##          If the data are to be plotted, the comps are grouped, with each
##          group corresponding to a panel in the plot, and groups are indexed
##          by plindx.
##    \item A single multiplier is calculated to apply to all the comps.
##  }
#'
##  @param fit Stock Synthesis output as read by r4SS function SS_output
##  @param type either 'len' (for length composition data), 'size' (for
##  generalized size composition data), 'age' (for age composition data),
##  or 'con' (for conditional age at length data)
##  @param fleet vector of one or more fleet numbers whose data are to
##  be analysed simultaneously (the output N multiplier applies
##  to all fleets combined)
##  @param fleetnames Vector of alternative fleet names to draw from for
##  plot titles and captions. It should have length equal to the number
##  of fleets in the model, not the number of fleets considered in this function.
##  @param part vector of one or more partition values; analysis is restricted
##  to composition data with one of these partition values.
##  Default is to include all partition values (0, 1, 2).
##  @param label.part Include labels indicating which partitions are included?
##  @param sexes vector of one or more values for Sexes; analysis is
##  restricted to composition data with one of these
##  Sexes values.  Ignored if type=='con'.
##  @param label.sex Include labels indicating which sexes are included?
##  @param seas string indicating how to treat data from multiple seasons
##  'comb' - combine seasonal data for each year and plot against Yr
##  'sep' - treat seasons separately, plotting against Yr.S
##  If is.null(seas) it is assumed that there is only one season in
##  the selected data (a warning is output if this is not true) and
##  option 'comb' is used.
##  @param method a vector of one or more size-frequency method numbers
##  (ignored unless type = 'size').
##  If !is.null(method), analysis is restricted to size-frequency
##  methods in this vector.  NB comps are separated by method
##  @param plotit if TRUE, make an illustrative plot like one or more
##  panels of Fig. 4 in Francis (2011).
##  @param printit if TRUE, print results to R console.
##  @param datonly if TRUE, don't show the model expectations
##  @param plotadj if TRUE, plot the confidence intervals associated with
##  the adjusted sample sizes (TRUE by default unless datonly = TRUE)
##  @param maxpanel maximum number of panels within a plot
##  @param set.pars Set the graphical parameters such as mar and mfrow.
##  Can be set to FALSE in order to add plots form multiple calls to
##  this function as separate panels in one larger figure.
##  @author Chris Francis, Andre Punt, Ian Taylor
##  @export
##  @seealso \code{\link{SSMethod.Cond.TA1.8}}
##  @references Francis, R.I.C.C. (2011). Data weighting in statistical
##  fisheries stock assessment models. Canadian Journal of
##  Fisheries and Aquatic Sciences 68: 1124-1138.
## -------------------------------------------r4ss
plotSS.francis <- function(
   fit, type, fleet, part=0:2, sexes=0:3, seas=NULL, method=NULL, 
   plotit=TRUE, printit=TRUE, datonly=FALSE, plotadj=!datonly,
   maxpanel=1000, fleetnames=NULL, label.part=TRUE, label.sex=TRUE, 
   set.pars=TRUE, col.obs=c("green3","green"), col.fit=lucent(c("blue2","cyan"),0.5),
   png=FALSE, pngres=400, PIN=c(8,9), outnam, lang=c("e","f") )
{
	oldpar = par(no.readonly=TRUE)
	fart = function(opar) { if (any("windows"%in%names(dev.list()))) par(opar) }
	on.exit(fart(oldpar))

	## Check the type is correct and the sexes is correct
	is.in <- function (x, y)!is.na(match(x, y))
	if(!is.in(type,c('age','len','size','con'))){
		stop('Illegal value for type (should be "age", "len", "size", or "con")')
	} else {
		if(sum(!is.in(sexes,c(0:3)))>0){
			stop('Unrecognised value for sexes')
		}
	}
	## Replace default fleetnames with user input if requested
	if(is.null(fleetnames)){
		# use fleetnames in the model
		fleetnames <- fit$FleetNames
	} else {
		# if custom names input, check length
		if(length(fleetnames) != fit$nfleets){
		stop('fleetnames needs to be NULL or have length = nfleets = ', fit$nfleets)
		}
	}
	## Select the type of datbase
	dbase <- fit[[paste(type,'dbase',sep='')]]
	## sel is vector of row indices selected for the plot/calculations
	## select row indices matching fleet and partition
	sel <- is.in(dbase$Fleet,fleet) & is.in(dbase$Part,part)
	## select row indices matching Sexes column
	if(type!='con'){
		## change column name on earlier SS versions to match change from
		## Pick_sex to Sexes in 3.30.12 (July 2018)
		names(dbase)[names(dbase)=='Pick_sex'] <- 'Sexes'
		#sel <- sel & is.in(dbase$'Sexes',sexes)
		sel <- sel & is.in(dbase$'Sex',sexes)
	}
	## For generalized size frequency comps, select chosen size method
	if(type=='size' & !is.null(method)){
		sel <- sel & is.in(dbase$method,method)
	}
	# If there are no rows selected, return empty result
	if(sum(sel)==0) return()
	## Subset comp database for selected rows
	dbase <- dbase[sel,]
	if(is.null(seas)){
		seas <- 'comb'
		if(length(unique(dbase$Seas))>1)
			cat('Warning: combining data from multiple seasons\n')
	}
	## if generalized size comp is used, check for mix of units
	if(type=='size'){
		if(length(unique(dbase$units))>1){
			cat('Warning: mix of units being compared:',unique(dbase$units),'\n')
		}
	}
	## create label for partitions
	partitions <- sort(unique(dbase$Part)) # values are 0, 1, or 2
	partition.labels <- c("whole","discarded","retained")[partitions+1]
	partition.labels <- paste0("(",paste(partition.labels,collapse="&")," catch)")
	## indx is string combining fleet, year, and potentially conditional bin
	indx <- paste(dbase$Fleet, dbase$Yr, if(type=='con') dbase$'Lbin_lo' else '', if(seas=='sep') dbase$Seas else '')
	## if subsetting by sex, add Sexes value to the indx strings
	sex.flag <- type!='con' & max(tapply(dbase$'Sexes', dbase$Fleet, function(x)length(unique(x))))>1
	if(sex.flag){
		indx <- paste(indx,dbase$'Sexes')
	}
	## if subsetting by generalized size-method, add that value to indx strings
	method.flag <- type=='size' && length(unique(dbase$method))>1
	if(method.flag){
		indx <- paste(indx,dbase$method)
	}
	## unique strings in indx vector
	uindx <- unique(indx)
	## test for length 1 results
	if(length(uindx)==1){
		## presumably the method is meaningless of there's only 1 point,
		## but it's good to be able to have the function play through
		cat('Warning: only one point to plot\n')
		return()
	}
	## create empty data.frame to store information on each observation
	pldat <- matrix(0,length(uindx),12, dimnames=list(uindx, c('Obsmn','Obslo','Obshi','semn','Expmn','Vexp','N','Std.res','ObsloAdj','ObshiAdj','Fleet','Yr')))
	## add columns of zeros to fill with values necessary for subsetting
	if(type=='con') pldat <- cbind(pldat,Lbin=0)
	if(sex.flag)    pldat <- cbind(pldat,sexes=0)
	if(type=='size'){
		pldat <- cbind(pldat,method=0)
		## vector to store units (which are strings and don't fit in pldat matrix)
		plunits <- rep(NA,nrow(pldat)) 
	}
	## Find the weighting factor for this combination of factors
	for(i in 1:length(uindx)){ ## each row of pldat is an individual comp
		subdbase <- dbase[indx==uindx[i],]
		xvar <- subdbase$Bin
		pldat[i,'Obsmn']   <- sum(subdbase$Obs*xvar)/sum(subdbase$Obs)
		pldat[i,'Expmn']   <- sum(subdbase$Exp*xvar)/sum(subdbase$Exp)
		pldat[i,'semn']    <- sqrt((sum(subdbase$Exp*xvar^2)/sum(subdbase$Exp) - pldat[i,'Expmn']^2)/mean(subdbase$Nsamp_adj))
		pldat[i,'Vexp']    <- sum(subdbase$Exp*xvar^2) - pldat[i,'Expmn']^2  ## Francis 2011, p.1137 (RH 210225)
		pldat[i,'N']       <- mean(subdbase$Nsamp_in,na.rm=TRUE) ## (RH 210225)
		pldat[i,'Obslo']   <- pldat[i,'Obsmn'] - 2 * pldat[i,'semn']
		pldat[i,'Obshi']   <- pldat[i,'Obsmn'] + 2 * pldat[i,'semn']
		pldat[i,'Std.res'] <- (pldat[i,'Obsmn'] - pldat[i,'Expmn'])/pldat[i,'semn']
		pldat[i,'Fleet']   <- mean(subdbase$Fleet)
		pldat[i,'Yr']      <- mean(if(seas=='comb')subdbase$Yr else subdbase$Yr.S)
		if(type=='con')
			pldat[i,'Lbin'] <- mean(subdbase$'Lbin_lo')
		if(sex.flag)
			pldat[i,'sexes'] <- mean(subdbase$'Sexes')
		if(type=='size'){
			pldat[i,'method'] <- mean(subdbase$method)
			plunits[i] <- subdbase$units[1] # units of size comps
		}
#if ("2 1994  " %in% uindx[i]) {browser();return()}
	}
	pldat = as.data.frame(pldat)
	lldat = split(pldat, pldat$Fleet)
	Nmult = sapply(lldat, function(x) { 1/var(x[,'Std.res'],na.rm=TRUE) })
	#Nmult <- 1/var(pldat[,'Std.res'],na.rm=TRUE)
#browser();return()
	wj   = sapply(lldat, function(x) { 1 / var((x[,'Obsmn']-x[,'Expmn'])/((x[,'Vexp']/x[,'N'])^0.5), na.rm=TRUE) })

	## Find the adjusted confidence intervals
	for(i in 1:length(uindx)){
		imult = Nmult[as.character(pldat[i,'Fleet'])]
		pldat[i,'ObsloAdj'] <- pldat[i,'Obsmn'] - 2 * pldat[i,'semn']/sqrt(imult)
		pldat[i,'ObshiAdj'] <- pldat[i,'Obsmn'] + 2 * pldat[i,'semn']/sqrt(imult)
		#pldat[i,'ObsloAdj'] <- pldat[i,'Obsmn'] - 2 * pldat[i,'semn']/sqrt(Nmult)
		#pldat[i,'ObshiAdj'] <- pldat[i,'Obsmn'] + 2 * pldat[i,'semn']/sqrt(Nmult)
	}
	Nfleet <- length(unique(pldat[,'Fleet']))

	## make plot if requested
	if(plotit){
		createFdir(lang)
		fout = fout.e = outnam
		for (l in lang) {
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )

			col.obs = rep(col.obs,2)[1:2]
			col.fit = rep(col.fit,2)[1:2]
			plindx <- if(type=='con'){
				paste(pldat[,'Fleet'], pldat[,'Yr'])
			} else {
				pldat[,'Fleet']
			}
			if(sex.flag)
				plindx <- paste(plindx,pldat[,'sexes'])
			if(method.flag)
				plindx <- paste(plindx,pldat[,'method'])
			uplindx <- unique(plindx)
	
				if (png) {
					clearFiles(paste0(fout,".png"))
					png(filename=paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				}
	
			## Select number of panels
			Npanel <- length(uplindx)
			## Ian T. 9/25/14: changing from having at least 4 panels to no minimum
			#NpanelSet <- max(4,min(length(uplindx),maxpanel))
			NpanelSet <- min(length(uplindx),maxpanel)
			#Nr <- ceiling(sqrt(NpanelSet)); Nc <- ceiling(NpanelSet/Nr)
			Nr = Npanel; Nc = 1
			if(set.pars){
				# save current graphical parameters
				par_current <- par()
				# set new parameters
				par(mfrow=c(Nr,Nc),mar=c(1.5,2,1,1)+0.1,mgp=c(0,0.5,0),oma=c(1.2,1.2,0,0), las=1)
				par(cex=1)
			}
			for(i in 1:Npanel){
				## loop over panels
				subpldat <- pldat[plindx==uplindx[i],,drop=FALSE]
				x <- subpldat[,ifelse(type=='con','Lbin','Yr')]
				# calculate ylim, including removing Inf values
				plot(x,subpldat[,'Obsmn'], type="n", xaxt="n", #pch='-', 
					xlim = if(length(x)>2) range(x) else c(min(x)-(0.5/length(x)),max(x)+(0.5/length(x))),
					ylim = range(subpldat[,c('Obslo','Obshi','ObsloAdj','ObshiAdj','Expmn')], na.rm=TRUE),
					xlab='',ylab='')
				bigtck = intersect(seq(1900, 3000, ifelse((max(x)-min(x))<=5,1,5)), floor(min(x)):ceiling(max(x)) )
				liltck = intersect(seq(1900, 3000, 1), floor(min(x)):ceiling(max(x)) )
				axis(1, at=liltck, labels=FALSE, tcl=-0.2)
				axis(1, at=bigtck, labels=bigtck, tcl=-0.4)
				segments(x, subpldat[,'Obslo'], x, subpldat[,'Obshi'], lwd=2, lend=3, col=col.obs[1])
				if(plotadj){
					arrows(x,subpldat[,'ObsloAdj'],x,subpldat[,'ObshiAdj'],lwd=2, length=0.04, angle=90, code=3, col=col.obs[1])
				}
				points(x, subpldat[,'Obsmn'], pch=21,col=col.obs[1], bg=col.obs[2], cex=1)
				#points(x,subpldat[,'Obsmn'],pch=21,bg='grey80')
				ord <- order(x)
				if(!datonly){
					if(length(x)>1){
						lines(x[ord], subpldat[ord,'Expmn'], lwd=2, col=col.fit[1])
						#points(x[ord], subpldat[ord,'Expmn'], pch=22, col=col.fit[1], bg=col.fit[2], cex=0.8)
					} else {
						lines(c(x-0.5,x+0.5), rep(subpldat[,'Expmn'],2), col=col.fit[1])
					}
				}
	#browser();return()
				## Lines
				fl <- fleetnames[subpldat[1,'Fleet']]
				yr <- paste(subpldat[1,'Yr'])
				lab <- if(type=='con')ifelse(Nfleet>1,paste(yr,fl),yr) else fl
				if(sex.flag & label.sex){
					lab <- paste(lab,ifelse(subpldat[1,'sexes']==0,'comb','sex'))
				}
				if(method.flag){
					lab <- paste(lab,'meth',subpldat[1,'method'])
				}
				if(label.part){
					lab <- paste(lab,partition.labels)
				}
				mtext(linguaFranca(lab,l),side=3) #,at=mean(x))
			}
			## define y-axis label
			ylab <- 'Mean age' # default as age unless replaced below
			if(type=="len"){
				ylab <- 'Mean length'
			}
			if(type=="size"){
				## probably more efficient ways to sort out these labels,
				## but lots of if-statements make logic easier to follow
				units <- unique(plunits[plindx %in% uplindx])
				if(length(units)==1){ # not sure if this will always be true or not
					if(units %in% c('kg','lb')){
						ylab <- paste0('Mean weight (',units,')')
					}
					if(units %in% c('cm','in')){
						ylab <- paste0('Mean length (',units,')')
					}
				} else {
				# just in case it's possible to have multiple units in one panel
				ylab <- paste0('Mean value (',paste(units, collapse=' or '),')')
				}
			}
			mtext(linguaFranca(ylab,l), side=2,las=0,outer=TRUE)
			mtext(linguaFranca(ifelse(type=='con','Length','Year'),l),side=1,outer=TRUE)
			if (png) dev.off()
	
			## restore previous graphics parameters (if changed to begin with)
			if(set.pars){
				par(mfrow=par_current$mfrow, mar=par_current$mar, mgp=par_current$mgp,
				oma=par_current$oma, las=par_current$las)
			}
		}; eop()
	}
	if(!datonly) {
		lldat = split(pldat, pldat$Fleet)
		lltmp = lapply(lldat, function(x){
			matrix(sample(x[,'Std.res'],1000*nrow(x),replace=TRUE),nrow(x))
		})
		confint = sapply(lltmp, function(x){
			ivar = apply(x,2,function(xx) {
				xxx = 1/var(xx,na.rm=TRUE)
				out = ifelse(is.finite(xxx),xxx,NA)
				return(out)
			})
			## Note: if there are only two indices (with different values), there are only 3 possible random samples, and only one with non-zero variance.
			qnts = quantile(ivar, c(0.025,0.975), na.rm=TRUE)
			pooh = as.vector(qnts)
			return(pooh)
		})
#browser();return()
		Output <- list(w=Nmult, lo=confint[1,], hi=confint[2,], w.francis=wj)
		Outs <- paste0("Francis Weights - ", type, ": ", fleetnames[fleet], ": ", round(Output[["w"]],4), " (", round(Output[["lo"]],4), "-", round(Output[["hi"]],4), ")")
		
		
		#tmp <- matrix(sample(pldat[,'Std.res'],1000*nrow(pldat),replace=TRUE),nrow(pldat))
#browser();return()
		#confint <- as.vector(quantile(apply(tmp,2,function(x)1/var(x,na.rm=TRUE)), c(0.025,0.975), na.rm=TRUE))
		#Output <- c(w=Nmult,lo=confint[1],hi=confint[2])
		#Outs <- paste("Francis Weights - ", type, ": ", fleetnames[fleet],": ", round(Nmult,4), " (",round(confint[1],4),"-",round(confint[2],4),")", sep="")
		if(printit){
			print(Outs)
		}
		w.francis = Output; ttput(w.francis)
		return(Output)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.francis


## plotSS.dmcmc-------------------------2021-05-25
##  Plot MCMC diagnostics (traces, split chains, ACFs).
## -----------------------------------PBSawatea|RH
plotSS.dmcmc = function(mcmcObj, mpdObj, ptypes, lang, pngres=400, PIN=c(9,9))
{
	## MCMC diagnostics
	## ----------------
	fout = fout.e = "traceBiomass"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
			bmcmc = mcmcObj$B[getYrIdx(colnames(mcmcObj$B))]/1000
			bmpd  = mpdObj$B[getYrIdx(names(mpdObj$B))]/1000
#browser();return()
			panelTraces(mcmc=bmcmc, mpd=bmpd, xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, same.limits=F, lang=l)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()

	fout = fout.e = "traceRecruits"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
			rmcmc = mcmcObj$R[getYrIdx(colnames(mcmcObj$R))]/1000
			rmpd  = mpdObj$R[getYrIdx(names(mpdObj$R))]/1000
			panelTraces(mcmc=rmcmc, mpd=rmpd, xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, same.limits=F, lang=l)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()

	fout = fout.e = "traceParams"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
			pmcmc = mcmcObj$P
			pmpd  = mpdObj$P
			panelTraces(mcmc=pmcmc, mpd=pmpd, xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, same.limits=F, lang=l)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
#browser();return()

	fout = fout.e = "splitChain"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			par(mfrow=c(1,1), mar=c(3,3,0.5,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))  ## mar and oma ignored, fixed in call to `mochaLatte'
			panelChains(mcmc=mcmcObj$P, axes=TRUE, pdisc=0, between=list(x=0, y=0), col.trace=c("red","blue","black"), xlab="Parameter Value", ylab="Cumulative Frequency", cex.axis=1.2, cex.lab=1.4, yaxt="n", lang=l)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()

	fout = fout.e = "paramACFs"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
			plotACFs(mcmc=mcmcObj$P, lag.max=60, lang=l)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
}


## plotSS.index-------------------------2021-02-12
##  Plot SS model fit to abundance index series
## ----------------------------------------r4ss|RH
plotSS.index = function (replist, subplots=c(1:10, 12), plot=TRUE, print=FALSE, 
   fleets="all", fleetnames="default", smooth=TRUE, add=FALSE, 
   datplot=TRUE, labels=list("Year", "Index", "Observed index", 
   "Expected index", "Log index", "Log observed index", 
   "Log expected index", "Standardized index", "Catchability (Q)", 
   "Time-varying catchability", "Vulnerable biomass", "Catchability vs. vulnerable biomass", 
   "Residual", "Deviation"), col1="default", col2="default", 
   col3="blue", col4="red", pch1=21, pch2=16, cex=1, 
   bg="white", legend=TRUE, legendloc="topright", seasnames=NULL, 
   pwidth=9, pheight=7, punits="in", res=400, ptsize=10, PIN=c(9,9),
   cex.main=1, mainTitle=FALSE, plotdir="default", minyr=NULL, 
   maxyr=NULL, maximum_ymax_ratio=Inf, show_input_uncertainty=TRUE, 
   verbose=TRUE, onepage=FALSE, outnam, lang="e", ...) 
{
	oldpar=par(no.readonly=TRUE)
	fart=function(opar) { if (any("windows"%in%names(dev.list()))) par(opar) }
	on.exit(fart(oldpar))

	cpue <- replist$cpue
	SS_versionNumeric <- replist$SS_versionNumeric
	if (is.null(dim(cpue))) {
		message("skipping index plots: no index data in this model")
		return()
	}
	pngfun <- function(file, caption=NA, lang="e") {
		createFdir(lang, dir=plotdir)
		changeLangOpts(L=lang)
		fout = switch(lang, 'e' = file.path(plotdir, file), 'f' = file.path(plotdir,"french", file) )
		clearFiles(fout)
		png(filename=fout, width=pwidth, height=pheight, units=punits, res=res, pointsize=ptsize)
		if (onepage)
			expandGraph(mfrow=rc, mar=c(3,3.5,0.5,0.5), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		else 
			expandGraph(mfrow=c(1,1), mar=c(3,3,1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		plotinfo <- rbind(plotinfo, data.frame(file=file, caption=caption))
		return(plotinfo)
	}
	## Deal with messy label details neglected by original program
	labels.default=list("Year", "Index", "Observed index", "Expected index", "Log index", "Log observed index", "Log expected index", "Standardized index", "Catchability (Q)", "Time-varying catchability", "Vulnerable biomass", "Catchability vs. vulnerable biomass", "Residual", "Deviation")
	names(labels.default)=1:12
	if (is.null(names(labels)))
		names(labels)=1:length(labels)
	labels.new=labels.default
	labels.new[names(labels)]=labels
	labels=labels.new
	nlabs =sapply(labels,countVec) ## how many labels are available
#browser();return()
	get.lab=function(i, j, k=nlabs)
	{
		labels[[i]][min(j,k[i])]
	}
	plotinfo <- NULL

	index.fn <- function(addexpected=TRUE, log=FALSE, l="e", ...) {
		if (error == -1 & log == TRUE) {
			return()
		}
		if (error == 0) {
			if (!log) {
				lower_total <- qlnorm(0.025, meanlog=log(y[include]), 
					sdlog=cpueuse$SE[include])
				upper_total <- qlnorm(0.975, meanlog=log(y[include]), 
					sdlog=cpueuse$SE[include])
			}
			else {
				lower_total <- qnorm(0.025, mean=log(y[include]), 
					sd=cpueuse$SE[include])
				upper_total <- qnorm(0.975, mean=log(y[include]), 
					sd=cpueuse$SE[include])
			}
		}
		if (error == -1) {
			lower_total <- qnorm(0.025, mean=y[include], sd=cpueuse$SE[include])
			upper_total <- qnorm(0.975, mean=y[include], sd=cpueuse$SE[include])
		}
		if (error > 0) {
			lower_total <- log(y[include]) + qt(0.025, df=error) * 
				cpueuse$SE[include]
			upper_total <- log(y[include]) + qt(0.975, df=error) * 
				cpueuse$SE[include]
			if (!log) {
				lower_total <- exp(lower_total)
				upper_total <- exp(upper_total)
			}
		}
		if (max(upper_total) == Inf) {
			warning("Removing upper interval on indices with infinite upper quantile values.\n", 
				"Check the uncertainty inputs for the indices.")
			upper_total[upper_total == Inf] <- 100 * max(cpueuse$Obs[upper_total == 
				Inf])
		}
		main <- paste0(get.lab(2,ifleet), Fleet)
		if (log) {
			main <- paste0(get.lab(5,ifleet), Fleet)
		}
		if (!mainTitle) {
			main <- ""
		}
		xlim <- c(max(minyr, min(x)), min(maxyr, max(x)))
		if (legend & length(colvec1) > 1) {
			xlim[2] <- xlim[2] + 0.25 * diff(xlim)
		}
		if (!add) {
			zmax <- NULL
			if (addexpected) {
				zmin <- min(z, na.rm=TRUE)
			}
			logzrange <- range(log(z))
			if (!log) {
				ylim <- c(0, 1.05 * min(max(upper_total, zmax, 
					na.rm=TRUE), max(maximum_ymax_ratio * y)))
			}
			if (log) {
				ylim <- range(c(lower_total, upper_total), na.rm=TRUE)
			}
			plot(x=x[include], y=y[include], type="n", xlim=xlim, ylim=ylim, yaxs=ifelse(log, "r", "i"), 
				xlab=linguaFranca(get.lab(1,ifleet),l), ylab=linguaFranca("Abundance Index",l), main=linguaFranca(main,l), cex.main=cex.main, cex.axis=1.2, cex.lab=1.5, ...)
			panlab = ifelse(!log, get.lab(2,ifleet), get.lab(5,ifleet))
			addLabel(0.925, 0.975, linguaFranca(panlab,l), adj=c(1,1), cex=1.5, col="blue")
#browser();return()
		}
		if (show_input_uncertainty && any(!is.null(cpueuse$SE_input[include]))) {
			if (error == 0) {
				if (!log) {
					lower_input <- qlnorm(0.025, meanlog=log(y[include]), 
					sdlog=cpueuse$SE_input[include])
					upper_input <- qlnorm(0.975, meanlog=log(y[include]), 
					sdlog=cpueuse$SE_input[include])
				}
				else {
					lower_input <- qnorm(0.025, mean=log(y[include]), 
					sd=cpueuse$SE_input[include])
					upper_input <- qnorm(0.975, mean=log(y[include]), 
					sd=cpueuse$SE_input[include])
				}
			}
			if (error == -1) {
				lower_input <- qnorm(0.025, mean=y[include], 
					sd=cpueuse$SE_input[include])
				upper_input <- qnorm(0.975, mean=y[include], 
					sd=cpueuse$SE_input[include])
			}
			if (error > 0) {
				lower_total <- log(y[include]) + qt(0.025, df=error) * 
					cpueuse$SE_input[include]
				upper_total <- log(y[include]) + qt(0.975, df=error) * 
					cpueuse$SE_input[include]
				if (!log) {
					lower_total <- exp(lower_total)
					upper_total <- exp(upper_total)
				}
			}
			segments(x[include], lower_input, x[include], upper_input, 
				col=colvec1[s], lwd=3, lend=1)
		}
		#arrows(x0=x[include], y0=lower_total, x1=x[include], y1=upper_total, length=0.03, angle=90, code=3, col=colvec1[s])
		arrows(x0=x[include], y0=lower_total, x1=x[include], y1=upper_total, length=0.03, angle=90, code=3, col="black", lwd=2)
		if (!log) {
			#points(x=x[include], y=y[include], pch=pch1, cex=cex, bg=bg, col=colvec1[s])
			points(x=x[include], y=y[include], pch=15, cex=1.2, col="red")
#browser();return()
			if (addexpected) {
				lines(x, z, lwd=2, col=col3)
			}
		}
		else {
			points(x=x[include], y=log(y[include]), pch=pch1, 
				cex=cex, bg=bg, col=colvec1[s])
			if (addexpected) {
				lines(x, log(z), lwd=2, col=col3)
			}
		}
		if (legend & length(colvec1) > 1) {
			legend(x=legendloc, legend=seasnames, pch=pch1, 
				col=colvec1, cex=cex)
		}
	}

	index_resids.fn <- function(option=1, l="e", ...) {
		if (option == 1) {
			ylab <- get.lab(13,ifleet)
			y <- (log(cpueuse$Obs) - log(cpueuse$Exp))/cpueuse$SE
		}
		if (error == 0 & option == 2) {
			ylab <- get.lab(13,ifleet)
			y <- (log(cpueuse$Obs) - log(cpueuse$Exp))/cpueuse$SE_input
		}
		if (option == 3) {
			ylab <- get.lab(14,ifleet)
			y <- cpueuse$Dev
		}
		main <- paste(ylab, Fleet)
		if (!mainTitle) {
			main <- ""
		}
		xlim <- c(max(minyr, min(x)), min(maxyr, max(x)))
		if (legend & length(colvec1) > 1) {
			xlim[2] <- xlim[2] + 0.25 * diff(xlim)
		}
		ylim <- c(-1.05, 1.05) * max(abs(y[include]))
		if (!add) {
			plot(x=x[include], y=y[include], type="n", xlim=xlim, ylim=ylim, yaxs="i", xlab=linguaFranca(get.lab(1,ifleet),l), ylab=linguaFranca(ylab,l), main=linguaFranca(main,l), cex.main=cex.main, ...)
		}
		points(x=x[include], y=y[include], pch=pch1, cex=cex, bg=adjustcolor(colvec1[s], alpha.f=0.7), col=adjustcolor(colvec1[s], alpha.f=0.7))
		abline(h=0, lty=3)
		if (legend & length(colvec1) > 1) {
			legend(x=legendloc, legend=linguaFranca(seasnames,l), pch=pch1, pt.bg=colvec1, col=colvec1, cex=cex)
		}
	}

	obs_vs_exp.fn <- function(log=FALSE, l="e", ...) {
		main <- paste(get.lab(2,ifleet), Fleet, sep=" ")
		if (!mainTitle) {
			main <- ""
		}
		if (!add) {
			if (!log) {
				plot(y[include], z[include], type="n", linguaFranca(xlab=get.lab(3,ifleet),l), 
					ylab=linguaFranca(get.lab(4,ifleet),l), main=linguaFranca(main,l), cex.main=cex.main, 
					ylim=c(0, 1.05 * max(z)), xlim=c(0, 1.05 * max(y)), xaxs="i", yaxs="i", ...)
			}
			else {
				plot(log(y[include]), log(z[include]), type="n", cex.main=cex.main, 
					xlab=linguaFranca(get.lab(6,ifleet),l), ylab=linguaFranca(get.lab(7,ifleet),l), main=linguaFranca(main,l) )
			}
		}
		if (!log) {
			points(y[include], z[include], col=colvec2[s], pch=pch2, cex=cex)
		}
		else {
			points(log(y[include]), log(z[include]), col=colvec2[s], pch=pch2, cex=cex)
		}
		abline(a=0, b=1, lty=3)
		if (smooth && npoints > 6 && diff(range(y)) > 0) {
			if (!log) {
				psmooth <- loess(z[include] ~ y[include], degree=1)
				lines(psmooth$x[order(psmooth$x)], psmooth$fit[order(psmooth$x)], lwd=1.2, col=col4, lty="dashed")
			}
			else {
				psmooth <- loess(log(z[include]) ~ log(y[include]), degree=1)
				lines(psmooth$x[order(psmooth$x)], psmooth$fit[order(psmooth$x)], lwd=1.2, col=col4, lty="dashed")
			}
		}
		if (legend & length(colvec2) > 1) {
			legend(x=legendloc, legend=linguaFranca(seasnames,l), pch=pch2, col=colvec2, cex=cex)
		}
	}

	timevarying_q.fn <- function(l="e") {
		main <- paste(get.lab(10,ifleet), Fleet, sep=" ")
		if (!mainTitle) 
			main <- ""
		q <- cpueuse$Calc_Q
		if (!add) 
			plot(x, q, type="o", cex.main=cex.main, col=colvec2[1], pch=pch2,
				xlab=linguaFranca(get.lab(1,ifleet),l), ylab=linguaFranca(get.lab(9,ifleet),l), main=linguaFranca(main,l) )
	}

	q_vs_vuln_bio.fn <- function(l="e") {
		main <- paste(get.lab(12,ifleet), Fleet, sep=" ")
		if (!mainTitle) 
			main <- ""
		v <- cpueuse$Vuln_bio
		q1 <- cpueuse$Calc_Q
		q2 <- cpueuse$Eff_Q
		if (all(q1 == q2)) 
			ylab <- get.lab(9,ifleet)
		else ylab <- "Effective catchability"
		if (!add) 
			plot(v, q2, type="o", cex.main=cex.main, col=colvec2[1], pch=pch2,
				xlab=linguaFranca(get.lab(11,ifleet),l), ylab=linguaFranca(ylab,l), main=linguaFranca(main,l) )
	}

	if (length(grep("supr_per", cpue$Supr_Per))) {
		warning("Some indices have superperiods. Values will be plotted\n", 
			"in year/season associated with data in report file.")
		cpue <- cpue[!is.na(cpue$Dev), ]
	}
	FleetNames <- replist$FleetNames
	nfleets <- replist$nfleets
	nseasons <- replist$nseasons
	parameters <- replist$parameters
	Q_extraSD_info <- parameters[grep("Q_extraSD", parameters$Label), ]
	nSDpars <- nrow(Q_extraSD_info)
	if (nSDpars > 0) {
		Q_extraSD_info$Fleet <- NA
		for (ipar in 1:nSDpars) {
			if (SS_versionNumeric >= 3.3) {
				num <- strsplit(Q_extraSD_info$Label[ipar], split="[()]", fixed=FALSE)[[1]][2]
			}
			else {
				num <- strsplit(substring(Q_extraSD_info$Label[ipar], nchar("Q_extraSD_") + 1), split="_", fixed=TRUE)[[1]][1]
			}
			Q_extraSD_info$Fleet[ipar] <- as.numeric(num)
		}
	}
	if (nseasons > 1) {
		cpue$YrSeas <- cpue$Yr + (cpue$Seas - 0.5)/nseasons
	}
	else {
		cpue$YrSeas <- cpue$Yr
	}
	if (plotdir == "default") 
		plotdir <- replist$inputs$dir
	if (fleetnames[1] == "default") 
		fleetnames <- FleetNames
	if (fleets[1] == "all") {
		fleets <- 1:nfleets
	}
	else {
		if (length(intersect(fleets, 1:nfleets)) != length(fleets)) {
			return("Input 'fleets' should be 'all' or a vector of values between 1 and nfleets.")
		}
	}
	fleetvec <- intersect(fleets, unique(as.numeric(cpue$Fleet)))
	allcpue <- data.frame()
	any_negative <- FALSE
	
	for (j in subplots) {
		jj  =switch(j,
			"cpuedata", "cpuefit", "obs.vs.exp", "log.cpuedata", "log.cpuefit",
			"log.obs.vs.exp", "time.varying.q", "q.vs.vuln.bio", "stand.cpue.all", "resids.se.total",
			"resids.se.input", "resids.se.total" )
		if (missing(outnam))
			onefile=paste0("index", j, ".", jj, ".", ifelse(j==9,"","fleets"), ".png")
		else
			onefile = paste0(outnam, ".png")
		if (onepage) {
			if (print) {
				createFdir(lang, dir=plotdir)
				changeLangOpts(L=lang)
				fout = switch(lang, 'e' = file.path(plotdir, onefile), 'f' = file.path(plotdir,"french", onefile) )
				clearFiles(fout)
				png(filename=fout, units="in", res=res, width=PIN[1], height=PIN[2])
			}
			rc=.findSquare(length(fleetvec))
			expandGraph(mfrow=rc, mar=c(3,3.5,0.5,0.5), oma=c(0,0,0,0), mgp = c(2,0.5,0))
			#expandGraph(mfrow=c(length(fleetvec),1), mar=c(3,3,1,1), oma=c(0,0,0,0))
		}
		## Loops through by fleet making it tricky to have multipanel plot
		for (ifleet in fleetvec) {
			usecol <- FALSE
			if (length(unique(cpue$Seas[cpue$Fleet == ifleet])) > 1) {
				usecol <- TRUE
			}
			if (!usecol) {
				legend <- FALSE
			}
			if (col1[1] == "default") {
				colvec1 <- "black"
				if (usecol & nseasons == 4) {
					colvec1 <- c("blue4", "green3", "orange2", "red3")
				}
				if (usecol & !nseasons %in% c(1, 4)) {
					colvec1 <- rich.colors.short(nseasons)
				}
			}
			else {
				colvec1 <- col1
				if (length(colvec1) < nseasons) {
					colvec1 <- rep(col1, nseasons)
				}
			}
			if (col2[1] == "default") {
				colvec2 <- "blue"
				if (usecol & nseasons == 4) {
					colvec2 <- c("blue4", "green3", "orange2", "red3")
				}
				if (usecol & !nseasons %in% c(1, 4)) {
					colvec2 <- rich.colors.short(nseasons)
				}
			}
			else {
				colvec2 <- col2
				if (length(colvec1) < nseasons) {
					colvec1 <- rep(col1, nseasons)
				}
			}
			if (is.null(seasnames)) 
				seasnames <- paste("Season", 1:nseasons, sep="")
			Fleet <- fleetnames[ifleet]
			error <- replist$survey_error[ifleet]
			if (error == 0) {
				error_caption <- "lognormal error"
			}
			if (error == -1) {
				error_caption <- "normal error"
			}
			if (error == 1) {
				error_caption <- paste0("T-distributed error with ", 
					error, " degree of freedom")
			}
			if (error > 1) {
				error_caption <- paste0("T-distributed error with ", 
					error, " degrees of freedom")
			}
			cpueuse <- cpue[cpue$Fleet == ifleet, ]
			cpueuse <- cpueuse[order(cpueuse$YrSeas), ]
			time <- diff(range(cpueuse$Calc_Q)) > 0
			time2 <- diff(range(cpueuse$Eff_Q)) > 0
			if (is.na(time2)) {
				time2 <- FALSE
			}
			if (exists("Q_extraSD_info") && ifleet %in% Q_extraSD_info$Fleet) {
				cpueuse$SE_input <- cpueuse$SE - Q_extraSD_info$Value[Q_extraSD_info$Fleet == 
					ifleet]
			}
			else {
				cpueuse$SE_input <- NULL
			}
			x <- cpueuse$YrSeas
			y <- cpueuse$Obs
			z <- cpueuse$Exp
			npoints <- length(z)
			include <- !is.na(cpueuse$Like)
	#browser();return()
			if (any(include)) {
				if (usecol) {
					s <- cpueuse$Seas[which(include)]
				}
				else {
					s <- 1
				}
				if (datplot) {
					if (min(cpueuse$Obs >= 0)) {
						cpueuse$Index <- rep(ifleet, length(cpueuse$YrSeas))
						cpueuse$stdvalue <- cpueuse$Obs/mean(cpueuse$Obs)
						tempcpue <- cbind(cpueuse$Index, cpueuse$YrSeas, 
						cpueuse$Obs, cpueuse$stdvalue)
						colnames(tempcpue) <- c("Index", "year", "value", "stdvalue")
						allcpue <- rbind(allcpue, tempcpue)
					}
					else {
						if (verbose) {
						message("Excluding fleet ", ifleet, " from index comparison figure because it has negative values")
						}
						any_negative <- TRUE
					}
				}
				addlegend <- function(pch, colvec) {
					names <- paste(seasnames, "observations")
				}
				if (plot) {
					if (1 %in% subplots & datplot)
						index.fn(addexpected=FALSE, l=lang)
					if (2 %in% subplots) 
						index.fn(l=lang)
					if (3 %in% subplots) 
						obs_vs_exp.fn(l=lang)
				}
				if (print && !onepage) {
					if (1 %in% subplots & datplot) {
						file=sub("fleets", tolower(gsub("[[:punct:]]+",".",Fleet)), onefile)
						#file <- paste0("index1_cpuedata_", Fleet, ".png")
						caption <- paste0("Index data for ", Fleet, 
						". ", "Lines indicate 95% uncertainty interval around index values ", 
						"based on the model assumption of ", error_caption, 
						". ", "Thicker lines (if present) indicate input uncertainty before addition of ", 
						"estimated additional uncertainty parameter.")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						index.fn(addexpected=FALSE, l=lang)
						dev.off(); eop()
					}
					if (2 %in% subplots) {
						file <- paste0("index2_cpuefit_", Fleet, ".png")
						caption <- paste0("Fit to index data for ", 
						Fleet, ". ", "Lines indicate 95% uncertainty interval around index values ", 
						"based on the model assumption of ", error_caption, 
						". ", "Thicker lines (if present) indicate input uncertainty before addition of ", 
						"estimated additional uncertainty parameter.")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						index.fn(l=lang)
						dev.off(); eop()
					}
					if (3 %in% subplots) {
						file <- paste0("index3_obs_vs_exp_", Fleet, 
						".png")
						caption <- paste("Observed vs. expected index values with smoother for", 
						Fleet)
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						obs_vs_exp.fn(l=lang)
						dev.off(); eop()
					}
				}
				if (error != -1) {
					if (plot) {
						if (4 %in% subplots & datplot) {
						index.fn(log=TRUE, addexpected=FALSE, l=lang)
						}
						if (5 %in% subplots) {
						index.fn(log=TRUE, l=lang)
						}
						if (6 %in% subplots) {
						obs_vs_exp.fn(log=TRUE, l=lang)
						}
					}
					if (print) {
						if (4 %in% subplots & datplot) {
						file <- paste0("index4_logcpuedata_", Fleet, 
							".png")
						caption <- paste0("Log index data for ", 
							Fleet, ". ", "Lines indicate 95% uncertainty interval around index values ", 
							"based on the model assumption of ", error_caption, 
							". ", "Thicker lines (if present) indicate input uncertainty before addition of ", 
							"estimated additional uncertainty parameter.")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						index.fn(log=TRUE, addexpected=FALSE, l=lang)
						dev.off(); eop()
						}
						if (5 %in% subplots) {
						file <- paste0("index5_logcpuefit_", Fleet, 
							".png")
						caption <- paste0("Fit to log index data on log scale for ", 
							Fleet, ". ", "Lines indicate 95% uncertainty interval around index values ", 
							"based on the model assumption of ", error_caption, 
							". ", "Thicker lines (if present) indicate input uncertainty before addition of ", 
							"estimated additional uncertainty parameter.")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						index.fn(log=TRUE, l=lang)
						dev.off(); eop()
						}
						if (6 %in% subplots) {
						file <- paste0("index6_log_obs_vs_exp_", 
							Fleet, ".png")
						caption <- paste("log(observed) vs. log(expected) index values with smoother for", 
							Fleet)
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						obs_vs_exp.fn(log=TRUE, l=lang)
						dev.off(); eop()
						}
					}
				}
				if (plot) {
					if (7 %in% subplots & time) {
						timevarying_q.fn(l=lang)
					}
					if (8 %in% subplots & time2) {
						q_vs_vuln_bio.fn(l=lang)
					}
				}
				if (print) {
					if (7 %in% subplots & time) {
						file <- paste0("index7_timevarying_q_", Fleet, 
						".png")
						caption <- paste("Timeseries of catchability for", 
						Fleet)
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						timevarying_q.fn(l=lang)
						dev.off(); eop()
					}
					if (8 %in% subplots & time2) {
						file <- paste0("index8_q_vs_vuln_bio_", Fleet, 
						".png")
						caption <- paste0("Catchability vs. vulnerable biomass for fleet ", 
						Fleet, "<br> \n", "This plot should illustrate curvature of nonlinear catchability relationship<br> \n", 
						"or reveal patterns associated with random-walk catchability.")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						q_vs_vuln_bio.fn(l=lang)
						dev.off(); eop()
					}
				}
				if (plot) {
					if (10 %in% subplots) {
						index_resids.fn(option=1, l=lang)
					}
					if (11 %in% subplots) {
						index_resids.fn(option=2, l=lang)
					}
					if (12 %in% subplots) {
						index_resids.fn(option=3, l=lang)
					}
				}
				if (print) {
					if (10 %in% subplots) {
						file <- paste0("index10_resids_SE_total_", 
						Fleet, ".png")
						caption <- paste0("Residuals of fit to index for ", 
						Fleet, ".")
						if (error == 0) {
						caption <- paste0(caption, "<br>Values are (log(Obs) - log(Exp))/SE ", 
							"where SE is the total standard error including any ", 
							"estimated additional uncertainty.")
						}
						else {
						caption <- paste0(caption, "<br>Values are based on the total standard error ", 
							"including any estimated additional uncertainty.")
						}
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						index_resids.fn(option=1, l=lang)
						dev.off(); eop()
					}
					if (11 %in% subplots & show_input_uncertainty && 
						any(!is.null(cpueuse$SE_input[include])) && 
						any(cpueuse$SE_input > cpueuse$SE)) {
						file <- paste0("index11_resids_SE_input_", 
						Fleet, ".png")
						caption <- paste0("Residuals for fit to index for ", 
						Fleet, ".")
						if (error == 0) {
						caption <- paste0(caption, "<br>Values are (log(Obs) - log(Exp))/SE_input ", 
							"where SE_input is the input standard error", 
							"excluding any estimated additional uncertainty.")
						}
						else {
						caption <- paste0(caption, "<br>Values are based on the input standard error ", 
							"excluding any estimated additional uncertainty.")
						}
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						index_resids.fn(option=2, l=lang)
						dev.off(); eop()
					}
					if (12 %in% subplots) {
						file <- paste0("index12_resids_SE_total_", 
						Fleet, ".png")
						caption <- paste0("Deviations for fit to index for ", 
						Fleet, ".")
						if (error != -1) {
						caption <- paste0(caption, "<br>Values are log(Obs) - log(Exp) ", 
							"and thus independent of index uncertainty.")
						}
						if (error == -1) {
						caption <- paste0(caption, "<br>Values are Obs - Exp ", 
							"and thus independent of index uncertainty.")
						}
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						index_resids.fn(option=3, l=lang)
						dev.off(); eop()
					}
				}
			}
		}
		if (onepage && print) {dev.off(); eop()}
	}
	if (datplot == TRUE & nrow(allcpue) > 0) {
		all_index.fn <- function(l="e") {
			main <- "All index plot"
			if (!mainTitle) {
				main <- ""
			}
			xlim <- c(min(allcpue$year, na.rm=TRUE) - 1, max(allcpue$year, na.rm=TRUE) + 1)
			xlim[1] <- max(xlim[1], minyr)
			xlim[2] <- min(xlim[2], maxyr)
			ylim <- c(range(allcpue$stdvalue, na.rm=TRUE))
			usecols <- rich.colors.short(max(allcpue$Index, na.rm=TRUE), alpha=0.7)
			if (max(allcpue$Index, na.rm=TRUE) >= 2) {
				usecols <- rich.colors.short(max(allcpue$Index, na.rm=TRUE) + 1, alpha=0.7)[-1]
			}
			if (!add) 
				plot(0, type="n", xlab=get.lab(1,ifleet), main=main, 
					cex.main=cex.main, col=usecols[1], ylab=get.lab(8,ifleet), xlim=xlim, ylim=ylim)
			for (ifleet in fleetvec) {
				points(x=allcpue$year[allcpue$Index == ifleet], 
					y=allcpue$stdvalue[allcpue$Index == ifleet], 
					pch=pch2, col=usecols[ifleet], cex=cex, 
					lwd=1, lty="dashed", type="o")
			}
		}
		if (plot & (9 %in% subplots)) {
			all_index.fn(l=lang)
		}
		if (print & (9 %in% subplots)) {
			file <- paste0("index9_standcpueall", ".png")
			caption <- paste("Standardized indices overlaid.", 
				"Each index is rescaled to have mean observation=1.0.")
			if (any_negative) {
				caption <- paste(caption, "Indices with negative observations have been excluded.")
			}
			plotinfo <- pngfun(file=file, caption=caption, lang=lang)
			all_index.fn(l=lang)
			dev.off(); eop()
		}
	}
	if (!is.null(plotinfo)) 
		plotinfo$category <- "Index"
	return(invisible(plotinfo))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.index


## plotSS.pairs-------------------------2021-05-26
##  Pairs|density plot comparison among parameters
## ---------------------------------------------RH
plotSS.pairs = function(P.mpd, P.mcmc, type="image", ptypes, lang=c("e","f"), pngres=400, PIN=c(10,10))
{
	panel.cor <- function(x, y, digits=2, prefix="", ...)
	{
		usr <- par("usr"); on.exit(par(usr))
		par(usr=c(0, 1, 0, 1))
		r <- abs(cor(x, y))
		txt <- format(c(r, 0.123456789), digits=digits)[1]
		txt <- paste(prefix, txt, sep="")
		text(0.5, 0.5, txt, cex=max(cex.min,((r)^0.1)*cex.fac))
	}
	panel.cor.small = eval(parse(text=sub("1\\.75", "1.25", deparse(panel.cor))))

	panel.dens <- function(x, y, type="image", ncol=100, ...)
	{
		cr = colorRampPalette(c("aliceblue","blue","green","yellow","orange","red"))
		if (type %in% "dots"){
			points (x,y,...)
		}
		if (type %in% c("contour","image")){
			xy = cbind(x, y)
			bw = apply(xy,2,function(x){quantile(abs(diff(x))*0.25,0.5)})
			est <- KernSmooth::bkde2D(xy, bandwidth=bw)#c(0.25,0.25))
			#contour(est$x1, est$x2, est$fhat, add=T)
			#image(est$x1,est$x2,est$fhat,col=gray(200:100/200), add=T)
			image(est$x1,est$x2,est$fhat,col=cr(ncol), add=TRUE)
			box()
		}
	}
	panel.dens.type = eval(parse(text=sub("type = \"image\"", paste0("type = \"",type,"\""), deparse(panel.dens))))

	panel.text = function(x, y, labels, cex, font, ...)
	{
		polygon(c(0,0,1,1), c(0,1,1,0), col="gainsboro")
		text(0.5,0.5,labels,cex=cex.max,font=2)
	}
#browser();return()

	## Doing 1 pairs plot with all parameters
	npp = 1
	nuP = length(P.mpd)
	npr = ceiling(nuP/npp)
	P.use = P.mcmc
	N.use = gsub("_"," ",colnames(P.use))
	colnames(P.use) = sapply(lapply(N.use,strwrap,width=10),paste0,collapse="\n")
	colnames(P.use) = gsub("_","\n",colnames(P.use))

	## Pairs plot
	for (i in 1:npp) {
		if (i<npp) ii = (1:npr)+(i-1)*npr
		else       ii = (nuP-npr+1):nuP
		#ii = ii[1:10] ## testing only
		if (length(ii)>20) {
			cex.min=0.6; cex.max=0.8; cex.fac=1.25
		} else if (length(ii)>10) {
			cex.min=0.65; cex.max=0.9; cex.fac=1.75
		} else {
			cex.min=0.7; cex.max=1; cex.fac=2.5
		}
		fout = fout.e = paste0("pairsPars", ifelse(npp>1, pad0(i,2), ""))
		for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
			createFdir(lang, dir=".")
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
				else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				par(mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(2,0.75,0))
				pairs(P.use[, ii], col=lucent("black",0.1), pch=20, cex=0.2, gap=0, lower.panel=panel.cor.small, upper.panel=panel.dens.type, text.panel=panel.text, cex.axis=1, cex.labels=1)
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.pairs


## plotSS.pars--------------------------2021-03-30
##  Plot parameter fits and priors.
## ----------------------------------------r4ss|RH
plotSS.pars = function (replist, plotdir=NULL, xlab="Parameter value", 
	ylab="Density", showmle=TRUE, showpost=TRUE, showprior=TRUE, 
	showinit=TRUE, showdev=FALSE, showlegend=TRUE, fitrange=FALSE, fitnudge=0,
	xaxs="i", xlim=NULL, ylim=NULL, verbose=TRUE, debug=FALSE, 
	nrows=3, ncols=3, ltyvec=c(1, 1, 3, 4), 
	colvec=c("blue", "red", "black", "green", rgb(0,0,0,0.5)), add=FALSE, 
	plot=TRUE, print=FALSE, punits="in", 
	ptsize=10, strings=NULL, exact=FALSE, newheaders=NULL,
	outnam, res=400, PIN=c(8,8), lang="e") 
{
	GetPrior <- function(Ptype, Pmin, Pmax, Pr, Psd, Pval) {
		Prior_Like <- NULL
		if (is.na(Ptype)) {
			warning("problem with prior type interpretation. Ptype:", Ptype)
		}
		Pconst <- 0.0001
		if (Ptype %in% c("No_prior", "")) {
			Prior_Like <- rep(0, length(Pval))
		}
		if (Ptype == "Normal") {
			Prior_Like <- 0.5 * ((Pval - Pr)/Psd)^2
		}
		if (Ptype == "Sym_Beta") {
			mu <- -(Psd * (log((Pmax + Pmin) * 0.5 - Pmin))) - (Psd * (log(0.5)))
			Prior_Like <- -(mu + (Psd * (log(Pval - Pmin + Pconst))) + (Psd * (log(1 - ((Pval - Pmin - Pconst)/(Pmax - Pmin))))))
		}
		if (Ptype == "Full_Beta") {
			mu <- (Pr - Pmin)/(Pmax - Pmin)
			tau <- (Pr - Pmin) * (Pmax - Pr)/(Psd^2) - 1
			Bprior <- tau * mu
			Aprior <- tau * (1 - mu)
			if (Bprior <= 1 | Aprior <= 1) {
				warning("bad Beta prior")
			}
			Prior_Like <- (1 - Bprior) * log(Pconst + Pval - Pmin) + (1 - Aprior) * log(Pconst + Pmax - Pval) - (1 - Bprior) * log(Pconst + Pr - Pmin) - (1 - Aprior) * log(Pconst + Pmax - Pr)
		}
		if (Ptype == "Log_Norm") {
			Prior_Like <- 0.5 * ((log(Pval) - Pr)/Psd)^2
		}
		if (Ptype == "Log_Norm_w/biasadj") {
			if (Pmin > 0) {
				Prior_Like <- 0.5 * ((log(Pval) - Pr + 0.5 * Psd^2)/Psd)^2
			}
			else {
				warning("cannot do prior in log space for parm with min <=0.0")
			}
		}
		if (Ptype == "Gamma") {
			scale <- (Psd^2)/Pr
			shape <- Pr/scale
			Prior_Like <- -1 * (-shape * log(scale) - lgamma(shape) + (shape - 1) * log(Pval) - Pval/scale)
		}
		if (Ptype == "F") {
			Prior_Like <- rep(0, length(Pval))
		}
		if (is.null(Prior_Like)) {
			warning("Problem calculating prior. The prior type doesn't match ", 
				"any of the options in the SSplotPars function.\n", 
				"Ptype: ", Ptype)
		}
		return(Prior_Like)
	}
	pngfun <- function(file, caption=NA, lang="e") {
		createFdir(lang, dir=plotdir)
		changeLangOpts(L=lang)
		fout = switch(lang, 'e' = file.path(plotdir, file), 'f' = file.path(plotdir,"french", file) )
		clearFiles(fout)
		png(filename=fout, width=PIN[1], height=PIN[2], units=punits, res=res, pointsize=ptsize)
		plotinfo <- rbind(plotinfo, data.frame(file=file, caption=caption))
		return(plotinfo)
	}
	plotinfo <- NULL
	if (!"parameters" %in% names(replist)) {
		stop("'replist' input needs to be a list created by the SS_output function")
	}
	if (is.null(plotdir)) {
		plotdir <- replist$inputs$dir
	}
	if (print & add) {
		stop("Inputs 'print' and 'add' can't both be TRUE")
	}
	if (print & plot) {
		warning("Inputs 'print' and 'plot' can't both be TRUE\n", 
			"changing to 'plot=FALSE'")
	}
	parameters <- replist$parameters
	parameters$Label = convPN(parameters$Label)
	allnames <- parameters$Label[!is.na(parameters$Active_Cnt)]
	if (!is.null(strings)) {
		goodnames <- NULL
		if (exact) {
			goodnames <- allnames[allnames %in% strings]
#browser();return()
		}
		else {
			for (i in 1:length(strings)) {
				goodnames <- c(goodnames, grep(strings[i], allnames, fixed=TRUE, value=TRUE))
			}
		}
		goodnames <- unique(goodnames)
		if (verbose) {
			message("Active parameters matching input vector 'strings':")
			print(goodnames)
		}
		if (length(goodnames) == 0) {
			warning("No active parameters match input vector 'strings'.")
			return()
		}
	}
	else {
		goodnames <- allnames
		if (length(goodnames) == 0) {
			warning("No active parameters.")
			return()
		}
	}
	skip <- grep("Impl_err_", goodnames)
	if (length(skip) > 0) {
		goodnames <- goodnames[-skip]
		message("Skipping 'Impl_err_' parameters which don't have bounds reported")
	}
	skip <- grep("F_fleet_", goodnames)
	if (length(skip) > 0) {
		goodnames <- goodnames[-skip]
		message("Skipping 'F_fleet_' parameters which aren't yet supported by this function")
	}
	if (!showdev) {
		devnames <- c("RecrDev", "InitAge", "ForeRecr", "DEVadd", "DEVmult", "DEVrwalk", "DEV_MR_rwalk", "ARDEV")
		devrows <- NULL
		for (iname in 1:length(devnames)) {
			devrows <- unique(c(devrows, grep(devnames[iname], goodnames)))
		}
		if (length(devrows) > 0) {
			goodnames <- goodnames[-devrows]
			if (verbose) {
				message("Excluding ", length(devrows), " deviation parameters because input 'showdev'=FALSE")
			}
			if (length(goodnames) == 0) {
				message("no parameters to plot")
				return()
			}
		}
	}
	else {
		if (length(grep("rwalk", x=goodnames)) > 0 | length(grep("DEVadd", x=goodnames)) > 0 |
			length(grep("DEVmult", x=goodnames)) > 0 | length(grep("ARDEV", x=goodnames)) > 0) {
			warning("Parameter deviates are not fully implemented in this function.\n", 
				"Prior and bounds unavailable so these are skipped and\n", 
				"fitrange is set to TRUE for those parameters.")
		}
	}
	stds <- parameters$Parm_StDev[parameters$Label %in% goodnames]
	if (showmle & (all(is.na(stds)) || min(stds, na.rm=TRUE) <= 0)) {
		message("Some parameters have std. dev. values in Report.sso equal to 0.\n", 
			"  Asymptotic uncertainty estimates will not be shown.\n", 
			"  Try re-running the model with the Hessian but no MCMC.")
	}
	npars <- length(goodnames)
	if (is.null(strings) & verbose) {
		messagetext <- paste0("Plotting distributions for ", npars, " estimated parameters.")
		if (!showdev) {
			messagetext <- gsub(pattern=".", replacement=" (deviations not included).", x=messagetext, fixed=TRUE)
		}
		message(messagetext)
	}
	npages <- ceiling(npars/(nrows * ncols))
	plotPars.fn <- function(l="e") {
		if (!add) {
			plot(0, type="n", xlim=xlim2, ylim=ylim2, xaxs=xaxs, yaxs="i", xlab="", ylab="", main="", cex.main=1, axes=FALSE)
			axis(1)
			addLabel(0.975,0.975, linguaFranca(gsub("_"," ",header),l), adj=c(1,1), cex=1.2, col="darkgreen", font=2)
		}
		colval <- colvec[4]
		if (showpost & goodpost) {
			plot(posthist, add=TRUE, freq=FALSE, col=colval, border="green3", lwd=0.5)
			abline(v=postmedian, col=colvec[5], lwd=2, lty=ltyvec[3])
		}
		if (!isdev & showprior) {
			lines(x, prior, lwd=2, lty=ltyvec[2])
		}
		if (showmle) {
			if (!is.na(parsd) && parsd > 0) {
				lines(x, mle, col=colvec[1], lwd=1, lty=ltyvec[1])
				lines(rep(finalval, 2), c(0, dnorm(finalval, finalval, parsd) * mlescale), col=colvec[1], lty=ltyvec[1])
			}
			else {
				abline(v=finalval, col=colvec[1], lty=ltyvec[1])
			}
		}
		if (showinit) {
			par(xpd=NA)
			points(initval, -0.02 * ymax, col=colvec[2], pch=17, cex=1.2)
			par(xpd=FALSE)
		}
		box()
		if (max(par("mfg")[1:2]) == 1) {
			mtext(linguaFranca(xlab,l), side=1, line=0.5, outer=TRUE)
			mtext(linguaFranca(ylab,l), side=2, line=0.5, outer=TRUE)
			if (showlegend) {
				showvec <- c(showprior, showmle, showpost, showpost, showinit)
				legend("topleft", cex=0.9, bty="n", pch=c(NA,NA,15,NA,17)[showvec], lty=c(ltyvec[2], 
					ltyvec[1], NA, ltyvec[3], NA)[showvec], lwd=c(2,1,NA,2,NA)[showvec], col=c(colvec[3], 
					colvec[1], colvec[4], colvec[5], colvec[2])[showvec], pt.cex=c(1,1,2,1,1)[showvec],
					legend=linguaFranca(c("prior", "max. likelihood", "posterior", "posterior median", "initial value")[showvec],l)
				)
			}
		}
	}
	if (debug) {
		message("Making plots of parameters:")
	}
	if (plot & !add) {
		par(mfrow=c(nrows, ncols), mar=c(1,0.5,1,0), oma=c(2,2,0,1), mgp=c(2,0.5,0))
	}
	for (ipar in 1:npars) {
		ipage <- floor(1 + (ipar - 1)/(nrows * ncols - 1))
		parname <- goodnames[ipar]
		if (debug) {
			message("	", parname)
		}
		parline <- parameters[parameters$Label == parname, ]
		initval <- parline$Init
		finalval <- parline$Value
		parsd <- parline$Parm_StDev
		Pmin <- parline$Min
		Pmax <- parline$Max
		Ptype <- parline$Pr_type
		Psd <- parline$Pr_SD
		Pr <- parline$Prior
#browser();return()
		if (is.na(Ptype) || Ptype == "dev") {
			Ptype <- "Normal"
			Pr <- 0
		}
		if (any(sapply(X=c("RecrDev", "InitAge", "ForeRecr"), 
			FUN=grepl, parname))) {
			Psd <- parameters$Value[parameters$Label == "SR_sigmaR"]
		}
		isdev <- FALSE
		if (length(grep("DEVrwalk", x=parname)) > 0 | length(grep("DEVadd", x=parname)) > 0 |
			length(grep("DEVmult", x=parname)) > 0 | length(grep("ARDEV", x=parname)) > 0) {
			initval <- 0
			isdev <- TRUE
		}
		ymax <- 0
		xmin <- NULL
		xmax <- NULL
		if (!isdev) {
			x <- seq(Pmin, Pmax, length=5000)
			negL_prior <- GetPrior(Ptype=Ptype, Pmin=Pmin, Pmax=Pmax, Pr=Pr, Psd=Psd, Pval=x)
			prior <- exp(-1 * negL_prior)
			if (length(prior) == 0) {
				prior <- rep(NA, length(x))
			}
		}
		else {
			x <- finalval + seq(-4 * parsd, 4 * parsd, length=5000)
		}
		if (!isdev & showprior) {
			prior <- prior/(sum(prior) * mean(diff(x)))
			ymax <- max(ymax, max(prior), na.rm=TRUE)
		}
		if (showmle) {
			if (!is.na(parsd) && parsd > 0) {
				mle <- dnorm(x, finalval, parsd)
				mlescale <- 1/(sum(mle) * mean(diff(x)))
				mle <- mle * mlescale
				ymax <- max(ymax, max(mle), na.rm=TRUE)
				xmin <- qnorm(0.001, finalval, parsd)
				xmax <- qnorm(0.999, finalval, parsd)
			}
			else {
				xmin <- xmax <- finalval
			}
		}
		mcmc <- replist$mcmc
		if (showpost && is.null(mcmc)) {
			message("$mcmc not found in input 'replist', changing input to 'showpost=FALSE'")
			showpost <- FALSE
		}
		if (showpost && length(mcmc) < 20) {
			message("mcmc output has fewer than 20 rows, changing input to 'showpost=FALSE'")
			showpost <- FALSE
		}
		goodpost <- FALSE
		if (showpost) {
			postparname <- parname
			if (substring(parname, 1, 1) == "_") {
				postparname <- paste0("X", postparname)
			}
			jpar <- (1:ncol(mcmc))[names(mcmc) == postparname]
			if (length(jpar) == 1) {
				post <- mcmc[, jpar]
				xmin <- min(xmin, quantile(post, 0.001))
				xmax <- max(xmax, quantile(post, 0.999))
				goodpost <- TRUE
			}
			else {
				warning("parameter '", postparname, "', not found in posteriors.")
			}
		}
		if (is.null(xlim)) {
			if (fitrange & ((!is.na(parsd) && parsd != 0) | showpost)) {
				if (fitnudge<=0) {
					xmin <- max(Pmin, xmin, na.rm=TRUE)
					xmax <- min(Pmax, xmax, na.rm=TRUE)
				} else {
					xspan = abs(diff(c(xmin,xmax)))
					xmin  = xmin - (fitnudge * xspan)
					xmax  = xmax + (fitnudge * xspan)
#browser();return()
				}
			}
			else {
				if (!isdev) 
				  xmin <- Pmin
				if (!isdev) 
				  xmax <- Pmax
			}
			xlim2 <- c(xmin, xmax)
		}
		else {
			xlim2 <- xlim
		}
		if (showpost & goodpost) {
			jpar <- (1:ncol(mcmc))[names(mcmc) == postparname]
			post <- mcmc[, jpar]
			breakvec <- seq(xmin, xmax, length=50)
			if (min(breakvec) > min(post)) 
				breakvec <- c(min(post), breakvec)
			if (max(breakvec) < max(post)) 
				breakvec <- c(breakvec, max(post))
			posthist <- hist(post, plot=FALSE, breaks=breakvec)
			postmedian <- median(post)
			ymax <- max(ymax, max(posthist$density), na.rm=FALSE)
		}
		if (is.null(newheaders)) {
			header <- parname
		}
		else {
			header <- newheaders[ipar]
		}
		if (is.null(ylim)) {
			ylim2 <- c(0, 1.1 * ymax)
		}
		else {
			ylim2 <- ylim
		}
		if (print && (ipar%%(nrows * ncols) == 1 | nrows*ncols==1)) {
			caption <- "Parameter distribution plots"
			pagetext <- ""
			if (npages > 1) {
				pagetext <- paste("_page", ipage, sep="")
				caption <- paste(caption, " (plot ", ipage, " of ", npages, ").", sep="")
			}
			if (ipar == 1) {
				if (!showdev) {
					caption <- paste(caption, "<br>Deviation parameters are not included.")
				}
				if (length(grep("F_fleet_", allnames)) > 0) {
					caption <- paste(caption, "<br>F parameters are not included.")
				}
				if (fitrange) {
					caption <- paste(caption, "<br>Plotting range is scaled to fit parameter estimates.", 
					"Use fitrange=FALSE to use parameter bounds instead.")
				}
				else {
					caption <- paste(caption, "<br>Plotting range is equal to input limits on parameters.", 
					"Use fitrange=TRUE to scale the range to the estimates.")
				}
			}
			if (missing(outnam)) 
				outnam = "parameter_distributions"
			file <- paste0(outnam, pagetext, ".png", sep="")
			plotinfo <- pngfun(file=file, caption=caption, lang=lang)
#browser();return()
			par(mfrow=c(nrows, ncols), mar=c(1,0.5,1,0), oma=c(2,2,0,1), mgp=c(2,0.5,0))
		}
		if (print | plot) {
#browser();return()
			plotPars.fn(l=lang)
		}
		if (print && (ipar%%(nrows * ncols) == 0 | ipar == npars)) {
			dev.off(); eop()
		}
	}
	if (!is.null(plotinfo)) {
		plotinfo$category <- "Pars"
	}
	return(invisible(plotinfo))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.pars


## plotSS.pmcmc-------------------------2020-11-04
##  Plot parameter MCMC quantile boxplots.
##  Modified from PBSawatea's 'plotBmcmcPOP' function.
##  Input comes from the output file 'derived.parameters.sso' (use r4ss::SSgetMCMC)
## -----------------------------------------AME|RH
plotSS.pmcmc=function(obj, p=tcall(quants5), xyType="quantBox",
   lineType=c(3,2,1,2,3), refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, tcl.val=-0.2, yrs, y0=TRUE,
   pyrs=NULL, LRP=NULL, USR=NULL, catpol=NULL,
   yaxis.by, yLab="Recruitment", outnam, lang=c("e","f"), 
   ptypes="win", pngres=400, PIN=c(8,6), ...)
{
	# See plt.quantBio if want other xyTypes, as took out here:
	plt.qB <- function(obj, xyType="lines", new=TRUE, xLim, yLim, yrs, pyrs, LRP, USR, pvec, ...) {
		if ( new ) {
			plot(xLim, yLim, type="n", xlab=linguaFranca("Year",l), ylab=linguaFranca(yLab,l), ...)
			if (!y0) abline(h=0,col="gainsboro")
		}
		#yrs <- as.numeric(dimnames(obj)[[2]])
		#yrs <- as.numeric(gsub("[^[:digit:]]","",dimnames(result1)[[2]]))

		# Quantile boxplots - assumes five quantiles.
		if ( xyType=="quantBox" ) {
			if(!is.null(USR))
				abline(h=c(USR), col=c("slategray"), lwd=1, lty=5)
			delta <- 0.25 ## width of half-box
			# Draw the outer whiskers.
			segments(yrs, obj[1,], yrs,obj[5,], lty=1,col=1 )
			# Overlay the box.
			for ( i in 1:length(yrs) )
				#rect( yrs[i]-delta,obj[2,i], yrs[i]+delta, obj[4,i],... ) ## AME
				polygon(x=c(rep(yrs[i]-delta,2),rep(yrs[i]+delta,2)), y=obj[c(2,4,4,2),i], border="black", col="gainsboro", ...) ## RH (190620)
			# Add the median.
			segments( yrs-delta,obj[3,],yrs+delta,obj[3,],lty=1,col=1 )
			if (!is.null(pyrs)){
				pp = as.character(pyrs)
				segments(pyrs, obj[1,pp], pyrs, obj[5,pp], lty=1,col="red" )
				for ( i in (length(yrs)-length(pyrs)+1):length(yrs) )
					polygon(x=c(rep(yrs[i]-delta,2),rep(yrs[i]+delta,2)), y=obj[c(2,4,4,2),i], border="red", col="pink", ...) ## RH (190620)
				segments( pyrs-delta, obj[3,pp], pyrs+delta, obj[3,pp], lty=1, col="red" )
#browser();return()
			}
		}
		# Uncertainty envelope - assumes five quantiles.
		if ( xyType=="envelope" ) {
			if(!is.null(LRP) && !is.null(USR))
				abline(h=c(LRP,USR), col=c("red","green4"), lwd=1, lty=5)
			x  = setdiff(yrs,pyrs)
			xx = as.character(x)
			polygon(c(x,rev(x)), c(obj[1,xx],rev(obj[5,xx])), col=lucent("blue",0.05), border=FALSE)
			lines(x, obj[3,xx], lty=1, lwd=3, col="blue")
			lines(c(x,NA,x), c(obj[1,xx],NA,obj[5,xx]), lty=2, lwd=2, col="blue")
			lines(c(x,NA,x), c(obj[2,xx],NA,obj[4,xx]), lty=3, lwd=1, col="blue")

			if (!is.null(pyrs)){
				kcol = if (length(pvec)==1) "red" else c("green3","orange2","red")
				for (k in 1:length(pvec)) {
					x  = c(pyrs[1]-1,pyrs)
					xx = as.character(x)
					yy  = obj[,xx]
					pp  = pvec[[k]]
					yy[rownames(pp),colnames(pp)] = pp
					polygon(c(x,rev(x)), c(yy[1,xx],rev(yy[5,xx])), col=lucent(kcol[k],0.1), border=FALSE)
#browser();return()
					lines(x, yy[3,xx], lty=1, lwd=3, col=kcol[k])
					lines(c(x,NA,x), c(yy[1,xx],NA,yy[5,xx]), lty=2, lwd=2, col=kcol[k])
					lines(c(x,NA,x), c(yy[2,xx],NA,yy[4,xx]), lty=3, lwd=1, col=kcol[k])
				}
			}
		}
	}
	# Plot quantiles of biomass using the posterior densities.
	yrs1 = yrs2 = result1 = result2 = NULL
	if(is.null(p)) p = c(0.05, 0.25, 0.50, 0.75, 0.95)

	# Calculate the quantiles of the reconstructed biomass.
	result1 <- apply( obj, 2, quantile, probs=p )
	yrs1 <- as.numeric(gsub("[^[:digit:]]","",dimnames(result1)[[2]]))
	if (!missing(yrs)) {
		zyrs = is.element(yrs1,c(yrs,pyrs))
		yrs1 = yrs1[zyrs]
		result1 = result1[,zyrs]
	}
	if (!is.null(catpol)) {
		cpolnam = dimnames(catpol)[[3]]
		cpolcat = avgCP[1,1,cpolnam]
		projmat = apply( catpol, 2:3, quantile, probs=p )
		projvec = lapply(1:dim(projmat)[3], function(k) {projmat[,,k]})
		names(projvec) = cpolcat
	} else {
		projvec = list('AC'=result1[,intersect(colnames(result1),as.character(pyrs))])
	}
#browser();return()
	if ( is.null(yLim) )
		yLim <- c(ifelse(y0,0,min(result1)), max(result1))
	if ( is.null(xLim) )
		xLim=range(yrs1)
	colnames(result1) = yrs1

	createFdir(lang)
	if (missing(outnam))
		outnam=gsub("[[:space:]]+",".",tolower(yLab))
	fout = fout.e = outnam
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps")      postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			expandGraph(mar=c(3.5,3.75,1.5,0.75), mgp=c(2,0.5,0))
			plt.qB(result1, xLim=xLim, yLim=yLim, xyType=xyType, yrs=yrs1, pyrs=pyrs, LRP=LRP, USR=USR, pvec=projvec, ...)
#browser();return()
			axis(1, at=intersect(seq(1900,3000,5), xLim[1]:xLim[2]), tcl=tcl.val, labels=FALSE)
			if (missing(yaxis.by)) {
				yint = diff(seq(par()$yaxp[1], par()$yaxp[2], len=par()$yaxp[3]+1))[1]
				ytck = seq(par()$yaxp[1]-yint, par()$yaxp[2]+yint, yint/ifelse(par()$yaxp[3]>6,2,4))
				axis(2, at=ytck, tcl=tcl.val, labels=FALSE)
			} else 
				axis(2, at=seq(par()$yaxp[1], par()$yaxp[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
			#axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
			if (p %in% c("eps","png")) dev.off()
		}
	}; eop()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.pmcmc


## plotSS.profile-----------------------2021-04-28
##  Plot SS likelihood profiles
## ----------------------------------------r4ss|RH
plotSS.profile = function(mydir, string="NatM", p.vec=seq(0.040,0.065,0.005))
{
	## note: don't run this in your main directory; make a copy in case something goes wrong
	#mydir <- "C:/ss/Simple - Copy"
	
	## The following commands related to starter.ss could be done by hand:
	## Read starter file
	starter <- SS_readstarter(file.path(mydir, 'starter.ss'))
	## Change control file name in the starter file
	starter$ctlfile <- "control_modified.ss"
	## Make sure the prior likelihood is calculated for non-estimated quantities
	starter$prior_like <- 1
	# Write modified starter file
	SS_writestarter(starter, dir=mydir, overwrite=TRUE)
	
	## Vector of values to profile over
	#h.vec <- seq(0.3,0.9,.1)
	Nprofile <- length(p.vec)
	
	# Run SS_profile command
	file.copy(from="C:/Users/haighr/Files/Archive/Bat/ss.exe", to=mydir, overwrite=FALSE, copy.date=TRUE)
	profile <- SS_profile(dir=mydir, model="ss", masterctlfile="control.ss_new", newctlfile="control_modified.ss", string=string, profilevec=p.vec)

	# Read the output files (with names like Report1.sso, Report2.sso, etc.)
	profilemodels <- SSgetoutput(dirvec=mydir, keyvec=1:Nprofile)
	ttput(profilemodels)
	# Summarize output
	profilesummary <- SSsummarize(profilemodels)
	ttput(profilesummary)
	
	# OPTIONAL COMMANDS TO ADD MODEL WITH PROFILE PARAMETER ESTIMATED
	#MLEmodel <- SS_output("C:/ss/SSv3.24l_Dec5/Simple")
	#profilemodels$MLE <- MLEmodel
	#profilesummary <- SSsummarize(profilemodels)
	# END OPTIONAL COMMANDS
	
	# plot profile using summary created above
	SSplotProfile(profilesummary, profile.string="NatM_p_1_Fem_GP_1", profile.label="Natural mortality (M)", plot=F, print=T, plotdir=mydir)
	
	# make timeseries plots comparing models in profile
	SSplotComparisons(profilesummary,legendlabels=paste("M =",p.vec), plot=F, print=T, plotdir=mydir)
	return(list(profilesummary=profilesummary, profilemodels=profilemodels))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.profile


## plotSS.rdevs-------------------------2021-04-01
##  Plot reruitment deviations
## ----------------------------------------r4ss|RH
plotSS.rdevs = function (replist, subplots=1:3, plot=TRUE, print=FALSE, 
   add=FALSE, uncertainty=TRUE, minyr=-Inf, maxyr=Inf, 
   forecastplot=FALSE, col1="black", col2="blue", col3="green3", 
   col4="red", legendloc="topleft", labels=c("Year", "Asymptotic standard error estimate", 
   "Log recruitment deviation", "Bias adjustment fraction, 1 - stddev^2 / sigmaR^2"), 
   pwidth=8, pheight=6, punits="in", res=400, ptsize=10, 
   cex.main=1, plotdir="default", verbose=TRUE, outnam, lang="e") 
{
	oldpar = par(no.readonly=TRUE)
	fart = function(opar) { if (any("windows"%in%names(dev.list()))) par(opar); eop() }
	on.exit(fart(oldpar))
	changeLangOpts(L=lang)
	if (missing(outnam))
		outnam = NULL
	ttput(outnam)

	pngfun <- function(file, caption=NA, lang="e") {
		createFdir(lang, dir=plotdir)
		changeLangOpts(L=lang)
		fout = switch(lang, 'e' = file.path(plotdir, file), 'f' = file.path(plotdir,"french", file) )
		clearFiles(fout)
		png(filename=fout, width=pwidth, height=pheight, units=punits, res=res, pointsize=ptsize)
		plotinfo <- rbind(plotinfo, data.frame(file=file, caption=caption))
		return(plotinfo)
	}
	plotinfo <- NULL
	if (plotdir == "default") 
		plotdir <- replist$inputs$dir
	parameters <- replist$parameters
	recruit <- replist$recruit
	startyr <- replist$startyr
	endyr <- replist$endyr
	sigma_R_in <- replist$sigma_R_in
	recdevEarly <- parameters[substring(parameters$Label, 1, 13) %in% c("Early_RecrDev"), ]
	early_initage <- parameters[substring(parameters$Label, 1, 13) %in% c("Early_InitAge"), ]
	main_initage <- parameters[substring(parameters$Label, 1, 12) %in% c("Main_InitAge"), ]
	recdev <- parameters[substring(parameters$Label, 1, 12) %in% c("Main_RecrDev"), ]
	recdevFore <- parameters[substring(parameters$Label, 1, 8) == "ForeRecr", ]
	recdevLate <- parameters[substring(parameters$Label, 1, 12) == "Late_RecrDev", ]
	if (nrow(recdev) == 0 || max(recdev$Value) == 0) {
		if (verbose) 
			cat("Skipped SSplotrecdevs - no rec devs estimated\n")
	}
	else {
		if (nrow(recdev) > 0) {
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
			}
			else {
				recdevFore$Yr <- NULL
			}
			if (nrow(recdevLate) > 0) {
				recdevLate$Yr <- as.numeric(substring(recdevLate$Label, 14))
				recdevFore <- rbind(recdevLate, recdevFore)
			}
			Yr <- c(recdevEarly$Yr, recdev$Yr, recdevFore$Yr)
			if (forecastplot) {
				goodyrs <- ifelse(Yr >= minyr & Yr <= maxyr, TRUE, FALSE)
			}
			else {
				goodyrs <- Yr <= endyr + 1 & Yr >= minyr & Yr <= maxyr
			}
			xlim <- range(Yr[goodyrs], na.rm=TRUE)
			ylim <- range(c(recdevEarly$Value, recdev$Value, 
				recdevFore$Value)[goodyrs], na.rm=TRUE)
			recdevfunc <- function(uncertainty, l="e") {
				alldevs <- rbind(recdevEarly, recdev, recdevFore)[goodyrs,]
				colvec <- c(rep(col2, nrow(recdevEarly)), rep(col1,nrow(recdev)), rep(col2, nrow(recdevFore)))[goodyrs]
				val <- alldevs$Value
				Yr <- alldevs$Yr
				if (uncertainty) {
					std <- alldevs$Parm_StDev
					recdev_hi <- val + 1.96 * std
					recdev_lo <- val - 1.96 * std
					ylim <- range(recdev_hi, recdev_lo, na.rm=TRUE)
				}
				else {
					ylim <- range(val, na.rm=TRUE)
				}
				expandGraph(mfrow=c(1,1), mar=c(3.5,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
				plot(Yr, Yr, type="n", xlab=linguaFranca(labels[1],l), ylab=linguaFranca(labels[3],l), ylim=ylim, cex.axis=1.2, cex.lab=1.5)
				abline(h=0, col="grey")
				if (uncertainty)
					arrows(Yr, recdev_lo, Yr, recdev_hi, length=0.03, code=3, angle=90, lwd=1.2, col=colvec)
				lines(Yr, val, lty=3)
				points(Yr, val, pch=16, col=colvec)
				R0  = replist$parameters[grep("R0",rownames(replist$parameters)),"Value"]
				R0a = show0(round(R0,3),3)
				R0b = formatC(exp(R0),format="d",big.mark=switch(l,'e'=",",'f'=" "))
				addLabel(0.98,0.96,bquote(LN~italic(R)[0]== .(R0a)), cex=1, adj=c(1,0), col="blue")
				addLabel(0.98,0.93,bquote(italic(R)[0]== .(R0b)), cex=1, adj=c(1,0), col="blue")
				recdevs =  alldevs[,c("Yr","Value","Parm_StDev")]
				ttput(recdevs)
#browser();return()
			} ## end recdevfun
			if (uncertainty) {
				recdevfunc3 <- function(l="e") {
					#par(mar=par("mar")[c(1:3, 2)])
					ymax <- 1.1 * max(recdev$Parm_StDev, recdevEarly$Parm_StDev, recdevFore$Parm_StDev, sigma_R_in, na.rm=TRUE)
					expandGraph(mfrow=c(1,1), mar=c(3.5,3.5,2,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
					plot(recdev$Yr, recdev$Parm_StDev, xlab=linguaFranca(labels[1],l), main=linguaFranca("Recruitment deviation variance",l), cex.main=cex.main, ylab=linguaFranca(labels[2],l), xlim=xlim, ylim=c(0, ymax), type="b")
					if (nrow(recdevEarly) > 0)
						lines(recdevEarly$Yr, recdevEarly$Parm_StDev, type="b", col=col2)
					if (forecastplot & nrow(recdevFore) > 0)
						lines(recdevFore$Yr, recdevFore$Parm_StDev, type="b", col=col2)
					abline(h=0, col="grey")
					abline(h=sigma_R_in, col=col4)
				}
			}
			if (plot) {
				if (1 %in% subplots) 
					recdevfunc(uncertainty=FALSE, l=lang)
				if (uncertainty) {
					if (2 %in% subplots) 
						recdevfunc(uncertainty=TRUE, l=lang)
					if (3 %in% subplots) 
						recdevfunc3(l=lang)
				}
			}
			if (print) {
				if (1 %in% subplots) {
					if (!is.null(ttcall(outnam)))
						file = paste0(sub("\\.png$","",outnam), ".png")
					else
						file <- "recdevs1_points.png"
					caption <- "Recruitment deviations"
					plotinfo <- pngfun(file=file, caption=caption, lang=lang)
					recdevfunc(uncertainty=FALSE, l=lang)
					dev.off(); eop()
				}
				if (uncertainty) {
					if (2 %in% subplots) {
						if (!is.null(ttcall(outnam)))
							file = paste0(sub("\\.png$","",outnam), ".png")
						else
							file <- "recdevs2_withbars.png"
						caption <- "Recruitment deviations with 95% intervals"
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						recdevfunc(uncertainty=TRUE, l=lang)
						dev.off(); eop()
					}
					if (3 %in% subplots) {
						if (!is.null(ttcall(outnam)))
							file = paste0(sub("\\.png$","",outnam), ".png")
						else
							file <- "recdevs3_varcheck.png"
						caption <- paste("Recruitment deviations variance check.<br>", 
						"See later figure of transformed variance values for comparison", 
						"with bias adjustment settings in the model.")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						recdevfunc3(l=lang)
						dev.off(); eop()
					}
				}
			}
		}
	}
	if (!is.null(plotinfo)) 
		plotinfo$category <- "RecDev"
	return(invisible(plotinfo))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.rdevs


## plotSS.rmcmc-------------------------2021-09-03
##  Plot routine output for MCMCs.
## ---------------------------------------------RH
plotSS.rmcmc = function(mcmcObj, mpdObj, ptypes, lang, pngres=400, PIN=c(9,9))
{
	#so("linguaFranca.r")
	unpackList(mpdObj); unpackList(mcmcObj)

	## Pairs|density plot comparison among parameters
	plotSS.pairs(P.mpd, P.mcmc, ptypes=ptypes, lang=lang, pngres=pngres, PIN=c(10,10), type="image")

#return(); # browser();return()

	## Parameter posteriors and priors
	for (l in lang) {
		plotSS.pars(replist, nrows=P.rc[1], ncols=P.rc[2], plot=F, print=T, fitrange=T, fitnudge=0.5, showpost=T, strings=names(P.mpd), exact=T, plotdir=getwd(), lang=l, outnam="pdfParameters")
	}

	## Snail plots of Bt/Bmsy and ut/umsy
	for (p in ptypes) {
		plotSS.snail(BoverBmsy, UoverUmsy, yrs=modYrs, p=tcall(quants3)[c(1,3)], ngear=ngear, Cnames="Trawl Plus", assYrs=assYrs, currYear=currYr, ptypes=p, lang=lang, outnam="snail.defunct")
		colnames(BoverBmsy) = sub("^SSB_","",colnames(BoverBmsy))
		colnames(UoverUmsy) = sub("^SSB_","",colnames(UoverUmsy))
		plotSnail(BoverBmsy, UoverUmsy, yrs=modYrs, p=tcall(quants3)[c(1,3)], xLim=NULL, yLim=NULL, ngear=ngear, assYrs=assYrs, outs=F, Cnames="Trawl+", ptypes=ptypes, outnam="snail", lang=l)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.rmcmc


## plotSS.selex-------------------------2020-08-25
##  Plot selectivity curves with maturity ogive.
##  Source: R package 'r4ss' v.1.39.1
## ----------------------------------------r4ss|RH
plotSS.selex=function (replist, infotable=NULL, fleets="all", fleetnames="default", 
   sizefactors=c("Lsel"), agefactors=c("Asel", "Asel2"), 
   years="endyr", minyr=-Inf, maxyr=Inf, maxage=50, season=1, sexes="all", 
   selexlines=1:6, subplot=101, skipAgeSelex10=TRUE, xlim,
   plot=TRUE, print=FALSE, add=FALSE,
   labels=c("Length (cm)", "Age (yr)", "Year", "Selectivity", "Retention", "Discard mortality"), 
   col1="red", col2="blue", lwd=2, spacepoints=5, staggerpoints=1, 
   legendloc="bottomright", pwidth=7, pheight=7, punits="in", 
   res=400, ptsize=12, cex.main=1, showmain=TRUE, plotdir="default", 
   verbose=TRUE, debug=FALSE, sobj=NULL, lang=c("e","f"))
{
	oldpar = par(no.readonly=TRUE)
	fart = function(opar) { if (any("windows"%in%names(dev.list()))) par(opar) }
	on.exit(fart(oldpar))

	## drag race herstory for browsing wayward data
	sumtingwong = function(mess, debug, elvira) {
		.flush.cat(mess, "\n")
		if (debug) {
			bozeng = lenv()
			for (i in ls(elvira))
				eval(parse(text=paste0("tget(",i,",tenv=elvira)")))
			browser();return()
		} else
			return(invisible(mess))
	}
	if (!is.element(subplot, c(1:9,11:15,21:22,101:102))) {
		mess = c(" \n","Choose another subplot from:",
			"  1 = Selectivity at length for multiple fleets",
			"  2 = Selectivity at age for multiple fleets",
			"  3 = Surface plot of time-varing length by fleet",
			"  4 = Countour plot of time-varing length by fleet",
			"  5 = Surface plot of time-varying retention by fleet",
			"  6 = Countour plot of time-varying retention by fleet",
			"  7 = Surface plot of time-varying mortality by fleet",
			"  8 = Contour plot of time-varying mortality by fleet",
			"  9 = Length by fleet?",
			" 11 = Surface plot of time-varying ??? by fleet",
			" 12 = Surface plot of time-varying ??? by fleet",
			" 13 = Age ??? by fleet",
			" 14 = Age ??? by fleet",
			" 15 = Semi-parametric (2D-AR1) selectivity",
			" 21 = Contour plot of age-length by fleet",
			" 22 = Uncertainty unschertainty ???",
			"101 = Awatea-type selectivity plot with maturity curve underlay",
			"102 = Like 101 but add Awatea's estimated selectivity for comparison"
		)
		stop(paste0(mess, collapse="\n"))
	}
	infotable2 <- NULL
	nsexes <- replist$nsexes
	nseasons <- replist$nseasons
	nfleets <- replist$nfleets
	lbinspop <- replist$lbinspop
	nlbinspop <- replist$nlbinspop
	sizeselex <- replist$sizeselex
	ageselex <- replist$ageselex
	accuage <- replist$accuage
	startyr <- replist$startyr
	endyr <- replist$endyr
	FleetNames <- replist$FleetNames
	growdat <- replist$endgrowth
	growthCVtype <- replist$growthCVtype
	mainmorphs <- replist$mainmorphs
	nareas <- replist$nareas
	ngpatterns <- replist$ngpatterns
	derived_quants <- replist$derived_quants
	if (is.null(ageselex)) {
		message("Skipping age-based selectivity plots: no output available")
	}
	if (is.null(sizeselex)) {
		message("Skipping length-based selectivity plots: no output available")
	}
	pngfun <- function(file, caption=NA) {
		metas = c("\\", "/", ":", "*", "?", "\"", "<", ">", "|")
		meta  = paste0(c("[",metas,"]"),collapse="")
		fnam  = gsub("[_ ]+", "_", gsub(meta,"",file))
		.flush.cat("Figure file:", fnam, "\nsaved to:", plotdir, "\n")
		png(filename=file.path(plotdir, fnam), width=pwidth, 
			height=pheight, units=punits, res=res, pointsize=ptsize)
		plotinfo <- rbind(plotinfo, data.frame(file=fnam, caption=caption))
		return(plotinfo)
	}
	plotinfo <- NULL
	if (plotdir == "default") 
		plotdir <- replist$inputs$dir
	ians_blues <- c("white", "grey", "lightblue", "skyblue", "steelblue1", "slateblue", topo.colors(6), "blue", "blue2", "blue3", "blue4", "black")
	if (fleets[1] == "all") {
		fleets <- 1:nfleets
	}
	else {
		if (length(intersect(fleets, 1:nfleets)) != length(fleets)) {
			return("Input 'fleets' should be 'all' or a vector of values between 1 and nfleets.")
		}
	}
	if (fleetnames[1] == "default") 
		fleetnames <- FleetNames
	if (sexes[1] == "all") {
		sexes <- 1:nsexes
	}
	else {
		if (length(intersect(sexes, 1:nsexes)) != length(sexes)) {
			return("Input 'sexes' should be 'all' or a vector of values between 1 and nsexes.")
		}
	}
	if (years[1] == "endyr") 
		years <- endyr
	plotAllSel <- function(factor="Lsel",xlim) {
		if (factor %in% unique(sizeselex$Factor)) {
			agebased <- FALSE
			allselex <- sizeselex[sizeselex$Factor == factor & 
				sizeselex$Fleet %in% fleets & sizeselex$Sex %in% 
				sexes, ]
		}
		if (factor %in% unique(ageselex$Factor)) {
			agebased <- TRUE
			allselex <- ageselex[ageselex$Factor == factor & 
				ageselex$Seas == season & ageselex$Fleet %in% 
				fleets & ageselex$Sex %in% sexes, ]
		}
		if (!factor %in% unique(c(sizeselex$Factor, ageselex$Factor))) {
			cat("  Factor '", factor, "' not found in age- or length-based selectivity.\n", 
				"  This may be due to having 'detailed age-structured reports'\n", 
				"  turned off in the starter file.\n", sep="")
			return()
		}
		if (nrow(allselex) == 0) {
			cat("  combination of season, fleets, & sexes didn't produce any results\n")
			return()
		}
		time <- rep(FALSE, nfleets)
		for (ifleet in fleets) time[ifleet] <- any(apply(allselex[allselex$Fleet == 
			ifleet & allselex$Yr %in% (startyr:endyr), ], 2, 
			function(x) {
				any(x != x[1])
			}))
		if (any(time)) {
			if (length(years) > 1 & length(fleets) > 1) 
				cat("plot not yet configured to work well with multiple years and multiple fleets\n")
			inputyears <- years
			years <- NULL
			years2 <- NULL
			year_ranges <- NULL
			for (i in 1:length(inputyears)) {
				if (inputyears[i] >= startyr) {
					newyear <- min(endyr, allselex$Yr[allselex$Yr >= inputyears[i]])
					newyear2 <- max(startyr, allselex$Yr[allselex$Yr <= inputyears[i]])
					if (newyear2 <= newyear) {
						newyear_range <- paste(newyear2, "-", newyear, sep="")
					if (newyear == newyear2 & newyear > startyr - 3) 
						newyear_range <- newyear
					if (!newyear_range %in% year_ranges) {
						years <- c(years, newyear)
						years2 <- c(years2, newyear2)
						year_ranges <- c(year_ranges, newyear_range)
					}
					}
				}
			}
			if (all(years2 == startyr & years == endyr)) {
				years <- endyr
				years2 <- startyr
				year_ranges <- paste(startyr, "-", endyr, sep="")
			}
			bad <- rep(FALSE, length(years))
			for (i in 1:length(years)) {
				y <- years[i]
				y2 <- years2[i]
				if (sum(years == y) > 1) 
					bad[years == y & years2 == y] <- TRUE
				if (sum(years2 == y2) > 1) 
					bad[years == y2 & years2 == y2] <- TRUE
			}
			years <- years[!bad]
			years2 <- years2[!bad]
			year_ranges <- year_ranges[!bad]
			if ((startyr - 3) %in% inputyears) {
				years <- c(years, startyr - 3)
				year_ranges <- c(year_ranges, "Benchmarks")
			}
		}
		else {
			year_ranges <- ""
		}
#browser();return()
		allselex <- allselex2 <- allselex[allselex$Yr %in% years, ]
		if (nrow(allselex) == 0) {
			cat("No values found for this combination of years and factor\n")
			return()
		}
		Sex <- allselex$Sex
		if (!agebased) {
			allselex <- allselex[, -(1:5)]
			xlab <- labels[1]
		}
		if (agebased) {
			allselex <- allselex[, -(1:7)]
			xlab <- labels[2]
		}
		if (!is.null(infotable)) {
			infotable2 <- infotable
			good <- Sex %in% infotable$Sex
			allselex <- allselex[good, ]
			allselex2 <- allselex2[good, ]
			if (nrow(infotable2) != nrow(allselex)) {
				stop("Problem with input 'infotable'. Number of rows doesn't match.")
			}
		}
		else {
			infotable2 <- allselex2[c("Fleet", "Sex", "Yr")]
			infotable2$ifleet <- NA
			infotable2$FleetName <- fleetnames[infotable2$Fleet]
			infotable2$longname <- infotable2$FleetName
			for (i in 1:nrow(infotable2)) {
				infotable2$Yr_range[i] <- year_ranges[years == 
					infotable2$Yr[i]]
			}
			if (length(unique(infotable2$Yr)) > 1) {
				infotable2$longname <- paste(infotable2$FleetName, 
					infotable2$Yr_range)
			}
			twosex <- all(1:2 %in% infotable2$Sex) && any(allselex[infotable2$Sex == 
				1, ] != allselex[infotable2$Sex == 2, ])
			if (!twosex) {
				good <- infotable2$Sex == min(infotable2$Sex)
				allselex <- allselex[good, ]
				allselex2 <- allselex2[good, ]
				infotable2 <- infotable2[good, ]
			}
			else {
				infotable2$longname <- paste(infotable2$longname, 
					c("(f)", "(m)")[infotable2$Sex])
			}
			allfleets <- sort(unique(infotable2$Fleet))
			for (ifleet in 1:length(allfleets)) infotable2$ifleet[infotable2$Fleet == 
				allfleets[ifleet]] <- ifleet
			colvec <- rich.colors.short(length(allfleets))
			infotable2$col <- colvec[infotable2$ifleet]
			infotable2$lty <- 1
			infotable2$lwd <- lwd
			if (twosex) 
				infotable2$lty <- infotable2$Sex
			allyears <- sort(unique(infotable2$Yr))
			if (length(allyears) > 1) {
				for (iyear in 1:length(allyears)) infotable2$lty[infotable2$Yr == 
					allyears[iyear]] <- iyear
				if (twosex) 
					infotable2$lwd[infotable2$Sex == 2] <- lwd/2
			}
			infotable2$pch <- infotable2$ifleet%%25
		}
		main <- factor
		if (factor == "Lsel") 
			main <- paste("Length-based selectivity")
		if (factor == "Asel") 
			main <- paste("Age-based selectivity")
		if (factor == "Asel2") 
			main <- paste("Derived age-based from length-based selectivity")
		if (factor == "Ret") 
			main <- paste("Retention")
		if (length(fleets) > 1) 
			main <- paste(main, "by fleet")
		if (length(fleets) == 1) 
			main <- paste(main, "for", fleetnames[fleets])
		if (length(unique(infotable2$Yr)) == 1) {
			main <- paste(main, "in", unique(infotable2$Yr))
		}
		if (!showmain) 
			main <- NULL
		bins <- as.numeric(names(allselex))

#browser();return()
		if(missing(xlim)) xlim = range(bins)
		if (!add) 
			plot(0, xlim=xlim, ylim=c(0, 1), type="n", main=main, cex.main=cex.main, xlab=xlab, ylab=labels[4])
		abline(h=0, col="grey")
		abline(h=1, col="grey")
		matplot(x=bins, y=t(allselex), col=infotable2$col, lty=infotable2$lty, lwd=infotable2$lwd, type="l", add=TRUE)
		allselex2 <- allselex
		if (spacepoints > 0) {
			for (iline in 1:nrow(allselex)) allselex2[iline, 
				(1:ncol(allselex))%%spacepoints != (staggerpoints * 
					iline)%%spacepoints] <- NA
			matplot(x=bins, y=t(allselex2), col=infotable2$col, 
				lwd=infotable2$lwd, pch=infotable2$pch, type="p", 
				add=TRUE)
		}
		else {
			infotable2$pch <- NA
		}
		if (nrow(infotable2) > 1) 
			legend(legendloc, inset=c(0, 0.05), legend=infotable2$longname, 
				col=infotable2$col, seg.len=4, lty=infotable2$lty, 
				pch=infotable2$pch, lwd=infotable2$lwd, bty="n")
		return(infotable2)
	}
	if (1 %in% subplot & !is.null(sizeselex)) {
		for (ifactor in 1:length(sizefactors)) {
			if (plot) 
				infotable2 <- plotAllSel(factor=sizefactors[ifactor])
			if (print) {
				file <- paste("sel01_multiple_fleets_length", ifactor, ".png", sep="")
				caption <- "Selectivity at length for multiple fleets."
				plotinfo <- pngfun(file=file, caption=caption)
				infotable2 <- plotAllSel(factor="Lsel")
				dev.off()
			}
		}
	}
	if (2 %in% subplot & !is.null(ageselex)) {
		for (ifactor in 1:length(agefactors)) {
			factor <- agefactors[ifactor]
#browser();return()
			if (plot) 
				infotable2 <- plotAllSel(factor=factor,xlim=xlim)
			if (print) {
				file <- paste("sel02_multiple_fleets_age", ifactor, ".png", sep="")
				caption <- "Selectivity at age for multiple fleets."
				if (factor == "Asel2") 
					caption <- paste("Selectivity at age derived from selectivity at length for multiple fleets.")
				plotinfo <- pngfun(file=file, caption=caption)
				infotable2 <- plotAllSel(factor=factor)
				dev.off()
			}
		}
	}
	if (any(3:9 %in% subplot) & !is.null(sizeselex)) {
		for (i in fleets) {
			for (m in sexes) {
				if (m == 1 & nsexes == 1) 
					sextitle1 <- "Time-"
				if (m == 1 & nsexes == 2) 
					sextitle1 <- "Female time-"
				if (m == 2) 
					sextitle1 <- "Male time-"
				if (m == 1 & nsexes == 1) 
					sextitle2 <- "Ending"
				if (m == 1 & nsexes == 2) 
					sextitle2 <- "Female ending"
				if (m == 2) 
					sextitle2 <- "Male ending"
				intret <- sizeselex[sizeselex$Factor == "Ret" & sizeselex$Yr != startyr - 3 & sizeselex$Sex == m, ]
				intmort <- sizeselex[sizeselex$Factor == "Mort" & sizeselex$Yr != startyr - 3 & sizeselex$Sex == m, ]
				intkeep <- sizeselex[sizeselex$Factor == "Keep" & sizeselex$Yr != startyr - 3 & sizeselex$Sex == m, ]
				intdead <- sizeselex[sizeselex$Factor == "Dead" & sizeselex$Yr != startyr - 3 & sizeselex$Sex == m, ]
				intselex <- sizeselex[sizeselex$Factor == "Lsel" & sizeselex$Yr != startyr - 3 & sizeselex$Sex == m, ]
				plotselex <- intselex[intselex$Fleet == i, ]
				plotret <- intret[intret$Fleet == i, ]
				plotmort <- intmort[intmort$Fleet == i, ]
				## Selectivity
				time <- any(apply(plotselex[-c(1, nrow(plotselex)), -(1:5)], 2, function(x) { any(x != x[1]) }))
				if (time) {
					x <- lbinspop
					subset <- plotselex$Yr >= minyr & plotselex$Yr <= 
					maxyr
					y <- plotselex$Yr[subset]
					z <- plotselex[subset, -(1:5)]
					z <- matrix(as.numeric(as.matrix(z)), ncol=ncol(z))
					z <- t(z)
					main <- paste(sextitle1, "varying selectivity for ", 
					fleetnames[i], sep="")
					if (plot) {
					if (3 %in% subplot) 
						persp(x, y, z, col="white", xlab=labels[1], ylab=labels[3], zlab=labels[4], expand=0.5, box=TRUE, main=main, cex.main=cex.main, ticktype="detailed", phi=35, theta=-10)
					if (4 %in% subplot) 
						contour(x, y, z, nlevels=5, xlab=labels[1], 
						ylab=labels[3], main=main, cex.main=cex.main, 
						col=ians_blues, lwd=lwd)
					}
					if (print) {
					if (3 %in% subplot) {
						file <- paste("sel03_len_timevary_surf_flt", i, "sex", m, ".png", sep="")
						caption <- paste("Surface plot of", main)
						plotinfo <- pngfun(file=file, caption=caption)
						persp(x, y, z, col="white", xlab=labels[1], ylab=labels[3], zlab=labels[4], expand=0.5, box=TRUE, main=main, cex.main=cex.main, ticktype="detailed", phi=35, theta=-10)
						dev.off()
					}
					if (4 %in% subplot) {
						file <- paste("sel04_len_timevary_contour_flt", 
						i, "sex", m, ".png", sep="")
						caption <- paste("Countour plot of", main)
						plotinfo <- pngfun(file=file, caption=caption)
						contour(x, y, z, nlevels=5, xlab=labels[1], 
						ylab=labels[3], main=main, cex.main=cex.main, 
						col=ians_blues, lwd=lwd)
						dev.off()
					}
					}
				}
				else if (any(c(3:4)%in%subplot))
					sumtingwong("Selectivity not available for this model (not length-based)", debug, penv())
				## Retention
				time2 <- any(apply(plotret[-nrow(plotret), -(1:5)], 2, function(x) { any(x != x[1]) }))
				if (time2) {
					x <- lbinspop
					subset <- intret$Yr >= minyr & intret$Yr <= 
					maxyr
					y <- intret$Yr[subset & intret$Fleet == i]
					z <- intret[subset & intret$Fleet == i, -(1:5)]
					z <- matrix(as.numeric(as.matrix(z)), ncol=ncol(z))
					z <- t(z)
					main <- paste(sextitle1, "varying retention for ", 
					fleetnames[i], sep="")
					if (plot) {
					if (5 %in% subplot) 
						persp(x, y, z, col="white", xlab=labels[1], 
						ylab=labels[3], zlab=labels[5], expand=0.5, 
						box=TRUE, main=main, cex.main=cex.main, 
						ticktype="detailed", phi=35, theta=-10)
					if (6 %in% subplot) 
						contour(x, y, z, nlevels=5, xlab=labels[1], 
						ylab=labels[3], main=main, cex.main=cex.main, 
						col=ians_blues, lwd=lwd)
					}
					if (print) {
					if (5 %in% subplot) {
						file <- paste("sel05_timevary_ret_surf_flt", 
						i, "sex", m, ".png", sep="")
						caption <- paste("Surface plot of", main)
						plotinfo <- pngfun(file=file, caption=caption)
						persp(x, y, z, col="white", xlab=labels[1], 
						ylab=labels[3], zlab=labels[5], expand=0.5, 
						box=TRUE, main=main, cex.main=cex.main, 
						ticktype="detailed", phi=35, theta=-10)
						dev.off()
					}
					if (6 %in% subplot) {
						file <- paste("sel06_timevary_ret_contour_flt", 
						i, "sex", m, ".png", sep="")
						caption <- paste("Countour plot of", main)
						plotinfo <- pngfun(file=file, caption=caption)
						contour(x, y, z, nlevels=5, xlab=labels[1], 
						ylab=labels[3], main=main, cex.main=cex.main, 
						col=ians_blues, lwd=lwd)
						dev.off()
					}
					}
				}
				else if (any(c(5:6)%in%subplot))
					sumtingwong("Retention not available for this model (not length-based)", debug, penv())
				## Discard mortality
				time3 <- any(apply(plotmort[-nrow(plotmort), -(1:5)], 2, function(x) { any(x != x[1]) }))
				if (time3) {
					x <- lbinspop
					subset <- intmort$Yr >= minyr & intmort$Yr <= 
					maxyr
					y <- intmort$Yr[subset & intmort$Fleet == i]
					z <- intmort[subset & intmort$Fleet == i, -(1:5)]
					z <- matrix(as.numeric(as.matrix(z)), ncol=ncol(z))
					z <- t(z)
					main <- paste(sextitle1, "varying discard mortality for ", fleetnames[i], sep="")
					if (plot) {
					if (7 %in% subplot) 
						persp(x, y, z, col="white", xlab=labels[1], 
						ylab=labels[3], zlab=labels[6], expand=0.5, 
						box=TRUE, main=main, cex.main=cex.main, 
						ticktype="detailed", phi=35, theta=-10, 
						zlim=c(0, max(z)))
					if (8 %in% subplot) 
						contour(x, y, z, nlevels=5, xlab=labels[1], 
						ylab=labels[3], main=main, cex.main=cex.main, 
						col=ians_blues, lwd=lwd)
					}
					if (print) {
					if (7 %in% subplot) {
						file <- paste("sel07_timevary_mort_surf_flt", 
						i, "sex", m, ".png", sep="")
						caption <- paste("Surface plot of", main)
						plotinfo <- pngfun(file=file, caption=caption)
						persp(x, y, z, col="white", xlab=labels[1], 
						ylab=labels[3], zlab=labels[6], expand=0.5, 
						box=TRUE, main=main, cex.main=cex.main, 
						ticktype="detailed", phi=35, theta=-10)
						dev.off()
					}
					if (8 %in% subplot) {
						file <- paste("sel08_timevary_mort_contour_flt", 
						i, "sex", m, ".png", sep="")
						caption <- paste("Surface plot of", main)
						plotinfo <- pngfun(file=file, caption=caption)
						contour(x, y, z, nlevels=5, xlab=labels[1], 
						ylab=labels[3], main=main, cex.main=cex.main, 
						col=ians_blues, lwd=lwd)
						dev.off()
					}
					}
				}
				else if (any(c(7:8)%in%subplot))
					sumtingwong("Discard mortality not available for this model (not length-based)", debug, penv())
				## Ending-year selectivity
				endselex <- plotselex[plotselex$Yr == endyr, -(1:5)]
				plotret <- plotret[nrow(plotret), -(1:5)]
				ylab <- labels[4]
				bins <- as.numeric(names(endselex))
				vals <- as.numeric(paste(endselex))
				retvals <- as.numeric(plotret)
				main <- paste(sextitle2, " year selectivity for ", fleetnames[i], sep="")
				selfunc <- function() {
					intret2 <- intret[intret$Fleet == i, ]
					retchecktemp <- as.vector(unlist(intret2[1, ]))
					retcheck <- as.numeric(retchecktemp[6:length(retchecktemp)])
					if (is.na(sum(retcheck))) 
					retcheckuse <- 0
					if (!is.na(sum(retcheck))) 
					retcheckuse <- 1 - min(retcheck)
					if (!add) 
					plot(bins, vals, xlab=labels[1], ylim=c(0, 1), main=main, cex.main=cex.main, ylab="", type="n")
					abline(h=0, col="grey")
					abline(h=1, col="grey")
					if (1 %in% selexlines) 
					lines(bins, vals, type="o", col=col2, cex=1.1)
					if (retcheckuse > 0) {
						useret <- intret[intret$Fleet == i, ]
						usekeep <- intkeep[intkeep$Fleet == i, ]
						usemort <- intmort[intmort$Fleet == i, ]
						usedead <- intdead[intdead$Fleet == i, ]
						if (endyr %in% as.numeric(useret$Yr)) {
							useyr <- endyr
						}
						else {
							useyr <- max(as.numeric(useret$Yr))
						}
						plotret <- useret[useret$Yr == useyr, ]
						plotkeep <- usekeep[usekeep$Yr == useyr, ]
						plotmort <- usemort[usemort$Yr == useyr, ]
						plotdead <- usedead[usedead$Yr == useyr, ]
						plotdisc <- plotret
						plotdisc[-(1:5)] <- vals * (1 - plotret[, -(1:5)])
						if (2 %in% selexlines) {
							lines((as.numeric(as.vector(names(plotret)[-(1:5)]))), (as.numeric(as.character(plotret[1, -(1:5)]))), col="red", type="o", pch=3, cex=0.9)
							ylab <- paste(ylab, ", Retention", sep="")
						}
						if (3 %in% selexlines) {
							lines((as.numeric(as.vector(names(plotmort)[-(1:5)]))), (as.numeric(as.character(plotmort[1, -(1:5)]))), col="orange", type="o", 
							pch=4, cex=0.9)
							ylab <- paste(ylab, ", Mortality", sep="")
						}
						if (4 %in% selexlines) 
							lines((as.numeric(as.vector(names(plotkeep)[-(1:5)]))), (as.numeric(as.character(plotkeep[1, -(1:5)]))), col="purple", type="o", 
							pch=2, cex=0.9)
						if (5 %in% selexlines) 
							lines((as.numeric(as.vector(names(plotdead)[-(1:5)]))), (as.numeric(as.character(plotdead[1, -(1:5)]))), col="green3", type="o", 
							pch=5, cex=0.9)
						if (6 %in% selexlines) 
							lines((as.numeric(as.vector(names(plotdead)[-(1:5)]))), (as.numeric(as.character(plotdisc[1, -(1:5)]))), col="grey50", type="o", 
							pch=6, cex=0.9)
						legend(legendloc, inset=c(0, 0.05), bty="n", c(labels[4], labels[5], labels[6], "Keep=Sel*Ret", "Dead=Sel*(Ret+(1-Ret)*Mort)", "Discard=Sel*(1-Ret)")[selexlines], lty=1, col=c("blue", "red", "orange", "purple", "green3", "grey50")[selexlines], pch=c(1, 3, 4, 2, 5, 6)[selexlines], pt.cex=c(1.1, 0.9, 0.9, 0.9, 0.9, 0.9)[selexlines])
					}
					mtext(ylab, side=2, line=3)
				}
				if ((min(vals) < 1 & max(vals) > 0) | (!is.na(diff(range(retvals))) && diff(range(retvals)) != 0)) {
					if (9 %in% subplot) {
						if (plot) 
							selfunc()
						if (print) {
							file <- paste("sel09_len_flt", i, "sex", m, ".png", sep="")
							caption <- main
							plotinfo <- pngfun(file=file, caption=caption)
							selfunc()
							dev.off()
						}
					}
				}
				else if (any(c(9)%in%subplot))
					sumtingwong("Selectivity not available for this model (not length-based)", debug, penv())
			}
		}
	}
	if (any(11:14 %in% subplot) & !is.null(ageselex)) {
		ylab <- labels[4]
		for (facnum in 1) {
			factor <- c("Asel", "Asel2")[facnum]
			for (i in fleets) {
				for (m in sexes) {
					if (m == 1 & nsexes == 1) 
					sextitle1 <- "Time-"
					if (m == 1 & nsexes == 2) 
					sextitle1 <- "Female time-"
					if (m == 2) 
					sextitle1 <- "Male time-"
					if (m == 1 & nsexes == 1) 
					sextitle2 <- "Ending"
					if (m == 1 & nsexes == 2) 
					sextitle2 <- "Female ending"
					if (m == 2) 
					sextitle2 <- "Male ending"
					ageselexcols <- (1:ncol(ageselex))[names(ageselex) %in% 
					as.character(0:accuage)]
					plotageselex <- ageselex[ageselex$Factor == 
					factor & ageselex$Fleet == i & ageselex$Yr != 
					startyr - 3 & ageselex$Sex == m, ]
					time <- any(apply(plotageselex[-c(1, nrow(plotageselex)), 
					ageselexcols], 2, function(x) {
					any(x != x[1])
					}))
					if (time) {
					if ((min(as.numeric(as.vector(t(plotageselex[, 
						-(1:7)])))) < 1)) {
						subset <- as.numeric(plotageselex$Yr) >= 
						minyr & as.numeric(plotageselex$Yr) <= 
						maxyr
						x <- seq(0, accuage, by=1)
						y <- as.numeric(plotageselex$Yr)[subset]
						z <- plotageselex[subset, -(1:7)]
						z <- matrix(as.numeric(as.matrix(z)), ncol=ncol(z))
						z <- t(z)
						main <- paste(sextitle1, "varying selectivity for ", 
						fleetnames[i], sep="")
						if (plot) {
						if (11 %in% subplot) 
							persp(x, y, z, col="white", xlab=labels[2], 
							ylab=labels[3], zlab=ylab, expand=0.5, 
							box=TRUE, main=main, cex.main=cex.main, 
							ticktype="detailed", phi=35, theta=-10)
						if (12 %in% subplot) 
							contour(x, y, z, nlevels=5, xlab=labels[2], 
							main=main, cex.main=cex.main, 
							col=ians_blues, lwd=lwd)
						}
						if (print) {
						if (11 %in% subplot) {
							file <- paste("sel11_timevary_surf_flt", 
							i, "sex", m, ".png", sep="")
							caption <- main
							plotinfo <- pngfun(file=file, caption=caption)
							persp(x, y, z, col="white", xlab=labels[2], 
							ylab=labels[3], zlab=ylab, expand=0.5, 
							box=TRUE, main=main, cex.main=cex.main, 
							ticktype="detailed", phi=35, theta=-10)
							dev.off()
						}
						if (12 %in% subplot) {
							file <- paste("sel12_timevary_contour_flt", 
							i, "sex", m, ".png", sep="")
							caption <- main
							plotinfo <- pngfun(file=file, caption=caption)
							contour(x, y, z, nlevels=5, xlab=labels[2], 
							main=main, cex.main=cex.main, 
							col=ians_blues, lwd=lwd)
							dev.off()
						}
						}
						plotageselex2 <- plotageselex[plotageselex$Yr %in% 
						c(max(as.numeric(plotageselex$Yr))), 
						]
						plotageselex2 <- plotageselex2[, -(1:7)]
						main <- paste(sextitle2, " year selectivity for ", 
						fleetnames[i], sep="")
						endselfunc <- function() {
						if (!add) 
							plot((as.numeric(names(plotageselex2))), 
							(as.numeric(paste(c(plotageselex2)))), 
							xlab=labels[2], ylim=c(0, 1), 
							main=main, cex.main=cex.main, 
							ylab=ylab, type="n", col=col2, 
							cex=1.1)
						lines((as.numeric(names(plotageselex2))), 
							(as.numeric(paste(c(plotageselex2)))), 
							type="o", col=col2, cex=1.1)
						abline(h=0, col="grey")
						}
						if (13 %in% subplot) {
						if (plot) 
							endselfunc()
						if (print) {
							file <- paste("sel13_age_flt", i, "sex", 
							m, ".png", sep="")
							caption <- main
							plotinfo <- pngfun(file=file, caption=caption)
							endselfunc()
							dev.off()
						}
						}
					}
					}
					if (!time) {
					plotageselex <- plotageselex[plotageselex$Yr == 
						max(plotageselex$Yr), ]
					plotageselex <- plotageselex[, -(1:7)]
					vals <- as.numeric(paste(c(plotageselex)))
					doplot <- nrow(plotageselex) > 0 && diff(range(vals)) != 
						0
					if (doplot & skipAgeSelex10) {
						doplot <- !(vals[1] == 0 & all(vals[-1] == 
						1))
					}
					if (doplot) {
						main <- paste(sextitle2, " year selectivity for ", 
						fleetnames[i], sep="")
						endselfunc2 <- function() {
						if (!add) 
							plot((as.numeric(names(plotageselex))), 
							vals, xlab=labels[2], ylim=c(0, 
								1), main=main, cex.main=cex.main, 
							ylab=ylab, type="n")
						lines((as.numeric(names(plotageselex))), 
							vals, type="o", col=col2, cex=1.1)
						abline(h=0, col="grey")
						}
						if (14 %in% subplot) {
						if (plot) 
							endselfunc2()
						if (print) {
							file <- paste("sel14_age_flt", i, "sex", 
							m, ".png", sep="")
							caption <- main
							plotinfo <- pngfun(file=file, caption=caption)
							endselfunc2()
							dev.off()
						}
						}
					}
					}
				}
			}
			flush.console()
		}
	}
	if (15 %in% subplot & !is.null(replist$seldev_matrix)) {
		seldev_pars <- replist$seldev_pars
		seldev_matrix <- replist$seldev_matrix
		devcol.fn <- colorRampPalette(colors=c("red", "white", 
			"blue"))
		seldev_func <- function(m, mar=c(4.1, 4.1, 1, 1)) {
			bins <- as.numeric(colnames(m))
			years <- as.numeric(rownames(m))
			par(mar=mar)
			image(x=bins, y=years, z=t(m), col=devcol.fn(10), 
				xlab=names(dimnames(m))[2], ylab=names(dimnames(m))[1], 
				axes=FALSE, ylim=rev(range(years) + c(-0.5, 
					0.5)))
			axis(1, at=bins)
			axis(2, at=years, las=1)
			box()
		}
		for (imatrix in 1:length(seldev_matrix)) {
			label <- names(seldev_matrix)[imatrix]
			main <- gsub(pattern="_", replacement=" ", x=label)
			main <- gsub(pattern="seldevs", replacement="selectivity deviations", 
				x=main)
			if (plot) {
				seldev_func(m=seldev_matrix[[imatrix]], mar=c(5, 
					4, 4, 1) + 0.1)
				title(main=main)
			}
			if (print) {
				file=paste("sel15_", label, ".png", sep="")
				caption <- gsub(pattern="selectivity ", replacement="", 
					x=main)
				caption <- paste0(caption, " for semi-parametric (2D-AR1) selectivity. ", 
					"Blue value are positive deviations and red values negative. ", 
					"The matrix of values is available in the list created by ", 
					"<code>SS_output()</code> as <code>$seldev_matrix</code> which ", 
					"is a list with an element for each combination of fleet and length or ", 
					"age which uses the semi-parametric selectivity.")
				plotinfo <- pngfun(file=file, caption=caption)
				seldev_func(m=seldev_matrix[[imatrix]])
				dev.off()
			}
		}
	}
	if (21 %in% subplot & !is.null(ngpatterns) && ngpatterns == 1 & !is.null(growdat) & !is.null(sizeselex) & !is.null(ageselex) & all(!is.na(lbinspop))) {
		growdat <- growdat[growdat$Seas == season, ]
		if (nseasons > 1) 
			cat("Warning: plots showing growth curve with selectivity are using season", 
				season, "growth,\nwhich may not match the timing of the fishery.\n")
		growdatF <- growdat[growdat$Sex == 1 & growdat$Morph == mainmorphs[1], ]
		growdatF$Sd_Size <- growdatF$SD_Mid
		if (growthCVtype == "logSD=f(A)") {
			growdatF$high <- qlnorm(0.975, meanlog=log(growdatF$Len_Mid), sdlog=growdatF$Sd_Size)
			growdatF$low <- qlnorm(0.025, meanlog=log(growdatF$Len_Mid), sdlog=growdatF$Sd_Size)
		}
		else {
			growdatF$high <- qnorm(0.975, mean=growdatF$Len_Mid, sd=growdatF$Sd_Size)
			growdatF$low <- qnorm(0.025, mean=growdatF$Len_Mid, sd=growdatF$Sd_Size)
		}
		if (nsexes > 1) {
			growdatM <- growdat[growdat$Sex == 2 & growdat$Morph == mainmorphs[2], ]
			growdatM$Sd_Size <- growdatM$SD_Mid
			if (growthCVtype == "logSD=f(A)") {
				growdatM$high <- qlnorm(0.975, meanlog=log(growdatM$Len_Mid), sdlog=growdatM$Sd_Size)
				growdatM$low <- qlnorm(0.025, meanlog=log(growdatM$Len_Mid), sdlog=growdatM$Sd_Size)
			}
			else {
				growdatM$high <- qnorm(0.975, mean=growdatM$Len_Mid, sd=growdatM$Sd_Size)
				growdatM$low <- qnorm(0.025, mean=growdatM$Len_Mid, sd=growdatM$Sd_Size)
			}
		}
		xlab <- labels[2]
		ylab <- labels[1]
		zlab <- labels[4]
		for (i in fleets) {
			for (m in sexes) {
				if (m == 1 & nsexes == 1) 
					sextitle2 <- "Ending"
				if (m == 1 & nsexes == 2) 
					sextitle2 <- "Female ending"
				if (m == 2) 
					sextitle2 <- "Male ending"
				plotlenselex <- as.numeric(sizeselex[sizeselex$Factor == "Lsel" & sizeselex$Yr == endyr & sizeselex$Fleet == i & sizeselex$Sex == m, -(1:5)])
				if (any(plotlenselex != 1)) {
					plotageselex <- as.numeric(ageselex[ageselex$Factor == 
					"Asel" & ageselex$Yr == endyr & ageselex$Fleet == 
					i & ageselex$Sex == m, -(1:7)])
					x <- seq(0, accuage, by=1)
					y <- lbinspop
					z <- plotageselex %o% plotlenselex
					main <- paste(sextitle2, " year selectivity and growth for ", 
					fleetnames[i], sep="")
					agelenselcontour <- function() {
					contour(x, y, z, nlevels=5, xlab=xlab, ylab=ylab, main=main, cex.main=cex.main, col=ians_blues, lwd=lwd)
					if (m == 1) {
						lines(x, growdatF$Len_Mid, col="white", lwd=5)
						lines(x, growdatF$Len_Mid, col=col1, lwd=3)
						lines(x, growdatF$high, col="white", lwd=1, lty=1)
						lines(x, growdatF$high, col=col1, lwd=1, lty="dashed")
						lines(x, growdatF$low, col="white", lwd=1, lty=1)
						lines(x, growdatF$low, col=col1, lwd=1, lty="dashed")
					}
					if (m == 2) {
						lines(x, growdatM$Len_Mid, col="white", lwd=5)
						lines(x, growdatM$Len_Mid, col=col2, lwd=3)
						lines(x, growdatM$high, col="white", lwd=1, lty=1)
						lines(x, growdatM$high, col=col2, lwd=1, lty="dashed")
						lines(x, growdatM$low, col="white", lwd=1, lty=1)
						lines(x, growdatM$low, col=col2, lwd=1, lty="dashed")
					}
					}
					if (plot) {
					if (21 %in% subplot) 
						agelenselcontour()
					}
					if (print) {
					if (21 %in% subplot) {
						file=paste("sel21_agelen_contour_flt", 
						i, "sex", m, ".png", sep="")
						caption <- main
						plotinfo <- pngfun(file=file, caption=caption)
						agelenselcontour()
						dev.off()
					}
					}
				}
				else
					sumtingwong("Size selectivity was not used in the model", debug, penv())
			}
		}
	}
	if (22 %in% subplot) {
		rows <- grep("Selex_std", derived_quants$Label)
#browser();return()
		if (length(rows) > 0) {
			sel <- derived_quants[rows, ]
			names <- sel$Label
			splitnames <- strsplit(names, "_")
			namesDF <- as.data.frame(matrix(unlist(strsplit(names, "_")), ncol=6, byrow=T))
			sel$Fleet <- as.numeric(as.character(namesDF$V3))
			sel$Sex <- as.character(namesDF$V4)
			sel$agelen <- as.character(namesDF$V5)
			sel$bin <- as.numeric(as.character(namesDF$V6))
			sel$lower <- pmax(qnorm(0.025, mean=sel$Value, sd=sel$StdDev), 0)
			sel$upper <- pmin(qnorm(0.975, mean=sel$Value, sd=sel$StdDev), 1)
			i <- sel$Fleet[1]
			agelen <- sel$agelen[1]
			xlab <- labels[1:2][1 + (sel$agelen[1] == "A")]
			for (m in intersect(unique(sel$Sex), c("Fem", "Mal")[sexes])) {
				seltemp <- sel[sel$Sex == m, ]
				if (m == "Fem" & nsexes == 1) 
					sextitle3 <- ""
				if (m == "Fem" & nsexes == 2) 
					sextitle3 <- "females"
				if (m == "Mal") 
					sextitle3 <- "males"
				main <- paste("Uncertainty in selectivity for", 
					fleetnames[i], sextitle3)
				no0 <- seltemp$StdDev > 0.001
				if (FALSE) {
					if (agelen == "L") 
					plotselex <- sizeselex[sizeselex$Factor == "Lsel" & ageselex$Fleet == i & sizeselex$Sex == m, ]
					if (agelen == "A") 
					plotselex <- ageselex[ageselex$Factor == "Asel" & ageselex$Fleet == i & ageselex$Sex == m, ]
				}
				plot_extra_selex_SD <- function() {
					if (!add) 
					plot(seltemp$bin, seltemp$Value, xlab=xlab, 
						ylim=c(0, 1), main=main, cex.main=cex.main, 
						ylab=labels[4], type="n", col=col2, 
						cex=1.1, xlim=c(0, max(seltemp$bin)))
					lines(seltemp$bin, seltemp$Value, xlab=xlab, 
					ylim=c(0, 1), main=main, cex.main=cex.main, 
					ylab=labels[4], type="o", col=col2, 
					cex=1.1, xlim=c(0, max(seltemp$bin)))
					arrows(x0=seltemp$bin[no0], y0=seltemp$lower[no0], 
					x1=seltemp$bin[no0], y1=seltemp$upper[no0], 
					length=0.01, angle=90, code=3, col=col2)
					abline(h=0, col="grey")
				}
				if (plot) 
					plot_extra_selex_SD()
				if (print) {
					file <- paste("sel22_uncertainty", "sex", m, ".png", sep="")
					caption <- main
					plotinfo <- pngfun(file=file, caption=caption)
					plot_extra_selex_SD()
					dev.off()
				}
			}
		}
		else 
			sumtingwong("'Selex_std' not found in object 'derived_quants'", debug, penv())
	}
	if (any(c(101:102) %in% subplot) & !is.null(ageselex)) {
		control.file=replist$Control_File; tput(control.file)
		if (!exists("species.name")) species.name=tcall(species.name)
		if (is.null(species.name))   species.name="Anonymous Fish"
		ptypes = c("win","png")[c(plot,print)]
		pngres = res
		#tget(ptypes); tget(pngres); tget(lang)
		adbase   = ageselex[is.element(ageselex$Factor,"Asel") & is.element(ageselex$Yr,endyr) & is.element(ageselex$Sex,sexes) & is.element(ageselex$Fleet,fleets),]
		SexNames = if (nsexes==1) "Unisex" else c("Female","Male")
		adbase$Index = paste0(FleetNames[adbase$Fleet],": ",SexNames[adbase$Sex])
		plt.selectivity(adbase, mainTitle=species.name, ptypes=ptypes, pngres=pngres, lang=lang, sobj=sobj, maxage=maxage)
#browser();return()
	}
	if (!is.null(plotinfo)) 
		plotinfo$category <- "Sel"
	return(invisible(list(infotable=infotable2, plotinfo=plotinfo)))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.selex


## plotSS.stdres------------------------2020-08-26
## Plot standardised residuals -- three plots on one page.
## Modified from PBSawatea code in 'PBSscape.r'.
## -----------------------------------------------
plotSS.stdres = function(replist, kind="AGE", fleets="all",
   fleetnames="default", sexes="all", datonly=FALSE, aggregates_by_mkt = FALSE,
   labels=c("Length (cm)", "Age (yr)", "Year", "Observed sample size",
   "Effective sample size", "Proportion", "cm", "Frequency", "Weight", "Length",
   "(mt)", "(numbers x1000)", "Stdev (Age)", "Conditional AAL plot, ", "Size bin"),
   plot=TRUE, print=FALSE, type="Multinomial",
   ptypes="png", pngres=400, PIN=c(7,9), outnam, lang="e", ...)
{
	oldpar = par(no.readonly=TRUE)
	fart = function(opar) { if (any("windows"%in%names(dev.list()))) par(opar); eop() }
	on.exit(fart(oldpar))
#	changeLangOpts(L=lang)
	if (missing(outnam))
		outnam = NULL
	ttput(outnam)

#	pngfun <- function(file, caption=NA) {
#		png(filename=file.path(plotdir, file), width=pwidth, 
#			height=pheight, units=punits, res=res, pointsize=ptsize)
#		plotinfo <- rbind(plotinfo, data.frame(file=file, caption=caption))
#		return(plotinfo)
#	}
	plotinfo <- NULL
	SS_versionNumeric <- replist$SS_versionNumeric
	agedbase          <- replist$agedbase
	nfleets           <- replist$nfleets
	nseasons          <- replist$nseasons
	seasfracs         <- replist$seasfracs
	FleetNames        <- replist$FleetNames
	nsexes            <- replist$nsexes
	accuage           <- replist$accuage
	titles            <- NULL
	titlemkt          <- ""
	if (fleets[1] == "all") {
		fleets <- 1:nfleets
	}
	else {
		if (length(intersect(fleets, 1:nfleets)) != length(fleets)) {
			stop("Input 'fleets' should be 'all' or a vector of values between 1 and nfleets.")
		}
	}
	if (fleetnames[1] == "default") {
		fleetnames <- FleetNames
	}
#browser();return()
	if (sexes[1] == "all") {
		sexes <- 0:nsexes
		titlesex <- "(sexes combined)"
	}
	if (nsexes == 1) {
		sexes <- 0:nsexes
	}
	if (nsexes == 1 && length(sexes) > 1) {
		titlesex <- ""
		filesex <- ""
	}
	if (nsexes > 1 && length(sexes) == 1) {
		if (sexes == 0) {
			titlesex <- "(sexes combined)"
			filesex <- "sex0"
		}
		if (sexes == 1) {
			titlesex <- "(female)"
			filesex <- "sex1"
		}
		if (sexes == 2) {
			titlesex <- "(male)"
			filesex <- "sex2"
		}
	}
	if (length(sexes)==2 && all(sexes %in% c(1:2))) {
		titlesex = "(M+F)"
		filesex <- "sex12"
	}
#browser();return()
	if (!(kind %in% c("LEN", "SIZE", "AGE", "cond", "GSTAGE", "GSTLEN", "L@A", "W@A"))) {
		stop("Input 'kind' to SSplotComps needs to be one of the following:\n  ", "'LEN','SIZE','AGE','cond','GSTAGE','GSTLEN','L@A','W@A'.")
	}
	#titlesex <- ifelse(printsex, titlesex, "")
	if (kind == "LEN") {
		dbase_kind <- lendbase
		kindlab=labels[1]
		if (datonly) {
			filenamestart <- "comp_lendat_"
			titledata <- "Length comp data, "
		}
		else {
			filenamestart <- "comp_lenfit_"
			titledata <- "Length comps, "
		}
	}
	if (kind == "GSTLEN") {
		dbase_kind       <- ghostlendbase
		kindlab=labels[1]
		if (datonly) {
			filenamestart <- "comp_gstlendat_"
			titledata     <- "Ghost length comp data, "
		}
		else {
			filenamestart <- "comp_gstlenfit_"
			titledata     <- "Ghost length comps, "
		}
	}
	if (kind == "SIZE") {
		dbase_kind       <- sizedbase[sizedbase$method == sizemethod,]
		if (!is.null(sizebinlabs)) {
			kindlab <- labels[15]
			axis1 <- sort(unique(dbase_kind$Bin))
			if (length(sizebinlabs) == length(axis1)) {
				axis1labs <- sizebinlabs
			}
			else {
				axis1labs <- axis1
				warning("Input 'sizebinlabs' differs in length from the unique Bin\n", "  values associated with sizemethod=", sizemethod, ". Using bin values instead.")
			}
		}
		else {
			sizeunits <- unique(dbase_kind$units)
			if (length(sizeunits) > 1) {
				stop("!error with size units in generalized size comp plots:\n", "   more than one unit value per method.\n")
			}
			if (sizeunits %in% c("in", "cm")) {
				kindlab <- paste(labels[10], " (", sizeunits, ")", sep="")
			}
			if (sizeunits %in% c("lb", "kg")) {
				kindlab <- paste(labels[9], " (", sizeunits, ")", sep="")
			}
		}
		if (datonly) {
			filenamestart <- "comp_sizedat_"
			titledata <- "Size comp data, "
		}
		else {
			filenamestart <- "comp_sizefit_"
			titledata <- "Size comps, "
		}
		if (length(unique(sizedbase$method)) > 1) {
			filenamestart <- paste0(filenamestart, "method", sizemethod, "_")
			titledata <- paste0(titledata, " size method ", sizemethod, ", ")
		}
	}
	if (kind == "AGE") {
		dbase_kind <- agedbase
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_agedat_"
			titledata <- "Age comp data, "
		}
		else {
			filenamestart <- "comp_agefit_"
			titledata <- "Age comps, "
		}
	}
	if (kind == "cond") {
		dbase_kind <- condbase
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_condAALdat_"
			titledata <- "Conditional age-at-length data, "
		}
		else {
			filenamestart <- "comp_condAALfit_"
			titledata <- "Conditional age-at-length, "
		}
	}
	if (kind == "GSTAGE") {
		dbase_kind <- ghostagedbase
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_gstagedat_"
			titledata <- "Ghost age comp data, "
		}
		else {
			filenamestart <- "comp_gstagefit_"
			titledata <- "Ghost age comps, "
		}
	}
	if (kind == "GSTcond") {
		dbase_kind <- ghostagedbase
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_gstCAALdat_"
			titledata <- "Ghost conditional age-at-length data, "
		}
		else {
			filenamestart <- "comp_gstCAALfit_"
			titledata <- "Ghost conditional age-at-length comps, "
		}
	}
	if (kind == "L@A") {
		dbase_kind <- ladbase[ladbase$Nsamp_adj != 0, ]
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_LAAdat_"
			titledata <- "Mean length at age data, "
		}
		else {
			filenamestart <- "comp_LAAfit_"
			titledata <- "Mean length at age fit, "
		}
		dbase_kind$SD <- dbase_kind$Lbin_lo/dbase_kind$Nsamp_adj
	}
	if (kind == "W@A") {
		dbase_kind <- wadbase[wadbase$Nsamp_adj != 0, ]
		kindlab=labels[2]
		if (datonly) {
			filenamestart <- "comp_WAAdat_"
			titledata <- "Mean weight at age data, "
		}
		else {
			filenamestart <- "comp_WAAfit_"
			titledata <- "Mean weight at age fit, "
		}
		dbase_kind$SD <- dbase_kind$Lbin_lo/dbase_kind$Nsamp_adj
	}
	if (nrow(dbase_kind) > 0) {
		if (aggregates_by_mkt) {
			dbase_kind$Part_group <- dbase_kind$Part
		}
		else {
			dbase_kind$Part_group <- -1
		}
	}
	if (any(dbase_kind$SuprPer == "Sup" & dbase_kind$Used == "skip")) {
		cat("Note: removing super-period composition values labeled 'skip'\n", "   and designating super-period values with a '*'\n")
		dbase_kind <- dbase_kind[dbase_kind$SuprPer == "No" | dbase_kind$Used != "skip", ]
		dbase_kind$YrSeasName <- paste(dbase_kind$YrSeasName, ifelse(dbase_kind$SuprPer == "Sup", "*", ""), sep="")
	}
	ageerr_warning <- TRUE
	dbase_kind <- dbase_kind[dbase_kind$Fleet %in% fleets & dbase_kind$sex %in% sexes, ]
	if (nrow(dbase_kind)==0)
		stop ("Sum Ting Wong")

	if (type %in% c("Multinomial","Fournier","Coleraine")) {
		dbase_kind$Pearson_orig = dbase_kind$Pearson
		dbase_temp = calcStdRes(dbase_kind,type=type)
#browser();return()
		dbase_kind$Pearson = dbase_temp$stdRes
	}

	for (f in fleets) {
		fdbase = dbase_kind[is.element(dbase_kind$Fleet,f),]
		## Already has Pearson residuals computed
		if (is.null(ttcall(outnam)))
			outnam = "ageresFleet"
		fnam = paste0(outnam,f,sub("mf$","",gsub("\\+","",tolower(extract.between(titlesex,"(",")")))))
		fout = fout.e = fnam
		for (l in lang) {
			changeLangOpts(L=l)
			fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
#browser();return()
			for (p in ptypes) {
				if (print && p=="eps") {
					clearFiles(paste0(fout,".eps"))
					postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
				} else if (print && p=="png"){
					clearFiles(paste0(fout,".png"))
					png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				}
				par(mfrow=c(3,1), mai=c(0.45,0.3,0.1,0.1), omi=c(0,0.25,0.4,0), mgp=c(2,0.75,0))
				plt.ageResids(fdbase, main="", lang=l, ...)   ## by age class
				plt.yearResids(fdbase, lang=l, ...)           ## by year
				plt.cohortResids(fdbase, lang=l, ...)         ## by cohort (year of birth)
				#title = paste(toUpper(tolower(gsub("[[:punct:]]+", " ", FleetNames[f]))),titlesex)
				#title = paste(gsub("URVEY", "urvey", gsub("ISHERY", "ishery", gsub("YNOPTIC", "ynoptic", gsub("[[:punct:]]+", " ", FleetNames[f])))),titlesex)
				title = paste(gsub("_", " ", FleetNames[f]),titlesex)
				mtext(linguaFranca(title,l), side=3, outer=TRUE, line=0.25, cex=1.5)
				mtext(linguaFranca("Standardised Residuals",l), side=2, outer=TRUE, line=0, cex=1.5)
				if (print && p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.stdres


## plotSS.stock.recruit-----------------2021-08-26
## Plot stock-recruitment function (based on MPDs)
## ----------------------------------------r4ss|RH
## xLimSR and yLimSR fixed here for YMR to have Run 26 and 27 figs
##  on same scales. Use these first two values to scale to data:
# xLimSR =c(0, max(obj$B$SB))
# yLimSR=c(0, max(obj$B$R, na.rm=TRUE))
#xLimSR=c(0, max(c(max(obj$B$SB),45000)))   # so it draw bigger if necessary
#yLimSR=c(0, max(c(max(obj$B$R, na.rm=TRUE),55000)))

plotSS.stock.recruit = function (replist, subplot = 1:3, add = FALSE, plot = TRUE, print = FALSE, 
   xlim = NULL, ylim = NULL, labels = c("Spawning biomass (mt)", 
   "Recruitment (1,000s)", "Spawning output", expression(paste("Spawning output (relative to ", 
   italic(B)[0], ")")), expression(paste("Recruitment (relative to  ", 
   italic(R)[0], ")")), "Log recruitment deviation"), 
   bioscale = "default", plotdir = "default", pwidth = 6.5, 
   pheight = 6.5, punits = "in", res = 300, ptsize = 10, verbose = TRUE, 
   colvec = c("blue", "black", "black", gray(0, 0.7)), ltyvec = c(1,2,1,NA),
   ptcol = "default", legend = TRUE, legendloc = NULL, 
   minyr = "default", textmindev = 0.5, relative = FALSE, expected = TRUE, 
   estimated = TRUE, bias_adjusted = TRUE, show_env = TRUE, 
   virg = TRUE, init = TRUE, forecast = FALSE,
   lang=c("f","e"), ptypes="win", pngres=400, PIN=c(8,7), outnam="stockRecruit") 
{
	## Calculating projected R gets tricky
	## Stock-recruitment function
	srFun=function(spawners, h=h.mpd, R0=R0.mpd, B0=B0.mpd) {
		# to input a vector of spawners in year t-1 and calculate recruits in year t 
		4 * h * R0 * spawners / ( ( 1 - h) * B0 + (5 * h - 1) * spawners)
	}
	parameters = replist$parameters
	h.mpd = parameters[grep("steep",rownames(parameters)),"Value"]

	recruit  <- replist[["recruit"]]
	if (minyr == "default") 
		minyr <- min(recruit[["Yr"]])
	recruit    <- recruit[recruit[["era"]] %in% c("Early", "Main", "Fixed", "Late", ifelse(forecast, "Forecast", NA)) &  recruit[["Yr"]] >= minyr, ]
	timeseries <- replist[["timeseries"]]
	nsexes     <- replist[["nsexes"]]
	if (bioscale == "default") {
		if (nsexes == 1) 
			bioscale <- 0.5
		else bioscale <- 1
	}
	recruit[["spawn_bio"]] <- bioscale * recruit[["SpawnBio"]]
	timeseries[["SpawnBio"]] <- bioscale * timeseries[["SpawnBio"]]
	
	if (is.null(ylim)) {
		ylim <- c(0, 1.1 * max(recruit[["pred_recr"]], recruit[["exp_recr"]], recruit[["bias_adjusted"]]))
	}
	x <- recruit[["spawn_bio"]]
	if (is.null(xlim)) {
		xlim <- c(0, 1.1 * max(x))
	}
	show_env <- show_env & any(recruit[["with_env"]] != recruit[["exp_recr"]])
	B0 <- sum(timeseries[["SpawnBio"]][timeseries[["Era"]] == "VIRG"], na.rm = TRUE)
	B1 <- sum(timeseries[["SpawnBio"]][timeseries[["Era"]] == "INIT"], na.rm = TRUE)
	R0 <- sum(timeseries[["Recruit_0"]][timeseries[["Era"]] == "VIRG"], na.rm = TRUE)
	R1 <- sum(timeseries[["Recruit_0"]][timeseries[["Era"]] == "INIT"], na.rm = TRUE)
	if (B0 == 0) {
		B0 <- head(recruit[["spawn_bio"]][recruit[["spawn_bio"]] != 0], 1)
	}
	if (R0 == 0) {
		R0 <- head(recruit[["exp_recr"]][recruit[["exp_recr"]] != 0], 1)
	}
	if (B0 == B1 & R0 == R1) {
		init <- FALSE
	}
	if (relative) {
		x.mult <- 1/B0
		y.mult <- 1/R0
	}
	else {
		x.mult <- 1
		y.mult <- 1
	}
	B0.mpd = B0
	R0.mpd = R0
#browser();return()
	years    <- recruit[["Yr"]]
	B = data.frame(year=recruit[["Yr"]], SB = x * x.mult, R  = recruit[["pred_recr"]] * y.mult)
	xLimSR = c(0, 1.5*max(B$SB,na.rm=TRUE))   # so it draw bigger if necessary
	xxx    = (seq(0, xLimSR[2], length.out=100))
	yyy    = srFun(xxx)
	yLimSR=c(0, 1.1*max(c(yyy,B$R),na.rm=TRUE))

	fout = fout.e = outnam
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=6.5, height=4, horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=6.5, height=4)
			par(mfrow=c(1,1), mar=c(3.25,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			ylab = linguaFranca("Recruitment",l)
			plot(xxx, yyy, lwd=2, xlim=xLimSR, ylim=yLimSR, type="l",
				xlab=bquote(.(linguaFranca("Spawning biomass",l)) ~ italic(B)[italic(t)] ~ "(tonnes)" ~ .(linguaFranca("in year ",l)) ~ italic(t)),
				#ylab=bquote(.(linguaFranca("Recruitment",l)) ~ " " ~ italic(R)[italic(t)] ~ " (1000s)" ~ .(linguaFranca("in year ",l)) ~ italic(t) ~ "+1") ) 
				ylab=bquote(.(linguaFranca("Recruitment",l)) ~ "" ~ italic(R)[italic(t)] ~ " (1000s)" ~ .(linguaFranca("in year ",l)) ~ italic(t) ~ "")
			)
			#text(B[-length(years), "SB"], B[-1, "R"], labels=substring(as.character(years[-length(years)]), 3), cex=0.6, col="blue")
			## R = age-0 fish
			text(B[,"SB"], B[,"R"], labels=substring(as.character(years), 3), cex=0.6, col="blue")
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.stock.recruit


## plotSS.ts----------------------------2021-03-30
##  Plot SS time series (Bt, BtB0, Age-0 recruits)
## ----------------------------------------r4ss|RH
plotSS.ts = function (replist, subplot, add=FALSE, areas="all", areacols="default", 
   areanames="default", forecastplot=TRUE, uncertainty=TRUE, 
   bioscale=1, minyr=-Inf, maxyr=Inf, plot=TRUE, print=FALSE, 
   plotdir="default", verbose=TRUE, btarg="default", minbthresh="default", 
   xlab="Year", labels=NULL, punits="in", ptsize=10, cex.main=1, sobj=NULL,
   res=400, outnam, PIN=c(9,9), lang="e") 
 
{
	oldpar = par(no.readonly=TRUE)
	fart = function(opar) { if (any("windows"%in%names(dev.list()))) par(opar) }
	on.exit(fart(oldpar))

	if (missing(subplot)) 
		stop("'subplot' input required")
	if (length(subplot) > 1) 
		stop("function can only do 1 subplot at a time")
	pngfun <- function(file, caption=NA, lang="e") {
		metas = c("\\", "/", ":", "*", "?", "\"", "<", ">", "|")
		meta  = paste0(c("[",metas,"]"),collapse="")
		fnam  = gsub("[_ ]+", "_", gsub(meta,"",file))
		createFdir(lang, dir=plotdir)
		changeLangOpts(L=lang)
		fout = switch(lang, 'e' = file.path(plotdir, file), 'f' = file.path(plotdir,"french", file) )
		clearFiles(fout)
		##.flush.cat("Figure file:", fout, "\nsaved to:", plotdir, "\n")
		png(filename=fout, width=PIN[1], height=PIN[2], units=punits, res=res, pointsize=ptsize)
		plotinfo <- rbind(plotinfo, data.frame(file=fout, caption=caption))
		return(plotinfo)
	}
	plotinfo <- NULL
	if (is.null(labels)) {
		labels <- c("Total biomass (t)", "Total biomass (t) at beginning of season", "Summary biomass (t)", "Summary biomass (t) at beginning of season", "Spawning biomass (t)", "Bt / B0", "Spawning output", "Age-0 recruits (1000s)", "Fraction of total age-0 recruits", "Management target", "Minimum stock size threshold")
	}
	SS_versionshort <- replist$SS_versionshort
	timeseries      <- replist$timeseries
	nseasons        <- replist$nseasons
	spawnseas       <- replist$spawnseas
	birthseas       <- replist$birthseas
	startyr         <- replist$startyr
	endyr           <- replist$endyr
	nsexes          <- replist$nsexes
	nareas          <- replist$nareas
	derived_quants  <- replist$derived_quants
	seasfracs       <- replist$seasfracs
	B_ratio_denominator <- replist$B_ratio_denominator
	recruitment_dist    <- replist$recruitment_dist
	if (btarg == "default") 
		btarg <- replist$btarg
	if (minbthresh == "default") 
		minbthresh <- replist$minbthresh
	if (areacols[1] == "default") {
		areacols <- rich.colors.short(nareas)
		if (nareas == 3) {
			areacols <- c("blue", "red", "green3")
		}
		if (nareas > 3) {
			areacols <- rich.colors.short(nareas + 1)[-1]
		}
	}
	if (!is.null(birthseas)) {
		nbirthseas <- length(birthseas)
		seascols <- rich.colors.short(nbirthseas)
		if (nbirthseas > 2) 
			seascols <- rich.colors.short(nbirthseas + 1)[-1]
	}
	if (is.null(B_ratio_denominator)) 
		B_ratio_denominator <- 1
	if (plotdir == "default") {
		plotdir <- replist$inputs$dir
	}
	if (is.null(replist$SpawnOutputUnits) || is.na(replist$SpawnOutputUnits) || 
		replist$SpawnOutputUnits == "numbers") {
		labels[5] <- labels[7]
		labels[6] <- gsub("biomass", "output", labels[6])
	}
	if (areas[1] == "all") {
		areas <- 1:nareas
	}
	else {
		if (length(intersect(areas, 1:nareas)) != length(areas)) 
			stop("Input 'areas' should be 'all' or a vector of values between 1 and nareas.")
	}
	if (nareas > 1 & areanames[1] == "default") {
		areanames <- paste("area", 1:nareas)
	}
	ts <- timeseries
	if (nseasons > 1) {
		if (SS_versionshort == "SS-V3.11") {
			ts$YrSeas <- ts$Yr + (ts$Seas - 1)/nseasons
		}
		else {
			ts$YrSeas <- ts$Yr + seasfracs
		}
	}
	else {
		ts$YrSeas <- ts$Yr
	}
	ts <- ts[ts$YrSeas >= minyr & ts$YrSeas <= maxyr, ]
#browser();return()
	ts$period <- "time"
	ts$period[ts$Yr < startyr] <- "equilibria"
	ts$period[ts$Yr > endyr + 1] <- "fore"
	if (!forecastplot) 
		ts$period[ts$Yr > endyr + 1] <- "exclude"
	biofunc <- function(subplot, outnam, l="e") {
		plot1 <- ts$Area == 1 & ts$Era == "VIRG"
		plot2 <- ts$Area == 1 & ts$period == "time" & ts$Era != "VIRG"
		plot3 <- ts$Area == 1 & ts$period == "fore" & ts$Era != "VIRG"
		if (subplot %in% c(3, 6, 7, 9)) {
			plot1 <- ts$Area == 1 & ts$Era == "VIRG" & ts$Seas == spawnseas
			plot2 <- ts$Area == 1 & ts$period == "time" & ts$Era != "VIRG" & ts$Seas == spawnseas
			plot3 <- ts$Area == 1 & ts$period == "fore" & ts$Era != "VIRG" & ts$Seas == spawnseas
		}
		if (subplot %in% 1:3) {
			yvals <- ts[,"Bio_all",drop=FALSE]  ##ts$Bio_all  ## includes the forecast if SS has been run properly as an MPD
			ylab <- labels[1]
			if (subplot == 3) {
				ylab <- paste(labels[2], spawnseas)
			}
		}
		if (subplot %in% 4:6) {
			yvals <- ts[,"Bio_smry",drop=FALSE] ##ts$Bio_smry
			ylab <- labels[3]
			if (subplot == 6) {
				ylab <- paste(labels[4], spawnseas)
			}
		}
		if (subplot %in% 7:8) {
			yvals <- bioscale * ts[,"SpawnBio",drop=FALSE] ##ts$SpawnBio
			ylab <- labels[5]
		}
		if (subplot %in% 9:10) {
			B0 = ts$SpawnBio[!is.na(ts$SpawnBio)][1]
			yvals <- ts[,"SpawnBio",drop=FALSE] / B0 # $SpawnBio/ts$SpawnBio[!is.na(ts$SpawnBio)][1]
			ylab <- labels[6]
		}
		if (subplot %in% 11:15) {
			yvals <- ts[,"Recruit_0",drop=FALSE] ##$Recruit_0
			ylab <- labels[8]
			if (all(yvals[ts$Era == "VIRG",] == 0 & max(ts$Seas == 1))) {
				yvals[ts$Era == "VIRG",] <- derived_quants["Recr_Virgin", "Value"]
			}
			if (all(yvals[ts$Era == "INIT",] == 0 & max(ts$Seas == 1))) {
				yvals[ts$Era == "INIT",] <- derived_quants["Recr_Unfished", "Value"]
			}
		}
		if (subplot %in% c(13, 15)) 
			ylab <- labels[9]
		if (subplot %in% c(101)){
			ylab = "Biomass comparisons"
			yvals = ts[,c("Bio_all","SmryBio_SX:1_GP:1","SmryBio_SX:2_GP:1","SpawnBio")] * bioscale
#browser();return()
		}
		if (subplot %in% c(102)){
			ylab = "Vulnerable biomass"
			nfleets = length(.su(grep("F:",colnames(ts),value=T)))
			for (j in 1:nfleets) {
				ts[,paste0("u",j)]  = 1-exp(-ts[,paste0("F:_",j)])
				ts[,paste0("V",j)]  = ts[,paste0("sel(B):_",j)] / ts[,paste0("u",j)]
				ts[1,paste0("V",j)] = ts[2,paste0("V",j)]
			}
			yvals = ts[,c("Bio_all",paste0("V",1:nfleets))] * bioscale
#browser();return()
		}
		if (!is.element(subplot, c(1:15,101:102))) {
			stop("subplot should be a value from 1 to 15 (r4ss) or 101 to 102 (PBSsynth)")
		}
		main=ylab
		yrshift <- 0
		if (!is.null(birthseas) && max(birthseas) < spawnseas) {
			yrshift <- 1
		}
		if (!is.null(replist$recruitment_dist$recruit_dist) && "Age" %in% names(replist$recruitment_dist$recruit_dist)) {
			yrshift <- min(as.numeric(replist$recruitment_dist$recruit_dist$Age, na.rm=TRUE))
		}
		if (!is.null(birthseas) && nbirthseas > 1) {
			if (subplot == 11) {
				for (y in ts$Yr) {
					yvals[ts$Yr == y & ts$Seas == 1,] <- sum(yvals[ts$Yr == y,], na.rm=TRUE)
					yvals[ts$Yr == y & ts$Seas > 1,]  <- 0
				}
			}
			if (subplot == 15) {
				for (y in ts$Yr) {
					yvals[ts$Yr == y,] <- yvals[ts$Yr == y,]/sum(yvals[ts$Yr == y,], na.rm=TRUE)
				}
			}
			if (subplot %in% c(14, 15)) 
				main=paste(main, "by birth season")
		}
		if (nareas > 1) {
			if (subplot %in% c(2, 3, 5, 6, 8, 10, 12, 13)) {
				main=paste(main, "by area")
			}
			if (subplot %in% c(1, 4, 7, 11, 13)) {
				#yvals2 <- rep(NA, length(ts$YrSeas))
				yvals2 <- as.data.frame(array(NA, dim=c(length(ts$YrSeas),1), dimnames=list(rownames(ts),"sumting")))
				for (iyr in 1:nrow(yvals)) {
					y <- ts$YrSeas[iyr]
					yvals2[iyr,] <- sum(yvals[ts$YrSeas == y,])
				}
				if (subplot == 13) {
					yvals <- yvals/yvals2
				}
				else {
					yvals <- yvals2
				}
			}
			if (subplot == 9) {
				#yvals2 <- rep(NA, length(ts$YrSeas))
				yvals2 <- as.data.frame(array(NA, dim=c(length(ts$YrSeas),1), dimnames=list(rownames(ts),"sumting")))
				for (iyr in 1:nrow(yvals)) {
					y <- ts$YrSeas[iyr]
					yvals[iyr,] <- sum(ts$SpawnBio[ts$YrSeas == y])
				}
				yvals <- yvals/yvals[!is.na(yvals)][1]
			}
			ymax <- max(yvals, 1, na.rm=TRUE)
			if (subplot == 10) {
				for (iarea in 1:nareas) {
					yvals <- ts[ts$Area == iarea,"SpawnBio",drop=FALSE]/(ts$SpawnBio[ts$Area == iarea & ts$Seas == spawnseas][1])
					ymax <- max(yvals, na.rm=TRUE)
				}
			}
		}
		if (subplot == 10) {
			yvals[1,] <- NA
		}
		if (forecastplot) 
			main <- paste(main, "[ trajectory + forecast ]")
		if (uncertainty & subplot %in% c(7, 9, 11, 101)) {
			main <- paste(main, "with ~95% C.I.") #"with ~95% asymptotic intervals")
			if (!"SSB_Virgin" %in% derived_quants$Label) {
				warning("Skipping spawning biomass with uncertainty plot because 'SSB_Virgin' not in derived quantites.\n", 
					"  Try changing 'min yr for Spbio_sdreport' in starter file to -1.\n")
				stdtable <- NULL
			}
			else {
				if (subplot %in% c(7,101)) {
					stdtable <- derived_quants[grep("SSB_Virgin", derived_quants[, 1]):(grep("Recr_Virgin", derived_quants[, 1]) - 1), 1:3]
					stdtable$Yr <- substring(stdtable$Label, 5)
					stdtable$Yr[1:2] <- as.numeric(stdtable$Yr[3]) - (2:1) - yrshift
					stdtable$Yr <- as.numeric(stdtable$Yr)
				}
				if (subplot == 9) {
					stdtable <- derived_quants[substring(derived_quants$Label, 1, 6) == "Bratio", ]
					stdtable$Yr <- as.numeric(substring(stdtable$Label, 8))
					bioscale <- B_ratio_denominator
				}
				if (subplot == 11) {
					stdtable <- derived_quants[substring(derived_quants$Label, 1, 5) == "Recr_", ]
					stdtable <- stdtable[tolower(stdtable$Label) != "recr_unfished", ]
					stdtable$Yr <- substring(stdtable$Label, 6)
					stdtable$Yr[1:2] <- as.numeric(stdtable$Yr[3]) - (2:1)
					stdtable$Yr <- as.numeric(stdtable$Yr) + yrshift
					bioscale <- 1
				}
				v <- stdtable$Value * bioscale
				std <- stdtable$StdDev * bioscale
				if (subplot == 11) {
					stdtable$logint <- sqrt(log(1 + (std/v)^2))
					stdtable$lower <- qlnorm(p=0.025, meanlog=log(v), sdlog=stdtable$logint)
					stdtable$upper <- qlnorm(p=0.975, meanlog=log(v), sdlog=stdtable$logint)
				}
				else {
					stdtable$upper <- v + 1.96 * std
					stdtable$lower <- pmax(v - 1.96 * std, 0)
				}
				if (max(stdtable$Yr) < max(floor(ts$YrSeas))) {
					warning(max(stdtable$Yr), " is the last year with uncertainty in Report file, but ", 
					max(ts$YrSeas), " is last year of time series. ", "Consider changing starter file input for ", "'max yr for sdreport outputs' to -2")
				}
				stdtable <- stdtable[stdtable$Yr >= minyr & stdtable$Yr <= maxyr, ]
			}
		}
#browser();return()

		if (nareas == 1) {
			ymax <- max(yvals[plot1 | plot2 | plot3,], na.rm=TRUE)
		}
		if (subplot %in% c(13, 15)) 
			ymax <- 1
		if (uncertainty & subplot %in% c(7, 9, 11, 101)) {
			ymax <- max(ymax, stdtable$upper, na.rm=TRUE)
		}
		if (print) {
			if (missing(outnam)){
				filename <- main
				filename <- gsub(",", "", filename, fixed=TRUE)
				filename <- gsub("~", "", filename, fixed=TRUE)
				filename <- gsub("%", "", filename, fixed=TRUE)
				if (forecastplot) 
					filename <- paste(filename, "forecast")
				if (uncertainty & subplot %in% c(5, 7, 9)) 
					filename <- paste(filename, "intervals")
				filename <- paste0("ts", subplot, "_", filename, ".png")
				filename <- gsub(pattern=" ", replacement="_", x=filename, fixed=TRUE)
			} else
				filename = paste0(outnam, ".png")
			plotinfo <- pngfun(file=filename, caption=main, lang=lang)
		}
		ts$Yr[ts$Era == "VIRG"] <- ts$Yr[ts$Era == "VIRG"] + 1 + yrshift
		ts$YrSeas[ts$Era == "VIRG"] <- ts$YrSeas[ts$Era == "VIRG"] + 1 + yrshift
		sp = as.character(subplot)
		if (!add) {
			yrvals <- ts$YrSeas[plot1 | plot2 | plot3]
			xlim <- range(yrvals)
			if (!is.null(sobj)) {
				x2 = sobj$Year
				y2 = switch( sp, '7'=sobj$SB, '9'=sobj$SB/sobj$SB[1], '11'=sobj$R, rep(0,length(yrvals)) )
				use.y2 = !all(y2==0)
				if (use.y2)
					ymax  = max(ymax, max(y2))
			}
			else use.y2 = FALSE
			expandGraph(mfrow=c(1,1), mar=c(3,3,1,0.5))
			#plot(yrvals, yvals[plot1 | plot2 | plot3], type="n", xlab=xlab, xlim=xlim, ylim=c(0, 1.05 * ymax), yaxs="i", ylab=ylab, main=main, cex.main=cex.main, font.main=1)
			plot(0, 0, type="n", xlim=xlim, ylim=c(0, 1.05 * ymax), yaxs="i", xlab=linguaFranca(xlab,l), ylab=linguaFranca(ylab,l), main=NULL, cex.main=cex.main, cex.lab=1.4, font.main=1)
#browser();return()
			axis(1, at=intersect(seq(1900,3000,5),xlim[1]:xlim[2]), tcl=-0.2, labels=FALSE)
		}
		if (subplot %in% c(9, 10)) {
			addtarg <- function() {
				if (btarg > 0 & btarg < 1) {
					abline(h=btarg, col="green4",lty=2)
					text(max(startyr, minyr) + 1, btarg + 0.02, linguaFranca(ifelse(subplot==9,paste0(btarg,"B0"),labels[10]),l), adj=0, col="green4")
				}
				if (minbthresh > 0 & minbthresh < 1) {
					abline(h=minbthresh, col="red", lty=3)
					text(max(startyr, minyr) + 1, minbthresh + 0.02, linguaFranca(ifelse(subplot==9,paste0(minbthresh,"B0"),labels[10]),l), adj=0, col="red")
				}
			}
			addtarg()
		}
		if (subplot %in% 7:8) {
			addtarg <- function() {
				if (btarg > 1) {
					abline(h=btarg, col="red")
					text(max(startyr, minyr) + 4, btarg + 0.02 * diff(par()$usr[3:4]), linguaFranca(labels[10],l), adj=0)
				}
				if (minbthresh > 1) {
					abline(h=minbthresh, col="red")
					text(max(startyr, minyr) + 4, minbthresh + 0.02 * diff(par()$usr[3:4]), linguaFranca(labels[11],l), adj=0)
				}
			}
			addtarg()
		}
		if (subplot %in% 14:15) {
			for (iseas in 1:nbirthseas) {
				s <- birthseas[iseas]
				mycol <- seascols[iseas]
				mytype <- "o"
				plot1 <- ts$Seas == s & ts$Era == "VIRG"
				plot2 <- ts$Seas == s & ts$period == "time" & 
					ts$Era != "VIRG"
				plot3 <- ts$Seas == s & ts$period == "fore" & 
					ts$Era != "VIRG"
				points(ts$Yr[plot1], yvals[plot1,], pch=19, col=mycol)
				lines(ts$Yr[plot2], yvals[plot2,], type=mytype, col=mycol)
				points(ts$Yr[plot3], yvals[plot3,], pch=19, col=mycol)
			}
			#legend("topright", legend=linguaFranca(paste("Season", birthseas),l), lty=1, pch=1, col=seascols, bty="n")
			addLegend(0.975, 0.975, legend=linguaFranca(paste("Season", birthseas),l), lty=1, pch=1, col=seascols, bty="n", xjust=1, yjust=1)
		}
		else {
			if (subplot %in% c(1,4,7,9,11,14,15, 101:102)) 
				myareas <- 1
			else myareas <- areas
			for (iarea in myareas) {
				if (subplot == 10) {
					yvals <- ts[,"SpawnBio",drop=FALSE]/(ts$SpawnBio[ts$Area == iarea & ts$Seas == spawnseas][1])
				}
				if (subplot %in% c(3, 6, 7, 8, 9, 10)) {
					plot1 <- ts$Area == iarea & ts$Era == "VIRG" & ts$Seas == spawnseas
					plot2 <- ts$Area == iarea & ts$period == "time" & ts$Era != "VIRG" & ts$Seas == spawnseas
					plot3 <- ts$Area == iarea & ts$period == "fore" & ts$Era != "VIRG" & ts$Seas == spawnseas
				}
				else {
					yvals.vec = apply(yvals,1,sum,na.rm=T)
					plot1 <- yvals.vec > 0 & ts$Area == iarea & ts$Era == "VIRG"
					plot2 <- yvals.vec > 0 & ts$Area == iarea & ts$period == "time" & ts$Era != "VIRG"
					plot3 <- yvals.vec > 0 & ts$Area == iarea & ts$period == "fore" & ts$Era != "VIRG"
#browser();return()

				}
				if (subplot %in% 9:10) {
					plot1 <- NULL
					plot2[3] <- FALSE
				}
				mycol <- areacols[iarea]
				mytype <- "o"
				if (subplot == 11 & uncertainty) 
					mytype <- "p"
				if (!uncertainty) {
					legtxt = legcol = leglty = NULL
					if (use.y2) {
						if (subplot %in% c(11)) {
							x2 = c(x2[1]-1, x2[1:(length(x2)-1)])  ## PJS requested shifting Awatea age-1 recruits back 1 year to match SS age-0 recruits
						}
						lines(x2, y2, lty=1, lwd=2, col="yellow")
						lines(x2, y2, lty=2, lwd=1, col="orange")
						legtxt = c(legtxt, paste0("Awatea ", switch(sp, '7'="Bt", '9'="Bt/B0", '11'="age-1 recruits shifted", "sumtingwong")))
						leglty = c(leglty, 2)
						legcol = c(legcol, "darkgoldenrod1")
					}
					points(ts$YrSeas[plot1], yvals[plot1,], pch=19, col=mycol)
					lines(ts$YrSeas[plot2], yvals[plot2,], type=mytype, col=mycol)
					points(ts$YrSeas[plot3], yvals[plot3,], pch=19, col=mycol)
					legtxt = c(legtxt, paste0("SS ", switch(sp, '7'="Bt", '9'="Bt/B0", "age-0 recruits")))
					leglty = c(leglty, 1)
					legcol = c(legcol, mycol)
					addLegend (0.975, 0.975, legend=linguaFranca(legtxt,l), lty=leglty, col=legcol, lwd=2, seg.len=3, bty="n", xjust=1, yjust=1)
#browser();return()
				}
				else {
					for (j in 1:ncol(yvals)) {
						jj = colnames(yvals)[j]
						if (grepl("SX:",jj)) {
							col.sex = ifelse(grepl("SX:1",jj),"orange",.colBlind["bluegreen"])
							points(ts$YrSeas[plot1], yvals[plot1,j], pch=22, col=col.sex, bg=colorspace::lighten(col.sex, amount=0.25)) #lucent(mycol,0.2))
							lines(ts$YrSeas[plot2|plot3], yvals[plot2|plot3,j], col=col.sex, lwd=2)
						}
						else {
							#col.pch = ifelse(grepl("Bio_all",jj),"black",mycol)
							col.fleet = c("blue","green3","red","purple") ## may need more colours
							col.pch = switch(jj, 'Bio_all'="black", 'V1'=col.fleet[1], 'V2'=col.fleet[2], 'V3'=col.fleet[3], 'V4'=col.fleet[4], mycol)
							points(ts$YrSeas[plot1], yvals[plot1,j], pch=22, col=col.pch, bg=colorspace::lighten(col.pch, amount=0.25)) #lucent(mycol,0.2))
							lines(ts$YrSeas[plot2], yvals[plot2,j], col="gainsboro", lwd=2)
							points(ts$YrSeas[plot2], yvals[plot2,j], pch=21, col=col.pch, bg=colorspace::lighten(col.pch, amount=0.25)) #lucent(mycol,0.2))
							points(ts$YrSeas[plot3], yvals[plot3,j], pch=24, cex=0.8, col=col.pch, bg=colorspace::lighten(col.pch, amount=0.25)) #lucent(mycol,0.2))
						}
						if (subplot %in% c(7,9,11) || colnames(yvals)[j] %in% c("SpawnBio")) {
							if (subplot == 7) {
								plot1 <- stdtable$Label == "SSB_Virgin"
								stdtable$Yr[plot1] <- stdtable$Yr[plot1] + yrshift
							}
							if (subplot == 9) {
								plot1 <- stdtable$Label == "Bratio_Virgin"
							}
							if (subplot == 11) {
								plot1 <- stdtable$Label == "Recr_Virgin"
								stdtable$Yr[plot1] <- stdtable$Yr[plot1] + 1
							}
							plot2 <- stdtable$Yr %in% ts$Yr[plot2]
							plot3 <- stdtable$Yr %in% ts$Yr[plot3]
							plotall <- plot1 | plot2 | plot3
						}
						if (subplot %in% c(7, 9) || colnames(yvals)[j] %in% c("SpawnBio")) {
							## virgin
							arrows(x0=rep(stdtable$Yr[plot1] + 1,2), x1=rep(stdtable$Yr[plot1] + 1,2),
								y0=rep(yvals[plot1,j],2), y1=c(stdtable$upper[plot1], stdtable$lower[plot1]), angle=90, length=0.05, col=mycol)
							## trajectory
							polygon(x=c(stdtable$Yr[plot2],rev(stdtable$Yr[plot2])), y=c(stdtable$lower[plot2],rev(stdtable$upper[plot2])), border=FALSE, col=lucent(mycol,0.2))
							lines(stdtable$Yr[plot2], stdtable$upper[plot2], lty=2, col=mycol)
							lines(stdtable$Yr[plot2], stdtable$lower[plot2], lty=2, col=mycol)
							## forecast
							if (any(plot3)) {
								plot33 = plot3 
								plot33[c(grep(TRUE,plot3)[1]-1,grep(TRUE,plot3))]=TRUE
								polygon(x=c(stdtable$Yr[plot33],rev(stdtable$Yr[plot33])), y=c(stdtable$lower[plot33],rev(stdtable$upper[plot33])), border=FALSE, col=lucent(mycol,0.10))
								lines(stdtable$Yr[plot33], stdtable$upper[plot33], lty=3, col=mycol)
								lines(stdtable$Yr[plot33], stdtable$lower[plot33], lty=3, col=mycol)
							}
						}
						if (subplot == 11) {
							old_warn <- options()$warn
							options(warn=-1)
							arrows(x0=stdtable$Yr[plotall], y0=stdtable$lower[plotall], y1=stdtable$upper[plotall], length=0.01, angle=90, code=3, col=mycol)
							options(warn=old_warn)
						}
					}
				}
			}
			if (nareas > 1 & subplot %in% c(2, 3, 5, 6, 8, 10, 12)) 
				addLegend(0.975, 0.975, legend=linguaFranca(areanames[areas],l), lty=1, pch=1, col=areacols[areas], bty="n", xjust=1, yjust=1)
			if (subplot %in% c(101:102)) {
				ts$catch = apply(ts[,grep("retain\\(B)",colnames(ts))],1,sum)
				drawBars(ts$Yr[plot2|plot3], ts$catch[plot2|plot3], col="gainsboro", fill="hotpink", width=1)
				legtxt = NULL
				if (subplot %in% 101) {
					legtxt = c("Total biomass","Female biomass", "Male biomass", "Spawning biomass", "Retained catch")
					legcol = col=c("black", "orange", .colBlind["bluegreen"], mycol, "hotpink")
				}
				else if (subplot %in% 102) {
					legtxt = c("Total biomass",paste0("Vulnerable biomass - fleet ",1:nfleets), "Retained catch")
					legcol = col=c("black", col.fleet[1:nfleets], "hotpink")
				}
				if (!is.null(legtxt))
					addLegend(0.975, 0.975, legend=linguaFranca(legtxt,l), pch=22, col=legcol, pt.bg=colorspace::lighten(legcol,0.2), bty="n", xjust=1, yjust=1)
				box()
			}
		}
#browser();return()
		if (verbose) {
			message("  finished time series subplot ", subplot, ": ", main)
		}
		if (print) 
			dev.off()
		attr(ts,"plotinfo") = plotinfo
		attr(ts,"category") = "Timeseries"
		return(ts)
	}
	skip <- FALSE
	if (nareas == 1 & subplot %in% c(2, 5, 8, 10, 12:13)) {
		skip <- TRUE; mess = paste0("subplot ", subplot, " not available for single area") }
	if (nseasons == 1 & subplot %in% c(3, 6)) {
		skip <- TRUE; mess = paste0("subplot ", subplot, " not available for single season") }
	if (subplot %in% c(14:15) & (is.null(birthseas) || nbirthseas == 1)) {
		skip <- TRUE; mess = paste0("subplot ", subplot, " not available for unspecified or single birth season") }
	if (!skip) {
		out <- biofunc(subplot=subplot, outnam=outnam, l=lang)
		return(invisible(out))
	} else {
		.flush.cat("Alert:", mess, "\n")
		return (mess)
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.ts


##====================================== Start suite of quantile boxes of age-fit residuals

## QUANTILE BOXES of AGE-FIT RESIDUALS by AGE, YEAR, COHORT

## plt.ageResids------------------------2020-11-04
## AME changing for POP, just plotting age class resids here,
## not qq-plot. Moving that to new function (so can do 4 on page).
## See popScape.r for original.
## -----------------------------------------AME|RH
plt.ageResids <- function( obj, ages=NULL, main=NULL, lang="e", yrs, ...)
{
	## par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )
	## Subset to required ages.
	if (!is.null(ages))
		obj <- obj[ (obj$Bin >= ages[1]) & (obj$Bin <=ages[2]), ]
	else
		ages = range(obj$Bin)
	if( max(diff(sort(unique(obj$Bin)))) > 1) {
		allAges=min(obj$Bin):max(obj$Bin)
		nodataAges=allAges[ !(allAges %in% obj$Bin)]
		xx=split(c(obj$Pearson, rep(NA, length(nodataAges))), c(obj$Bin, nodataAges))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantBox(xx, xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xage = setdiff(allAges, nodataAges)
		if (length(xage)>20)
			xage = sort(unique(ceiling(xage/5)*5))
		axis(1, at=match(xage,as.numeric(xpos$names)), labels=xage, ..., mgp=c(2,0.5,0))
	} else {
		#xpos <- boxplot( split( obj$Pearson, obj$Bin ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantBox(split(obj$Pearson,obj$Bin), xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xage = sort(unique(ceiling(obj$Bin/5)*5))
		axis(1, at=match(xage,as.numeric(xpos$names)), labels=xage, ..., mgp=c(2,0.5,0))
	}
	abline( h=0, lty=2, col="red" )
	axis(2, ...)
	mtext( side=1, line=2, ..., text=linguaFranca("Age class",lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.ageResids

## plt.yearResids-----------------------2019-07-18
## AME adding to plot age residuals by year. Is called for comm and survs.
##  fill.in=TRUE is to add the missing years for boxplot
##  ..POP does not do qq plot. See popScape.r for previous.
## -----------------------------------------AME|RH
plt.yearResids <- function(obj, ages=NULL, main=NULL, fill.in=TRUE, lang="e", ...)
{
	# Subset to required ages - still do as don't want age 1.
	if (!is.null(ages))
		obj <- obj[ (obj$Bin >= ages[1]) & (obj$Bin <=ages[2]), ]
	if(fill.in) {
		allYears = min(obj$Yr):max(obj$Yr)
		nodataYears = allYears[ !(allYears %in% obj$Yr)]
		xx = split(c(obj$Pearson, rep(NA, length(nodataYears))), c(obj$Yr, nodataYears))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE, ... )     #AME outline=FALSE removes outliers
		xpos = quantBox(xx, xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars), ...)
		xyrs = setdiff(allYears, nodataYears)
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, ..., mgp=c(2,0.5,0))
#browser();return()
	} else {
		#xpos <- boxplot( split( obj$Pearson, obj$Yr ), whisklty=1, xlab="", ylab="", outline=FALSE, ... ) #AME outline=FALSE removes outliers
		xpos = quantBox( split(obj$Pearson,obj$Yr), xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars), ...)
		xyrs = sort(unique(ceiling(obj$Yr/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, ..., mgp=c(2,0.5,0))
	}
	abline( h=0, lty=2, col="red" )
	axis(2, ...)
	mtext( side=1, line=2, ..., text=linguaFranca("Year",lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.yearResids

## plt.cohortResids---------------------2019-07-18
## Plot age residuals by cohort.
## -----------------------------------------AME|RH
plt.cohortResids <- function( obj, ages=NULL, main=NULL, lang="e", use.rdevs=T, ...)
{
	## Input is the CAc object from a Awatea res file. Ages to 59 as
	##  plus-age class will mess up year-of-birth calculation. Not automated.
	## par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )
	## Subset to required ages - still do as don't want age 1 or 60 for cohorts
	if (!is.null(ages))
		obj <- obj[ (obj$Bin >= ages[1]) & (obj$Bin <=ages[2]), ]
	## obj$Pearson has residuals for each age, year and both sexes. Need
	##  to assign a year of birth for each age as an extra column, then
	##  presumably just do the boxplot split using that.
	obj$birthyr = obj$Yr - obj$Bin
	if( max(diff(sort(unique(obj$birthyr)))) > 1) {
		allYears = min(obj$birthyr):max(obj$birthyr)
		nodataYears = allYears[ !(allYears %in% obj$birthyr)]
		xx = split(c(obj$Pearson, rep(NA, length(nodataYears))), c(obj$birthyr, nodataYears))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantBox(xx, xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xyrs = setdiff(allYears, nodataYears)
		if (length(xyrs)>20)
			xyrs = sort(unique(ceiling(xyrs/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, ..., mgp=c(2,0.5,0))
	} else {
		#xpos=boxplot( split( obj$Pearson, obj$birthyr ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantBox(split(obj$Pearson,obj$birthyr), xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xyrs = sort(unique(ceiling(obj$birthyr/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, ..., mgp=c(2,0.5,0))
	}
	if (use.rdevs) {
		plotSS.rdevs(replist, subplot=1, plot=F, print=F)
		if (!is.null(ttcall(recdevs))){ ## created by function 'plotSS.rdevs'
			ttget(recdevs)
			zyrs = recdevs$Yr[recdevs$Value>0]
			zpar = upar = tcall(boxpars)
			zpar$boxfill=lucent("green",0.25); zpar$medcol="green4"
			upar$boxfill=lucent("red",0.25); upar$medcol="red"
			#zpar$boxfill=lucent("red",0.5)
#browser();return()
			zbox = ubox = split(obj$Pearson,obj$birthyr)
			zbox[!is.element(names(zbox),zyrs)] = NA
			ubox[is.element(names(zbox),zyrs)] = NA
			quantBox(zbox, xaxt="n", yaxt="n", xlab="", ylab="", pars=zpar, add=TRUE)
			quantBox(ubox, xaxt="n", yaxt="n", xlab="", ylab="", pars=upar, add=TRUE)
		}
	}
	abline( h=0, lty=2, col="blue" )
	axis(2, ...)
	mtext( side=1, line=2, ..., text=linguaFranca("Year of birth",lang) )
	if ( !is.null(main) )
		mtext( side=3, line=-0.5, ..., outer=TRUE, text=linguaFranca(main,lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.cohortResids

##====================================== End suite of quantile boxes of age-fit residuals


## plt.selectivity----------------------2019-07-18
## Transferred selectivity code from PBSscape.r
## into plt function (RH 190718)
## sobj = second object (RH 201119)
## -----------------------------------------AME|RH
plt.selectivity <- function( obj, sobj=NULL, mainTitle="Rockfish", maxage,
   ptypes="win", pngres=400, PIN=c(7,8), lang=c("e","f"))
{
	## Plot the selectivity.
	panel.selex = function(x, y, ...){
		dots = list(...)[[1]]
		unpackList(dots)
		nsex = abs(getNpan()%%2-2)
		#unpackList(tcall(panel.dots))
		abline(h=seq(0.2,0.8,0.2), v=seq(10,max(x)-1,10), lty=3, col="seashell4")
		if (!is.null(sobj)) {
			z  = is.element(sobj$Index, sub("TRAWL\\+","TRAWL",index.e[getNpan()]))
			x2 = sobj$Age[z]; y2=sobj$Sel[z]
			lines(x, y, lty=1, col=scol[nsex], lwd=2)
			lines(x2, y2, lty=2, col=scol[nsex], lwd=2)
			m2 = sobj$Sel[grep("MATURITY: Female",sobj$Index)]
			lines(x2, m2, lty=2, col="grey20")
		}
		else
			lines(x, y, lty=1, col=scol[nsex], lwd=slwd)
		mvec = rep("m", length(m))
		if (length(m)>50)
			mvec[seq(2,length(m),2)]="\225"
		if (getNpan()%%2==1)
			text(x, m, mvec, col=lucent("black",0.75), cex=0.9)
#browser();return()
		#points(x, m, pch=109, col=lucent("black",0.75), cex=0.8) #, bg=lucent("gold",0.25)
		#lines(x, m, lty=1, col="pink", lwd=2); lines(x, m, lty=2, col="grey", lwd=2)
	}
	Pvec = intersect(as.character(0:500), colnames(obj))
	if (!missing(maxage)){
		realmaxage = as.numeric(rev(Pvec)[1])
		Pvec = as.character(0:min(maxage,realmaxage))
	}
	Xfile = NULL
	for (i in 1:nrow(obj)) {
		ifile = data.frame(Index=rep(obj[i,"Index"],length(Pvec)), Age=as.numeric(Pvec), Sel=unlist(obj[i,Pvec]))
		Xfile = rbind(Xfile,ifile)
	}
	## Alter xmax (needs changes in object names)
	#xmax = max(as.numeric(sapply(selP,function(x){names(x[is.element(x,1)])[1]})),na.rm=TRUE) #maximum minimum age when P first hits 1
	#if (is.na(xmax)) xmax = obj$extra$general$Nages
	## Check for dome-selectivity estimation
	#if ( any(unlist(sapply(currentRes$extra$priors[c("log_varRest_prior","log_surveyvarR_prior")],function(x){x[,1]}))>0) )
	#	xmax = obj$extra$general$Nages
	#obj$Sel = obj$Sel[obj$Sel$Age <= xmax,]

	## linguafranca takes too long on big vectors so shorten to unique
	Xfile.f  = Xfile
	index.e  = unique(Xfile$Index)
	index.f  = linguaFranca(index.e,"f")
	Xfile.f$Index = index.f[Xfile$Index]

	rc = c(ceiling(length(index.e)/2),2)

	## Get the maturity ogive
#browser();return()
	mats = getSS.control(tcall(control.file))$maturity_ogive[1:length(Pvec)]

	strip.list = list(col=lucent("black",0.5), bg=lucent("moccasin",0.5), height=0.1, cex=1.8)

	#panel.dots = list(scol=c("red","blue"), lwd=4)
	#tput(panel.dots)
	createFdir(lang)
	fout = fout.e = "selectivity"
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		if (is.null(ptypes)) ptypes="win"
		for (p in ptypes) {
			if (p=="eps")      postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			## regular par settings are ignored in lattice::xyplot -- need to fiddle with lattice settings
			#par(mfrow=c(1,1), mar=c(3.2,3.2,0.5,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			#par.sel = list(layout.widths = list(left.padding=-0.5, right.padding=-0.5),
			#	layout.heights=list(top.padding=0.5, main.key.padding=-0.75, bottom.padding=-2))
			if (l=="f")
				mochaLatte(Xfile.f,"Age","Sel","Index", panel=panel.selex, m=mats, strip=strip.list, scol=c("red","blue"), slwd=2, mar=c(0,0,0,0), oma=c(4,4,4,1), xlab=linguaFranca("Age (years)","f"), ylab=linguaFranca("Proportion","f"), mainTitle=paste0("s\u{00E9}lectivit\u{00E9} du ", linguaFranca(mainTitle,l)), cex.axis=1.2, cex.main=1.5, sobj=sobj, rc=rc, ylim=c(0,1))
			else{
				mochaLatte(Xfile,"Age","Sel","Index", panel=panel.selex, m=mats, strip=strip.list, scol=c("red","blue"), slwd=2, mar=c(0,0,0,0), oma=c(4,4,4,1), xlab="Age (years)", ylab="Proportion", mainTitle = paste0(mainTitle, " Selectivity"), cex.axis=1.2, sobj=sobj, rc=rc, ylim=c(0,1))
				#mtext(paste0(mainTitle, " Selectivity"), side=3, outer=FALSE, line=0, cex=1.5)
#browser();return()
			}
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.selectivity


#getYrIdx-------------------------------2011-08-31
# Purpose is to return selected years for plotting.
# Default is to select 5 year increments.
#----------------------------------------------AME
getYrIdx <- function( yrNames,mod=5 )
{
  # Coerce to numerice and select the years modulo "mod".
  yrVals <- as.numeric( gsub("[^[:digit:]]","",yrNames) )
  idx <- yrVals %% mod==0

  # Select years from character vector yrNames.
  result <- yrNames[ idx ]
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getYrIdx


