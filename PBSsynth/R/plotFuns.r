##================================================2026-02-23
## PBS Stock Synthesis plotting functions:
## ---------------------------------------
## calcRhat..............Plot the split-Rhat statistic and ESS
## compAE................Compare AE vectors (replaces Excel foofarall)
## compTS................Compare time series medians
## compWT................Compare weights or lengths among survey series
## mochaLatte............An alternative to lattice plots (mockLattice)
## panelBoxes............Plot quantile plots using 'nchains' to delimit separate boxes
## panelChains...........Plots cumulative fequency of 'nchains' by partitioning one trace
## panelTraces...........Plots sequential  trace of MCMC samples with running median and (0.05, 0.95) quantiles
## plot.cres.new.........Produce 5-panel overview of residuals (pkg=compResidual)
## plotACFs..............Plot ACFs for the estimated parameters
## plotBh................Plot biomass (B) versus steepness (h)
## plotSS.dmcmc..........Plot diagnostics (traces, split chains, ACFs) for MCMCs
## plotSS.pairs..........Pairs|density plot comparison among parameters
## plotSS.pmcmc..........Plot parameter quantile boxplots for MCMCs
## plotSS.rmcmc..........Plot routine output for MCMCs
## plotTraj..............Show all median trajectories (base+sens) in one figure
## plt.ageResids.........Plot age residuals by age class (mod.PBSawatea)
## plt.cohortResids......Plot age residuals by cohort (mod.PBSawatea)
## plt.yearResids........Plot age residuals by year (mod.PBSawatea)
## plt.selectivity.......Transferred selectivity code from PBSscape.r (mod.PBSawatea)
##==========================================================


## calcRhat ----------------------------2025-12-04
##  Calculate R-hat and ESS statistics for MCMC chains.
## ---------------------------------------AL|AV|RH
calcRhat <- function(dir=".", nchains=8, parpos, rhat.only=FALSE,
   only.bad=FALSE, badhat=1.01, offset=0.0025, recdevs=FALSE, 
   zero.bar=TRUE, barcols=.colBlind, xlim=NULL, exlax=NULL,
   htype=c("histogram", "density"), stacked=FALSE, smoother="density",
   png=FALSE, pngres=400, PIN=c(10,7.5),
   lang="e", outnam)
{
	## -----------------------------------------------
	## One way to monitor whether a chain has converged to the equilibrium
	##   distribution is to compare its behavior to other randomly initialized chains.
	##   This is the motivation for the potential scale reduction statistic, split-Rhat.
	## The split-Rhat statistic measures the ratio of the average variance of draws
	##   within each chain to the variance of the pooled draws across chains;
	##   if all chains are at equilibrium, these will be the same and Rhat=1.
	##   If the chains have not converged to a common distribution, the Rhat
	##   statistic will be greater than one.
	## Rhat is computed for each scalar quantity of interest, as the standard deviation
	##   of that quantity from all the chains included together, divided by the root mean
	##   square of the separate within-chain standard deviations.
	## Code sent by Adam Langley
	##   https://cran.r-project.org/web/packages/bayesplot/vignettes/visual-mcmc-diagnostics.html#general-mcmc-diagnostics
	## -----------------------------------------------

	## Subfunctions from rstan------------------------
	##  (because I cannot get rstan to install)
	is_constant <- function(x, tol = .Machine$double.eps) {
		abs(max(x) - min(x)) < tol
	}
	should_return_NA <- function(x) {
		# should NA be returned by a convergence diagnostic?
		anyNA(x) || any(!is.finite(x)) || is_constant(x)
	}
	autocovariance <- function(y) {
	  N <- length(y)
	  M <- fft_next_good_size(N)
	  Mt2 <- 2 * M
	  yc <- y - mean(y)
	  yc <- c(yc, rep.int(0, Mt2 - N))
	  transform <- fft(yc)
	  ac <- fft(Conj(transform) * transform, inverse = TRUE)
	  # use "biased" estimate as recommended by Geyer (1992)
	  ac <- Re(ac)[1:N] / (N^2 * 2)
	  ac
	}
	fft_next_good_size <- function(N) {
	  # Find the optimal next size for the FFT so that
	  # a minimum number of zeros are padded.
	  if (N <= 2)
	    return(2)
	  while (TRUE) {
	    m <- N
	    while ((m %% 2) == 0) m <- m / 2
	    while ((m %% 3) == 0) m <- m / 3
	    while ((m %% 5) == 0) m <- m / 5
	    if (m <= 1)
	      return(N)
	    N <- N + 1
	  }
	}
	rhat_rfun <- function(sims) {
		#' @references
		#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
		#' Paul-Christian B\"{u}rkner (2019). Rank-normalization, folding, and
		#' localization: An improved R-hat for assessing convergence of
		#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
		if (anyNA(sims)) {
			return(NA)
		}
		if (any(!is.finite(sims))) {
			return(NaN)
		}
		if (is_constant(sims)) {
			return(NA)
		}
		if (is.vector(sims)) {
			dim(sims) <- c(length(sims), 1)
		}
		chains <- ncol(sims)
		n_samples <- nrow(sims)
		chain_mean <- numeric(chains)
		chain_var <- numeric(chains)
		for (i in seq_len(chains)) {
			chain_mean[i] <- mean(sims[, i])
			chain_var[i] <- var(sims[, i])
		}
		var_between <- n_samples * var(chain_mean)
		var_within <- mean(chain_var)
		sqrt((var_between / var_within + n_samples - 1) / n_samples)
	}
	ess_rfun <- function(sims) {
		#' @references
		#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
		#' Paul-Christian B\"{u}rkner (2019). Rank-normalization, folding, and
		#' localization: An improved R-hat for assessing convergence of
		#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
		if (is.vector(sims)) {
			dim(sims) <- c(length(sims), 1)
		}
		chains <- ncol(sims)
		n_samples <- nrow(sims)
		if (n_samples < 3L || should_return_NA(sims)) {
			return(NA_real_)
		}
		acov <- lapply(seq_len(chains), function(i) autocovariance(sims[, i]))
		acov <- do.call(cbind, acov)
		chain_mean <- apply(sims, 2, mean)
		mean_var <- mean(acov[1, ]) * n_samples / (n_samples - 1)
		var_plus <- mean_var * (n_samples - 1) / n_samples
		if (chains > 1)
			var_plus <- var_plus + var(chain_mean)
	
		# Geyer's initial positive sequence
		rho_hat_t <- rep.int(0, n_samples)
		t <- 0
		rho_hat_even <- 1
		rho_hat_t[t + 1] <- rho_hat_even
		rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
		rho_hat_t[t + 2] <- rho_hat_odd
		while (t < nrow(acov) - 5 && !is.nan(rho_hat_even + rho_hat_odd) &&
					 (rho_hat_even + rho_hat_odd > 0)) {
			t <- t + 2
			rho_hat_even = 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
			rho_hat_odd = 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
			if ((rho_hat_even + rho_hat_odd) >= 0) {
				rho_hat_t[t + 1] <- rho_hat_even
				rho_hat_t[t + 2] <- rho_hat_odd
			}
		}
		max_t <- t
		# this is used in the improved estimate
		if (rho_hat_even>0)
				rho_hat_t[max_t + 1] <- rho_hat_even
		
		# Geyer's initial monotone sequence
		t <- 0
		while (t <= max_t - 4) {
			t <- t + 2
			if (rho_hat_t[t + 1] + rho_hat_t[t + 2] >
					rho_hat_t[t - 1] + rho_hat_t[t]) {
				rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2;
				rho_hat_t[t + 2] = rho_hat_t[t + 1];
			}
		}
		ess <- chains * n_samples
		# Geyer's truncated estimate
		# tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t])
		# Improved estimate reduces variance in antithetic case
		tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t+1]
		# Safety check for negative values and with max ess equal to ess*log10(ess)
		tau_hat <- max(tau_hat, 1/log10(ess))
		ess <- ess / tau_hat
		ess
	}
	addESS <- function (ess, zvec, xlim, bcol) {
		if (sum(zvec)==0)
			return()
		xoff = abs(diff(xlim)) * 0.1 #0.075
		yoff = 0.5
		if (zero.bar) ess = c(0,ess)
		sess = scaleVec(V=ess, Tmin=xlim[1]+ifelse(zero.bar,0,xoff), Tmax=xlim[2]-xoff)
		if (zero.bar) sess = sess[-1]
		zess = sess[zvec]
		x0   = rep(xlim[1], sum(zvec))
		x1   = zess
		y0   = y1 = (1:length(zvec))[zvec]
		xpol = rbind(x0, x0, x1, x1, NA)
		ypol = rbind(y0-yoff, y0+yoff, y1+yoff, y1-yoff, NA)
		polygon(as.vector(xpol), as.vector(ypol), col=lucent(bcol,0.25), border="grey", lwd=0.5)
	}
	dcurve <- function(x) {
	#dcurve <- function(x) {
		if (smoother %in% c("density")) {
			bw="nrd0"; adjust=1; kernel="gaussian"; gridsize=512; cut=3; na.rm=TRUE
			do.call(density, args=list(x=x, bw=bw, adjust=adjust, kernel=kernel, n=gridsize, cut=cut, na.rm=na.rm))
		} else if (smoother %in% c("bkde")) {
			require(KernSmooth)  ## (RH 251125) best to let algorithm choose bandwith (each data set seems to require a different bandwidth)
			bw=0.05; kernel="normal"; gridsize=401L
			do.call(bkde, args=list(x=x, kernel=kernel, gridsize=gridsize ))
		}
	}
	## End subfunctions ---------------------------

	## Need to find MCMC output when extra=TRUE
	if (dir==".") dir = mcmc.dir = getwd()
	if (dir.exists(file.path(dir,"sso")))
		mcmc.dir = file.path(dir,"sso")
	else
		mcmc.dir = dir

	filename = file.path(mcmc.dir,"posteriors.sso")
	if (!file.exists(filename))
		stop("No file name 'posteriors.sso' at\n\t", dir)
	dat <- read.table(filename, header=T)
	dat$chain = rep(1:8, each=nrow(dat)/nchains)
#browser();return()
	if (missing(parpos))
		parpos = 1:ncol(dat)
	cdat <- split(dat,dat$chain)
	## rhat is computed for each scalar quantity of interest, as
	## the standard deviation of that quantity from all the chains included together,
	## divided by the root mean square of the separate within-chain standard deviations. 
	rhat <- Rhat <- matrix(NA, nrow=length(parpos), ncol=1)
	Rhat.rstan   <- matrix(NA, nrow=length(parpos), ncol=2)
	## Gather parameters for comparing chain histograms
	phis = list()

	for (i in 1:length(parpos)) {
		ii    = parpos[i]
		iii   = colnames(dat)[ii]
#.flush.cat(iii,"\n")
		if (iii %in% c("Iter","Objective_function","chain")) next
		ilist = lapply(cdat, function(x) { return(x[,ii]) })
		## imat only for testing rstan's function 'rhat_rfun'
		imat  = do.call("cbind", lapply(ilist, data.frame, stringsAsFactors=FALSE))
		colnames(imat) = paste0("chain",1:ncol(imat))

		## Check : https://avehtari.github.io/masterclass/slides_rhat_neff.pdf
		nn = length(ilist[[1]])
		mm = nchains
		## Within-chains variance W
		Vj = sapply(ilist, function(x) { (1/(nn-1)) * sum((x - mean(x))^2) })
		W  = (1/mm) * sum(Vj)
		## Between-chains variance B
		psibar.j = sapply(ilist, function(x){ (1/nn) * sum(x) })        ## mm values
		psibar.. = (1/mm) * sum(sapply(ilist, function(x){ mean(x) }))  ## one value
		B  = (nn/(mm-1)) * sum((psibar.j -psibar..)^2)                  ## B/nn = variance of means of chains (psibar.j)
		## Estimate  total variance as a weighted mean of W
		Vhat = ((nn-1)/nn) * W + (1/nn) * B
		## As Vhat overestimates and W underestimates, compute Rhat:
		Rhat[i,1]  = sqrt(Vhat / W)

		## rstan calculation is same as Rhat because both use Aki Vehtari's 'improved' rhat
		Rhat.rstan[i,1] = rhat_rfun(as.matrix(imat))
		Rhat.rstan[i,2] = ess_rfun(as.matrix(imat))
		
#if (round(Rhat.rstan[i,1],7)!=round(Rhat[i,1],7)) {browser();return()}  ## test

		## Adam Langley's simple Rhat calculation
		sdall = sqrt(var(dat[,ii]))
		sdchn = sapply(ilist, function(x) { return(sqrt(var(x))) })
		rms   = sqrt(sum(sdchn^2)/nchains)
		rhat[i,1]  = sdall/rms

		## Gather a few parameters for histograms of chains
		if (grepl("R0|NatM|BH", iii)) {
			phis[[iii]] = ilist
		}
	} ## end i loop through parameters
	phis = phis[c(1,2,4,3)] ## re-order to NatM_F, NatM_M, BH_h, LN_R0
#browser();return()
	colnames(rhat) = "rhat"
	colnames(Rhat) = "Rhat"
	colnames(Rhat.rstan) = c("Rhat","ESS")
	rownames(rhat) = rownames(Rhat) = rownames(Rhat.rstan) = colnames(dat)[parpos]

	if (!recdevs && !rhat.only) {  ## extra plots
		comp.adam = FALSE
		if (comp.adam) {
			## Compare rhat (Adam Langley) vs. Rhat (Aki Vehtari)
			rhat.comp = cbind(rhat,Rhat)
			fout.e = paste0("rhat.comp.", tolower(basename(dir)))
			for (l in lang) {
				changeLangOpts(L=l)
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				if (png) {
					clearFiles(paste0(fout,".png"))
					png(file=paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
				}
				expandGraph(mfrow=c(1,1), mar=c(3.5,3.5,1,1), mgp=c(2,0.5,0))
				plot(0,0, xlim=range(rhat.comp, na.rm=TRUE), ylim=range(rhat.comp,na.rm=T), type="n", xlab=linguaFranca("simple rhat (Adam Langley)",l), ylab=linguaFranca("Rhat (Aki Vehtari)",l), cex.axis=1.2, cex.lab=1.5)
				abline(a=0,b=1,col="blue",lwd=2)
				points(rhat.comp[,"rhat"], rhat.comp[,"Rhat"], pch=21, bg="yellow", cex=1.1)
				if (png) dev.off()
			}; eop()  ## end lang loop
			if (rhat.only) {browser();return()}
		}

		## Plot histograms of chains for selectct parameters
		if ("histogram" %in% htype) {
			fout.e = paste0("rhat.", ifelse(stacked, "shist", "ohist"), ".", tolower(basename(dir)))
			if (!stacked)
				barcols = c("blue","red","green4","deepskyblue","orange","yellowgreen","purple","sienna") ## (RH 251125)
			barcols = rep(barcols,nchains)[1:nchains]
			for (l in lang) {
				changeLangOpts(L=l)
				#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				if (png) {
					clearFiles(paste0(fout,".png"))
					png(file=paste0(fout,".png"), width=10, height=9, units="in", res=pngres)
				}
				expandGraph(mfrow=.findSquare(length(phis)), mar=c(3.5,3.5,1,1), mgp=c(2,0.5,0))  ## for histograms
				for (j in  1:length(phis)) {
					jj    = names(phis)[j]
					jlist = phis[[jj]]
					jdat  = unlist(jlist)
					jlab = ifelse(grepl("R0",jj),"Log Recruitment (R0)", ifelse(grepl("NatM",jj), paste0("Natural Mortality (",ifelse(grepl("Fem",jj),"female","male"),")"), ifelse(grepl("BH",jj), "Beverton-Holt Steepness", jj)))
					brks = seq(min(jdat), max(jdat), length.out=50)
					if (stacked) {
						hist(unlist(jlist), breaks=brks, col="transparent", border=FALSE, main="", cex.axis=1.2, cex.lab=1.5, cex.main=1.5, las=1, xlab=linguaFranca(jlab,l), ylab=linguaFranca("Frequency",l), yaxs="i")
					} else {
						cdlim = sapply(jlist, function(x) {
							h  = hist(x=x, breaks=brks, plot=F)
							c(cmax=max(h$counts), dmax=max(h$density), xmin=min(h$breaks), xmax=max(h$breaks))
						})
						clim = c(0, max(cdlim[1,]) + 0.025 * max(cdlim[1,]))
						dlim = c(0, max(cdlim[2,]) + 0.025 * max(cdlim[2,]))
						xlim = range(cdlim[3:4,])
						hist(jdat, breaks=brks, col="transparent", border=FALSE, main="", cex.axis=1.2, cex.lab=1.5, cex.main=1.5, las=1, xlab=linguaFranca(jlab,l), ylab=linguaFranca("Frequency",l), yaxs="i", ylim=clim)
					}
					ncol = length(jlist)
					for (k in ncol:1) {  ## overlay successively smaller sets of samples from phis 
						kk = if (stacked) 1:k else (ncol + 1) - k
						eval(parse(text=paste0("kcol = switch(k ,\"", paste0(barcols, collapse="\", \""), "\")")))
						if (!stacked)
							kcol = lucent(kcol,0.20)
						kkk  = unlist(jlist[kk])
						hist(kkk, breaks=brks, col=kcol, border=ifelse(stacked,"gainsboro",kcol), lwd=0.5, add=T)
					}
#browser();return()
					if (all(par()$mfg[1:2]==1))
						addStrip(0.825, 0.95, col=barcols, lab=paste0(linguaFranca("chain",l,little=5)," ",1:ncol), xwidth=0.02, yheight=0.3)
				} ## end histograms of select parameters
				if (png) dev.off()
			}; eop()  ## end lang loop
		} ## end htype=histogram

		## Plot density overlays for chains of select parameters
		if ("density" %in% htype) {
			fout.e = paste0("rhat.", substring(smoother,1,4), ".", tolower(basename(dir)))
			barcols = rep(barcols,nchains)[1:nchains]
			for (l in lang) {
				changeLangOpts(L=l)
				fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
				if (png) {
					clearFiles(paste0(fout,".png"))
					png(file=paste0(fout,".png"), width=10, height=9, units="in", res=pngres)
				}
				expandGraph(mfrow=.findSquare(length(phis)), mar=c(3.5,3.5,1,1), mgp=c(2,0.5,0))  ## for density curves
				for (j in  1:length(phis)) {
					jj    = names(phis)[j]
					jlist = phis[[jj]]
					jdat  = unlist(jlist)
					jlab = ifelse(grepl("R0",jj),"Log Recruitment (R0)", ifelse(grepl("NatM",jj), paste0("Natural Mortality (",ifelse(grepl("Fem",jj),"female","male"),")"), ifelse(grepl("BH",jj), "Beverton-Holt Steepness", jj)))
					#xlim = range(unlist(lapply(jlist,range)))
					xylim = sapply(jlist, function(x) {
						d  = dcurve(x=x)
						c(xmin=min(d$x), dmax=max(d$x), ymax=max(d$y), zmin=min(x), zmax=max(x))
					})
					xlim = range(xylim[1:2,])
					ylim = c(0, max(xylim[3,]) + 0.025 * max(xylim[3,]))
					zlim = range(xylim[4:5,])
					bigden= dcurve(x=jdat)
					plot(bigden, col="transparent", main="", cex.axis=1.2, cex.lab=1.5, cex.main=1.5, las=1, xlab=linguaFranca(jlab,l), ylab=linguaFranca("Density",l), yaxs="i", xlim=xlim, ylim=ylim)
					abline(v=zlim, lty=2, col="black")
					ncol = length(jlist)
					#dencols = barcols; dencols[1]="gainsboro"; names(dencols)[1]="grey"
					#dencols = .colGnuplot
					dencols = c("blue","red","green4","deepskyblue","orange","yellowgreen","purple","sienna") ## (RH 251125)
#if(j==3) {browser();return()}
					for (k in ncol:1) {  ## overlay density curves
						kk = (ncol + 1) - k
						eval(parse(text=paste0("kcol = switch(k ,\"", paste0(dencols, collapse="\", \""), "\")")))
						kkk  = unlist(jlist[kk])
						kden = dcurve(x=kkk)
						xden = c(kden$x[1], kden$x, rev(kden$x)[1])
						yden = c(0, kden$y, 0)
						polygon(xden, yden, col=lucent(kcol,0.05), border=kcol, lwd=2)
#if(j==3 && k==6) {browser();return()}
					}
					box(lwd=2)
					if (all(par()$mfg[1:2]==1))
						addStrip(0.70, 0.95, col=dencols, lab=paste0(linguaFranca("chain",l,little=5)," ",1:ncol), xwidth=0.02, yheight=0.3)
				} ## end histograms of select parameters
				if (png) dev.off()
			}; eop()  ## end lang loop
		} ## end htype=density
	} ## end if not recdevs
#browser();return()

	rhat.langley = rhat; rhat = Rhat.rstan  ## switch over to using rstan's calculations (RH 240529)
	flotsam = c("rhat", "rhat.langley")
	## get rid of all the DEVadd parameters (for now)
	if ( any(grepl("DEVadd", rownames(rhat))) ){
		shunt = grep("DEVadd", rownames(rhat))
		keep  = grep("DEVadd", rownames(rhat), invert=TRUE)
		rhat.devadd = rhat[shunt,]
		rhat = rhat[keep,]
		flotsam = c(flotsam, "rhat.devadd")
	}
	keep  = if (recdevs) grep("RecrDev|ForeRecr", rownames(rhat))  ## need to exclude 'RecrDist'
		else grep("Iter|Objective_function|chain|^Early|^Main|^Late|^Fore", rownames(rhat), invert=TRUE)
	ess = if (ncol(rhat)>1) rhat[keep,2,drop=FALSE] else NULL  ## determine if ESS metrics exist
#browser();return()
	rhat  = rhat[keep,1,drop=FALSE] 
	nrhat = nrow(rhat)
	if (only.bad) {
		bad   = apply(rhat,1,function(x,bad) { x>bad }, bad=badhat)
		if (!any(bad)) {
			addLabel (0.5,0.5,"No bad pars") ; return()
		}
		rhat  = rhat[bad,1,drop=FALSE]
		if (!is.null(ess))
			ess   = ess[bad,1,drop=FALSE]
		nrhat = nrow(rhat)
	}
	bad  = apply(rhat,1,function(x,bad) { x>bad }, bad=badhat)
	good = !bad
	## Fix weird names in posteriors.sso
	onames = rownames(rhat)  ## original names
	rnames = convPN(onames)
	## get rid of dots etc until 'convPN' gets its shit together
	rnames = gsub("(\\_)?\\.","_",rnames)
	rnames = gsub("\\_base","",rnames)
	rnames = gsub("\\_GP1","",rnames)
#browser();return()
	rownames(rhat) = rnames
	if (!is.null(ess))
		rownames(ess) = rownames(rhat)
	flotsam = c(flotsam, "ess")
	save(list=flotsam, file="rhat.flotsam.rda")

	## Plot results
	fout.e = if (missing(outnam)) paste0("rhat.", tolower(basename(dir))) else outnam
	#xlim = c(min(0.98,min(rhat), na.rm=T), max(1.06,max(rhat), na.rm=T))
	if (is.null(xlim))
		xlim = c(min(rhat, na.rm=T) - offset, max(badhat + offset, max(rhat) + 2*offset, na.rm=T))
	xtck = pretty(xlim, n=10)
	ylim = c(dim(rhat)[1],1) + c(1,-2) * 0.05 * diff(c(1,dim(rhat)[1]))
	if (recdevs) {
		fout.e = if (missing(outnam)) paste0("rhat.recdev.", tolower(basename(dir))) else outnam
		clrs = rep("black", nrow(rhat))
		clrs[grep("^Late",onames)] = "blue"
		clrs[grep("^Fore",onames)] = "red"
		bgs = rep("gainsboro", nrow(rhat))
		bgs[grep("^Late",onames)] = "cyan"
		bgs[grep("^Fore",onames)] = "pink"
		yrs  = as.numeric(gsub("[^0-9.-]","",rownames(rhat)))
		ylim = xlim
		ytck = xtck
		xlim = range(yrs)
#browser();return()
	}
	## add in language loop
	for (l in lang) {
		changeLangOpts(L=l)
		#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			clearFiles(paste0(fout,".png"))
			png(file=paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
		}
		if (recdevs) {
			expandGraph(mfrow=c(1,1), mar=c(3,3.5,1,1), mgp=c(2,0.5,0))
			plot(0,0, xlim=xlim, ylim=ylim, type="n", xlab=linguaFranca("Year",l), ylab=linguaFranca("split R-hat",l), yaxs="i", yaxt="n", cex.lab=1.2, las=1)
			axis(2, at=ytck, tcl=-0.4)
			abline(h=c(1,badhat), lty=2:3, col=c("slategray","red"))
			x0 = yrs;  y0 = rep(0, length(yrs))  
			x1 = yrs;  y1 = rhat[,1]
			segments(x0,y0,x1,y1, lwd=2, col=clrs)
			points(x1,y1, pch=21, col=clrs, bg=bgs, cex=1)
		} else {
			ylab = rownames(rhat)
#browser();return()
			expandGraph(mfrow=c(1,1), mar=c(3, max(nchar(ylab))^ifelse(length(ylab)>20,0.7,0.75), 3,1), mgp=c(1.5,0.5,0))
			plot(0,0, xlim=xlim, ylim=ylim, type="n", xlab=linguaFranca("split R-hat (lines)",l), ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i", cex.lab=1.2)
			axis(1, at=xtck, tcl=-0.4)
			if (!is.null(ess)) {
				essvec = ess[,1]
				xoff   = abs(diff(xlim)) * 0.1 #0.075
				#if (zero.bar) xoff[1] = 0
				yoff   = 0.5
				mark   = (nchains * 50) + 0.001  ## add small fraction to fool axis labels
				xrng   = range(essvec)
				if (zero.bar)
					xrng[1] = 0
				xbig   = pretty(xrng,n=10)
				zuse   = xbig>min(xrng) & xbig<max(xrng)
				xtry   = .su(c(xrng, xbig[zuse]))
				if (mark>xtry[1] && mark<rev(xtry)[1]) {
					xtry = .su(c(mark,xtry))
					zmark  = match(mark,xtry)
				} else
					zmark = NULL
				stry   = scaleVec(V=xtry, Tmin=xlim[1]+ifelse(zero.bar,0,xoff), Tmax=xlim[2]-xoff)
#browser();return()
				excl   = c(zmark,length(stry))
				if (!zero.bar) excl = c(1, excl)
				axis(3, at=stry[-excl], labels=xtry[-excl], cex.lab=1.2)
				if (!is.null(zmark))
					segments(x0=stry[zmark], y0=par()$usr[4], x1=stry[zmark], y1=length(essvec)+yoff, lty=5, col="purple")
				mtext(linguaFranca("Effective Sample Size (bars)",l), side=3, line=1.5, cex=1.2)
			}
			#abline(v=c(1,badhat), lty=1:3, col=c("gainsboro","red"))
			segments(x0=c(1,badhat), y0=par()$usr[3], x1=c(1,badhat), y1=1-yoff, lty=c(1,2), col=c("gainsboro","red"))
			mtext(linguaFranca(ylab,l), side=2, at=1:dim(rhat)[1], las=1, line=0.25, cex=ifelse(length(ylab)>20,0.8,1))
			x0 = rep(0, nrhat); y0 = 1:nrhat
			x1 = rhat[,1]     ; y1 = 1:nrhat
			if(!is.null(ess)) {
				badess = essvec < mark
				goodess = !badess
				if (any(badess))
					addESS (ess=essvec, zvec=badess, xlim=xlim, bcol=.colBlind["vermillion"])
				if (any(goodess))
					addESS (ess=essvec, zvec=goodess, xlim=xlim, bcol=.colBlind["bluegreen"])
			}
			if (any(bad)) {
				x0b=x0[bad]; x1b=x1[bad]; y0b=y0[bad]; y1b=y1[bad]
				segments(x0b, y0b, x1b, y1b, lwd=2, col="red")
				points(x1b, y1b, pch=21, col="red", bg="pink", cex=1)
				text(par()$usr[2], y1b, show0(x1b,3,add2int=T,round2n=T), pos=2, col="red", cex=0.8)
			}
			if (any(good)) {
				x0g=x0[good]; x1g=x1[good]; y0g=y0[good]; y1g=y1[good]
				segments(x0[good],y0[good],x1[good],y1[good], lwd=2, col="blue")
				points(x1[good],y1[good], pch=21, col="blue", bg="cyan", cex=1)
				text(par()$usr[2], y1[good], show0(x1[good],3,add2int=T,round2n=T), pos=2, col="blue", cex=0.8)
			}
		}
		runlab = basename(dir)
		if (!is.null(exlax))
			runlab = paste0(exlax, runlab)
		## Try to map sensitivities  ## too complicated
		#if (exists("strSpp") && strSpp=="405") {
		#	runno = strsplit(runlab,split="\\.")[[1]][[2]]
		#	senno = switch(runno, '28'="B1", '29'="B1", '32'="S01", '33'="S01")
		#}
		addLabel(0.99,0.98, runlab, adj=c(1,0.5), cex=0.9, col="slategrey")
#browser();return()
		box()
		if (png) dev.off()
	}; eop()  ## end lang loop
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcRhat


## compAE ------------------------------2026-02-12
## Compare AE vectors (replaces Excel foofarall)
## ---------------------------------------------RH
compAE <- function(fnlen, fnage, fnmod=NULL, CASAL=0.10, maxage=45,
   strSpp="417", tabs=TRUE, png=FALSE, pngres=400, PIN=c(10,7), lang=c("f","e"))
{
	ici = lenv()
	data ("species", package="PBSdata", envir=ici)
	spcode = species[strSpp,"code3"]
	spname = species[strSpp,"name"]
	f1 = if (is.null(fnlen)) NULL else read.csv(fnlen)
	f2 = if (is.null(fnage)) NULL else read.csv(fnage)
	f3 = NULL
	if (!is.null(fnmod)) {
		f3 = read.csv(fnmod)
		t3 = test=t(f3[,-1])
		dimnames(t3) = list(age=sub("Age\\.","",rownames(t3)), vals=f3[,1])
		f3 = as.data.frame(t3)
		if (strSpp == "405") {  ## use old names
			d3 = read.csv(file.path(dirname(fnmod),"sigAge0.csv"))  ## calcAE: before filtering for records where n>5 (sigAge not used in AEfit anyway)
		} else {
			d3 = read.csv(file.path(dirname(fnmod),paste0("calcAE(", strSpp, ")-sigAge0.csv")))
		}
	}
	ages = (0:maxage) + 0.5
	nage = length(ages)
	for (i in 1:3) {
		fi = switch(i, f1, f2, f3)
		if (is.null(fi)) next
		if (nrow(fi) < nage) {
			message("Sumtingwong") ; browser(); return()
		}
		fi = fi[1:nage,]
		assign(switch(i, "f1", "f2", "f3"), fi)
	}
#browser();return()
	atab = data.frame(
		a   = ages, AE0=rep(0.001, nage) )
	if (!is.null(f1)) {
		atab = data.frame (atab, 
			nL  = f1$N, muL = f1$MU, sdL = f1$SD, cvL = f1$CV,
			AE1 = f1$SDa, AE2 =f1$SDsm )
	}
	if (!is.null(f2)){ 
		atab = data.frame (atab, 
			nA  = f2$N, muA = f2$MU, sdA = f2$SD, cvA = f2$CV,
			AE3 = f2$SDa, AE4 =f2$SDsm )
	}
	if (!is.null(CASAL)){ 
		atab = data.frame (atab, 
		AE5 = ages * CASAL )
	}
	if (!is.null(f3)) {
		atab = data.frame (atab, AE6 = f3$SD)
	}
	useAE=c(T,rep(!is.null(f1),2),rep(!is.null(f2),2),!is.null(CASAL),!is.null(f3))
	AElab = c("noAE", "CVlen", "smCVlen", "CVage", "smCVage", "CASAL", "Punt")[useAE]

	csvout = paste0("./tables/compAE(", strSpp, ")-", paste0(AElab,collapse="+"), ".csv")
	if (tabs) {
		clearFiles(csvout)
		write.csv(atab, file=csvout, row.names=FALSE)
	}
#browser();return()

	## Plot the results
	## --------------------------------------------
	xlim = c(0,maxage)
	ylim = range(atab[,grep("AE",colnames(atab))],na.rm=TRUE)
	if (useAE[7])
		ylim= range(ylim, d3[,"sigAtAge"], na.rm=TRUE)
	ylim[1] = 0
	fout.e = paste0("compAE(", strSpp, ")-", paste0(AElab,collapse="+"))
	for (l in lang) {
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png){ 
			createFdir(l)
			clearFiles(paste0(fout,".png"))
			png(file=paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		expandGraph(mfrow=c(1,1), mar=c(3.5,3.5,1,1), oma=c(0,0,0,0), mgp=c(2,0.5,0))
		plot(NA, xlim=xlim, ylim=ylim, xlab=linguaFranca("Expected age",l), ylab=linguaFranca("Standard deviation",l), cex.axis=1.25, cex.lab=1.5)
		abline(h=if (max(atab[,-1])<=10) seq(1,ylim[2],1) else seq(5,ylim[2],5), col="gainsboro")
		if (useAE[4]) {
			## Add polgons with missing age data
			xnot = rep(NA, length(ages))
			z0   = atab$sdA==0
			xnot[z0] = ages[z0]
			xgrp = split(na.omit(xnot),cumsum(is.na(xnot))[!is.na(xnot)])  ## https://stackoverflow.com/questions/74674411/split-vector-by-each-na-in-r
			xpol = lapply(xgrp, function(x) { c(min(x)-0.5, max(x)+0.5) })
			ypol = lapply(xpol, function(x) { par()$usr[3:4] })
			rubbish = sapply(1:length(xpol), function(i)  {
				x = xpol[[i]]; y = ypol[[i]]
				polygon(x[c(1,1,2,2)], y[c(1,2,2,1)], col=lucent("grey",0.5), border=FALSE)
			})
		}
#browser();return()
		## Add AE1 and AE3
		if (useAE[1]) lines(ages, atab$AE0, col="blue", lty=3, lwd=2)
		if (useAE[2]) lines(ages, atab$AE1, col="blue", lty=3)
		if (useAE[4]) lines(ages, atab$AE3, col="red", lty=3)
		if (useAE[2]) points(ages, atab$AE1, col="blue", pch=20)
		if (useAE[4]) points(ages, atab$AE3, col="red", pch=18)
		## Add AE2 and AE4
		if (useAE[3]) lines(ages, atab$AE2, col="blue", lty=1, lwd=2)
		if (useAE[5]) lines(ages, atab$AE4, col="red", lty=1, lwd=2)
		## Add AE5 and AE6
		if (useAE[6]) lines(ages, atab$AE5, col="purple", lty=2, lwd=3)
		if (useAE[7]) { ## AgeingError using Punt model
			xoff = 0
			points(d3[,"Primary_Age"] + xoff, d3[,"sigAtAge"], col="blue", pch=16)
			lines(ages, atab$AE6, col="forestgreen", lty=4, lwd=3)
		}
		if (all(useAE)) {
			legtxt = linguaFranca(c("AE1 : length-at-age CV", "AE2 : loess-smooth AE1", "AE3 : precision-reader CV", "AE4 : loess-smooth AE3", "AE5 : CASAL constant CV=0.1", "AE6 : AgeingError (Punt model)"), l)
			addLegend(0.20, 0.97, col=c("blue","blue","red","red","purple","forestgreen"), pch=c(20,NA,18,NA,NA,NA), lty=c(3,1,3,1,2,4), lwd=c(1,2,1,2,2,2), seg.len=4, legend=legtxt, bty="n", xjust=0, title="SD derivation", title.font=2, title.adj=0, title.cex=1.2)
		}
		if (sum(useAE)==2 && useAE[1] && useAE[7]) {
			legtxt = linguaFranca(c("no ageing error", "AgeingError (Punt) model input", "AgeingError (Punt) model fit"), l)
			addLegend(0.10, 0.9, col=c("blue","blue","forestgreen"), pch=c(NA,16,NA), lty=c(3,NA,4), lwd=c(2,1,2), seg.len=4, legend=legtxt, bty="n", xjust=0, title=linguaFranca(spname,l), title.font=2, title.adj=0, title.cex=1.2)
		}
#browser();return()
		box()
		if (png) dev.off()
	}; eop()
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compAE


## compTS-------------------------------2022-07-04
## Compare time series medians
## ---------------------------------------------RH
compTS <- function(x, runs=c(24,39), val="Rtdev", type="bars")
{
	yrs = as.numeric(dimnames(x)[[2]])
	tslist= list()
	for (i in runs) {
		ii = paste0("R", pad0(i,2))
		its = x[grep(paste0("^",i),dimnames(x)[[1]]),,val]
		tslist[[ii]] =  apply(its,2,median)
	}
#browser();return()
	tsdf = do.call("cbind", lapply(tslist, data.frame, stringsAsFactors=FALSE))
	colnames(tsdf) = names(tslist)
	barplot(t(as.matrix(tsdf)),beside=T,col=c("black","green"),space=c(0,0))
	addLegend(0.05, 0.95, fill=c("black","green"), legend=colnames(tsdf), title=val, bty="n")
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compTS


## compWT-------------------------------2026-02-11
##  Compare weights or lengths among survey series
## ---------------------------------------------RH
compWT <- function(dat, ssid, sex=2, wfld="len", ylim=c(0,60),
   Nmin=10, good=c(1,3,4,16,21,79), bad=c(22,36), outnam, strSpp="417",
   lang="e", png=FALSE, pngres=400, PIN=c(9,8))
{
	fnam = as.character(substitute(dat))
	spp  = sub("bio","",fnam)
	spp  = strSpp
	data("species", package="PBSdata", envir=.PBStoolEnv)
	spp3 = ttcall(species)[spp,"code3"]
	#ssid.all = c('1'="QCS synoptic", '2'="HS assemblage", '3'="HS synoptic", '4'="WCVI synoptic", '14'="IPHC longline", '16'="WCHG synoptic",
	#	'21'="GIG historical", '22'="HBLL north", '34'="SoG Hake", '36'="HBLL south", '39'="IRFLL north", '40'="IRFLL south",
	#	'48'="Dogfish", '68'="Hake acoustic", '79'="NMFS triennial", '670'="Shrimp trawl", '820'="Jig surveys")
	ssid.all = gfb.codes$ssid

	dat = dat[dat[,wfld] > 0 & !is.na(dat[,wfld]),]
	dat = dat[is.element(dat$sex, sex),]  ## qualify by sex
	if (missing(ssid)) {                  ## ssid = survey series ID
		ssid = .su(dat$SSID)
	}
	nfish = table(dat$SSID)
	zfish = nfish >= Nmin
	ussid = names(nfish)[zfish] ## use SSIDs with at least Nmin fish specimens
#browser();return()
	dat = dat[is.element(dat$SSID, ussid ),]
	wts = split(dat[,wfld], dat$SSID)
	z1  = is.element(names(wts), as.character(good))
	z2  = is.element(names(wts), as.character(bad))
	boxfill = rep("gainsboro", length(wts))
	boxfill[z1] = "green"
	boxfill[z2] = "coral"
	medcol = rep("slategray", length(wts))
	medcol[z1] = "blue"
	medcol[z2] = "black"
	names(wts) = ssid.all[names(wts)]
	pars = list(outpch="+", outcol="gray80", outcex=1.2, boxfill=boxfill, medcol=medcol)

	## PJS suggest adding # records
	nfish = sapply(wts,countVec)


	createFdir(lang)
	if(missing(outnam))
		outnam = paste0("compWT(", spp, ")-", switch(sex,'1'="Male-", '2'="Female-", ""), switch(wfld, 'len'="Length", 'wt'="Weight", 'age'="Age"), "-by-SSID")
#browser();return()
	fout.e = outnam
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		expandGraph(mfrow=c(1,1), mar=c(ifelse(l=="e",8,8),3,1,1))
		out = quantBox (wts, pars=pars, xaxt="n", las=1, ylim=ylim, boxwex=0.6)
		text(1:length(out$n), out$stats[3,], round(out$stats[3,],0), pos=3, offset=0.2, col=medcol, font=2, cex=0.8)
#browser();return()
		xlab = linguaFranca(names(wts),l)
		xlab = sub(" du merlu","\ndu merlu",sub("merlu du ","merlu du\n",sub("morue ","morue\n",xlab)))
		axis(side=1, at=1:length(wts), labels=xlab, las=2, tcl=-0.3)
		mtext(linguaFranca(switch(wfld, 'len'="Length (cm)", 'wt'="Weight (kg)", 'age'="Age (y)"),l), side=2, line=1.75, cex=1.5, las=3)
		text(1:length(wts), 0, labels=format(nfish, big.mark=ifelse(is.null(options()$big.mark), ",", options()$big.mark)), col=sub("black","red",medcol), font=2 )
		addLabel(0.975, 0.975, txt=switch(sex, '1'="Males", '2'="Females"), col=switch(sex, '1'="green4", '2'="darkorchid4"), cex=1.4, font=2, adj=c(1,1))
		addLabel(0.015, 0.029, txt="N=", col="black", cex=1.2, adj=c(0,0))
		if (png) dev.off()
	}; eop()
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~compWT


## mochaLatte---------------------------2024-07-05
##  An alternative to lattice plots (mockLattice)
## ---------------------------------------------RH
mochaLatte <- function(dat, xfld, yfld, ffld, panel, byrow=FALSE, 
   strip=list(col=lucent("black",0.5), bg=lucent("moccasin",0.5), 
   height=0.1, cex=1.4), rc, ...)
{
	opar = par(no.readonly=TRUE)
	on.exit(par(opar))
	options(scipen=5)

	getlab <- function(lim, p=0.1) {
		tck = pretty(lim, n=6)
		if (p<=0) return(as.character(tck))
		lim = lim + (diff(lim) * p * c(1,-1))
		pos = tck>=lim[1] & tck<=lim[2]
		sho = rep("",length(tck))
		sho[pos] = as.character(round(tck[pos],5))  ## RH 200528 -- fornow limit out to 5 decimal places
		return(list(tck=tck, sho=sho))
	}
	facs = as.character(unique(dat[,ffld]))
	nfac = length(facs)
	if (missing(rc))
		rc = .findSquare(nfac)
	#if (rc[1]>rc[2]) rc = rev(rc)  ## more visually pleasing when more columns than rows
	#if (rc[1]<rc[2]) rc = rev(rc)  ## more visually pleasing when more columns than rows
	#strip$cex = strip$cex * (rc[1]^(-0.2))
	dots = list(...)
	if (is.null(dots$mar)) mar=c(0,3.8,0,0)  else mar = dots$mar
	if (is.null(dots$oma)) oma=c(4,2,0.5,1)  else oma = dots$oma
	if (is.null(dots$mgp)) mgp=c(2,0.5,0)    else mgp = dots$mgp
	if (is.null(dots$fn.ylim)) fn.ylim=function(x){range(x,na.rm=T)}  else fn.ylim   = dots$fn.ylim    ## RH 220222
	if (is.null(dots$cex.strip)) strip$cex=strip$cex * (rc[1]^(-0.2)) else strip$cex = dots$cex.strip

	hzero = mar[2]==0  ## plots joined horizontally
	vzero = mar[1]==0  ## plots joined vertically
	if (byrow) {
		par(mfrow=rc, mar=mar, oma=oma, mgp=mgp)
	} else {
		par(mfcol=rc, mar=mar, oma=oma, mgp=mgp)
		## If plotting by column and facs include Males, need jiggery pokery
		if (any(grep("[Mm].le$",facs))) {
			tmat = cbind(matrix(seq(1,prod(rc),2), ncol=2), matrix(seq(2,prod(rc),2), ncol=2)) ## temporary matrix
			tmat = cbind(matrix(seq(1,prod(rc),2), ncol=rc[2]/2), matrix(seq(2,prod(rc),2), ncol=rc[2]/2)) ## temporary matrix
			omat = tmat[,c(seq(1,ncol(tmat),2),seq(2,ncol(tmat),2))] ## ordered matrix
			ovec = as.numeric(omat)
			facs = facs[ovec]
		}
	}
#browser();return()
	for (i in facs) {
		if (is.na(i)) {
			frame(); next } ## needed for mfcol when factors have separate female and male components
		#inum = grep(i,facs)  ## grep doesn't work if fac names contain control characters like '()'
		inum = (1:length(facs))[is.element(facs,i)]  ## RH 230321

		idat = dat[is.element(dat[,ffld],i),]
		if (is.null(dots$xlim)) 
			xlim = range(idat[,xfld], na.rm=TRUE) 
		else
			xlim = dots$xlim
		if (is.null(dots$ylim)){
			yval = idat[,yfld]
			names(yval) = idat[,xfld]
			ylim = round(fn.ylim(yval),5)
		} else {
			ylim = dots$ylim
		}
		yticks    = getlab(ylim, p=0.05)
		ylim[2]   = ylim[2] + diff(ylim)*strip$height ## add space for latte foam
		#do.call("plot", c(list(x=0, y=0, type="n", xlim=xlim, ylim=ylim, xaxt=ifelse(vzero,"n","n"), yaxt=ifelse(hzero,"n","n"), xlab="", ylab=""), dots[setdiff(names(dots), c("xlim","ylim","xlab","ylab","xfac","yfac"))]))
		exclude = c("xlim","ylim","xlab","ylab","xfac","yfac","outline")
		evalCall(plot, c(list(x=0, y=0, type="n", xlim=xlim, ylim=ylim, xaxt=ifelse(vzero,"n","n"), yaxt=ifelse(hzero,"n","n"), xlab="", ylab=""), dots[setdiff(names(dots), exclude)]), checkdef=T, checkpar=T)
		#dots[setdiff(names(dots), c("xlim","ylim","xlab","ylab","xfac","yfac"))]), checkdef=T, checkpar=T)
		#if ((!vzero&&!hzero) || vzero || (hzero && par()$mfg[2]==1)){
		if (byrow) {
			do.xlab = (!vzero&&!hzero) || (hzero && !vzero) || i%in%rev(facs)[1:rc[2]] 
			do.ylab = (!vzero&&!hzero) || (vzero&&!hzero) || (hzero && par()$mfg[2]==1)
		} else {
			do.xlab = (!vzero&&!hzero) || (hzero && !vzero) || inum%%par()$mfg[3]==0 || is.na(facs[inum+1])
			do.ylab = (!vzero&&!hzero) || (vzero&&!hzero) || (hzero && par()$mfg[2]==1)
		}
#browser();return()
		if ( do.xlab ) {
			if (!is.null(dots$xfac)) {
				xticks = unique(dat[,xfld])
				exclude = c("xlim","ylim","xfac","yfac", "outline","xaxt")
				#do.call("axis", c(list(side=1, at=xticks, labels=dots$xfac), dots[setdiff(names(dots),c("xlim","ylim","xfac","yfac"))]))
				evalCall(axis, c(list(side=1, at=xticks, labels=dots$xfac), dots[setdiff(names(dots),exclude)]), checkdef=T, checkpar=T)
			} else {
				xticks = getlab(xlim, p=0.025)
				exclude = c("xlim","ylim","xfac","yfac","exclude","xaxt")
				#do.call("axis", c(list(side=1, at=xticks[["tck"]], labels=xticks[["sho"]]), dots[setdiff(names(dots),c("xlim","ylim","xfac","yfac"))]))
				evalCall(axis, c(list(side=1, at=xticks[["tck"]], labels=xticks[["sho"]]), dots[setdiff(names(dots),exclude)]), checkdef=T, checkpar=T)
			}
		}
		if ( do.ylab ){
			exclude = c("xlim","ylim","xfac","yfac","outline","yaxt")
			#do.call("axis", c(list(side=2, at=yticks[["tck"]], labels=yticks[["sho"]]), dots[setdiff(names(dots),c("xlim","ylim","xfac","yfac"))]))
			evalCall(axis, c(list(side=2, at=yticks[["tck"]], labels=yticks[["sho"]]), dots[setdiff(names(dots),exclude)]), checkdef=T, checkpar=T)
#if(i=="q_2") {browser(); return()}
		}
#if (inum==19) {browser();return()}
		panel(x=idat[,xfld], y=idat[,yfld], dots[setdiff(names(dots),c("xfac","yfac"))])
		#legend("topleft", legend=i, x.intersp=0, box.col=strip$col, bg=strip$bg)
		strip$xbox = par()$usr[c(1,1,2,2)]
		strip$ybox = rep(par()$usr[4],4)
		strip$yoff = c(-diff(par()$usr[3:4])*strip$height,0,0,-diff(par()$usr[3:4])*strip$height)
		strip$ybox = strip$ybox + strip$yoff
		polygon(strip$xbox,strip$ybox,col=strip$bg, border=strip$col)
		text(x=mean(strip$xbox), y=mean(strip$ybox), labels=i, cex=strip$cex)
#if (i=="M_1") {browser();return()}
		box()
	}
	if (!is.null(dots$xlab))
		mtext(text=dots$xlab, side=1, line=par()$oma[1]*0.6, outer=TRUE, cex=ifelse(is.null(dots$cex.lab),1.5,dots$cex.lab))
	if (!is.null(dots$ylab))
		#mtext(text=dots$ylab, side=2, line=(par()$oma[2] + par()$mar[2])*0.6, outer=TRUE, cex=ifelse(is.null(dots$cex.lab),1.5,dots$cex.lab))
		mtext(text=dots$ylab, side=2, line=par()$oma[2]*ifelse(par()$mar[2]>=1.5,0.1,0.6), outer=TRUE, cex=ifelse(is.null(dots$cex.lab),1.5,dots$cex.lab))
	if (!is.null(dots$mainTitle))
		mtext(text=dots$mainTitle, side=3, line=par()$oma[3]*0.2, outer=TRUE, cex=ifelse(is.null(dots$cex.main),1.75,dots$cex.main))
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~mochaLatte


## panelBoxes---------------------------2019-11-24
##  Plot quantile plots using 'nchains' to delimit separate boxes.
##  mcmc=data.frame e.g, 'currentMCMC$P' from object created by 'importMCMC'.
##  Very difficult to manipulate trellis plots (RH)
## -----------------------------------------------
panelBoxes <- function (mcmc, nchains=9, pdisc=0, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, xlim=c(0.25,9.75), 
   boxfill=paste0(rep(c("cyan","green","coral"),each=3), rep(1:3,3)),
   cex.main=1.2, cex.lab=1.2, cex.strip=0.9, cex.axis=0.9, las=0, 
   tck=0.4, tick.number=5, xfac=paste0("B",1:nchains), outline=TRUE,
	lang="e", ...)
{
	panel.box <- function(x, y, ...) {
		dots = list(...)[[1]]
		unpackList(dots)
		if (is.null(dots$outline)) outline = TRUE
		## xlim and ylim determined by 'mochaLatte'
		#if (is.null(dots$xlim)) xlim = range(x,na.rm=TRUE)
		#if (is.null(dots$ylim)) ylim = if (outline) range(y,na.rm=TRUE) else quantile(y, quants5[c(1,5)])
		#if (is.null(dots$xfac)) xfac = unique(x)
		basecol = "slategray"
		#boxfill = paste0(rep(c("cyan","green","coral"),each=3)
		boxpars = list(boxwex=0.5, boxfill=boxfill, boxcol=basecol, outpch=3, outcex=0.3, outcol=lucent(basecol,0.25), medlwd=1, whisklty=1, whiskcol=basecol)
		chainlink = rep(1:nchains,ff)
		chainbox  = split(y,x)
		#quantbox(chainbox, add=T, xaxt="n", yaxt="n", pars=boxpars, ...)
		mess =  sapply(dots,deparse)
		mess = paste0(paste0(names(mess),"=",mess),collapse=",") 
		messy = paste0("list(add=TRUE, xaxt=\"n\", yaxt=\"n\", pars=boxpars, ", mess,")")
		argos = eval(parse(text=messy))
		#do.call("quantbox", args=list(x=chainbox, add=TRUE, xaxt="n", yaxt="n", pars=boxpars, deparse(mess)) )
		do.call("quantbox", args=c(x=list(chainbox),argos) )
		#evalCall(quantbox, args=argos, checkdef=T, checkpar=T )
#browser(); return()
	}

	if (pdisc>0 && pdisc<1)
		mcmc = mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	ff = rep(round(n/nchains),nchains-1)
	ff = c(ff,n-sum(ff))
	p <- ncol(mcmc)
	dat <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Chain=rep(rep(1:nchains,ff),p), Value=as.vector(as.matrix(mcmc)))
	dat$Index = paste(dat$Factor,dat$Chain,sep="-")

	dots=list(...); fn.ylim=dots$fn.ylim
	if (is.null(fn.ylim)){
		if (outline) fn.ylim <- function(x){range(x, na.rm=TRUE)} 
		else         fn.ylim <- function(x){extendrange(sapply(split(x,names(x)), quantile,quants5[c(1,5)], na.rm=TRUE))}
	}
#browser();return()
	mochaLatte(dat, xfld="Chain", yfld="Value", ffld="Factor", panel=panel.box, xlim=xlim, tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang), xfac=xfac, fn.ylim=fn.ylim , outline=outline, byrow=T, ...) #, mar=c(0,3.8,0,0), oma=c(4,3,0.5,1)
	gc(verbose=FALSE)
	invisible(dat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~panelBoxes


## panelChains---------------------------2023-08-22
##  Plots cumulative fequency of 'nchains' by partitioning one trace.
##  Revised from 'plotTracePOP'
##  mcmc=data.frame e.g, 'currentMCMC$P' from object created by 'importMCMC'.
##  Very difficult to manipulate trellis plots (RH)
## -----------------------------------------------
panelChains <- function (mcmc, nchains=3, pdisc=0.1, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1, span=1/4,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
   cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.8, las=0, 
   tck=0.4, tick.number=5, lty.trace=1, lwd.trace=1, col.trace="grey", 
   lty.median=1, lwd.median=1, col.median="black", lty.quant=2, lwd.quant=1, 
   col.quant="black", plot=TRUE, probs=tcall(quants3), lang="e", ...)
{
	panel.chain <- function(x, y, ...) {
		dots = list(...)
		unpackList(dots)
		if (is.null(dots$xlim)) xlim = range(x,na.rm=TRUE)
		if (is.null(dots$ylim)) ylim = range(y,na.rm=TRUE)
		abline (h=0.5, lty=3, lwd=1, col="grey")
		chainlink = rep(1:nchains,ff)
		for (i in 1:nchains) {
			z = is.element(chainlink,i)
			lines(x[z], y[z], lty=rep(lty.trace,nchains)[i], lwd=lwd.trace, col=rep(col.trace,nchains)[i])
			#lines(x[z], y[z], lty=1, lwd=2, col=c("red","green4","blue"))
		}
	}

	if (pdisc>0 && pdisc<1)
		mcmc = mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
	relation <- if (same.limits) "same" else "free"
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	n <- nrow(mcmc)
	ff = rep(round(n/nchains),nchains-1)
	ff = c(ff,n-sum(ff))
	p <- ncol(mcmc)
	dat <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Chain=rep(rep(1:nchains,ff),p), Value=as.vector(as.matrix(mcmc)))

	dat$Index = paste(dat$Factor,dat$Chain,sep="-")
	vList     = split(dat$Value,dat$Index)
	qList     = sapply(vList,function(x){
		xsort  = sort(x)
		xscal  = xsort - min(xsort)
		ycumu  = cumsum(xscal)/sum(xscal)
		out    = cbind(x=xsort,y=ycumu)
		return(out)
	}, simplify = FALSE )
	dat$CumFreq = dat$ValueSort = NA
	for (i in names(qList)) {
		z = is.element(dat$Index,i)
		dat$ValueSort[z] = qList[[i]][,"x"]
		dat$CumFreq[z]   = qList[[i]][,"y"]
	}
	#mar=c(0,3,0,0), oma=c(4,3,0.5,1)
	#mochaLatte(dat,xfld="ValueSort",yfld="CumFreq",ffld="Factor", panel=panel.chain, ylim=c(0,1), mar=c(2,0,0,0), oma=c(1.5,4.5,0.5,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang))
	#mochaLatte(dat,xfld="ValueSort",yfld="CumFreq",ffld="Factor", panel=panel.chain, ylim=c(0,1), mar=c(2,2,0,0), oma=c(4,4,0.5,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang))
	mochaLatte(dat, xfld="ValueSort", yfld="CumFreq", ffld="Factor", panel=panel.chain, ylim=c(0,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang), byrow=T, ...)
	invisible(dat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~panelChains


## panelTraces-------------------------2022-02-22
##  Plots the sequential  trace of MCMC samples 
##  with running median and (0.05, 0.95) quantiles.
## ---------------------------------------------RH
panelTraces <- function (mcmc, mpd=mcmc[1,], nchains=1, pdisc=0, 
   axes=FALSE, same.limits=FALSE, between=list(x=axes,y=axes), div=1,
   log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
   cex.main=1.2, cex.lab=1.2, cex.strip=0.9, cex.axis=0.9, las=0, 
   tck=0.4, tick.number=5, xfac=NULL, 
	lang="e", ...)
{
	panel.trace <- function(x, y, ...) {
		dots = list(...)
		unpackList(dots)
		if (is.null(dots$xlim)) xlim = range(x,na.rm=TRUE)
		if (is.null(dots$ylim)) ylim = range(y,na.rm=TRUE)
		lty.trace  = 1;  lwd.trace  = 1;  col.trace  = "grey"
		lty.quant  = 2;  lwd.quant  = 1;
		col.quant  = "black"
		lty.median = 1;  lwd.median = 1.5;  col.median = "blue"
		lines(x, y, lty=lty.trace, lwd=lwd.trace, col=col.trace)
		z = !is.na(y); xx=x[z]; yy=y[z]  ## RH 220222
#browser();return()
		if (any(is.finite(yy)) && var(yy) > 0) {
			lines(xx, cquantile.vec(yy, prob=tcall(quants3)[1]), lty=lty.quant,  lwd=lwd.quant,  col=col.quant)
			lines(xx, cquantile.vec(yy, prob=tcall(quants3)[2]), lty=lty.median, lwd=lwd.median, col=col.median)
			lines(xx, cquantile.vec(yy, prob=tcall(quants3)[3]), lty=lty.quant,  lwd=lwd.quant,  col=col.quant)
			points(xx[1], mpd[getNpan()], pch=21, col="black", bg="red", cex=2)
		}
	}

	if (pdisc>0 && pdisc<1)
		mcmc = mcmc[(round(pdisc*nrow(mcmc))+1):nrow(mcmc),]  # get rid of the first 'pdisc' (e.g., 10%)
	if (is.null(dim(mcmc))) {
		mcmc.name <- rev(as.character(substitute(mcmc)))[1]
		mcmc <- matrix(mcmc, dimnames=list(NULL, mcmc.name))
	}
	mcmc <- if (log) 
		log(mcmc/div, base=base)
	else mcmc/div
	mcmc <- as.data.frame(mcmc)
	ylim = NULL
	if (same.limits)
		ylim = range(mcmc,mpd,na.rm=TRUE)
	n <- nrow(mcmc)
	ff = rep(round(n/nchains),nchains-1)
	ff = c(ff,n-sum(ff))
	p <- ncol(mcmc)
	dat <- data.frame(Factor=ordered(rep(names(mcmc), each=n), 
		names(mcmc)), Draw=rep(1:n, p), Chain=rep(rep(1:nchains,ff),p), Value=as.vector(as.matrix(mcmc)))
	#dat$Index = paste(dat$Factor,dat$Chain,sep="-")
#browser();return()
	#mochaLatte(dat, xfld="Draw", yfld="Value", ffld="Factor", panel=panel.trace, xlim=c(0,nrow(mcmc)), ylim=ylim, mar=c(0,3,0,0), oma=c(4,4,0.5,1), tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang) )
	mochaLatte(dat, xfld="Draw", yfld="Value", ffld="Factor", panel=panel.trace, xlim=c(0,nrow(mcmc)), ylim=ylim, tcl=-0.3, las=1, cex.axis=cex.axis, cex.lab=cex.lab, cex.strip=cex.strip, xlab=linguaFranca(xlab,lang), ylab=linguaFranca(ylab,lang), byrow=T, ... )
#browser();return()
	invisible(dat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~panelTraces


## plot_cres ---------------------------2026-02-23
##  Produce 5-panel overview of residuals (pkg=compResidual) by default,
##  or one specific plot if the 'pick_one' argument is provided.
## Changed name from 'plot.cres.new' to 'plot_cres' to avoid methods crap (RH 260220)
## ---------------------------------------VT|AN|RH
plot_cres <- function (x, pick_one=NULL, maxLag=NULL, main="",
   png=FALSE, pngres=400, PIN=c(10,9), lang=c("f","e"), outnam=NULL, ...) 
{
	## ---------------------
	## Plot 1: bubble plot of the residuals
	## Plot 2: autocorrelation and cross-correlation plot for 
	##  all dimensions of the residual matrix (rows, columns and diagonal)
	## Plot 3: Q-Q plot of the residuals, 
	##  each color is a row of the residual matrix (i.e. compositional group) 
	## Plot 4 & 5: simple plot of the residuals,
	##  [deprecated] each color is a row of the residual matrix (i.e. compositional group)
	##  two residual plots by year and by age.
	## ---------------------
	## first adapted from compResidual package : 2024-10-02
	## ---------------------

	op = par(no.readonly=TRUE)
	on.exit(par(op))

	if (!is.null(rownames(x))) 
		yname <- rownames(x)
	else yname = 1:nrow(x)
	if (!is.null(colnames(x))) 
		xname <- colnames(x)
	else xname = 1:ncol(x)
	col <- rep(1:nrow(x), ncol(x))

	dots = list(...)
	probs = dots$probs
	if (is.null(probs))  probs = seq(0,1,0.1)
	brk = dots$brk
	if (is.null(brk))    brk = 0.2
	bscale = dots$bscale
	if (is.null(bscale)) bscale = 0.4
	udots = dots[setdiff(names(dots), c("oma","mar","probs","brk","bscale"))]

	createFdir(lang)
	if (is.null(outnam))
		fout.e = "OSA-residuals"
	else
		fout.e = outnam

	for (l in lang) {  ## to implement 'linguaFranca'.
		changeLangOpts(L=l)
		#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		if (missing(pick_one)) {
			#expandGraph(mfrow=c(2,2), mar=c(3.3, 3.5, 0.5, 0.75), oma=c(0.5, 0, 2, 0), cex.lab=1.5, mgp=c(2.2,0.5,0), las=1 )
			expandGraph(mfrow=c(1,1), mar=c(3, 3, 0.1, 0.5), oma=c(0,0,1.5,0), cex.axis=1, cex.lab=1.2, mgp=c(1.75,0.5,0), las=1 )
		}
		else {
			if (!pick_one %in% 1:4) 
				stop("pick_one should be a number between 1 and 4")
			par(mfrow=c(1,1))
		}

		## OSA residuals by age (row) and year (column)
		if (pick_one == 1 || missing(pick_one)) {
			if (missing(pick_one))
				par(fig=c(0,0.5,0.4,1), new=FALSE)
			#add_legend <- function(x, cex.text=1, ...) {
			add_legend <- function(x, cex.text=1, ...) {
				zscale <- pretty(x, min.n = 4)
				uu <- par("usr")
				#yy <- rep(uu[3] + 0.03 * (uu[4] - uu[3]), length(zscale))
				yy <- rep(uu[4] - 0.05 * (uu[4] - uu[3]), length(zscale))
				xx <- seq(uu[1] + 0.175 * (uu[2] - uu[1]), uu[1] + 0.475 * (uu[2] - uu[1]), length=length(zscale))
				yl <- uu[4] - 0.04 * (uu[4] - uu[3])
				xl <- uu[1] + 0.025 * (uu[2] - uu[1])
				text(xl, yl, labels=linguaFranca("Scale",l), cex=1.2*cex.text, font=2, adj=0)
				text(xx, yy, labels = zscale, cex = cex.text, font=2, pos=3)
				colb <- ifelse(zscale < 0, rgb(1, 0, 0, alpha = 0.5), rgb(0, 0, 1, alpha = 0.5))
				points(xx, yy, cex = sqrt(abs(zscale))/max(sqrt(abs(zscale)), na.rm = TRUE) * 5 * bscale, pch = 19, col = colb)
			}
			sample <- rep(1:ncol(x), each = nrow(x))
			composition <- rep(1:nrow(x), ncol(x))
			xlim = extendrange(sample, f=0.025)
			ylim = extendrange(composition, f=c(0.02,0.06))
			#plotby(sample, composition, x, bubblescale=bubblescale, xlab=linguaFranca("Year",l), ylab=linguaFranca("Age (y)",l), yaxt="n", xaxt="n", xlim=xlim, ylim=ylim, ...)
			do.call(plotby, args=list(x=sample, y=composition, z=x, bubblescale=bscale, xlab=linguaFranca("Year",l), ylab=linguaFranca("Age (y)",l), yaxt="n", xaxt="n", xlim=xlim, ylim=ylim, udots))
			xTicks <- pretty(1:ncol(x))
			xTicks <- xTicks[xTicks != 0]
			yTicks <- pretty(1:nrow(x))
			yTicks <- yTicks[yTicks != 0]
			axis(1, at = xTicks, labels = xname[xTicks])#, cex.axis=1.2)
			axis(2, at = yTicks, labels = yname[yTicks])#, cex.axis=1.2)
			add_legend(x, cex.text = 0.8, bubblescale=bscale)
			addLabel(0.95,0.95, linguaFranca("OSA residuals",l), cex=1, adj=c(1,0))
#browser();return()
		}
		## Autocorrelation plot (remove lag 0)
		if (ncol(x)>1 && (pick_one == 2 || missing(pick_one))) {
			if (missing(pick_one))
				par(fig=c(0,0.5,0,0.4), new=TRUE)
			mess = c(
				"acfr <- compResidual:::acf_res(x, \"row\")",
				"acfc <- compResidual:::acf_res(x, \"column\")",
				"acfd <- compResidual:::acf_res(x, \"diagonal\")"  )
			eval(parse(text=paste0(mess, collapse="; ")))
			ci <- 1.96/sqrt(length(x))
			## Get rid of lag 0
			acfr$acf = acfr$acf[-1,,,drop=F]; acfr$lag = acfr$lag[-1,,,drop=F]
			acfc$acf = acfc$acf[-1,,,drop=F]; acfc$lag = acfc$lag[-1,,,drop=F]
			acfd$acf = acfd$acf[-1,,,drop=F]; acfd$lag = acfd$lag[-1,,,drop=F]
			maxlag <- c(max(acfr$lag), max(acfc$lag), max(acfd$lag))
			if (is.null(maxLag)) 
				maxLag <- min( c(20, max(maxlag) ) )
			pmaxlag = pmin(maxLag,maxlag)
			low  <- as.numeric( min( c(acfr$acf[1:pmaxlag[1],,], acfc$acf[1:pmaxlag[2],,], acfd$acf[1:pmaxlag[3],,]), na.rm=T ) )
			ylow <- min(c(-ci, low))
			high <- as.numeric( max( c(acfr$acf[1:pmaxlag[1],,], acfc$acf[1:pmaxlag[2],,], acfd$acf[1:pmaxlag[3],,]), na.rm=T ) )
			yhigh <- max(c(ci, high))
			ylim = c(ylow, yhigh); ylim = ylim + c(-1,1) * 0.05*diff(ylim)
			#lcol = c("#7fc97f", "#beaed4", "#fdc086")
			lcol = .colBlind[c("redpurple","bluegreen","vermillion")]
			#plot(y=as.vector(acfr$acf), x=as.vector(acfr$lag), type="h", ylim=ylim, xlim=c(0, maxLag), col=lcol[1], lwd=3, ylab=linguaFranca("ACF",l), xlab=linguaFranca("Lag",l), cex.axis=1.2, ...)
#browser(); return()
			do.call(plot.default, args=list(y=as.vector(acfr$acf), x=as.vector(acfr$lag), type="h", ylim=ylim, xlim=c(0, maxLag), col=lcol[1], lwd=3, ylab=linguaFranca("ACF",l), xlab=linguaFranca("Lag",l), log="", udots)) #, cex.axis=1.2
			lines(y=as.vector(acfc$acf), x=as.vector(acfc$lag) + 0.15, type="h", col=lcol[2], lwd=3)
			lines(y=as.vector(acfd$acf), x=as.vector(acfd$lag) + 0.3, type="h", col=lcol[3], lwd=3)
			abline(h=0)
			abline(h=c(-ci, ci), lty=2)
			legend("topright", col=lcol, legend=linguaFranca(c("row","column","diagonal"),l,little=8), bty="n", lwd=2)
		} else {
			frame(); addLabel(0.5,0.5,"No ACF", col="slategrey", cex=3)
		}
		## Quantile-quantile plot
		if (pick_one == 3 || missing(pick_one)) {
			if (missing(pick_one))
				par(fig=c(0.5,1,0.67,1), new=TRUE)
			qqnorm(x, col=col, main="", xlab=linguaFranca("Theoretical quantiles",l), ylab=linguaFranca("Sample quantiles",l), pch=20, cex=1.2)#, cex.axis=1.2)
			abline(0, 1)
			#legend("topleft", col=col, legend=yname, pch=20, bty="n", cex=1, pt.cex=1.2, ncol=4)
			legend("topleft", legend=linguaFranca("Colours indicate age",l), bty="n", cex=1)
		}
		## Residuals by matrix element
		if (pick_one == 4 || missing(pick_one)) {
			#if (missing(pick_one))
			#	par(fig=c(0.5,1,0.33,0.67), new=TRUE)
			#plot(as.vector(x), col=col, ylab=linguaFranca("Residuals",l), xlab=linguaFranca("Matrix index",l), pch=20, cex=1.2, cex.axis=1.2)

			## Subfunctions--------------------------
			decimalplaces <- function(x) { ## https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r
#browser();return()
				if ((x %% 1) != 0) {
					period = ifelse(exists("l"), switch(l,'f'=",",'e'="."), ".")
					nchar(strsplit(sub('0+$', '', as.character(x)), period, fixed=TRUE)[[1]][[2]])
				} else {
					return(0)
				}
			}
			easyBake <- function(rval, bycol=TRUE, brk=0.2) {
				if (!bycol) rval = t(rval)
				ncol = ncol(rval)
				cnam = colnames(rval)
				if (is.null(cnam)) cnam = as.character(1:ncol)
				lval = list()
				for (i in 1:ncol(rval)) lval[[i]] = rval[,i]
				names(lval) = colnames(rval)
				lout = lapply(lval, function(z){
					xval   = table(ceiling(z/brk) * brk)
					xval   = xval/max(xval)
					## https://stackoverflow.com/questions/15236440/as-numeric-with-comma-decimal-separators
					yval   = as.numeric(sub(",",".",names(xval)))  ## unusual tweak for as.numeric
					yvalue = round(seq(yval[1], rev(yval)[1], brk), decimalplaces(brk))
					xvalue = rep(0,length(yvalue))
					names(xvalue) = yvalue
					xvalue[names(xval)] = xval
if(length(xvalue)!=length(yvalue)){browser(); return()}
#browser();return()
					ymed = median(z,na.rm=TRUE)
					return(list(xval=xvalue, yval=yvalue, ymed=ymed))
				})
			}
			##---------------------------------------
			xx    = as.matrix(x)
			#qages = apply(xx,1,quantile, probs=probs)
			#qyear = apply(xx,2,quantile, probs=probs)
			for (bycol in c(TRUE, FALSE)) {
				out = easyBake(xx, bycol=bycol, brk=brk)
				if (is.null(names(out))) names(out) = 1:length(out)  ## ages
#if (!bycol) {browser();return()}
				nout = as.numeric(names(out)[1]):as.numeric(rev(names(out))[1])
				allout = as.list(rep(NA,length(nout)))
				names(allout) = nout
				allout[names(out)] = out  ## just use out (not enough plotting space)
				## Plot this mess around
				if (missing(pick_one))
					par(fig=c(0.5,1,ifelse(bycol,0.33,0),ifelse(bycol,0.67,0.33)), new=TRUE)
				xlim=c(0.5,length(out)+ 0.5)
				ylim =  range(lapply(out, function(x) {range(x$yval)}))
				ymed =  lapply(out, function(x) {x$ymed})
				ymed = do.call("rbind", lapply(ymed, data.frame, stringsAsFactors=FALSE))
				plot(0,0, type="n", xlim=xlim, ylim=ylim, xlab=ifelse(bycol,"",linguaFranca("Age",l)), ylab=linguaFranca("Residuals",l), xaxt="n")
				abline(h=c(0,-2,2), col="black", lwd=c(1,1,1), lty=c(1,2,2))
				for(i in 1:length(out)){
					ixL = i - (out[[i]]$xval * 0.5)
					ixR = i + (out[[i]]$xval * 0.5)
					iy = out[[i]]$yval
					polygon(x=c(ixL,rev(ixR)), y=c(iy,rev(iy)), border=ifelse(bycol,"royalblue","sienna"), col=ifelse(bycol,"aquamarine","moccasin"), lwd=0.5)
				}
				points(1:nrow(ymed), ymed[,1], pch=21, cex=0.7, lwd=0.5, col="black", bg=ifelse(bycol,"blue","orangered"))
				nfac = length(out)
				zint = ifelse(bycol,ifelse(nfac>15,3,1),5)
				zlab = .su(c(1,seq(zint,nfac-1,zint),nfac))
				#zlab = floor(seq(1,nfac,length=12))
				axis(side=1, at=1:nfac, labels=FALSE, tick=T, tcl=-0.1)
				axis(side=1, at=zlab, labels=names(out)[zlab], tick=T, tcl=-0.3, las=2) ## rotate labels
#browser();return()
				addLabel(0.99,0.99, txt=paste0(linguaFranca("breaks",l,little=6), " = ", brk), cex=0.8, adj=c(1,1))
			}
		} ## end 4th and/or 5th plot
		mtext(linguaFranca(main,l), outer=TRUE, side=3, line=0, cex=1.5)
#browser(); return()
		if (png) dev.off()
	}; eop()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot_cres


## plotACFs-----------------------------2023-08-31
##  Plot ACFs for the estimated parameters.
##  Control eps and png from PBScape.r in plt.mcmcGraphs
##----------------------------------------------RH
plotACFs <- function(mcmc, lag.max=60, lang="e", ...) #, ptypes=tcall(PBSawatea)$ptype, pngres=400)
{
	#if (!is.null(dev.list())) on.exit(expandGraph(mfrow=c(1,1)))
	acfs  = apply(mcmc, 2, function(x){acf(x,lag.max=lag.max,plot=FALSE)$acf})
	ylim  = range(acfs[round(acfs,5)>-1 & round(acfs,5)<1])
	idx   = apply(mcmc, 2, allEqual)
	mcmcP = mcmc[,!idx,drop=FALSE]
	dots  = list(...)
	rc    = dots$rc
	if (is.null(rc))
		rc = .findSquare(ncol(mcmcP))
#browser();return()

	#for (p in ptypes) {
		#if (p=="eps") postscript("paramAcf.eps", width=8, height=8, horizontal=FALSE,  paper="special")
		#else if (p=="png") png("paramAcf.png", width=8, height=8, units="in", res=pngres)
		#expandGraph(mfrow=rc, mar=c(0.5,3.5,0,0), oma=c(4,0.5,0.5,0.5))
		expandGraph(mfrow=rc, mar=c(0.5,0.5,0,0), oma=c(3.7,4,0.5,0.5))
		sapply(1:ncol(mcmcP), function(i){
			ii  = colnames(mcmcP)[i]
			mcP = mcmcP[,i]
			xy = acf(mcP, lag.max=lag.max, ylim=ylim, xaxt="n", yaxt="n", xlab="", ylab="", plot=F)
			xy$lag = xy$lag[2:(lag.max+1),1,1,drop=FALSE]
			xy$acf = xy$acf[2:(lag.max+1),1,1,drop=FALSE]
			#ylim = extendrange(xy$acf)
			plot(xy, xlim=c(1,60), ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n", type="n")
			abline(h=setdiff(seq(-0.2,0.8,0.2),0), col="slategray3", lty=3)
			do.call(lines, args=c(list(x=1:lag.max, y=xy$acf, type="h"), dots[grep("rc",names(dots),invert=T)]))
#browser();return()
			#evalCall(lines, argu=c(list(x=1:lag.max, y=xy$acf, type="h"), dots[grep("rc",names(dots),invert=T)]), checkdef=T)  ## doesn't seem to work [needs debugging]
			#lines(1:lag.max, xy$acf, type="h", ...)
			xgloss = ifelse(i > (ncol(mcmcP)-rc[2]),TRUE,FALSE)
			axis(1, labels=xgloss, tick=xgloss, cex.axis=1.2)
			ygloss = ifelse(par()$mfg[2]==1,TRUE,FALSE)
			axis(2, labels=ygloss, tick=ygloss, cex.axis=1.2, las=1)
			addLabel(0.95, 0.975, ii, cex=ifelse(rc[2]<4,1.25,1.0), adj=c(1,1), col="navy")
			box(lwd=1.5)
		})
		mtext(linguaFranca("Lag",lang),side=1,outer=TRUE,line=2.25,cex=1.5)
		mtext(ifelse(lang=="f","FAC","ACF"),side=2,outer=TRUE,line=2.25,cex=1.5)
		#dev.off()
	#}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotACFs


## plotBh ------------------------------2025-11-21
##  Plot biomass (B) versus steepness (h)
## ---------------------------------------------RH
plotBh <- function(dmcmc, Bmcmc, Bmsy, xlim, xfld="BH_h", ytype="rel",
   currYr, base, top, hbin=0.015, png=FALSE, pngres=400, PIN=c(10,7.5))
{
	## Subroutine----------------------------------
	addBox <- function(x, y, q=c(0.05,0.25,0.5,0.75,0.95), outs=F, boxcol, lwd=2, width=0.0075, offset=0.0075) ## (RH 251121)
	{
		## Assume x is scalar and y is vector ( for now)
		yqnt = quantile(y, q)
		yrng = range(y, na.rm=TRUE)
		xdif = diff(par()$usr[1:2]) * width
		xadj = x + offset
		xpol = c(xadj[c(1,1)] - xdif, xadj[c(1,1)] + xdif)
		ypol = yqnt[c(2,4,4,2)]
		polygon(xpol, ypol, border=boxcol, col=darkenRGB(boxcol,-0.8), lwd=lwd)
		segments(x0=xadj[c(1,1)], y0=yqnt[c(2,4)], x1=xadj[c(1,1)], y1=yqnt[c(1,5)], col=boxcol, lwd=lwd)
		segments(x0=xadj-xdif, y0=yqnt[c(3)], x1=xadj+xdif, y1=yqnt[c(3)], col=boxcol, lwd=lwd)
		if (outs) {
			segments(x0=xadj[c(1,1)], y0=yqnt[c(1,5)], x1=xadj[c(1,1)], y1=yrng[c(1,2)], col=boxcol, lwd=1, lty=3)
		}
	}
	## ----------------------------------subroutine

	## Gather data
	df = data.frame(x=dmcmc[,xfld], B0=Bmcmc[,1], Bmsy=Bmsy, Bcurr=Bmcmc[,paste0("SSB_",as.character(currYr))])
	df$Bcurr_B0   = df$Bcurr/df$B0
	df$Bcurr_Bmsy = df$Bcurr/df$Bmsy
	df$xbin = ceiling(df$x / hbin) * hbin
	xbin = table(df$xbin)
	df$y1 = switch(ytype,'abs'=df$B0, 'rel'=df$Bcurr_B0)
	df$y2 = switch(ytype,'abs'=df$Bmsy, 'rel'=df$Bcurr_Bmsy)
	ylim  = range(c(df$y1, df$y2), na.rm=TRUE)
	if (missing(xlim))
		xlim = range(df$x, na.rm=TRUE)
	
	## Find minimum h where MSY value is estimable
	z = !is.na(df$y2) & is.finite(df$y2)
	hlow = df$x[z][findPV(min(df$x[z]),df$x[z])]

	if (missing(top))
		top = quantile(c(df$y1, df$y2), probs=0.99, na.rm=T)
	if (missing(base))
		base = -(top * 0.2)

	## Start plotting
	fout.e = paste0("plotBh-", switch(ytype, 'abs'="absolute", 'rel'="relative",""))
	for (l in lang) { 
		createFdir(lang=l)
		changeLangOpts(L=l)
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (png) {
			clearFiles(paste0(fout,".png"))
			png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		}
		expandGraph(mfrow=c(1,1), mar=c(3.0,3.5,0.75,1), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
		plot(0, 0, type="n", xlim=xlim, ylim=c(base,top), xlab="Steepness", ylab=paste(switch(ytype,'abs'="Absolute", 'rel'="Relative",""), "Biomass"), cex.axis=1.2, cex.lab=1.5)
		#abline(v=hlow, lty=2, col="slategray")
		grey = "darkslategray4"
		segments(x0=hlow, y0=0, x1=hlow, y1=par()$usr[4]-(0.2*diff(par()$usr[3:4])), lty=2, col=grey)
		text(hlow, 0, labels=show0(hlow,2,round2n=TRUE), col=grey, pos=1, cex=0.9)
		points(df$x,df$y1, pch=20, col="lightblue")
		points(df$x,df$y2, pch=20, col="pink")
		lines(loess.smooth(df$x,df$y1),col="blue",lwd=3)
		lines(loess.smooth(df$x,df$y2),col="red",lwd=3)
		addBox(x=max(df$x,na.rm=T), y=df$y1, boxcol="blue", offset=0.015)
		addBox(x=max(df$x,na.rm=T), y=df$y2, boxcol="red", offset=0.02)
		rout = scaleVec(xbin, Tmin=base, Tmax=0)
		drawBars(as.numeric(names(rout)), base=base, rout, col=grey, fill=lucent("green",0.5), width=hbin)
#browser();return()
		legtxt = switch(ytype,
			'abs' = c(expression(italic(B)[0]), expression(italic(B)[MSY])),
			'rel' = c(bquote(italic(B)[.(currYr)]/italic(B)[0]), bquote(italic(B)[.(currYr)]/italic(B)[MSY]))
		)
		htxt = sub("h","italic(h)",sub("_","~~",xfld))
		legtxt = c(parse(text=htxt), legtxt)
		legend("topleft", legend=legtxt, pch=c(22,NA,NA), col=c("black","blue","red"), pt.bg=c(lucent("green",0.25),NA,NA), pt.cex=c(3,NA,NA), lty=c(NA,1,1), lwd=c(1,3,3), bty="n", cex=1.5, x.intersp=c(1,1,1), inset=0.025)
		if (png) dev.off()
	}; eop()  ## end language loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotBh


## plotSS.dmcmc-------------------------2025-09-16
##  Plot MCMC diagnostics (traces, split chains, ACFs).
##  Functions in PBSawatea.
## -----------------------------------PBSawatea|RH
plotSS.dmcmc <- function(mcmcObj, mpdObj, ptypes, lang, pngres=400, PIN=c(9,9),
	BRP=rep(TRUE,3), rhat=TRUE, steep=TRUE)
{
	so("panelTraces.r","synth"); so("mochaLatte.r","synth"); so("calcRhat.r","synth")
	fout.e = "sumtingwong"

#sumting = T
#if (sumting) {

	## MCMC diagnostics
	## ----------------
	if (BRP[1]) {
		fout.e = "traceBiomass"
		for (l in lang) {
			changeLangOpts(L=l)
			#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
			bmcmc = mcmcObj$B[getYrIdx(colnames(mcmcObj$B))]/1000
			bmpd  = mpdObj$B[getYrIdx(names(mpdObj$B))]/1000
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
				else if (p=="png") {
					clearFiles(paste0(fout,".png"))
					png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				}
				panelTraces(mcmc=bmcmc, mpd=bmpd, xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, same.limits=F, lang=l)
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
#browser();return()

	if (BRP[2]) {
		fout.e = "traceRecruits"
		for (l in lang) {
			changeLangOpts(L=l)
			#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
			rmcmc = mcmcObj$R[getYrIdx(colnames(mcmcObj$R))]/1000
			rmpd  = mpdObj$R[getYrIdx(names(mpdObj$R))]/1000
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
				else if (p=="png") {
					clearFiles(paste0(fout,".png"))
					png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				}
				panelTraces(mcmc=rmcmc, mpd=rmpd, xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, same.limits=F, lang=l, mar=c(0,1.5,0,1), oma=c(3.5,2.5,0.5,0), mgp=c(1.75,0.5,0))
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}

#} ## end sumting

	if (BRP[3]) {
		fout.e = "traceParams"
		for (l in lang) {
			changeLangOpts(L=l)
			#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
			pmcmc = mcmcObj$P
			names(pmcmc) = linguaFranca(names(pmcmc), l)
			pmpd  = mpdObj$P
			for (p in ptypes) {
				if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
				else if (p=="png") {
					clearFiles(paste0(fout,".png"))
					png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				}
				panelTraces(mcmc=pmcmc, mpd=pmpd, xlab="Samples", ylab="Parameter value", cex.axis=1.2, cex.lab=1.5, same.limits=F, lang=l, mar=c(0,3,0,1), oma=c(3.5,2,0.5,0), mgp=c(1.75,0.5,0))
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
#browser();return()

	fout.e = "splitChain"
	nchains   = 8
	#col.trace = rep(rev(c("red", "blue", "black", "green", "orange", "purple", "cyan", "gold")), nchains)[1:nchains]
	col.trace = rep(rev(c("red", "blue", "black", "green", "orange", "purple", "brown", "pink")), nchains)[1:nchains]
	for (l in lang) {
		changeLangOpts(L=l)
		#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		mcmcP = mcmcObj$P
		names(mcmcP) = linguaFranca(names(mcmcP), l)
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") {
				clearFiles(paste0(fout,".png"))
				png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			}
			panelChains(mcmc=mcmcP, nchains=nchains, axes=TRUE, pdisc=0, between=list(x=0, y=0), col.trace=col.trace, xlab="Parameter Value", ylab="Cumulative Frequency", cex.axis=1.2, cex.lab=1.4, yaxt="n", lang=l, mar=c(1.5,2,0.5,1), oma=c(2,2,0,0), mgp=c(1.75,0.5,0), lwd.trace=1)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
#browser();return()

	fout.e = "paramACFs"
	for (l in lang) {
		changeLangOpts(L=l)
		#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		mcmcP = mcmcObj$P
		names(mcmcP) = linguaFranca(names(mcmcP), l)
		for (p in ptypes) {
			if (p=="eps") postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") {
				clearFiles(paste0(fout,".png"))
				png(paste0(fout,".png"), width=PIN[1], height=PIN[2], units="in", res=pngres)
			}
			plotACFs(mcmc=mcmcP, lag.max=60, lang=l)
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()

	## Plot Rhat ratios
	if (rhat) {
		for (i in 0:1) {
			ii = as.logical(i)
			png = "png" %in% ptypes
			theme="Hawaii"; barcols = colorspace::sequential_hcl(8,theme, rev=F)
			#calcRhat(dir=getwd(), recdevs=ii, png=png, outnam=ifelse(ii,"rhat.recdev","rhat"), lang=lang, pngres=pngres, PIN=c(10,6))
			calcRhat(dir=sub("/sso","",mcmc.dir), recdevs=ii, png=png, lang=lang, pngres=pngres, PIN=c(10,6))  ## let 'calcRhat' create output names
		}
	}

	## Plot B vs h (steepness)
	if (steep) {
		png = "png" %in% ptypes
		plotBh(d.mcmc, B.mcmc, Bmsy.mcmc, ytype="rel", currYr=currYr, png=png)
#browser();return()
		plotBh(d.mcmc, B.mcmc, Bmsy.mcmc, ytype="abs", currYr=currYr, png=png)
	}
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.dmcmc


## plotSS.pairs-------------------------2024-04-10
##  Pairs|density plot comparison among parameters
## ---------------------------------------------RH
plotSS.pairs <- function(P.mpd, P.mcmc, type="image", ptypes, 
   lang=c("e","f"), pngres=400, PIN=c(10,10))
{
	panel.cor <- function(x, y, digits=2, prefix="", ...)
	{
		usr <- par("usr"); on.exit(par(usr=usr)) ## (RH 240410)
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

	panel.text <- function(x, y, labels, cex, font, ...)
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
#browser();return()
				if (p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.pairs


## plotSS.pmcmc-------------------------2025-12-03
##  Plot boxplots of quantiles for various MCMC parameters.
##  Modified from PBSawatea's 'plotBmcmcPOP' function.
##  Input comes from the output file 'derived.parameters.sso' (use r4ss::SSgetMCMC)
## -----------------------------------------AME|RH
plotSS.pmcmc <- function(obj, pqs=NULL, xyType="quantBox",
   lineType=c(3,2,1,2,3), refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, tcl.val=-0.2, yrs, y0=TRUE,
   pyrs=NULL, refpts=list(LRP=NULL,USR=NULL,TRP=NULL,RRR=NULL), catpol=NULL,
   yaxis.by, yLab="Recruitment", outnam, lang=c("e","f"),
   ptypes="win", pngres=400, PIN=c(8,6), gmu=TRUE, ...)
{
#browser();return()
	# See plt.quantBio if want other xyTypes, as took out here:
	plt.qB <- function(qobj, xyType="lines", new=TRUE, xLim, yLim, yrs, pyrs, refpts, pvec, ...)
	{
		if ( new ) {
			unpackList(refpts, scope="L")
			col.lrp = .colBlind["redpurple"]
			col.usr = .colBlind["bluegreen"]
			col.trp = .colBlind["blue"]
			col.rrr = .colBlind["bluegreen"]
			expandGraph(mar=c(3,3.5,1.2,1.2), mgp = c(1.6,0.5,0))
			plot(xLim, yLim, type="n", xlab=linguaFranca("Year",l), ylab=linguaFranca(yLab,l), ...)
			if (exists("currYr"))
				abline(v=currYr, col="purple", lty=2)
			if (!y0) abline(h=0,col="gainsboro")
			if (!is.null(LRP)) abline(h=LRP, col=col.lrp, lty=4, lwd=1.25)  #"red", lty=5)  ## (RH 240912)
			if (!is.null(USR)) abline(h=USR, col=col.usr, lty=5, lwd=1.25)  #"green4", lty=5)  ## (RH 240912)
			if (!is.null(TRP)) abline(h=TRP, col=col.trp, lty=2, lwd=1.25)  ## (RH 250917)
			if (!is.null(RRR)) abline(h=RRR, col=col.rrr, lty=5, lwd=1.25)  ## (RH 250917)
		}
		#yrs <- as.numeric(dimnames(qobj)[[2]])
		#yrs <- as.numeric(gsub("[^[:digit:]]","",dimnames(result1)[[2]]))

		# Quantile boxplots - assumes five quantiles.
		if ( xyType=="quantBox" ) {
			#colvec      <- c(rep("green3",nrow(recdevEarly)), rep("black",nrow(recdevMain)), rep("blue",nrow(recdevLate)), rep("red", nrow(recdevFore)))
			#bgvec       <- c(rep("green",nrow(recdevEarly)), rep("gainsboro",nrow(recdevMain)), rep("cyan",nrow(recdevLate)), rep("pink", nrow(recdevFore)))
			## Above not robust to all situations (RH 250408)
			colvec      <- c(rep("green3",length(earlyYrs)), rep("black",length(mainYrs)), rep("blue",length(lateYrs)), rep("red", length(foreYrs)) )
			names(colvec) = c(earlyYrs, mainYrs, lateYrs, foreYrs)
			bgvec      <- c(rep("green",length(earlyYrs)), rep("gainsboro",length(mainYrs)), rep("cyan",length(lateYrs)), rep("pink", length(foreYrs)) )
			names(bgvec) = c(earlyYrs, mainYrs, lateYrs, foreYrs)

			if ("afR.mcmc" %in% fnam) { ## Rdist only occurs during recdevMain
				areacol = switch(area,'5ABC'="blue",'3CD'="green4",'5DE'="red","black")
				areabg  = switch(area,'5ABC'="cyan",'3CD'="green",'5DE'="pink","gainsboro")
				colvec  = rep(areacol,length(yrs))
				bgvec   = rep(areabg,length(yrs))
				names(colvec) = names(bgvec) = yrs  ## (RH 250425)
				medval  = mean(qobj[3,])
				abline(h=medval, col=areacol, lty=3)
				addLabel(0.02,0.02,paste0("mean of medians = ",show0(round(medval,3),3)), col=areacol, adj=c(0,0),cex=1)
			}
#browser();return()
			#if(!is.null(USR)) abline(h=c(USR), col=c("slategray"), lwd=1, lty=5)  ## already covered by common frame above
			delta <- 0.25 ## width of half-box
			# Draw the outer whiskers.
			segments(yrs, qobj[1,], yrs,qobj[5,], lty=1, col="black", lwd=ifelse(narea==1,1,0.5))
			# Overlay the box.
			for ( i in 1:length(yrs) ) {
				ii = yrs[i]; iii = as.character(ii)
				#rect( yrs[i]-delta,qobj[2,i], yrs[i]+delta, qobj[4,i],... ) ## AME
				#polygon(x=c(rep(yrs[i]-delta,2),rep(yrs[i]+delta,2)), y=qobj[c(2,4,4,2),i], border="black", col="gainsboro", ...) ## RH (190620)
				polygon(x=c(rep(yrs[i]-delta,2),rep(yrs[i]+delta,2)), y=qobj[c(2,4,4,2),iii], border=colvec[iii], col=bgvec[iii], lwd=ifelse(narea==1,1,0.5)) ## RH (230816)
#browser();return()
			}
			# Add the median.
			segments(yrs-delta, qobj[3,], yrs+delta, qobj[3,], lty=1, col=colvec, lwd=ifelse(narea==1,1,0.5))
#browser();return()

			## Revised handling of projections to fix Rdev plot (RH 230809)
			addyrs = setdiff(allyrs,yrs)
			prjyrs = as.character(pyrs)
#browser();return()
			if (!is.null(pyrs) && all(prjyrs %in% colnames(qobj)) ) {
				## Following now handled above (RH 230816)
				#pyrs = sort(unique(c(addyrs,pyrs)))
				#pp   = as.character(pyrs)
				#segments(pyrs, qobj[1,pp], pyrs, qobj[5,pp], lty=1,col="red" )
				#for ( i in (length(yrs)-length(pyrs)+1):length(yrs) )
				#	polygon(x=c(rep(yrs[i]-delta,2),rep(yrs[i]+delta,2)), y=qobj[c(2,4,4,2),i], border="red", col="pink", ...) ## RH (190620)
				#segments( pyrs-delta, qobj[3,pp], pyrs+delta, qobj[3,pp], lty=1, col="red" )
			} else {
				if (fnam %in% c("Rdev.mcmc","Rtdev.mcmc")) {  ## messy and may be prone to error
					#parameters <- replist$parameters
					#recdevEarly <- parameters[substring(parameters$Label, 1, 13) %in% c("Early_RecrDev"), ]
					#recdev <- parameters[substring(parameters$Label, 1, 12) %in% c("Main_RecrDev"), ]
					#recdevFore <- parameters[substring(parameters$Label, 1, 8) == "ForeRecr", ]
					#recdevLate <- parameters[substring(parameters$Label, 1, 12) == "Late_RecrDev", ]
					#colvec <- c(rep("blue",nrow(recdevLate)), rep("red", nrow(recdevFore)))
					nullyrs = setdiff(names(colvec),yrs)  ## (RH 250409)
					if (length(nullyrs)>0) {
						points(as.numeric(nullyrs), rep(0,length(nullyrs)), pch=20, col=colvec[nullyrs], cex=1)
					}
				}
			}
			box()
		}  ## end quantBox
		# Uncertainty envelope - assumes five quantiles.
		if ( xyType=="envelope" ) {
			x  = setdiff(yrs,pyrs)
			#x  = union(earlyYrs, mainYrs)
			x.late = NULL
			#if (nrow(recdevLate)>0) {  ## (RH 230821)
			if (length(lateYrs)>0) {  ## (RH 250409)
				#xx.late = c(sub("Late_RecrDev_","",rownames(recdevLate)), as.character(currYr))
				#x.late  = as.numeric(xx.late)
				x.late = c(lateYrs, currYr) ## include currYr in the late period (RH 251028) looks better
				#x.late = lateYrs  ## currYr goes into projections, extend main to start of late period
#browser();return()
				x = setdiff(x,.su(c(x.late,currYr)))
			}
			xx = as.character(x)
			## Main recruitment period
			xpol = x; xpol[1]=xpol[1]-0.5; xpol[length(xpol)]=xpol[length(xpol)]+0.5  ## buffering for polygon (RH 250410)
			col.main = "black"
			polygon(c(xpol,rev(xpol)), c(qobj[1,xx],rev(qobj[5,xx])), col=lucent(col.main,0.05), border=FALSE)
			lines(x, qobj[3,xx], lty=1, lwd=ifelse(narea==1,3,2), col=col.main)
			lines(c(x,NA,x), c(qobj[1,xx],NA,qobj[5,xx]), lty=2, lwd=ifelse(narea==1,2,1), col=col.main)
			lines(c(x,NA,x), c(qobj[2,xx],NA,qobj[4,xx]), lty=3, lwd=ifelse(narea==1,1,1), col=col.main)
			if (!is.null(x.late)) {  ## (RH 230821)
				xpol = x.late; xpol[1]=xpol[1]-0.5; xpol[length(xpol)]=xpol[length(xpol)]+0.5  ## buffering for polygon (RH 250410)
				ypol = as.character(x.late)
				xlin = xx.late = c(x.late[1] - 1, x.late)  ## extend x.late from last main year for discrete data
				ylin = as.character(xlin)
				col.late = "blue"
				polygon(c(xpol,rev(xpol)), c(qobj[1,ypol],rev(qobj[5,ypol])), col=lucent(col.late,0.05), border=FALSE)
				lines(xlin, qobj[3,ylin], lty=1, lwd=ifelse(narea==1,3,2), col=col.late)
				lines(c(xlin,NA,xlin), c(qobj[1,ylin],NA,qobj[5,ylin]), lty=2, lwd=ifelse(narea==1,2,1), col=col.late)
				lines(c(xlin,NA,xlin), c(qobj[2,ylin],NA,qobj[4,ylin]), lty=3, lwd=ifelse(narea==1,1,1), col=col.late)
			}
			#x.fore = .su(c(currYr, pyrs))
			x.fore = setdiff(pyrs, currYr)  ## exclude current year from projections (RH 251028)
			if (all(!is.na(pvec)) && all(x.fore %in% yrs)){
				kcol = if (length(pvec)==1) "red" else c("green3","orange2","red")
				xpol = x.fore
				if (length(x.fore)==1) {
					xpol = c(x.fore, x.fore)
				}
				ypol = as.character(xpol)
				#xpol[1]=xpol[1]-0.5; xpol[length(xpol)]=xpol[length(xpol)]+0.5  ## buffering for polygon (RH 250410)
				xpol[1]=xpol[1]-1; xpol[length(xpol)]=xpol[length(xpol)]+0  ## buffering for polygon (RH 251028)
				xlin = xx.fore = c(x.fore[1] - 1, x.fore)  ## extend x.fore from last main year for discrete data
				ylin = as.character(xlin)
				for (k in 1:length(pvec)) {
					pp  = pvec[[k]]
					yy  = cbind(qobj[,ylin[1]], pp); colnames(yy)[1] = ylin[1]
					polygon(c(xpol,rev(xpol)), c(yy[1,ypol],rev(yy[5,ypol])), col=lucent(kcol[k],0.1), border=FALSE)
					lines(xlin, yy[3,ylin], lty=1, lwd=ifelse(narea==1,3,2), col=kcol[k])
					lines(c(xlin,NA,xlin), c(yy[1,ylin],NA,yy[5,ylin]), lty=2, lwd=ifelse(narea==1,2,1), col=kcol[k])
					lines(c(xlin,NA,xlin), c(yy[2,ylin],NA,yy[4,ylin]), lty=3, lwd=ifelse(narea==1,1,1), col=kcol[k])
#abline(v=c(2015,2024)) ## just checking start years for late and proj
				}
#browser();return()
			}
		}  ## end envelope
		if (gmu && !all(sapply(refpts,is.null))) {
			xrfp = par()$usr[1] + 0.05*diff(par()$usr[1:2])
			yrfp = c(LRP,USR,TRP,RRR) + ifelse("png"%in%ptypes,0.02,0.02) * diff(par()$usr[3:4])
			arfp = !sapply(list(LRP,USR,TRP,RRR),is.null)  ## (RH 250917) active ref point
			text(x=xrfp, y=yrfp, labels=linguaFranca(c("LRP","USR","TRP","RRR")[arfp],l), col=c(col.lrp,col.usr,col.trp,col.rrr)[arfp], font=2, cex=ifelse(prod(par()$mfrow)==1,1,0.8))
		}
	} ## end subfunction plt.qB

	## Deal with refpts list object
	if (!is.list(refpts))
		stop ("Supply a list of reference points for 'refpts'")
	if (!all(sapply(refpts,is.null))) {
		refpts.in = refpts
		refpts = list(LRP=NULL, USR=NULL, TRP=NULL, RRR=NULL)
		refpts[names(refpts.in)] = refpts.in
	}
#browser();return()

	## Save all years to counteract r4ss' propensity to choose available data
	allyrs = c(yrs, pyrs)
	fnam   = as.character(substitute(obj))

	## Need to get the various classifications for recdevs (and other parameters) for PJS (RH 230816)
#browser();return()
	if (!exists("replist")) {
		if (exists("central.mpd.dir") && dir.exists(central.mpd.dir)) {
		#if (exists("central.mpd.dir", envir=sys.frame(-1))) {
		#	c.mpd.dir = get("central.mpd.dir", envir=sys.frame(-1))
		#	if (dir.exists(c.mpd.dir)) {
				assign("replist", SS_output(central.mpd.dir, verbose=FALSE, printstats=FALSE), envir=.GlobalEnv)
		#	}
		} else {
			message ("sumtingwong with replist")
			browser();return()
		}
	}
	parameters  <- replist$parameters
	recdevEarly <- parameters[substring(parameters$Label, 1, 13) %in% c("Early_RecrDev"), ]
	recdevMain  <- parameters[substring(parameters$Label, 1, 12) %in% c("Main_RecrDev"), ]
	recdevLate  <- parameters[substring(parameters$Label, 1, 12) == "Late_RecrDev", ]
	recdevFore  <- parameters[substring(parameters$Label, 1, 8) == "ForeRecr", ]

	## Switch to specific years rather than rely on recdev information only (RH 250409)
	mainYrs  <- as.numeric(sub("Main_RecrDev_","",rownames(recdevMain)))
	earlyYrs <- setdiff(startYear:(mainYrs[1]-1), mainYrs) ## this set is most sensitivr to inputs
	lateYrs  <- as.numeric(sub("Late_RecrDev_","",rownames(recdevLate)))
	foreYrs  <- as.numeric(sub("ForeRecr_","",rownames(recdevFore)))
#browser();return()

	if (inherits(obj,"list") &&  is.null(dim(obj)))
		narea = length(obj)
	else
		narea = 1
	createFdir(lang)
	if (missing(outnam))
		outnam = gsub("[[:space:]]+",".",tolower(yLab))
	fout.e = outnam
	for (l in lang) { 
		changeLangOpts(L=l)
		#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		for (p in ptypes) {
			if (p=="eps")      postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") {
				clearFiles(paste0(fout,".png"))
				png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			}
#browser();return()
			## narea is an internal object determined above
			rc = if(narea==3) c(3,1) else .findSquare(narea)
			expandGraph(mfrow=rc, mar=c(3.0,3.5,0.75,0.5), oma=c(0,0,0,0), mgp=c(1.75,0.5,0))
			for (i in 1:narea) {
				if (narea==10) {
					area = if (length(area.name)>1) switch(strSpp,'396'="CST","BC coast") else area.name  ## get area name from 'set.controls.r' (RH 251001)
					aobj = obj
					acatpol = catpol  ## could be NULL
				} else {
					area = names(obj)[i]
					aobj = obj[[i]]
					acatpol = catpol[[i]]
				}
				# Plot quantiles of biomass using the posterior densities.
				yrs1 = yrs2 = result1 = result2 = NULL
				if(is.null(pqs)) pqs = c(0.05, 0.25, 0.50, 0.75, 0.95)
		
				# Calculate the quantiles of the reconstructed biomass.
#browser();return()
				result1 <- apply( aobj, 2, quantile, probs=pqs, na.rm=T )
				colnames(result1) = gsub("SSB_|u_|F_","",colnames(result1))
				yrs1 <- as.numeric(gsub("[^[:digit:]]","",dimnames(result1)[[2]]))
				if (!missing(yrs)) {
					zyrs = is.element(yrs1,c(yrs,pyrs))
					zyrs = zyrs & apply(result1,2,function(x){!all(is.na(x))})  ## to emulate original input
					yrs1 = yrs1[zyrs]
					result1 = result1[,zyrs]
				}
				cpyrs = .su(c(currYr, pyrs))  ## pyrs must include at least the current year projection
				if (!is.null(cpyrs) && any(cpyrs%in%colnames(result1))) {
					if (!is.null(acatpol) && !all(dimnames(acatpol)$proj=="AC.00")) {
						cpolnam = dimnames(acatpol)$proj
						#if (is.null(cpolnam)) cpolnam = "AC.00"
						if (narea==1)
							cpolcat = avgCP[1,1,cpolnam]
						else
							cpolcat = xavgCP[1,sub("^BC$|^CST$","total",area),cpolnam]
						projmat = apply(acatpol, 2:3, quantile, probs=pqs, na.rm=T )
						projvec = lapply(1:dim(projmat)[3], function(k) {projmat[,,k]})
						names(projvec) = cpolcat
					} else {
						projvec = list('AC'=result1[,intersect(colnames(result1),as.character(cpyrs))])
					}
				} else {
					projvec = NA
				}
				if (!all(is.na(projvec)) && length(cpyrs)==1) ## adjust matrix 
					projvec = lapply(projvec, function(x,y){xx=t(t(x)); colnames(xx)=y; xx}, y=cpyrs)

				if ( is.null(yLim) || narea>1 ){
					yLim <- c(ifelse(y0, 0, min(unlist(result1),unlist(projvec),na.rm=T)), max(unlist(result1),unlist(projvec),na.rm=T))
#browser();return()
					if (length(fnam)==1 && fnam=="afR.mcmc")
						yLim = c(0,1)  ## (RH 250731)
				}
				if ( is.null(xLim) || narea>1 )
					xLim = range(allyrs) #range(yrs1)
				colnames(result1) = yrs1

#browser();return()
				plt.qB(qobj=result1, xLim=xLim, yLim=yLim, xyType=xyType, yrs=yrs1, pyrs=pyrs, refpts=refpts, pvec=projvec, ...)
				axis(1, at=intersect(seq(1900,3000,5), xLim[1]:xLim[2]), tcl=tcl.val, labels=FALSE)
				if (missing(yaxis.by)) {
					yint = diff(seq(par()$yaxp[1], par()$yaxp[2], len=par()$yaxp[3]+1))[1]
					yint.char = as.character(yint) ## becasue comparing numerics is flaky
					#ytck = seq(par()$yaxp[1]-yint, par()$yaxp[2]+yint, yint/ifelse(par()$yaxp[3]>6,2,4))
					if ( yint.char %in% as.character(5*10^(-10:10)) ) by = yint/5
					else if ( yint.char %in% as.character(4*10^(-10:10)) ) by = yint/4
					else if ( yint.char %in% as.character(2*10^(-10:10)) ) by = yint/2
					else if ( yint.char %in% as.character(1*10^(-10:10)) ) by = yint/5
					else { message("stupid tick checker has failed"); browser(); return() }
					ytck = seq(par()$yaxp[1]-yint, par()$yaxp[2]+yint, by=by)
					axis(2, at=ytck, tcl=tcl.val, labels=FALSE)
				} else 
					axis(2, at=seq(par()$yaxp[1], par()$yaxp[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
				#axis(2, at=seq(0, yLim[2], by=yaxis.by), tcl=tcl.val, labels=FALSE)
				if (length(projvec)==3)
					addLegend(0.5,0.975,legend=paste0(names(projvec)," t"), lty=1, lwd=ifelse(narea==1,3,2), seg.len=3, col=c("green3","orange2","red"), bty="n", xjust=0, yjust=1,title=linguaFranca("Projected catch",l), cex=ifelse(narea==1,1,0.8))
#browser();return()
				if (narea>1 || strSpp=="396") {## retrofit for POP 2023
					area = sub("CST|5ABC-5DE-3CD", "BC", area)  ## adjust for SGR 2025
					#addLabel(0.975, 0.975, linguaFranca(area,l), adj=c(1,1), font=2, cex=ifelse(narea==1,1.2,1.2), col=switch(area,'5ABC'="blue",'3CD'="green4",'5DE'="red","black") )
					addLabel(0.05, 0.975, linguaFranca(area,l), adj=c(0,1), font=2, cex=ifelse(narea==1,1.2,1.2), col=switch(area,'5ABC'="blue",'3CD'="green4",'5DE'="red","black") )
				}
			}  ## end i loop (area)
			if (p %in% c("eps","png")) dev.off()
		}
	}; eop()  ## end language loop
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.pmcmc


## plotSS.rmcmc-------------------------2025-09-15
##  Plot routine output for MCMCs.
## ---------------------------------------------RH
plotSS.rmcmc <- function(mcmcObj, mpdObj, ptypes, lang, pngres=400, PIN=c(9,9))
{
	so("linguaFranca.r")
	unpackList(mpdObj); unpackList(mcmcObj)

	## Pairs|density plot comparison among parameters
	plotSS.pairs(P.mpd, P.mcmc, ptypes=ptypes, lang=lang, pngres=pngres, PIN=c(10,10), type="image")

	## Parameter posteriors and priors
	for (l in lang) {
		for (p in ptypes) {
			if (p=="win") { plot=TRUE; print=FALSE }
			else          { plot=FALSE; print=TRUE }
			plotSS.pars(replist, nrows=P.rc[1], ncols=P.rc[2], plot=plot, print=print, fitrange=T, fitnudge=0.1, showpost=T, strings=names(P.mpd), exact=T, plotdir=getwd(), lang=l, outnam="pdfParameters")
		}
	}

	## Snail plots of Bt/Bmsy and ut/umsy
	colnames(BoverBmsy) = sub("^SSB_","",colnames(BoverBmsy))
	colnames(UoverUmsy) = sub("^u_","",colnames(UoverUmsy))
	for (p in ptypes) {  ## needs 'lang' not 'l'
		plotSnail(BoverBmsy, UoverUmsy, yrs=modYrs, p=tcall(quants3)[c(1,3)], xLim=NULL, yLim=NULL, ngear=ngear, assYrs=assYrs, outs=F, Cnames=fleets.lab[1:ngear], ptypes=p, outnam="snail", lang=lang, labYrs=c(1950,seq(1960,1975,5), assYrs))
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.rmcmc


## plotTraj-----------------------------2020-06-04
## Show all median trajectories (base+sens) in one figure.
## ---------------------------------------------RH
plotTraj <- function(sdat, index, traj="B", bdat=NULL, y0=FALSE, lab.stock, 
   startYear=1935, currYear=2025, sen.lab,
   col = c("black","green4","blue","red","purple","orange"), 
   lty = c(2:6), logR=FALSE, 
   png=FALSE, pngres=400, PIN=c(8,8), lang=c("e","f"), ...)
{
	opar = par(no.readonly=TRUE); on.exit(par(opar))
	Ntraj     = length(traj)
	trajYears = startYear:currYear
	Nyears    = length(trajYears)
	if (is.null(tcall(run.sens)))
		run.sens  = names(sdat$currentMCMC.sens)[index]
	else
		run.sens = tcall(run.sens)[index]
#browser();return()
	do.call("assign", args=list(x="run.sens", value=run.sens, envir=.GlobalEnv))
	Nruns     = length(run.sens)
	trajMat   = array(NA, dim=c(Nyears,Nruns,Ntraj), dimnames=list(year=trajYears, run=run.sens, traj=traj))
	#trajMat   = array(NA, dim=c(Nyears,Nruns,Ntraj+ifelse("BtB0"%in%traj,1,0)), dimnames=list(year=trajYears, run=run.sens, traj=if("BtB0"%in%traj) c(traj,"BmsyB0") else traj))
	Bmsy.q5   = list()  ## not always needed
	if (!is.null(bdat))
		col = setdiff(col,"black") ## save black for base case
	## ylabs needs to be a list if mixing character elements with expressions
	ylabs     = list("Spawning Biomass","Vulnerable Biomass","Recruitment","Exploitation Rate",expression(italic(B[t])/italic(B)[0]),"Unknown")
	names(ylabs) = c("B","VB","R","U","BtB0","NA")
	for (i in 1:Nruns) {
		ii = run.sens[i]
		for (j in 1:Ntraj) {
			jj = traj[j]
			if (jj=="BtB0") {
				jdat = sdat[["currentMCMC.sens"]][[ii]][["B"]]
				sdat[["currentMCMC.sens"]][[ii]][[jj]] = sweep(jdat,1,jdat[,1],"/")
				## Add Bmsy stuff (thanks PJS for the complication); and then he changed his mind...
				sdat[["currentMCMC.sens"]][[ii]][["BmsyB0"]] = sdat[["currentMSY.sens"]][[ii]][["B"]]/jdat[,1]
				Bmsy.q5[[ii]] = quantile(0.8*sdat[["currentMCMC.sens"]][[ii]][["BmsyB0"]],quants5)
			}
			jtmp = sdat[["currentMCMC.sens"]][[ii]][[jj]]
			if (jj %in% c("U","VB")) jtmp = splitGear(jtmp)[[1]]
			jval = apply(jtmp,2,median)
			if (jj=="R" && logR) jval = log10(jval)
			trajMat[substring(names(jval),1,4),ii,jj] = jval
		}
	}
	if (missing(lab.stock)) lab.stock = istock
	
	createFdir(lang)
	fout = fout.e = paste0(lab.stock,".sens.traj.",paste0(traj,collapse="+"))
	for (l in lang) {  ## could switch to other languages if available in 'linguaFranca'.
		changeLangOpts(L=l)
		fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )

		## tcol and legtxt need to be reset inside language loop because they are altered when bdat is supplied
		tcol    = rep(col,Nruns)[1:Nruns]
		legtxt  = sen.lab[index]  ## defined in global data object 'stock'
#browser();return()
		tlty    = rep(lty,Nruns)[1:Nruns]

		if (png) png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
		par(mfrow=.findSquare(Ntraj), mar=c(3,3.5,0.5,0), oma=c(0,0,0,1), mgp=c(2,0.5,0))
		x = trajYears; xlim = range(x)
		for (j in 1:Ntraj) {
			jj   = traj[j]
			jmat = trajMat[,,jj]
			ylim = range(jmat,na.rm=TRUE)
			#if (jj %in% c("BtB0")) ylim[1]=0
			if (y0) ylim[1] = 0
			if (jj=="R" && logR) ylim[1]=1
			plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="",log=ifelse(jj=="RRR","y",""))
			if (jj %in% c("BtB0")){
				#abline(h=seq(0.2,1,0.2),col="gainsboro")
				abline(h=c(0.2,0.4,1), lty=5, col="grey20") #col=c("salmon","darkorchid","navy"))
			}
			x  = as.numeric(rownames(jmat))
			for (k in 1:ncol(jmat)){
				kk = run.sens[k]
				y  = jmat[,k]
				lines(x,y, col=tcol[k], lwd=ifelse(i==1,2,2), lty=tlty[k]) ## no longer include base case in first position
			}
			if (!is.null(bdat)){  ## Assume for now that this is the Base Case (could be other runs)
				bline = bdat[[jj]]
#browser();return()
				lines(x=as.numeric(names(bline)), y=if (jj=="R" && logR) log10(bline) else bline, col="black", lty=1, lwd=3)
			}
			mtext(linguaFranca("Year",l), side=1, line=1.75, cex=ifelse(Ntraj==1,1.5,1))
			#ylab = ifelse(jj=="B","Spawning Biomass",ifelse(jj=="VB","Vulnerable Biomass",ifelse(jj=="R","Recruitment",ifelse(jj=="U","Exploitation Rate","Unknown"))))
			mtext(linguaFranca(ylabs[[jj]],l), side=2, line=1.8, cex=ifelse(Ntraj==1,1.5,1.2))
			if (j==1){
				legtxt = gsub("_"," ",legtxt)
				if (!is.null(bdat)){
					legtxt = c(bdat$label, legtxt[1:Nruns])
					tcol   = c("black",tcol)
					tlty    = c("solid", tlty)
				}
#browser():return()
				addLegend(ifelse(jj%in%c("R","U"),0.025,0.05), ifelse(jj%in%c("BtB0"),0.025,0.975), col=tcol, seg.len=5, legend=linguaFranca(legtxt,l), bty="o", box.col="grey", bg=ifelse(jj%in%c("BtB0"),"white","transparent"), xjust=ifelse(jj%in%c("R","U"),0,0), yjust=ifelse(jj%in%c("BtB0"),0,1), lwd=2, lty=tlty, ...)
			}
			if (Ntraj==1) {
				axis(1,at=intersect(seq(1900,2500,5),x),labels=FALSE,tcl=-0.2)
				axis(1,at=intersect(seq(1900,2500,10),x),labels=FALSE)
			}
			box()
		}
		if (png) dev.off()
	}; eop()
	do.call("assign", args=list(x="trajSens", value=trajMat, envir=.GlobalEnv))
	invisible()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotTraj


##====================================== Start suite of quantile boxes of age-fit residuals

## QUANTILE BOXES of AGE-FIT RESIDUALS by AGE, YEAR, COHORT

## plt.ageResids------------------------2024-04-08 (eclipse)
## AME changing for POP, just plotting age class resids here,
## not qq-plot. Moving that to new function (so can do 4 on page).
## See popScape.r for original.
## -----------------------------------------AME|RH
plt.ageResids <- function( obj, ages=NULL, main=NULL, lang="e", yrs, resfld="Pearson", ...)
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
		xx=split(c(obj[,resfld], rep(NA, length(nodataAges))), c(obj$Bin, nodataAges))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantBox(xx, xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xage = setdiff(allAges, nodataAges)
		if (length(xage)>20)
			xage = sort(unique(ceiling(xage/5)*5))
		axis(1, at=match(xage,as.numeric(xpos$names)), labels=xage, ..., mgp=c(2,0.5,0))
	} else {
		#xpos <- boxplot( split( obj[,resfld], obj$Bin ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		xpos = quantBox(split(obj[,resfld],obj$Bin), xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars))
		xage = sort(unique(ceiling(obj$Bin/5)*5))
		axis(1, at=match(xage,as.numeric(xpos$names)), labels=xage, ..., mgp=c(2,0.5,0))
	}
	abline( h=0, lty=2, col="red" )
	axis(2, ...)
	mtext( side=1, line=2, ..., text=linguaFranca("Age class",lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.ageResids

## plt.cohortResids---------------------2024-04-08 (eclipse)
## Plot age residuals by cohort.
## Removed multiple calls to quantBox when colouring by recdevs (RH 230502)
## -----------------------------------------AME|RH
plt.cohortResids <- function( obj, ages=NULL, main=NULL, lang="e", use.rdevs=T, resfld="Pearson", ...)
{
	## Input is the CAc object from a Awatea res file. Ages to 59 as
	##  plus-age class will mess up year-of-birth calculation. Not automated.
	## par( oma=c(2,1,1,1), mar=c(2,2,2,1), mfrow=c(2,1) )
	## Subset to required ages - still do as don't want age 1 or 60 for cohorts
	if (!is.null(ages))
		obj <- obj[ (obj$Bin >= ages[1]) & (obj$Bin <=ages[2]), ]
	## obj[,resfld] has residuals for each age, year and both sexes. Need
	##  to assign a year of birth for each age as an extra column, then
	##  presumably just do the boxplot split using that.
	obj$birthyr = obj$Yr - obj$Bin

	upar = tcall(boxpars)
	if (use.rdevs) {
		recdevs = getSS.rdevs(replist)
		#if (!is.null(ttcall(recdevs))){ ## created by function 'getSS.rdevs'
		if (!is.null(recdevs)){ ## created by function 'getSS.rdevs'
			#ttget(recdevs)
			ayrs = recdevs$Yr
			zyrs = recdevs$Yr[recdevs$Value<=0]
			medcol = rep("green4",length(ayrs)); names(medcol) = ayrs
			medcol[as.character(zyrs)] = "red"
			boxfill = rep(lucent("green",0.25),length(ayrs)); names(boxfill) = ayrs
			boxfill[as.character(zyrs)] = lucent("red",0.25)
			upar$medcol = medcol
			upar$boxfill = boxfill
		}
	}
	if( max(diff(sort(unique(obj$birthyr)))) > 1) {
		allYears = min(obj$birthyr):max(obj$birthyr)
		nodataYears = allYears[ !(allYears %in% obj$birthyr)]
		xx = split(c(obj[,resfld], rep(NA, length(nodataYears))), c(obj$birthyr, nodataYears))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		#xpos = quantBox(xx, xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars), add=FALSE)
		xpos = quantBox(xx, xaxt="n", yaxt="n", xlab="", ylab="", pars=upar, add=FALSE)
		xyrs = setdiff(allYears, nodataYears)
		if (length(xyrs)>20)
			xyrs = sort(unique(ceiling(xyrs/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, ..., mgp=c(2,0.5,0))
	} else {
		#xpos=boxplot( split( obj[,resfld], obj$birthyr ), whisklty=1, xlab="", ylab="", outline=FALSE )     #AME outline=FALSE removes outliers
		#xpos = quantBox(split(obj[,resfld],obj$birthyr), xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars), add=FALSE)
		xpos = quantBox(split(obj[,resfld],obj$birthyr), xaxt="n", yaxt="n", xlab="", ylab="", pars=upar, add=FALSE)
		xyrs = sort(unique(ceiling(obj$birthyr/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, ..., mgp=c(2,0.5,0))
	}
#browser();return()
	#if (use.rdevs) {
	#	plotSS.rdevs(replist, subplot=1, plot=F, print=F)
	#	if (!is.null(ttcall(recdevs))){ ## created by function 'plotSS.rdevs'
	#		ttget(recdevs)
	#		zyrs = recdevs$Yr[recdevs$Value>0]
	#		zpar = upar = tcall(boxpars)
	#		zpar$boxfill=lucent("green",0.25); zpar$medcol="green4"
	#		upar$boxfill=lucent("red",0.25); upar$medcol="red"
	#		#zpar$boxfill=lucent("red",0.5)
	#		zbox = ubox = split(obj[,resfld],obj$birthyr)
	#		zbox[!is.element(names(zbox),zyrs)] = NA
	#		ubox[is.element(names(zbox),zyrs)] = NA
	#		quantBox(zbox, xaxt="n", yaxt="n", xlab="", ylab="", pars=zpar, add=TRUE)  ## this only seems to work if there are no more frames available
	#		quantBox(ubox, xaxt="n", yaxt="n", xlab="", ylab="", pars=upar, add=TRUE)
	#	}
	#}
	abline( h=0, lty=2, col="blue" )
	axis(2, ...)
	mtext( side=1, line=2, ..., text=linguaFranca("Year of birth",lang) )
	if ( !is.null(main) )
		mtext( side=3, line=-0.5, ..., outer=TRUE, text=linguaFranca(main,lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.cohortResids

## plt.yearResids-----------------------2024-04-08 (eclipse)
## AME adding to plot age residuals by year. Is called for comm and survs.
##  fill.in=TRUE is to add the missing years for boxplot
##  ..POP does not do qq plot. See popScape.r for previous.
## -----------------------------------------AME|RH
plt.yearResids <- function(obj, ages=NULL, main=NULL, fill.in=TRUE, lang="e", resfld="Pearson", ...)
{
	# Subset to required ages - still do as don't want age 1.
	if (!is.null(ages))
		obj <- obj[ (obj$Bin >= ages[1]) & (obj$Bin <=ages[2]), ]
	if(fill.in) {
		allYears = min(obj$Yr):max(obj$Yr)
		nodataYears = allYears[ !(allYears %in% obj$Yr)]
		xx = split(c(obj[,resfld], rep(NA, length(nodataYears))), c(obj$Yr, nodataYears))
		#xpos <- boxplot( xx, whisklty=1, xlab="", ylab="", outline=FALSE, ... )     #AME outline=FALSE removes outliers
		xpos = quantBox(xx, xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars), ...)
		xyrs = setdiff(allYears, nodataYears)
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, ..., mgp=c(2,0.5,0))
#browser();return()
	} else {
		#xpos <- boxplot( split( obj[,resfld], obj$Yr ), whisklty=1, xlab="", ylab="", outline=FALSE, ... ) #AME outline=FALSE removes outliers
		xpos = quantBox( split(obj[,resfld],obj$Yr), xaxt="n", yaxt="n", xlab="", ylab="", pars=tcall(boxpars), ...)
		xyrs = sort(unique(ceiling(obj$Yr/5)*5))
		axis(1, at=match(xyrs,as.numeric(xpos$names)), labels=xyrs, ..., mgp=c(2,0.5,0))
	}
	abline( h=0, lty=2, col="red" )
	axis(2, ...)
	mtext( side=1, line=2, ..., text=linguaFranca("Year",lang) )
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.yearResids

##====================================== End suite of quantile boxes of age-fit residuals


## plt.selectivity----------------------2025-07-21
## Transferred selectivity code from PBSscape.r
## into plt function (RH 190718)
## sobj = second object (RH 201119)
## -----------------------------------------AME|RH
plt.selectivity <- function( obj, sobj=NULL, mainTitle="Rockfish", maxage,
   ptypes="win", pngres=400, PIN=c(7,8), lang=c("e","f"))
{
	## Plot the selectivity.
	panel.selex <- function(x, y, ...){
		dots = list(...)[[1]]
		unpackList(dots)
		nsex = abs(getNpan()%%2-2)
		#unpackList(tcall(panel.dots))
		agebin = if (maxage>25) 10 else 5
		abline(h=seq(0.2,0.8,0.2), v=seq(agebin,max(x)-1,agebin), lty=3, col="seashell4")
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
		mvec = NULL
		if (!is.null(mats)) {
			mvec = rep("m", length(m))
			if (length(m)>50)
				mvec[seq(2,length(m),2)] = convUTF("\\u{2022}") ## "\225" octals no longer supported well
			if (getNpan()%%2==1)
				text(x, m, mvec, col=lucent("black",0.75), cex=0.9)
		}
#browser();return()
		#points(x, m, pch=109, col=lucent("black",0.75), cex=0.8) #, bg=lucent("gold",0.25)
		#lines(x, m, lty=1, col="pink", lwd=2); lines(x, m, lty=2, col="grey", lwd=2)
	} ## end panel.selex

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
#browser();return()
	Xfile.f$Index = index.f[Xfile$Index]

	if (length(index.e)<12)
		rc = c(ceiling(length(index.e)/2),2)
	else 
		rc = c(ceiling(length(index.e)/4),4)
#browser();return()


	## Get the maturity ogive
	#mats = getSS.control(tcall(control.file))$maturity_ogive[1:length(Pvec)]
	mats = SS_readctl(tcall(control.file))$Age_Maturity[1:length(Pvec)]  ## (RH 240820)

	strip.list = list(col=lucent("black",0.5), bg=lucent("moccasin",0.5), height=0.1, cex=1.5)

	#panel.dots = list(scol=c("red","blue"), lwd=4)
	#tput(panel.dots)
	createFdir(lang)
	fout.e = "selectivity"
	for (l in lang) {
		changeLangOpts(L=l)
		#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
		fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
		if (is.null(ptypes)) ptypes="win"
		for (p in ptypes) {
			if (p=="eps")      postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
			else if (p=="png") {
				clearFiles(paste0(fout,".png"))
				png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
			}
			## regular par settings are ignored in lattice::xyplot -- need to fiddle with lattice settings
			#par(mfrow=c(1,1), mar=c(3.2,3.2,0.5,1), oma=c(0,0,0,0), mgp=c(2,0.75,0))
			#par.sel = list(layout.widths = list(left.padding=-0.5, right.padding=-0.5),
			#	layout.heights=list(top.padding=0.5, main.key.padding=-0.75, bottom.padding=-2))
			if (l=="f")
				mochaLatte(Xfile.f,"Age","Sel","Index", panel=panel.selex, m=mats, strip=strip.list, scol=c("red","blue"), slwd=2, mar=c(0,0,0,0), oma=c(4,4,4,1), xlab=linguaFranca("Age (years)","f"), ylab=linguaFranca("Proportion","f"), mainTitle=paste0("s\u{00E9}lectivit\u{00E9} du ", linguaFranca(mainTitle,l)), cex.axis=1.1, cex.main=1.4, sobj=sobj, rc=rc, ylim=c(0,1), byrow=FALSE)
			else{
#browser();return()
				mochaLatte(Xfile,"Age","Sel","Index", panel=panel.selex, m=mats, strip=strip.list, scol=c("red","blue"), slwd=2, mar=c(0,0,0,0), oma=c(4,4,4,1), xlab="Age (years)", ylab="Proportion", mainTitle = paste0(mainTitle, " Selectivity"), cex.axis=1.2, cex.lab=1.3, cex.main=1.4, sobj=sobj, rc=rc, ylim=c(0,1), byrow=F)
				#mtext(paste0(mainTitle, " Selectivity"), side=3, outer=FALSE, line=0, cex=1.5)
			}
			if (p %in% c("eps","png")) dev.off()
		} ## end p (ptypes) loop
	}; eop()
#browser();return()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plt.selectivity


