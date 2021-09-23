##==========================================================
## PBS Stock Synthesis:
##  Overwrite existing functions from other packages (only adnuts at the moment).
## ---------------------------------------------------------
## launch_shinyadmb......Launch shiny but useless GUI system to explore fits (mod.adnuts).
## pairs_admb............Plot pairwise parameter posteriors (mod.adnuts).
## sample_admb...........Bayesian inference of an ADMB model using the no-U-turn sampler (mod.adnuts)
## sample_admb_nuts......Sample ADMB using NUTS (No U-Turn Sampling) algorithm
## sample_admb_parallel..Sample ADMB MCMCs using parallel cores (mod.adnuts)
## sample_admb_rwm.......Sample ADMB using RWM (Random Walk Metropolis) algorithm
##==========================================================


## launch_shinyadmb---------------------2020-11-11
## adnuts Functions: launch_shinyadmb
## --------------------------------------adnuts|RH
launch_shinyadmb=function(fit)
{
	tmp_params=fit$sampler_params
	nms_params=lapply(tmp_params,function(x){x[,setdiff(colnames(x),"lp__")]})
	fit$sampler_params=nms_params
	ttput(fit)
	shinystan::launch_shinystan(.as.shinyadnuts(fit))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~launch_shinyadmb


## pairs_admb---------------------------2021-03-10
## adnuts function: 'pairs_admb'
## --------------------------------------adnuts|RH
pairs_admb = function (fit, diag=c("trace", "acf", "hist"), acf.ylim=c(-1,1), 
   ymult=NULL, axis.col=gray(0.5), pars=NULL, label.cex=0.5, limits=NULL, ...) 
{
	mle <- fit$mle
	posterior <- extract_samples(fit, inc_lp=TRUE)
	chains <- rep(1:dim(fit$samples)[2], each=dim(fit$samples)[1] - fit$warmup)
	divs <- if (fit$algorithm == "NUTS") 
		extract_sampler_params(fit)$divergent__
	else NULL
	ptcex <- 0.1
	divcex <- 0.5
	chaincols <- 1:length(unique(chains))
	wp <- function(par.name) {
		x <- which(mle$par.names == par.name)
		if (length(x) == 0) 
			return(NA)
		if (length(x) > 1) 
			stop("Par matched multiple??")
		return(x)
	}
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	diag <- match.arg(diag)
	par.names <- names(posterior)
	if (is.null(pars)) {
		if (NCOL(posterior) > 10) {
			warning("Only showing first 10 parameters, use 'pars' argument to adjust")
			pars <- par.names[1:10]
		}
		else {
			pars <- par.names[1:NCOL(posterior)]
		}
	}
	else if (is.numeric(pars[1])) {
		pars <- par.names[pars]
	}
	pars.bad <- match(x=pars, table=names(posterior))
	if (any(is.na(pars.bad))) {
		warning("Some par names did not match -- dropped")
		pars <- pars[!is.na(pars.bad)]
	}
	n <- length(pars)
	if (n == 1) 
		stop("This function is only meaningful for >1 parameter")
	posterior <- posterior[, pars]

	## Fix sizing of correlations to be relative to chosen parameters RH 201026
	post.par = posterior
	post.cor = cex.cor = round(cor(post.par),2)
	z.cor    = abs(post.cor) < 1
	cex.cor[z.cor] = scaleVec(abs(post.cor[z.cor]),1,3)

	if (is.null(ymult)) 
		ymult <- rep(1.3, n)
	if (is.null(limits)) {
		limits <- list()
		for (i in 1:n) {
			pp <- wp(pars[i])
			limit.temp <- if (is.na(pp)) NULL
			else mle$est[pp] + c(-1, 1) * 1.96 * mle$se[pp]
			min.temp <- min(posterior[, i], limit.temp[1])
			max.temp <- max(posterior[, i], limit.temp[2])
			margin <- 0.15 * (max.temp - min.temp)
			limits[[i]] <- c(min.temp - margin, max.temp + margin)
		}
	}
	N <- NROW(posterior)
	mypch <- 16
	mycol <- 1
	if (N >= 1000) {
		mycol <- rgb(0, 0, 0, 0.5)
	}
	else if (N >= 10000) {
		mycol <- rgb(0, 0, 0, 0.05)
	}
	if (is.null(divs))
		divs <- rep(0, N)
	par(mfrow=c(n, n), mar=0 * c(0.1, 0.1, 0.1, 0.1), yaxs="i", xaxs="i", mgp=c(0.25, 0.25, 0), tck=-0.02, cex.axis=1, col.axis=axis.col, oma=c(2, 2, 2, 2))
	temp.box <- function() box(col=axis.col, lwd=0.5)
	for (row in 1:n) {
		for (col in 1:n) {
			if (row == col) {
				if (diag == "hist") {
					h <- hist(posterior[, row], plot=F)
					if (is.null(limits)) {
					hist(posterior[, row], axes=F, freq=FALSE, 
						ann=F, ylim=c(0, ymult[row] * max(h$density)), 
						col=gray(0.8), border=gray(0.5))
					}
					else {
					hist(posterior[, row], axes=F, freq=FALSE, 
						ann=F, ylim=c(0, ymult[row] * max(h$density)), 
						col=gray(0.8), border=gray(0.5), xlim=limits[[row]])
					}
					temp.box()
				}
				else if (diag == "acf") {
					acf(posterior[, row], axes=F, ann=F, ylim=acf.ylim)
					temp.box()
				}
				else if (diag == "trace") {
					xlim <- c(1, length(chains[chains == 1]))
					plot(x=0, y=0, type="n", axes=FALSE, ann=FALSE, ylim=limits[[row]], xlim=xlim)
					for (ll in unique(chains)) {
						lines(posterior[chains == ll, row], col=chaincols[ll], lwd=0.1)
					}
					temp.box()
				}
			}
			if (row > col) {
				par(xaxs="r", yaxs="r")
				plot(x=posterior[, col], y=posterior[, row], axes=FALSE, ann=FALSE, pch=mypch, cex=ptcex, col=mycol, xlim=limits[[col]], ylim=limits[[row]], ...)
				points(x=posterior[which(divs == 1), col], y=posterior[which(divs == 1), row], pch=mypch, cex=divcex, col="red")
				p1 <- wp(pars[row])
				p2 <- wp(pars[col])
				if (!is.na(p1) & !is.na(p2)) {
					points(x=mle$est[p2], y=mle$est[p1], pch=16, 
					cex=0.5, col=2)
					if (!requireNamespace("ellipse", quietly=TRUE)) {
					warning("ellipse package needs to be installed to show ellipses")
					}
					else {
					ellipse.temp <- ellipse:::ellipse(x=mle$cor[p2, p1], 
						scale=mle$se[c(p2, p1)], centre=mle$est[c(p2, 
						p1)], npoints=1000, level=0.95)
					lines(ellipse.temp, lwd=0.5, lty=1, col="red")
					}
				}
				par(xaxs="i", yaxs="i")
				temp.box()
			}
			if (row < col) {
				plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1), axes=F, ann=F)
				temp.cor <- round(cor(posterior[, c(row, col)])[1, 2], 2)
#browser();return()
				#legend("center", legend=NA, title=temp.cor, cex=(3 * abs(temp.cor) + 0.25) * 0.5, bty="n")
				#addLegend(0.5,0.5, legend=post.cor[row,col], cex=cex.cor[row,col], adj=c(0.5,0.5), bty="n")
				text(0.5, 0.5, labels=post.cor[row,col], cex=cex.cor[row,col])
				temp.box()
			}
			if (row == n) {
				par(mgp=c(0.05, ifelse(col%%2 == 0, 0, 0.5), 0))
				axis(1, col=axis.col, lwd=0.5)
			}
			if (col == 1 & row > 1) {
				par(mgp=c(0.05, ifelse(row%%2 == 1, 0.15, 0.5), 0))
				axis(2, col=axis.col, lwd=0.5)
			}
			if (col == 1 & row == 1) {
				par(mgp=c(0.05, ifelse(row%%2 == 1, 0.15, 0.5), 0))
				axis(2, col=axis.col, lwd=0.5)
			}
			if (row == 1) 
				mtext(pars[col], line=ifelse(col%%2 == 1, 0.1, 1.1), cex=label.cex)
			if (col == n) 
				mtext(pars[row], side=4, line=ifelse(row%%2 == 1, 0, 1), cex=label.cex)
		}
	}
	return(list(post.cor=post.cor, cex.cor=cex.cor))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pairs_admb


## sample_admb--------------------------2021-03-25
## adnuts Function: sample_admb
## --------------------------------------adnuts|RH
sample_admb = function (model, path=getwd(), iter=2000, init=NULL, chains=3, 
   warmup=NULL, seeds=NULL, thin=1, mceval=FALSE, duration=NULL, 
   parallel=FALSE, cores=NULL, control=NULL, algorithm="NUTS", ...) 
{
	stopifnot(thin >= 1)
	stopifnot(chains >= 1)
	if (is.null(seeds)) 
		seeds <- sample.int(1e+07, size=chains)
	if (iter < 10 | !is.numeric(iter)) 
		stop("iter must be > 10")
	stopifnot(is.character(path))
	stopifnot(is.character(model))
	if (!dir.exists(path)) 
		stop(paste("Folder", path, "does not exist. Check argument 'path'"))
	if (.Platform$OS.type == "windows") {
		sspath="C:/Users/haighr/Files/Archive/Bat"            ## RH 201021
		ff <- file.path(sspath, paste(sub("\\.exe$","",model), ".exe", sep=""))  ## RH 201021
	}
	else {
		ff <- file.path(path, paste("./", model, sep=""))
	}
	if (!file.exists(ff)) 
		stop(paste("File", ff, "not found. Check 'path' and 'model' arguments"))
	control <- adnuts:::.update_control(control)
#browser();return()
	if (is.null(warmup)) 
		warmup <- floor(iter/2)
	if (!(algorithm %in% c("NUTS", "RWM"))) 
		stop("Invalid algorithm specified")
	if (is.null(init)) {
		warning("Using MLE inits for each chain -- strongly recommended to use dispersed inits")
	}
	else if (is.function(init)) {
		init <- lapply(1:chains, function(x) init())
	}
	else if (!is.list(init)) {
		stop("init must be NULL, a list, or a function")
	}
	if (!is.null(init) & length(init) != chains) {
		stop("Length of init does not equal number of chains.")
	}
	trash <- file.remove(list.files()[grep(".psv", x=list.files())])
#browser();return()
	if (!parallel) {
		if (algorithm == "NUTS") {
			#mcmc.out <- lapply(1:chains, function(i) adnuts:::sample_admb_nuts(path=path, model=model, iter=iter, init=init[[i]], chain=i, thin=thin, warmup=warmup, seed=seeds[i], duration=duration, control=control, parallel=FALSE, ...))
			mcmc.out <- lapply(1:chains, function(i) sample_admb_nuts(path=path, model=model, iter=iter, init=init[[i]], chain=i, thin=thin, warmup=warmup, seed=seeds[i], duration=duration, control=control, parallel=FALSE, ...))
		}
		else {
			mcmc.out <- lapply(1:chains, function(i) sample_admb_rwm(path=path, model=model, iter=iter, thin=thin, warmup=warmup, init=init[[i]], chain=i, seed=seeds[i], control=control, duration=duration, parallel=FALSE, ...))
		}
	}
	else {
		if (!requireNamespace("snowfall", quietly=TRUE)) 
			stop("snowfall package not found")
		snowfall::sfInit(parallel=TRUE, cpus=cores)
		snowfall::sfExportAll()  ## export R objects like function 'sample_admb_parallel'
		on.exit(snowfall::sfStop())
		mcmc.out <- snowfall::sfLapply(1:chains, function(i){
			sample_admb_parallel(parallel_number=i, path=path, model=model, duration=duration, 
			algorithm=algorithm, iter=iter, init=init[[i]], warmup=warmup, seed=seeds[i],
			thin=thin, control=control, ...)})
	}
#browser();return()

	## Build output list
	warmup <- mcmc.out[[1]]$warmup
	mle <- adnuts:::.read_mle_fit(model=model, path=path)
	if (is.null(mle)) {
		par.names <- dimnames(mcmc.out[[1]]$samples)[[2]]
		par.names <- par.names[-length(par.names)]
	}
	else {
		par.names <- mle$par.names
	}
	iters <- unlist(lapply(mcmc.out, function(x) dim(x$samples)[1]))
	if (any(iters != iter/thin)) {
		N <- min(iters)
		warning(paste("Variable chain lengths, truncating to minimum=", N))
	}
	else {
		N <- iter/thin
	}
	samples <- array(NA, dim=c(N, chains, 1 + length(par.names)), dimnames=list(NULL, NULL, c(par.names, "lp__")))
	for (i in 1:chains) samples[, i, ] <- mcmc.out[[i]]$samples[1:N,]
	if (algorithm == "NUTS") 
		sampler_params <- lapply(mcmc.out, function(x) x$sampler_params[1:N,])
	else sampler_params <- NULL
	time.warmup <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.warmup)))
	time.total <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.total)))
	cmd <- unlist(lapply(mcmc.out, function(x) x$cmd))
	if (N < warmup) 
		warning("Duration too short to finish warmup period")
	message(paste("... Merging post-warmup chains into main folder:", path))
	samples2 <- do.call(rbind, lapply(1:chains, function(i) samples[-(1:warmup), i, -dim(samples)[3]]))
	adnuts:::.write_psv(fn=model, samples=samples2, model.path=path)
	unbounded <- do.call(rbind, lapply(mcmc.out, function(x) x$unbounded))
#browser();return()
	oldwd <- getwd()
	on.exit(setwd(oldwd))
	setwd(path)
	write.table(unbounded, file="unbounded.csv", sep=",", col.names=FALSE, row.names=FALSE)
	if (mceval) {
		message("... Running -mceval on merged chains")
		#system(paste(model, "-mceval -maxI 0 -nohess"), ignore.stdout=FALSE)  ## this only reports one sample of derived_paramters
		system(paste(model, "-mceval"), ignore.stdout=FALSE)
	}
	covar.est <- cov(unbounded)
	result <- list(samples=samples, sampler_params=sampler_params, time.warmup=time.warmup, time.total=time.total, algorithm=algorithm, warmup=warmup, model=model, max_treedepth=mcmc.out[[1]]$max_treedepth, cmd=cmd, covar.est=covar.est, mle=mle)
	attr(result,"class")=c("adfit","list")
	return(invisible(result))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~sample_admb


## sample_admb_nuts---------------------2021-07-06
##  Sample ADMB using NUTS (No U-Turn Sampling) algorithm
## --------------------------------------adnuts|RH
sample_admb_nuts = function (path, model, iter=2000, init=NULL, chain=1, thin=1, 
   warmup=ceiling(iter/2), seed=NULL, duration=NULL, control=NULL, 
   skip_optimization=TRUE, verbose=TRUE, admb_args=NULL, parallel=FALSE) 
{
	wd.old <- getwd()
	on.exit(setwd(wd.old))
	setwd(path)
	control <- adnuts:::.update_control(control)
	eps <- control$stepsize
	stopifnot(iter >= 1)
	stopifnot(warmup <= iter)
	stopifnot(duration > 0)
	stopifnot(thin >= 1)
	if (is.null(warmup)) 
		stop("Must provide warmup")
	if (thin < 1 | thin > iter) 
		stop("Thin must be >1 and < iter")
	max_td <- control$max_treedepth
	adapt_delta <- control$adapt_delta
	model2 <- adnuts:::.update_model(model)
	if (skip_optimization) {
		cmd <- paste(model2, "-nox -nohess -maxfn 0 -phase 1000 -nuts -mcmc ", iter)
	}
	else {
		cmd <- paste(model2, "-hbf -nuts -mcmc ", iter)
	}
	cmd <- paste(cmd, "-warmup", warmup, "-chain", chain)
	if (!is.null(seed)) 
		cmd <- paste(cmd, "-mcseed", seed)
	if (!is.null(duration)) 
		cmd <- paste(cmd, "-duration", duration)
	cmd <- paste(cmd, "-max_treedepth", max_td, "-adapt_delta", adapt_delta)
	if (!is.null(eps)) 
		cmd <- paste(cmd, "-hyeps", eps)
	if (!is.null(control$adapt_init_buffer)) 
		cmd <- paste(cmd, "-adapt_init_buffer", control$adapt_init_buffer)
	if (!is.null(control$adapt_term_buffer)) 
		cmd <- paste(cmd, "-adapt_term_buffer", control$adapt_term_buffer)
	if (!is.null(control$adapt_window)) 
		cmd <- paste(cmd, "-adapt_window", control$adapt_window)
	if (!is.null(control$refresh)) 
		cmd <- paste(cmd, "-refresh", control$refresh)
	if (control$adapt_mass) 
		cmd <- paste(cmd, "-adapt_mass")
	if (control$adapt_mass_dense) 
		cmd <- paste(cmd, "-adapt_mass_dense")
	metric <- control$metric
	stopifnot(!is.null(metric))
	if (is.matrix(metric)) {
		if (!requireNamespace("matrixcalc", quietly = TRUE)) 
			stop("Package 'matrixcalc' is required to pass a matrix.\n Install it and try again.")
		cor.user <- metric/sqrt(diag(metric) %o% diag(metric))
		if (!matrixcalc::is.positive.definite(x = cor.user)) 
			stop("Invalid mass matrix passed: it is not positive definite.\n Check 'metric' argument or use different option.")
		adnuts:::.write.admb.cov(metric, hbf = 1)
		warning("admodel.cov overwritten, revert admodel_original.cov if needed")
	}
	else if (is.character(metric) && metric == "unit") {
		cmd <- paste(cmd, "-mcdiag")
	}
	else if (is.character(metric) && metric == "mle") {
	}
	else {
		stop("Invalid metric option")
	}
	if (!is.null(init)) {
		cmd <- paste(cmd, "-mcpin init.pin")
		write.table(file = "init.pin", x = unlist(init), row.names = F, col.names = F)
	}
	else {
	}
	if (!is.null(admb_args)) 
		cmd <- paste(cmd, admb_args)
	model2 <- adnuts:::.update_model(model)
	console <- adnuts:::.check_console_printing(parallel)
	progress <- NULL
#browser();return()
	if (console) {
		time <- system.time(system2(model2, cmd, stdout = ifelse(verbose, "", FALSE)))[3]
	}
	else {
		fn <- "mcmc_progress.txt"
		if (file.exists(fn)) 
			file.remove(fn)
		time <- system.time(system2(model2, cmd, stdout = ifelse(verbose, fn, FALSE)))[3]
		if (file.exists(fn)) {
			progress <- readLines("mcmc_progress.txt")
		}
		else {
			warning("Progress output file not found. Try troubleshooting in serial model")
		}
	}
	if (!file.exists("adaptation.csv") | !file.exists("unbounded.csv")) 
		stop(paste0("NUTS failed to run. Command attempted was:\n", cmd))
	sampler_params <- as.matrix(read.csv("adaptation.csv"))
	unbounded <- as.matrix(read.csv("unbounded.csv", header = FALSE))
	dimnames(unbounded) <- NULL
	pars <- adnuts:::.get_psv(model)
	par.names <- names(pars)
	if (!"lp__" %in% dimnames(sampler_params)[[2]]) {
		pars[, "log-posterior"] <- sampler_params[, "energy__"]
	}
	else {
		pars[, "log-posterior"] <- sampler_params[, "lp__"]
		sampler_params <- sampler_params[, -7]
	}
	pars <- as.matrix(pars)
	pars <- pars[seq(1, nrow(pars), by = thin), ]
	unbounded <- unbounded[seq(1, nrow(unbounded), by = thin), ]
	sampler_params <- sampler_params[seq(1, nrow(sampler_params), by = thin), ]
	time.total <- time
	time.warmup <- NA
	warmup <- warmup/thin
	gdump = gc(verbose=FALSE)
	return(list(samples = pars, sampler_params = sampler_params, 
		time.total = time.total, time.warmup = time.warmup, warmup = warmup, 
		max_treedepth = max_td, model = model, par.names = par.names, 
		cmd = cmd, unbounded = unbounded, progress = progress))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~sample_admb_nuts


## sample_admb_parallel-----------------2021-03-26
## adnuts Function: sample_admb_parallel
## --------------------------------------adnuts|RH
sample_admb_parallel = function (parallel_number, path, algorithm, ...) 
{
	olddir <- getwd()
	#newdir <- paste0(file.path(getwd(), path), "_chain_", parallel_number)  ## WTF
	newdir <- paste0(sub("/$","",file.path(path)), "/chain_", parallel_number) ## RH 201021
	adieu = function() {
		unlink(newdir, recursive=TRUE, force=TRUE)
		setwd(olddir)
	}
	on.exit(adieu())
	if (dir.exists(newdir)) {
		unlink(newdir, recursive=TRUE, force=TRUE)
		if (dir.exists(newdir)) 
			stop(paste("Could not remove folder:", newdir))
	}
	dir.create(newdir)
	f.need <- list.files(path, full.names=TRUE, pattern="\\.ss$|\\.cov$|\\.hes$")
	f.need <- grep("runnumber",f.need,invert=TRUE,value=TRUE)
	trash  <- file.copy(from=f.need, to=newdir, overwrite=TRUE, copy.date=TRUE)
	#trash  <- file.copy(from=list.files(path, full.names=TRUE), to=newdir, overwrite=TRUE, copy.date=TRUE)
	trash2 <- file.copy(from="C:/Users/haighr/Files/Archive/Bat/ss.exe", to=newdir, overwrite=TRUE, copy.date=TRUE) ## copy ss executable to temporary chain directories
#return(trash)
#browser();return()
	if (algorithm == "RWM") 
		fit <- sample_admb_rwm(path=newdir, chain=parallel_number, parallel=FALSE, ...)
	if (algorithm == "NUTS") 
		fit <- sample_admb_nuts(path=newdir, chain=parallel_number, parallel=FALSE, ...)
	#unlink(newdir, recursive=TRUE, force=TRUE)
	return(fit)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~sample_admb_parallel


## sample_admb_rwm----------------------2021-07-06
##  Sample ADMB using RWM (Random Walk Metropolis) algorithm
## --------------------------------------adnuts|RH
sample_admb_rwm = function (path, model, iter=2000, thin=1, warmup=ceiling(iter/2), 
   init=NULL, chain=1, seed=NULL, control=NULL, verbose=TRUE, 
   duration=NULL, admb_args=NULL, skip_optimization=TRUE, parallel=FALSE) 
{
	wd.old <- getwd()
	on.exit(setwd(wd.old))
	setwd(path)
	if (any(names(control) != "refresh")) 
		warning("Only refresh control argument is used with RWM, ignoring: ", 
			paste(names(control)[names(control) != "refresh"], 
				collapse = ", "), call. = FALSE)
	refresh <- control$refresh
	if (!is.null(refresh) & !is.numeric(refresh)) 
		stop("Invalid refresh value ", refresh)
	metric <- "mle"
	stopifnot(iter >= 1)
	stopifnot(warmup <= iter)
	stopifnot(duration > 0)
	stopifnot(thin >= 1)
	if (is.null(warmup)) 
		stop("Must provide warmup")
	if (thin < 1 | thin > iter) 
		stop("Thin must be >1 and < iter")
	if (skip_optimization) {
		cmd <- paste("-nox -nohess -maxfn 0 -phase 1000 -rwm -mcmc ", iter)
	}
	else {
		cmd <- paste("-rwm -mcmc ", iter)
	}
	cmd <- paste(cmd, "-mcscale", warmup, "-chain", chain)
	if (!is.null(seed)) 
		cmd <- paste(cmd, "-mcseed", seed)
	if (!is.null(duration)) 
		cmd <- paste(cmd, "-duration", duration)
	cmd <- paste(cmd, "-mcsave", thin)
	if (is.matrix(metric)) {
		cor.user <- metric/sqrt(diag(metric) %o% diag(metric))
		if (!matrixcalc::is.positive.definite(x = cor.user)) 
			stop("Invalid mass matrix, not positive definite")
		adnuts:::.write.admb.cov(metric)
	}
	else if (is.null(metric)) {
	}
	else if (metric == "mle") {
	}
	else if (metric == "unit") {
		cmd <- paste(cmd, "-mcdiag")
	}
	else {
		stop("Invalid metric option")
	}
	if (!is.null(init)) {
		cmd <- paste(cmd, "-mcpin init.pin")
		write.table(file = "init.pin", x = unlist(init), row.names = F, 
			col.names = F)
	}
	if (!is.null(refresh)) 
		cmd <- paste(cmd, "-refresh", refresh)
	if (!is.null(admb_args)) 
		cmd <- paste(cmd, admb_args)
#browser();return()
	model2 <- adnuts:::.update_model(model)
	console <- adnuts:::.check_console_printing(parallel)
	progress <- NULL
	if (console) {
		time <- system.time(system2(model2, cmd, stdout = ifelse(verbose, "", FALSE)))[3]
	}
	else {
		fn <- "mcmc_progress.txt"
		if (file.exists(fn)) 
			file.remove(fn)
		time <- system.time(system2(model2, cmd, stdout = ifelse(verbose, 
			fn, FALSE)))[3]
		if (file.exists(fn)) {
			progress <- readLines("mcmc_progress.txt")
		}
		else {
			warning("Progress output file not found. Try troubleshooting in serial model")
		}
	}
	if (!file.exists("unbounded.csv")) 
		stop(paste0("RWM failed to run. Command attempted was:\n", cmd))
	unbounded <- as.matrix(read.csv("unbounded.csv", header = FALSE))
	dimnames(unbounded) <- NULL
	pars <- adnuts:::.get_psv(model)
	par.names <- names(pars)
	lp <- as.vector(read.table("rwm_lp.txt", header = TRUE)[, 1])
	pars[, "log-posterior"] <- lp
	pars <- as.matrix(pars)
	time.total <- time
	time.warmup <- NA
	warmup <- warmup/thin
	gdump = gc(verbose=FALSE)
	return(list(samples = pars, sampler_params = NULL, time.total = time.total, 
		time.warmup = time.warmup, warmup = warmup, model = model, 
		par.names = par.names, cmd = cmd, unbounded = unbounded, 
		progress = progress))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~sample_admb_rwm

