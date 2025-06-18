##==========================================================
## PBS Stock Synthesis small helper functions:
## ------------------------------------------
## allEqual.........Determine if all elements in vector equal the first element
## cquantile.vec....Calculate cumulative quantile as a vector
## getNpan..........Get panel number when inside a multi-panel plot
## getYrIdx.........Select years for plotting
## is.numStr........Check if strings can be converted to numerics
## med5.95..........Print median (0.05, 0.95) to text
## ptab.............Prior tabulation (LaTeX) for lines in table of priors
## qtab.............Quantile tabulation summary using decimal places
## stab.............Quantile tabulation summary using significant digits
##==========================================================

allEqual <- function(x)
{
  result <- all( x==x[1] )
  result
}

## cquantile.vec------------------------2010-10-20
##  Calculate cumulative quantile as a vector
##  AME doing this, just do one prob at a time 
##  (so it returns a vector not a matrix)
## --------------------------------------------AME
cquantile.vec <- function(z, prob)  # cumulative quantile of vector
{                                   #  prob is a single number
  cquant <- rep(NA, length(z))
  if(length(prob) != 1) stop("length prob should be 1")
  for (i in 1:length(z)) {
    cquant[i] <- quantile(z[1:i], probs = prob, names = FALSE)
  }
  return(cquant)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~cquantile.vec

## getNpan------------------------------2019-05-10
##  Get panel number when inside a multi-panel plot.
## ---------------------------------------------RH
getNpan <- function()
{
	mfg=par()$mfg
	mfg[2]+(mfg[1]-1)*mfg[4]
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getNpan

## getYrIdx-------------------------------2011-08-31
##  Purpose is to return selected years for plotting.
##  Default is to select 5 year increments.
##---------------------------------------------AME
getYrIdx <- function(yrNames, mod=5)
{
  ## Coerce to numeric and select the years modulo "mod".
  yrVals <- as.numeric( gsub("[^[:digit:]]","",yrNames) )
  idx <- yrVals %% mod==0

  ## Select years from character vector yrNames.
  result <- yrNames[ idx ]
  result
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getYrIdx

## is.numStr----------------------------2020-09-20
##  Check if strings can be converted to numerics.
##----------------------------------------------RH
is.numStr <- function(x)
{
	out = sapply(x, function(xx) {
		xx = as.character(xx)
		all(grepl("[[:digit:]]|\\-|\\.|[eE]|[[:space:]]", strsplit(xx,"")[[1]]))
	})
	return(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~is.numStr

## print median (0.05, 0.95) to text
med5.95 <- function(xx.MCMC, dig=0, quants3=tcall(quants3))
{  ## dig is number of dec places
	if (is.null(quants3)) quants3=c(0.05,0.5,0.95)
	big.mark = options()$big.mark; if (is.null(big.mark)) big.mark=","
	mess = paste0( c( prettyNum(round(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=big.mark),
		"~(", prettyNum(round(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=big.mark), 
		",\\,", prettyNum(round(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=big.mark), ")"), collapse="")
	print(mess)
}

## ptab---------------------------------2025-04-08
##  Function to use for priors in table (adapted from PBSawatea).
## ---------------------------------------------RH
ptab <- function(xx)
{
	xx = sub("^\\s+", "", xx)  ## remove leading and trailing whitespace
	xlab = gsub("\\_+"," ",xx[1])
#browser();return()
	xnum =xx[-1]
	xnum[4] = switch(xnum[4], 'Normal'=6, 'No_prior'=0, 'Full_Beta'=2, 'Sym_Beta'=1)
	xnum = lapply(xnum,function(x){
		if(is.numStr(x))   as.numeric(x)
		else if (is.na(x)) "--"
		else x
	})
	xout = paste0(c(xlab, " & ", xnum[[1]], " & (", xnum[[2]], ", ", xnum[[3]], ") & ", xnum[[4]], " & (", xnum[[5]], ", ", xnum[[6]], ") & ", xnum[[7]], " & ", show0(round(xnum[[8]],3),3), " \\\\\\\\"), collapse="")
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ptab

## qtab---------------------------------2021-04-19
## Quantile tabulation summary using decimal places
## --------------------------------------------AME
qtab <- function(xx.MCMC, dig=0, quants3=tcall(quants3))
{  ## dig is number of dec places
	if (is.null(quants3)) quants3=c(0.05,0.5,0.95)
	big.mark = options()$big.mark; if (is.null(big.mark)) big.mark=","
	print(paste0( c( prettyNum(round(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=big.mark),
		" & ", prettyNum(round(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=big.mark),
		" & ", prettyNum(round(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=big.mark)), collapse=""))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~qtab

## stab---------------------------------2021-04-19
## Quantile tabulation summary using significant digits
## --------------------------------------------AME
stab <- function(xx.MCMC, dig=3, quants3=tcall(quants3), print=TRUE)
{  ## dig is number sig digits
	if (is.null(quants3)) quants3=c(0.05,0.5,0.95)
	big.mark = options()$big.mark; if (is.null(big.mark)) big.mark=","
	out = paste0( c( prettyNum(signif(quantile(xx.MCMC, quants3[1]), digits=dig), big.mark=big.mark), 
		" & ", prettyNum(signif(quantile(xx.MCMC, quants3[2]), digits=dig), big.mark=big.mark),
		" & ", prettyNum(signif(quantile(xx.MCMC, quants3[3]), digits=dig), big.mark=big.mark)), collapse="")
	if (print) print(out)
	invisible(out)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~stab
