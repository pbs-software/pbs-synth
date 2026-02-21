##==========================================================
## PBS Stock Synthesis:
##  Borrowed (commandeered) r4ss functions and modified them
##  for use by Offshore Rockfish stock assessments.
## ---------------------------------------------------------
## getSS.output..........SS_output modified to work with ss3.exe v.3.30.23.00
## make.multifig.........Make a multi-panel figure (mod.r4ss)
## plotSS.comparisons....Compare model output from a summary of multiple models
## plotSS.comps..........Plot age proportions and model fit (mod.r4ss)
## plotSS.francis........Plot mean age fits using Francis (2011) methodology
## plotSS.index..........Plot SS model fit to abundance index series
## plotSS.pars...........Plot parameter fits and priors
## plotSS.profile........Plot SS likelihood profiles
## plotSS.rdevs..........Plot reruitment deviations
## plotSS.selex..........Plot selectivity curves with maturity ogive (mod.r4ss)
## plotSS.stdres.........Plot standardised residuals (mod.PBSawatea)
## plotSS.stock.recruit..Plot stock-recruitment function (based on MPDs)
## plotSS.ts.............Plot time series from MPD output from SS (mod.r4ss)
## plotSS.yield..........Plot the equilibrium yield curve for an MPD run
##==========================================================


## getSS.output ------------------------2025-01-02
##  r4ss function 'SS_output' modified to work with
##  ss3.exe v.3.30.23.00, compiled 4 Nov 2024 by
##  Richard Methot (NOAA) using ADMB 13.2
## ----------------------------------------r4ss|RH
getSS.output <- function (dir="C:/myfiles/mymodels/myrun/", dir.mcmc=NULL, 
	repfile="Report.sso", compfile="CompReport.sso", covarfile="covar.sso", 
	forefile="Forecast-report.sso", wtfile="wtatage.ss_new", 
	warnfile="warning.sso", ncols=lifecycle::deprecated(), 
	forecast=TRUE, warn=TRUE, covar=TRUE, readwt=TRUE, 
	verbose=TRUE, printstats=TRUE, hidewarn=FALSE, NoCompOK=TRUE, 
	aalmaxbinrange=4) 
{
	flush.console()
	emptytest <- function(x) {
		sum(!is.na(x) & x == "")/length(x)
	}
	match_report_line <- function(string, obj=rawrep[, 1], substr1=TRUE) {
		match(string, if (substr1) {
			substring(obj, 1, nchar(string))
		}
		else {
			obj
		})
	}
	match_report_table <- function(string1, adjust1, string2=NULL, 
		adjust2=-1, which_blank=1, cols="nonblank", matchcol1=1, 
		matchcol2=1, obj=rawrep, blank_lines=rep_blank_or_hash_lines, 
		substr1=TRUE, substr2=TRUE, header=FALSE, type.convert=FALSE) {
		line1 <- match(string1, if (substr1) {
			substring(obj[, matchcol1], 1, nchar(string1))
		}
		else {
			obj[, matchcol1]
		})
		if (is.null(string2)) {
			line2 <- blank_lines[blank_lines > line1][which_blank]
			if (is.na(line2)) {
				line2 <- nrow(obj)
			}
		}
		else {
			line2 <- match(string2, if (substr2) {
				substring(obj[, matchcol2], 1, nchar(string2))
			}
			else {
				obj[, matchcol2]
			})
		}
		if (is.na(line1) | is.na(line2)) {
			return(NULL)
		}
		if (is.numeric(cols)) {
			out <- obj[(line1 + adjust1):(line2 + adjust2), cols]
		}
		if (cols[1] == "all") {
			out <- obj[(line1 + adjust1):(line2 + adjust2), ]
		}
		if (cols[1] == "nonblank") {
			out <- obj[(line1 + adjust1):(line2 + adjust2), ]
			out <- out[, apply(out, 2, emptytest) < 1]
		}
		if (header && nrow(out) > 0) {
			out[1, out[1, ] == ""] <- "NoName"
			names(out) <- out[1, ]
			out <- out[-1, ]
		}
		if (type.convert) {
			out <- type.convert(out, as.is=TRUE)
		}
		return(out)
	}
	df.rename <- function(df, oldnames, newnames) {
		if (!is.null(df)) {
			for (iname in seq_along(oldnames)) {
				names(df)[names(df) == oldnames[iname]] <- newnames[iname]
			}
		}
		return(df)
	}
	if (lifecycle::is_present(ncols)) {
		lifecycle::deprecate_warn(when="1.46.0", what="getSS.output(ncols)", 
			details="Input 'ncols' no longer needed.")
	}
	if (!is.character(dir) | length(dir) != 1) {
		stop("Input 'dir' should be a character string for a directory")
	}
	shortrepfile <- repfile
	repfile <- file.path(dir, repfile)
	parfile <- r4ss:::get_par_name(dir)
	if (is.na(parfile)) {
		if (!hidewarn) {
			message("Some stats skipped because the .par file not found.")
		}
	}
	if (file.exists(repfile)) {
		if (file.info(repfile)$size > 0) {
			if (verbose) {
				message("Getting header info from:\n  ", repfile)
			}
		}
		else {
			stop("report file is empty: ", repfile)
		}
	}
	else {
		stop("can't find report file: ", repfile)
	}
	rephead <- readLines(con=repfile, n=50)
	SS_versionCode <- rephead[grep("#V", rephead)]
	SS_version <- rephead[grep("Stock_Synthesis", rephead)]
	SS_version <- SS_version[substring(SS_version, 1, 2) != "#C"]
	SS_version <- SS_version[1]
	if (substring(SS_version, 1, 2) == "#V") {
		SS_version <- substring(SS_version, 3)
	}
	if (substring(SS_version, 1, 4) == "3.30") {
		SS_versionshort <- "3.30"
		SS_versionNumeric <- as.numeric(SS_versionshort)
	}
	else {
		SS_versionshort <- toupper(substr(SS_version, 1, 8))
		SS_versionNumeric <- as.numeric(substring(SS_versionshort, 5))
	}
	SS_versionMax <- 3.3
	SS_versionMin <- 3.24
	if (SS_versionNumeric < SS_versionMin | SS_versionNumeric > 
		SS_versionMax) {
		warning("This function tested on SS versions 3.24 and 3.30.\n", 
			"  You are using ", strsplit(SS_version, split=";")[[1]][1], 
			" which MIGHT NOT WORK with this package.")
	}
	else {
		if (verbose) {
			message("This function tested on SS versions 3.24 and 3.30.\n",
				"  You are using ", strsplit(SS_version, split=";")[[1]][1],
				" which SHOULD work with this package.")
		}
	}
	SS_versionFull=strsplit(SS_version,split=";")[[1]][1]  ## (RH 241217)

	findtime <- function(lines) {
		time <- strsplit(lines[grep("ime", lines)], "ime: ")[[1]]
		if (length(time) < 2) {
			return()
		}
		else {
			return(time[2])
		}
	}
	repfiletime <- findtime(rephead)
	if (verbose) {
		message("Report file time:", repfiletime)
	}
	comp <- FALSE
	if (is.null(compfile)) {
		if (verbose) {
			message("Skipping CompReport because 'compfile=NULL'")
		}
	}
	else {
		if (file.exists(file.path(dir, compfile))) {
			compfile <- file.path(dir, compfile)
			comphead <- readLines(con=compfile, n=30)
			compskip <- grep("Composition_Database", comphead)
			if (length(compskip) == 0) {
				if (verbose) {
					message("No composition data, possibly because detailed output", " is turned off in the starter file.")
				}
			}
			else {				compend <- grep(" end ", comphead)
				if (length(compend) == 0) {
					compend <- 999
				}
				comptime <- findtime(comphead)
				if (is.null(comptime) || is.null(repfiletime)) {
					message("problem comparing the file creation times:\n", 
					"  Report.sso:", repfiletime, "\n", "  CompReport.sso:", comptime, "\n")
				}
				else {
					if (comptime != repfiletime) {
						message("CompReport time:", comptime, "\n")
						stop(shortrepfile, " and ", compfile, " were from different model runs.")
					}
				}
				comp <- TRUE
			}
		}
		else {
			if (!is.null(compfile)) {
				if (!NoCompOK) {
					stop("Missing ", compfile, ". Change the 'compfile' input, rerun model to get the file,", 
					" or change input to 'NoCompOK=TRUE'")
				}
				else {
					message("Composition file not found: ", compfile)
				}
			}
		}
	}
	if (verbose) {
		message("Reading full report file")
	}
	flush.console()
	ncols  <- r4ss:::get_ncol(repfile)
	rawrep <- read.table(file=repfile, col.names=1:ncols, 
		fill=TRUE, quote="", colClasses="character", nrows=-1, 
		comment.char="", blank.lines.skip=FALSE)
	rep_blank_lines <- which(apply(rawrep, 1, emptytest) == 1)
	rep_hash_lines <- which(rawrep[, 1] == "#" & apply(rawrep[, -1], 1, emptytest) == 1)
	rep_blank_or_hash_lines <- sort(unique(c(rep_blank_lines, rep_hash_lines)))
	nonblanks <- apply(rawrep, 2, emptytest) < 1
	maxnonblank <- max(0, (1:ncols)[nonblanks == TRUE])
	if (maxnonblank == ncols) {
		stop("all columns are used and some data may been missed,\n", 
			"  increase 'ncols' input above current value (ncols=", ncols, ")")
	}
	custom <- !is.na(match_report_line(string="report:1", obj=rawrep[, 2]))
	if (verbose) {
		if ((maxnonblank + 1) == ncols) {
			message("Got all columns using ncols=", ncols)
		}
		if ((maxnonblank + 1) < ncols) {
			message("Got all columns. To speed code, use ncols=", maxnonblank + 1, " in the future.")
		}
		message("Got Report file")
	}
	flush.console()
	if (forecast) {
		forecastname <- file.path(dir, forefile)
		temp <- file.info(forecastname)$size
		if (is.na(temp) | temp == 0) {
			if (verbose) {
				message("Forecast-report.sso file is missing or empty.")
			}
		}
		else {
			rawforecast1 <- read.table(file=forecastname, col.names=1:ncols, 
				fill=TRUE, quote="", colClasses="character", nrows=-1)
			grab <- rawforecast1[, 1]
			nforecastyears <- as.numeric(rawforecast1[grab %in% c("N_forecast_yrs:"), 2])
			nforecastyears <- nforecastyears[1]
			sprtarg <- as.numeric(rawforecast1[match_report_line("SPR_target", rawforecast1[, 1]), 2])
			target_definitions <- grep("_as_target", rawforecast1[, 1], value=TRUE)
			if (length(target_definitions) == 0) {
				btarg <- as.numeric(rawforecast1[match_report_line("Btarget", rawforecast1[, 1]), 2])
			}
			else {
				if ("Ratio_SSB/B0_as_target" %in% target_definitions) {
					btarg <- as.numeric(rawforecast1[match_report_line("Ratio_target", rawforecast1[, 1]), 2])
				}
				if ("F0.1_as_target" %in% target_definitions) {
					btarg <- -999
				}
			}
		}
	}
	else {
		if (verbose) {
			message("You skipped the forecast file.")
		}
	}
	if (!exists("btarg")) {
		nforecastyears <- NA
		sprtarg <- -999
		btarg <- -999
		if (verbose) {
			message("  setting SPR target and Biomass target to -999.", 
				"  Lines won't be drawn for these targets by SS_plots unless", 
				"  'sprtarg' and 'btarg' are provided as inputs.")
		}
	}
	minbthresh <- -999
	if (!is.na(btarg) & btarg == 0.4) {
		if (verbose) {
			message("Setting minimum biomass threshhold to 0.25", 
				"  based on US west coast assumption associated with biomass target of 0.4.", 
				"  (can replace or override in SS_plots by setting 'minbthresh')")
		}
		minbthresh <- 0.25
	}
	if (!is.na(btarg) & btarg == 0.25) {
		if (verbose) {
			message("Setting minimum biomass threshhold to 0.125", 
				"  based on US west coast assumption associated with flatfish target of 0.25.", 
				"  (can replace or override in SS_plots by setting 'minbthresh')")
		}
		minbthresh <- 0.125
	}
	flush.console()
	logfile_name <- dir(dir, pattern=".log$")
	logfile_name <- logfile_name[logfile_name != "fmin.log"]
	if (length(logfile_name) > 1) {
		filetimes <- file.info(file.path(dir, logfile_name))$mtime
		logfile_name <- logfile_name[filetimes == max(filetimes)]
		if (verbose) {
			message("Multiple files in directory match pattern *.log\n", 
				"choosing most recently modified file:", logfile_name, 
				"\n")
		}
	}
	if (length(logfile_name) == 1 && file.info(file.path(dir, 
		logfile_name))$size > 0) {
		logfile <- readLines(file.path(dir, logfile_name))
		logfile <- grep("^size", logfile, value=TRUE)
		if (length(logfile) == 0) {
			warning(logfile_name, " does not contain information on the size of temporary files.")
			logfile <- NA
		}
		else {
			logfile <- tidyr::separate(as.data.frame(logfile), 
				col=1, into=c("File", "Size"), sep="=")
			names(logfile) <- c("TempFile", "Size")
			logfile[["Size"]] <- as.numeric(logfile[["Size"]])
			maxtemp <- max(logfile[["Size"]])
			if (verbose) {
				if (maxtemp == 0) {
					message("Got log file. There were NO temporary files were written in this run.")
				}
				else {
					message("Temporary files were written in this run.")
				}
			}
		}
	}
	else {
		logfile <- NA
		if (verbose) {
			message("No non-empty log file in directory or too many files ", 
				" matching pattern *.log")
		}
	}
	if (warn) {
		warnname <- file.path(dir, warnfile)
		if (!file.exists(warnname)) {
			message(warnfile, " file not found")
			warnrows <- NA
			warnlines <- NA
		}
		else {
			warnlines <- readLines(warnname, warn=FALSE)
			warnrows <- length(warnlines)
			if (verbose && warnrows > 0) {
				message("Got warning file. Final line:", tail(warnlines, 1))
			}
		}
	}
	else {
		if (verbose) {
			message("You skipped the warnings file")
		}
		warnrows <- NA
		warnlines <- NA
	}
	if (verbose) {
		message("Finished reading files")
	}
	flush.console()
	sizeselex <- match_report_table("LEN_SELEX", 6, header=TRUE, 
		type.convert=TRUE)
	sizeselex <- df.rename(sizeselex, oldnames=c("fleet", "year", 
		"seas", "gender", "morph", "label"), newnames=c("Fleet", 
		"Yr", "Seas", "Sex", "Morph", "Label"))
	rawdefs <- match_report_table("DEFINITIONS", 1, which_blank=1, 
		blank_lines=rep_blank_lines)
	Length_comp_error_controls <- NULL
	Age_comp_error_controls <- NULL
	Size_comp_error_controls <- NULL
	if ("Jitter:" %in% rawdefs[["X1"]]) {
		get.def <- function(string) {
			row <- grep(string, rawdefs[["X1"]])[1]
			if (length(row) > 0) {
				return(as.numeric(rawdefs[row, 2]))
			}
			else {
				return(NULL)
			}
		}
		N_seasons <- nseasons <- get.def("N_seasons")
		N_sub_seasons <- get.def("N_sub_seasons")
		Season_Durations <- seasdurations <- as.numeric(rawdefs[grep("Season_Durations", 
			rawdefs[["X1"]]), 1 + 1:nseasons])
		Spawn_month <- spawnmonth <- get.def("Spawn_month")
		Spawn_seas <- spawnseas <- get.def("Spawn_seas")
		Spawn_timing_in_season <- get.def("Spawn_timing_in_season")
		N_areas <- nareas <- get.def("N_areas")
		Start_year <- startyr <- get.def("Start_year")
		End_year <- endyr <- get.def("End_year")
		Retro_year <- get.def("Retro_year")
		N_forecast_yrs <- get.def("N_forecast_yrs")
		N_sexes <- nsexes <- get.def("N_sexes")
		Max_age <- accuage <- get.def("Max_age")
		Empirical_wt_at_age <- get.def("Empirical_wt_at_age")
		N_bio_patterns <- get.def("N_bio_patterns")
		N_platoons <- get.def("N_platoons")
		NatMort_option <- get.def("NatMort")
		GrowthModel_option <- get.def("GrowthModel")
		Maturity_option <- get.def("Maturity")
		Fecundity_option <- get.def("Fecundity")
		Start_from_par <- get.def("Start_from_par")
		Do_all_priors <- get.def("Do_all_priors")
		Use_softbound <- get.def("Use_softbound")
		N_nudata <- get.def("N_nudata")
		Max_phase <- get.def("Max_phase")
		Current_phase <- get.def("Current_phase")
		Jitter <- get.def("Jitter")
		ALK_tolerance <- get.def("ALK_tolerance")
		fleetdefs <- rawdefs[tail(grep("Fleet", rawdefs[["X1"]]), 
			1):nrow(rawdefs), ]
		names(fleetdefs) <- fleetdefs[1, ]
		fleetdefs <- fleetdefs[-1, ]
		fleetdefs <- fleetdefs[, 1:grep("fleet_name", tolower(names(fleetdefs)))]
		fleetdefs <- type.convert(fleetdefs, as.is=TRUE)
		fleetdefs <- df.rename(fleetdefs, oldnames=c("fleet_name"), 
			newnames=c("Fleet_name"))
		fleet_type <- fleetdefs[["fleet_type"]]
		fleet_timing <- fleetdefs[["timing"]]
		fleet_area <- fleetdefs[["area"]]
		catch_units <- fleetdefs[["catch_units"]]
		survey_units <- fleetdefs[["survey_units"]]
		survey_error <- fleetdefs[["survey_error"]]
		fleet_ID <- fleetdefs[["Fleet"]]
		IsFishFleet <- fleet_type <= 2
		nfishfleets <- sum(IsFishFleet)
		FleetNames <- fleetdefs[["Fleet_name"]]
		nfleets <- max(fleet_ID)
		seasfracs <- round(12 * cumsum(seasdurations))/12
		seasfracs <- seasfracs - seasdurations/2
		if ("Length_comp_error_controls" %in% rawdefs[["X1"]]) {
			Length_comp_error_controls <- match_report_table("Length_comp_error_controls", 
				adjust1=1, header=TRUE, type.convert=TRUE)
			if (nrow(Length_comp_error_controls) > 0) {
				present_Length_comp_error_controls <- TRUE
			}
		}
		if (exists("Length_comp_error_controls") & exists("present_Length_comp_error_controls")) {
			names(Length_comp_error_controls)[names(Length_comp_error_controls) == 
				"NoName"] <- c("NoName", "Fleet_name")
			Length_comp_error_controls <- dplyr::select(Length_comp_error_controls, 
				-NoName)
		}
		if ("Age_comp_error_controls" %in% rawdefs[["X1"]]) {
			Age_comp_error_controls <- match_report_table("Age_comp_error_controls", 
				adjust1=1, header=TRUE, type.convert=TRUE)
			if (nrow(Age_comp_error_controls) > 0) {
				present_Age_comp_error_controls <- TRUE
			}
		}
		if (exists("Age_comp_error_controls") & exists("present_Age_comp_error_controls") > 0) {
			names(Age_comp_error_controls)[names(Age_comp_error_controls) == 
				"NoName"] <- c("NoName", "Fleet_name")
			Age_comp_error_controls <- dplyr::select(Age_comp_error_controls, 
				-NoName)
		}
		if ("Size_comp_error_controls" %in% rawdefs[["X1"]]) {
			Size_comp_error_controls <- dplyr::rename(match_report_table("Size_comp_error_controls", 
				adjust1=1, header=TRUE, type.convert=TRUE), 
				Sz_method="#_Sz_method")
		}
	}
	else {
		nseasons <- as.numeric(rawdefs[grep("N_seasons", rawdefs[, 1]), 2])
		seasdurations <- as.numeric(rawdefs[grep("Season_Durations", 
			rawdefs[, 1]), 1 + 1:nseasons])
		seasfracs <- round(12 * cumsum(seasdurations))/12
		seasfracs <- seasfracs - seasdurations/2
		if (SS_versionNumeric >= 3.3) {
			FleetNames <- as.character(rawdefs[grep("fleet_names", 
				rawdefs[["X1"]]), -1])
			FleetNames <- FleetNames[!is.na(FleetNames) & FleetNames != 
				""]
			nfleets <- length(FleetNames)
			fleet_ID <- 1:nfleets
			fleetdefs <- tail(rawdefs, nfleets + 1)
			fleetdefs <- fleetdefs[, apply(rawdefs[-(1:3), ], 
				2, emptytest) < 1]
			fleetdefs[fleetdefs == ""] <- NA
			if (fleetdefs[1, 1] == "#_rows") {
				fleetdefs <- fleetdefs[-1, 1:7]
				names(fleetdefs) <- c("fleet_type", "timing", "area",
					"catch_units", "catch_mult", "survey_units", "survey_error")
			}
			else {
				names(fleetdefs) <- fleetdefs[1, ]
				names(fleetdefs)[1] <- "fleet"
				fleetdefs <- fleetdefs[-1, ]
			}
			fleetdefs <- type.convert(fleetdefs, as.is=TRUE)
			fleet_type <- fleetdefs[["fleet_type"]]
			fleet_timing <- fleetdefs[["timing"]]
			fleet_area <- fleetdefs[["area"]]
			catch_units <- fleetdefs[["catch_units"]]
			equ_catch_se <- fleetdefs[["equ_catch_se"]]
			catch_se <- fleetdefs[["catch_se"]]
			survey_units <- fleetdefs[["survey_units"]]
			survey_error <- fleetdefs[["survey_error"]]
			IsFishFleet <- fleet_type <= 2
		}
		else {
			fleetdefs <- rawdefs[-(1:3), apply(rawdefs[-(1:3), ], 2, emptytest) < 1]
			fleetdefs[fleetdefs == ""] <- NA
			lab <- fleetdefs[["X1"]]
			fleet_ID <- as.numeric(fleetdefs[grep("fleet_ID", lab), -1])
			names(fleetdefs) <- c("Label", paste("Fleet", fleet_ID, sep=""))
			FleetNames <- as.character(fleetdefs[grep("fleet_names", lab), -1])
			fleet_area <- as.numeric(fleetdefs[grep("fleet_area", lab), -1])
			catch_units <- as.numeric(fleetdefs[grep("Catch_units", lab), -1])
			catch_error <- as.numeric(fleetdefs[grep("Catch_error", lab), -1])
			survey_units <- as.numeric(fleetdefs[grep("Survey_units", lab), -1])
			survey_error <- as.numeric(fleetdefs[grep("Survey_error", lab), -1])
			IsFishFleet <- !is.na(catch_units)
			nfleets <- length(FleetNames)
		}
		begin <- match_report_line("TIME_SERIES") + 2
		end <- match_report_line("SPR_series") - 2
		nfishfleets <- sum(IsFishFleet)
		nsexes <- length(unique(as.numeric(sizeselex[["Sex"]])))
		nareas <- max(as.numeric(rawrep[begin:end, 1]))
		startyr <- min(as.numeric(rawrep[begin:end, 2])) + 2
		temptime <- rawrep[begin:end, 2:3]
		endyr <- max(as.numeric(temptime[temptime[, 2] == "TIME", 1]))
		tempaccu <- as.character(rawrep[match_report_line("Natural_Mortality") + 1, -(1:5)])
		accuage <- max(as.numeric(tempaccu[tempaccu != ""]))
	}
	if (comp) {
		ncols.compfile <- r4ss:::get_ncol(compfile, skip=3)
		allbins <- read.table(file=compfile, col.names=1:ncols.compfile, fill=TRUE, colClasses="character", skip=3, nrows=25)
		lbins <- as.numeric(allbins[grep("Size_Bins_dat", allbins[, 1]) + 2, -1])
		lbins <- lbins[!is.na(lbins)]
		nlbins <- length(lbins)
		lbinspop <- as.numeric(allbins[grep("Size_Bins_pop", allbins[, 1]) + 2, -1])
		lbinspop <- lbinspop[!is.na(lbinspop)]
		nlbinspop <- length(lbinspop)
		Lbin_method <- as.numeric(allbins[match_report_line("Method_for_Lbin_definition", allbins[, 1]), 2])
		if (compend == compskip + 2) {
			message("It appears that there is no composition data in CompReport.sso")
			comp <- FALSE
			agebins <- NA
			sizebinlist <- NA
			nagebins <- length(agebins)
		}
		else {
			col.names <- as.character(read.table(file=compfile, skip=compskip, nrows=1, colClasses="character"))
			rawcompdbase <- read.table(file=compfile, col.names=col.names, fill=TRUE, colClasses="character", skip=compskip, nrows=-1)
			names(rawcompdbase) <- rawcompdbase[1, ]
			names(rawcompdbase)[names(rawcompdbase) == "Used?"] <- "Used"
			endfile <- grep("End_comp_data", rawcompdbase[, 1])
			compdbase <- rawcompdbase[2:(endfile - 2), ]
			compdbase <- df.rename(compdbase, oldnames=c("Pick_sex", "Pick_gender", "Gender", "N", "Rep"), 
				newnames=c("Sexes", "Sexes", "Sex", "Nsamp_adj", "Repl."))
			duplicates <- duplicated(dplyr::select(compdbase, -Cum_obs, -Cum_exp))
			if (verbose) {
				message("Removing ", sum(duplicates), " out of ", nrow(compdbase),
					" rows in CompReport.sso which are duplicates.")
			}
			compdbase <- compdbase[!duplicates, ]
			compdbase[["sex"]] <- compdbase[["Sexes"]]
			compdbase[["sex"]][compdbase[["Sexes"]] == 3] <- compdbase[["Sex"]][compdbase[["Sexes"]] == 3]
			if (substr(SS_version, 1, 9) == "SS-V3.24f") {
				if (!hidewarn) {
					message("Correcting for bug in tag data output associated with SSv3.24f\n")
				}
				tag1rows <- compdbase[["Sexes"]] == "TAG1"
				if (any(tag1rows)) {
					tag1 <- compdbase[tag1rows, ]
					tag1new <- tag1
					tag1new[, 4:23] <- tag1new[, 3:22]
					tag1new[["Yr.S"]] <- tag1new[["Yr"]]
					tag1new[["Yr"]] <- floor(as.numeric(tag1new[["Yr"]]))
					compdbase[tag1rows, ] <- tag1new
				}
			}
			compdbase <- compdbase[compdbase[["Obs"]] != "", ]
			compdbase[compdbase == "_"] <- NA
			compdbase[["Used"]][is.na(compdbase[["Used"]])] <- "yes"
			if (!("SuprPer" %in% names(compdbase))) {
				compdbase[["SuprPer"]] <- "No"
			}
			compdbase[["SuprPer"]][is.na(compdbase[["SuprPer"]])] <- "No"
			n <- sum(is.na(compdbase[["Nsamp_adj"]]) & compdbase[["Used"]] != 
				"skip" & compdbase[["Kind"]] != "TAG2")
			if (n > 0) {
				warning(n, " rows from composition database have NA sample size\n", 
					"but are not part of a super-period. (Maybe input as N=0?)\n")
			}
			compdbase <- type.convert(compdbase, as.is=TRUE)
			if (nseasons > 1) {
				compdbase[["YrSeasName"]] <- paste(floor(compdbase[["Yr"]]), "s", compdbase[["Seas"]], sep="")
			}
			else {
				compdbase[["YrSeasName"]] <- compdbase[["Yr"]]
			}
			if (!"Yr.S" %in% names(compdbase)) {
				if (any(floor(compdbase[["Yr"]]) != compdbase[["Yr"]])) {
					compdbase[["Yr.S"]] <- compdbase[["Yr"]]
					compdbase[["Yr"]] <- floor(compdbase[["Yr"]])
				}
				else {
					compdbase[["Yr.S"]] <- compdbase[["Yr"]] + (0.5/nseasons) * compdbase[["Seas"]]
				}
			}
			compdbase[["Lbin_range"]] <- compdbase[["Lbin_hi"]] - compdbase[["Lbin_lo"]]
			compdbase[["Lbin_mid"]] <- 0.5 * (compdbase[["Lbin_lo"]] + compdbase[["Lbin_hi"]])
			Lbin_range <- compdbase[["Lbin_range"]]
			if (is.null(Lbin_range)) {
				notconditional <- TRUE
				conditional <- FALSE
			}
			else {
				notconditional <- !is.na(Lbin_range) & Lbin_range > aalmaxbinrange
				conditional <- !is.na(Lbin_range) & Lbin_range <= aalmaxbinrange
			}
			if ("skip" %in% compdbase[["SuprPer"]]) {
				compdbase[["Used"]][compdbase[["SuprPer"]] == "skip"] <- "skip"
				compdbase[["SuprPer"]][compdbase[["SuprPer"]] == "No"]
			}
			if (SS_versionNumeric >= 3.22) {
				lendbase <- compdbase[compdbase[["Kind"]]  == "LEN" & compdbase[["Used"]] != "skip", ]
				sizedbase <- compdbase[compdbase[["Kind"]] == "SIZE" & compdbase[["Used"]] != "skip", ]
				agedbase <- compdbase[compdbase[["Kind"]]  == "AGE" & compdbase[["Used"]] != "skip" & notconditional, ]
				condbase <- compdbase[compdbase[["Kind"]]  == "AGE" & compdbase[["Used"]] != "skip" & conditional, ]
				morphcompdbase <- compdbase[compdbase[["Kind"]] ==  "GP%" & compdbase[["Used"]] != "skip", ]
			}
			else {
				lendbase <- compdbase[compdbase[["Kind"]] == "LEN" & (compdbase[["SuprPer"]] == "Sup" | 
					(!is.na(compdbase[["Nsamp_adj"]]) & compdbase[["Nsamp_adj"]] > 0)), ]
				sizedbase <- compdbase[compdbase[["Kind"]] == "SIZE" & (compdbase[["SuprPer"]] == "Sup" | 
					(!is.na(compdbase[["Nsamp_adj"]]) & compdbase[["Nsamp_adj"]] > 0)), ]
				agedbase <- compdbase[compdbase[["Kind"]] == "AGE" & (compdbase[["SuprPer"]] == "Sup" | 
					(!is.na(compdbase[["Nsamp_adj"]]) & compdbase[["Nsamp_adj"]] > 0)) & notconditional, ]
				condbase <- compdbase[compdbase[["Kind"]] == "AGE" & (compdbase[["SuprPer"]] == "Sup" | 
					(!is.na(compdbase[["Nsamp_adj"]]) & compdbase[["Nsamp_adj"]] > 0)) & conditional, ]
			}
			ghostagedbase <- compdbase[compdbase[["Kind"]] == "AGE" & 
				compdbase[["Used"]] == "skip" & compdbase[["SuprPer"]] == "No" & notconditional, ]
			ghostcondbase <- compdbase[compdbase[["Kind"]] == "AGE" & 
				compdbase[["Used"]] == "skip" & compdbase[["SuprPer"]] == "No" & conditional, ]
			ghostlendbase <- compdbase[compdbase[["Kind"]] == "LEN" & 
				compdbase[["Used"]] == "skip" & compdbase[["SuprPer"]] == "No", ]
			compdbase[["Kind"]][compdbase[["Kind"]] == "L@A" & 
				compdbase[["Ageerr"]] < 0] <- "W@A"
			if (!is.null(sizedbase) && nrow(sizedbase) > 0) {
				sizedbase[["bio.or.num"]] <- c("bio", "num")[sizedbase[["Lbin_lo"]]]
				sizedbase[["units"]] <- c("kg", "lb", "cm", "in")[sizedbase[["Lbin_hi"]]]
				sizedbase[["method"]] <- sizedbase[["Ageerr"]]
				if (any(sizedbase[["units"]] %in% c("lb", "in"))) {
					if (verbose) {
						message("Note: converting bins in generalized size comp data ", 
							" in sizedbase back to the original units of lbs or inches.")
					}
				}
				sizedbase[["Bin"]][sizedbase[["units"]] == "lb"] <- sizedbase[["Bin"]][sizedbase[["units"]] == "lb"]/0.4536
				sizedbase[["Bin"]][sizedbase[["units"]] == "in"] <- sizedbase[["Bin"]][sizedbase[["units"]] == "in"]/2.54
				sizebinlist <- list()
				for (imethod in 1:max(sizedbase[["method"]])) {
					tmp <- sort(unique(sizedbase[["Bin"]][sizedbase[["method"]] == imethod]))
					if (length(tmp) == 0) 
						tmp <- NULL
					sizebinlist[[paste("size_method_", imethod, sep="")]] <- tmp
				}
			}
			else {
				sizebinlist <- NA
			}
			if (is.null(compdbase[["Nsamp_adj"]])) {
				good <- TRUE
			}
			else {
				good <- !is.na(compdbase[["Nsamp_adj"]])
			}
			ladbase <- compdbase[compdbase[["Kind"]] == "L@A" & good, ]
			wadbase <- compdbase[compdbase[["Kind"]] == "W@A" & good, ]
			tagdbase1 <- compdbase[compdbase[["Kind"]] == "TAG1", ]
			tagdbase2 <- compdbase[compdbase[["Kind"]] == "TAG2", ]
			if (verbose) {
				message("CompReport file separated by this code as follows", " (rows=Ncomps*Nbins):\n", 
					if (nrow(lendbase) > 0) {
						paste0("  ", nrow(lendbase), " rows of length comp data\n")
					}, if (nrow(sizedbase) > 0) {
						paste0("  ", nrow(sizedbase), " rows of generalized size comp data\n")
					}, if (nrow(agedbase) > 0) {
						paste0("  ", nrow(agedbase), " rows of age comp data\n")
					}, if (nrow(condbase) > 0) {
						paste0("  ", nrow(condbase), " rows of conditional age-at-length data\n")
					}, if (nrow(ghostagedbase) > 0) {
						paste0("  ", nrow(ghostagedbase), " rows of ghost fleet age comp data\n")
					}, if (nrow(ghostcondbase) > 0) {
						paste0("  ", nrow(ghostcondbase), " rows of ghost fleet conditional age-at-length data\n")
					}, if (nrow(ghostlendbase) > 0) {
						paste0("  ", nrow(ghostlendbase), " rows of ghost fleet length comp data\n")
					}, if (nrow(ladbase) > 0) {
						paste0("  ", nrow(ladbase), " rows of mean length at age data\n")
					}, if (nrow(wadbase) > 0) {
						paste0("  ", nrow(wadbase), " rows of mean weight at age data\n")
					}, if (nrow(tagdbase1) > 0) {
						paste0("  ", nrow(tagdbase1), " rows of 'TAG1' comp data\n")
					}, if (nrow(tagdbase2) > 0) {
						paste0("  ", nrow(tagdbase2), " rows of 'TAG2' comp data")
					}, if (nrow(morphcompdbase) > 0) {
						paste0("  ", nrow(morphcompdbase), " rows of morph comp data")
					}
				)
			}
			if (nrow(agedbase) > 0) {
				Lbin_ranges <- as.data.frame(table(agedbase[["Lbin_range"]]))
				names(Lbin_ranges)[1] <- "Lbin_hi-Lbin_lo"
				if (length(unique(agedbase[["Lbin_range"]])) > 1) {
					warning("different ranges of Lbin_lo to Lbin_hi found in age comps.\n", 
						paste(utils::capture.output(print(Lbin_ranges)), collapse="\n"),
						"\n consider increasing 'aalmaxbinrange' to designate\n", 
						"some of these data as conditional age-at-length.")
				}
				agebins <- sort(unique(agedbase[["Bin"]][!is.na(agedbase[["Bin"]])]))
			}
			else {
				if (nrow(condbase) > 0) {
					agebins <- sort(unique(condbase[["Bin"]][!is.na(condbase[["Bin"]])]))
				}
				else {
					agebins <- NA
				}
			}
			nagebins <- length(agebins)
		}
	}
	else {
		lbins     <- NA
		nlbins    <- NA
		lbinspop  <- NA
		nlbinspop <- ncol(sizeselex) - 5
		agebins   <- NA
		nagebins  <- NA
		Lbin_method <- 2
		sizebinlist <- NA
	}
	morph_indexing <- match_report_table("MORPH_INDEXING", 1, header=TRUE, type.convert=TRUE)
	morph_indexing <- df.rename(morph_indexing, oldnames=c("Gpattern", 
		"Bseas", "BirthSeason", "Gender"), newnames=c("GP", 
		"BirthSeas", "BirthSeas", "Sex"))
	if (!is.null(morph_indexing)) {
		ngpatterns <- max(morph_indexing[["GP"]])
	}
	else {
		ngpatterns <- NULL
	}
	if (verbose) {
		message("Finished dimensioning")
	}
	flush.console()
	stats <- list()
	stats[["SS_version"]] <- SS_version
	stats[["SS_versionshort"]] <- SS_versionshort
	stats[["SS_versionNumeric"]] <- SS_versionNumeric
	stats[["StartTime"]] <- paste(as.character(match_report_table("StartTime", 0, "StartTime", 0, cols=1:6)), collapse=" ")
	stats[["RunTime"]] <- paste(as.character(match_report_table("StartTime", 2, "StartTime", 2, cols=4:9)), collapse=" ")

	returndat <- list()
	returndat[["SS_Version"]]=SS_versionFull

	tempfiles <- match_report_table("Data_File", 0, "Control_File", 0, cols=1:2)
	stats[["Files_used"]] <- paste(c(tempfiles[1, ], tempfiles[2, ]), collapse=" ")
	returndat[["Data_File"]] <- tempfiles[1, 2]
	returndat[["Control_File"]] <- tempfiles[2, 2]

	log_det_hessian <- match_report_table("Hessian", 0, "Hessian", 0, cols=2)
	if (log_det_hessian == "Not") {
		covar <- FALSE
		log_det_hessian <- NA
	}
	stats[["log_det_hessian"]] <- as.numeric(log_det_hessian)
	Final_phase <- match_report_table("Final_phase", 0, "Final_phase", 0, cols=2)
	if (!is.null(Final_phase)) {
		stats[["Final_phase"]] <- as.numeric(Final_phase)
	}
	N_iterations <- match_report_table("N_iterations", 0, "N_iterations", 0, cols=2)
	if (!is.null(N_iterations)) {
		stats[["N_iterations"]] <- as.numeric(N_iterations)
	}
	stats[["Nwarnings"]] <- warnrows
	if (length(warn) > 20) {
		warn <- c(warn[1:20], paste("Note:", length(warn) - 20, 
			"additional lines truncated. Look in", warnfile, "file to see full list."))
	}
	stats[["warnings"]] <- warnlines
	rawlike <- match_report_table("LIKELIHOOD", 2, "Fleet:", -2)
	laplace_line <- which(rawlike[, 1] == "#_info_for_Laplace_calculations")
	if (length(laplace_line) > 0) {
		rawlike <- rawlike[-laplace_line, ]
	}
	like <- data.frame(signif(as.numeric(rawlike[, 2]), digits=7))
	names(like) <- "values"
	rownames(like) <- rawlike[, 1]
	lambdas <- rawlike[, 3]
	lambdas[lambdas == ""] <- NA
	lambdas <- as.numeric(lambdas)
	like[["lambdas"]] <- lambdas
	if (length(laplace_line) > 0) {
		stats[["likelihoods_used"]] <- like[1:(laplace_line - 1), ]
		stats[["likelihoods_laplace"]] <- like[laplace_line:nrow(like), ]
	}
	else {
		stats[["likelihoods_used"]] <- like
		stats[["likelihoods_laplace"]] <- NULL
	}
	likelihoods_by_fleet <- match_report_table("Fleet:", 0, header=TRUE)
	if (!is.null(likelihoods_by_fleet) && "Parm_devs_detail" %in% likelihoods_by_fleet[, 1]) {
		likelihoods_by_fleet <- match_report_table("Fleet:", 0, "Parm_devs_detail", -1, header=TRUE)
	}
	likelihoods_by_fleet[likelihoods_by_fleet == "_"] <- NA
	likelihoods_by_fleet <- type.convert(likelihoods_by_fleet, as.is=TRUE)
	names(likelihoods_by_fleet) <- c("Label", "ALL", FleetNames)
	labs <- likelihoods_by_fleet[["Label"]]
	for (irow in seq_along(labs)) {
		labs[irow] <- substr(labs[irow], 1, nchar(labs[irow]) - 1)
	}
	likelihoods_by_fleet[["Label"]] <- labs
	stats[["likelihoods_by_fleet"]] <- likelihoods_by_fleet
	likelihoods_by_tag_group <- match_report_table("Tag_Group:", 0, header=TRUE)
	if (!is.null(likelihoods_by_tag_group)) {
		likelihoods_by_tag_group[likelihoods_by_tag_group == "_"] <- NA
		likelihoods_by_tag_group <- type.convert(likelihoods_by_tag_group, as.is=TRUE)
		names(likelihoods_by_tag_group) <- c("Label", "ALL", 
			paste0("TagGroup_", names(likelihoods_by_tag_group)[-(1:2)]))
		likelihoods_by_tag_group[["Label"]][1] <- "Tag_Group"
		stats[["likelihoods_by_tag_group"]] <- likelihoods_by_tag_group
	}
	Parm_devs_detail <- match_report_table("Parm_devs_detail", 
		1, header=TRUE, type.convert=TRUE)
	stats[["Parm_devs_detail"]] <- Parm_devs_detail
	parameters <- match_report_table("PARAMETERS", 1, header=TRUE)
	parameters <- df.rename(parameters, oldnames=c("PR_type", 
		"Prior_Like"), newnames=c("Pr_type", "Pr_Like"))
	parameters[parameters == "_"] <- NA
	parameters[parameters == " "] <- NA
	parameters[parameters == "1.#INF"] <- Inf
	if (SS_versionNumeric == 3.21) {
		temp <- names(parameters)
		message("Inserting new 13th column heading in parameters section", 
			"due to error in Report.sso in SSv3.21f")
		temp <- c(temp[1:12], "PR_type_code", temp[-(1:12)])
		temp <- temp[-length(temp)]
		names(parameters) <- temp
	}
	if ("Gradient" %in% names(parameters) && any(parameters[["Gradient"]] %in% c("dev", "F"))) {
		bad <- parameters[["Gradient"]] %in% c("dev", "F")
		parameters[["Pr_type"]][bad] <- parameters[["Gradient"]][bad]
		parameters[["Gradient"]][bad] <- NA
	}
	parameters <- type.convert(parameters, as.is=TRUE)
	if (SS_versionNumeric < 3.21) {
		parameters[["Pr_type_numeric"]] <- parameters[["Pr_type"]]
		parameters[["Pr_type"]][parameters[["Pr_type_numeric"]] == -1] <- "No_prior"
		parameters[["Pr_type"]][parameters[["Pr_type_numeric"]] ==  0] <- "Normal"
		parameters[["Pr_type"]][parameters[["Pr_type_numeric"]] ==  1] <- "Sym_Beta"
		parameters[["Pr_type"]][parameters[["Pr_type_numeric"]] ==  2] <- "Full_Beta"
		parameters[["Pr_type"]][parameters[["Pr_type_numeric"]] ==  3] <- "Log_Norm"
		parameters[["Pr_type"]][parameters[["Pr_type_numeric"]] ==  4] <- "Log_Norm_adjusted"
	}
	ParmLabels <- parameters[["Label"]]
	ParmLabels[duplicated(ParmLabels)] <- paste0(ParmLabels[duplicated(ParmLabels)], "_2")
	rownames(parameters) <- ParmLabels
	if (!is.na(parfile)) {
		parline <- read.table(file.path(dir, parfile), fill=TRUE, comment.char="", nrows=1)
	}
	else {
		parline <- matrix(NA, 1, 16)
	}
	stats[["N_estimated_parameters"]] <- parline[1, 6]
	pars <- parameters[!is.na(parameters[["Active_Cnt"]]), ]
	if (nrow(pars) > 0) {
		pars[["Afterbound"]] <- ""
		pars[["checkdiff"]]  <- pars[["Value"]] - pars[["Min"]]
		pars[["checkdiff2"]] <- pars[["Max"]] - pars[["Value"]]
		pars[["checkdiff3"]] <- abs(pars[["Value"]] - (pars[["Max"]] - 
			(pars[["Max"]] - pars[["Min"]])/2))
		pars[["Afterbound"]][pars[["checkdiff"]] < 0.001 | pars[["checkdiff2"]] < 
			0.001 | pars[["checkdiff2"]] < 0.001] <- "CHECK"
		pars[["Afterbound"]][!pars[["Afterbound"]] %in% "CHECK"] <- "OK"
	}
	stats[["table_of_phases"]] <- table(parameters[["Phase"]])
	estimated_non_dev_parameters <- pars[, names(pars) %in% c("Value", 
		"Phase", "Min", "Max", "Init", "Prior", "Gradient", "Pr_type", 
		"Pr_SD", "Pr_Like", "Parm_StDev", "Status", "Afterbound")]
	devnames <- c("RecrDev", "InitAge", "ForeRecr", "DEVadd", 
		"DEVmult", "DEVrwalk", "DEV_MR_rwalk", "ARDEV")
	devrows <- NULL
	for (iname in seq_along(devnames)) {
		devrows <- unique(c(devrows, grep(devnames[iname], rownames(estimated_non_dev_parameters))))
	}
	if (!is.null(devrows) & length(devrows) > 0) {
		estimated_non_dev_parameters <- estimated_non_dev_parameters[-devrows, ]
	}
	stats[["estimated_non_dev_parameters"]] <- estimated_non_dev_parameters
	seldev_pars <- parameters[grep("ARDEV", parameters[["Label"]], 
		fixed=TRUE), names(parameters) %in% c("Label", "Value")]
	if (nrow(seldev_pars) == 0) {
		seldev_pars <- NULL
		seldev_matrix <- NULL
	}
	else {
		if (any(duplicated(FleetNames))) {
			warning("Duplicated fleet names will cause only the semi-parametric", 
				" selectivity to be available for the first of the duplicates.")
		}
		seldev_label_info <- strsplit(seldev_pars[["Label"]], split="_")
		seldev_label_info <- data.frame(do.call(rbind, lapply(seldev_label_info, rbind)))
		seldev_pars[["Fleet"]] <- seldev_label_info[["X1"]]
		yr_col <- grep("^y\\d\\d\\d\\d$", seldev_label_info[1, ])
		type_bin_col <- grep("^[aAlL][[:alpha:]]{0,3}\\d$", seldev_label_info[1, ])
		seldev_pars[["Year"]] <- as.numeric(substring(seldev_label_info[[yr_col]], 2))
		seldev_pars[["Type"]] <- ifelse(substring(seldev_label_info[[type_bin_col]], 
			1, 1) %in% c("A", "a"), yes="age", no="length")
		first_bin_digit <- ifelse(seldev_pars[["Type"]] == "age", 2, 5)
		seldev_pars[["Bin"]] <- as.numeric(substring(seldev_label_info[[type_bin_col]], first_bin_digit))
		seldev_pars <- seldev_pars[, -1]
		seldev_matrix <- list()
		for (fleet in sort(unique(seldev_pars[["Fleet"]]))) {
			seldev_pars_f <- seldev_pars[seldev_pars[["Fleet"]] == fleet, ]
			for (type in unique(seldev_pars_f[["Type"]])) {
				seldev_pars_sub <- seldev_pars_f[seldev_pars_f[["Type"]] == type, ]
				seldev_label <- paste0(fleet, "_", type, "_seldevs")
				seldev_yrs <- sort(unique(seldev_pars_sub[["Year"]]))
				seldev_bins <- sort(unique(seldev_pars_sub[["Bin"]]))
				if (type == "length") {
					seldev_matrix[[seldev_label]] <- matrix(nrow=length(seldev_yrs), 
						ncol=length(seldev_bins), dimnames=list(Year=seldev_yrs, Lbin=seldev_bins))
				}
				if (type == "age") {
					seldev_matrix[[seldev_label]] <- matrix(nrow=length(seldev_yrs), 
						ncol=length(seldev_bins), dimnames=list(Year=seldev_yrs, Age=seldev_bins))
				}
				for (y in seldev_yrs) {
					for (bin in seldev_bins) {
						seldev_matrix[[seldev_label]][paste(y), paste(bin)] <- 
							seldev_pars_sub[["Value"]][seldev_pars_sub[["Year"]] == y & 
							seldev_pars_sub[["Bin"]] == bin][1]
					}
				}
			}
		}
	}
	DM_pars <- parameters[grep("ln\\((EffN_mult)|(DM_theta)\\)", 
		parameters[["Label"]]), names(parameters) %in% c("Value", "Phase", "Min", "Max")]
	DM_pars[["Theta"]] <- exp(DM_pars[["Value"]])
	DM_pars$"Theta/(1+Theta)" <- DM_pars[["Theta"]]/(1 + DM_pars[["Theta"]])
	if (covar) {
		covarfile <- file.path(dir, covarfile)
		if (!file.exists(covarfile)) {
			message("covar file not found, input 'covar' changed to FALSE")
			covar <- FALSE
		}
		else {
			covarhead <- readLines(con=covarfile, n=10)
			covarskip <- grep("active-i", covarhead) - 1
			covartime <- findtime(covarhead)
			if (is.null(covartime) || is.null(repfiletime)) {
				message("problem comparing the file creation times:\n", 
				  "  Report.sso:", repfiletime, "\n", "  covar.sso:", covartime)
			}
			else {
				if (covartime != repfiletime) {
					message("covar time:", covartime)
					stop(shortrepfile, " and ", covarfile, " were from different model runs. Change input to covar=FALSE")
				}
			}
			nowrite <- grep("do not write", covarhead)
			if (length(nowrite) > 0) {
				warning("covar file contains the warning\n", 
					" '", covarhead[nowrite], "'\n", "  input 'covar' changed to FALSE.\n")
				covar <- FALSE
			}
		}
	}
	if (covar) {
		CoVar <- read.table(covarfile, header=TRUE, colClasses=c(rep("numeric", 
			4), rep("character", 4), "numeric"), skip=covarskip)
		if (verbose) {
			message("Got covar file.")
		}
		stdtable <- CoVar[CoVar[["Par..j"]] == "Std", c(7, 9, 5)]
		names(stdtable) <- c("name", "std", "type")
		N_estimated_parameters2 <- sum(stdtable[["type"]] == "Par")
		if (is.na(stats[["N_estimated_parameters"]])) {
			stats[["N_estimated_parameters"]] <- N_estimated_parameters2
		}
		else {
			if (stats[["N_estimated_parameters"]] != N_estimated_parameters2) {
				warning(stats[["N_estimated_parameters"]], " estimated parameters indicated by the par file\n  ", 
					N_estimated_parameters2, " estimated parameters shown in the covar file\n  ", 
					"Returning the par file value: ", stats[["N_estimated_parameters"]])
			}
		}
		if (any(is.na(stdtable[["std"]]))) {
			warning("NA value for parameter uncertainty found in ", 
				sum(is.na(stdtable[["std"]])), " rows of covar.sso file. ", 
				"First par with NA: ", stdtable[["name"]][is.na(stdtable[["std"]])])
		}
		Nstd <- sum(stdtable[["std"]] > 0, na.rm=TRUE)
		checkbadrun <- unique(stdtable[["std"]])
		if (length(checkbadrun) == 1) {
			if (checkbadrun %in% c(NA, "NaN", "na")) {
				stop(paste0("No quantities were estimated in the covar file \nand all", 
					"estimates of standard deviation are ", checkbadrun, 
					". \nTry re-running", "stock synthesis."))
			}
		}
		if (Nstd <= 1) {
			stop("Too few estimated quantities in covar file (n=", 
				Nstd, "). Change input to covar=FALSE.")
		}
	}
	else {
		if (verbose) {
			message("You skipped the covar file")
		}
	}
	flush.console()
	wtatage <- NULL
	if (readwt) {
		wtfile <- file.path(dir, wtfile)
		wtatage <- r4ss:::SS_readwtatage(file=wtfile, verbose=verbose)
	}
	if (is.null(dir.mcmc)) {
		mcmc <- NULL
	}
	else {
		dir.mcmc.full <- NULL
		if (dir.exists(dir.mcmc)) {
			dir.mcmc.full <- dir.mcmc
		}
		if (dir.exists(file.path(dir, dir.mcmc))) {
			dir.mcmc.full <- file.path(dir, dir.mcmc)
		}
		if (is.null(dir.mcmc.full)) {
			warning("'dir.mcmc' directory not found either as an absolute path ", 
				"or relative to the 'dir' input")
			mcmc <- NULL
		}
		else {
			if ("posteriors.sso" %in% dir(dir.mcmc.full)) {
				if (verbose) {
					message("Running 'SSgetMCMC' to get MCMC output")
				}
				mcmc <- SSgetMCMC(dir=dir.mcmc.full)
			}
			else {
				warning("skipping reading MCMC output because posterior.sso file", 
				  " not found in \n", dir.mcmc.full)
				mcmc <- NULL
			}
		}
	}
	der <- match_report_table("DERIVED_QUANTITIES", 4, header=TRUE)
	der <- df.rename(der, oldnames="LABEL", newnames="Label")
	der <- der[der[["Label"]] != "Bzero_again", ]
	der[der == "_"] <- NA
	der[der == ""]  <- NA
	test <- grep("Parm_dev_details", der[["Label"]])
	if (length(test) > 0) {
		der <- der[1:(min(test) - 1), ]
	}
	der <- type.convert(der, as.is=TRUE)
	der[["Label"]] <- gsub("SPB_", "SSB_", der[["Label"]], fixed=TRUE)
	rownames(der)[!duplicated(der[["Label"]])] <- der[["Label"]][!duplicated(der[["Label"]])]
	managementratiolabels <- match_report_table("DERIVED_QUANTITIES", 
		1, "DERIVED_QUANTITIES", 3, cols=1:2)
	names(managementratiolabels) <- c("Ratio", "Label")
	forecast_selectivity <- grep("forecast_selectivity", rawrep[, 1], value=TRUE)
	if (length(forecast_selectivity) == 0) {
		forecast_selectivity <- NA
		offset <- -1
	}
	else {
		offset <- -2
	}
	MGparmAdj <- match_report_table("MGparm_By_Year_after_adjustments", 
		1, header=TRUE, type.convert=TRUE)
	MGparmAdj <- df.rename(MGparmAdj, oldnames="Year", newnames="Yr")
	SelSizeAdj <- match_report_table("selparm(Size)_By_Year_after_adjustments", 2)
	if (is.null(SelSizeAdj) || nrow(SelSizeAdj) <= 2) {
		SelSizeAdj <- NULL
	}
	else {
		SelSizeAdj <- SelSizeAdj[, apply(SelSizeAdj, 2, emptytest) < 1]
		SelSizeAdj[SelSizeAdj == ""] <- NA
		SelSizeAdj <- type.convert(SelSizeAdj, as.is=TRUE)
		if (rawrep[match_report_line("selparm(Size)_By_Year_after_adjustments") + 1, 3] == "Change?") {
			names(SelSizeAdj) <- c("Fleet", "Yr", "Change?", paste0("Par", 1:(ncol(SelSizeAdj) - 3)))
		}
		else {
			names(SelSizeAdj) <- c("Fleet", "Yr", paste0("Par", 1:(ncol(SelSizeAdj) - 2)))
		}
	}
	SelAgeAdj <- match_report_table("selparm(Age)_By_Year_after_adjustments", 2)
	if (!is.null(SelAgeAdj) && nrow(SelAgeAdj) > 2) {
		SelAgeAdj <- SelAgeAdj[, apply(SelAgeAdj, 2, emptytest) < 1]
		SelAgeAdj[SelAgeAdj == ""] <- NA
		if (SelAgeAdj[1, 1] == "RECRUITMENT_DIST") {
			SelAgeAdj <- NA
		}
		else {
			SelAgeAdj <- type.convert(SelAgeAdj, as.is=TRUE)
			names(SelAgeAdj) <- c("Flt", "Yr", paste0("Par", 1:(ncol(SelAgeAdj) - 2)))
			if (rawrep[match_report_line("selparm(Age)_By_Year_after_adjustments") + 1, 3] == "Change?") {
				names(SelAgeAdj) <- c("Fleet", "Yr", "Change?", paste0("Par", 1:(ncol(SelAgeAdj) - 3)))
			}
			else {
				names(SelAgeAdj) <- c("Fleet", "Yr", paste0("Par", 1:(ncol(SelAgeAdj) - 2)))
			}
		}
	}
	else {
		SelAgeAdj <- NULL
	}
	recruitment_dist <- match_report_table("RECRUITMENT_DIST", 1, header=TRUE, type.convert=TRUE)
	if (!is.null(recruitment_dist)) {
		if ("Frac/sex" %in% names(recruitment_dist)) {
			first_seas_with_recruits <- min(recruitment_dist[["Seas"]][recruitment_dist$"Frac/sex" > 0])
		}
		else {
			#first_seas_with_recruits <- min(recruitment_dist[["Seas"]][recruitment_dist[["Value"]] >  0])  ## buggy
			first_seas_with_recruits <- min(recruitment_dist[["Seas"]][recruitment_dist[["recr_dist_F"]] >  0])  ## (RH 241217)
		}
		recruit_dist_Bmark <- match_report_table("RECRUITMENT_DIST_B", 1, header=TRUE, type.convert=TRUE)
		if (!is.null(recruit_dist_Bmark)) {
			if (SS_versionNumeric < 3.3) {
				recruit_dist_endyr <- match_report_table("RECRUITMENT_DIST_FORECAST", 1, header=TRUE, type.convert=TRUE)
			}
			else {
				recruit_dist_endyr <- match_report_table("RECRUITMENT_DIST_endyr", 1, header=TRUE, type.convert=TRUE)
				if (length(grep("RECRUITMENT_DIST_TIMESERIES", recruit_dist_endyr[["Settle#"]])) == 1) {
					tmp_brk_line <- grep("RECRUITMENT_DIST_TIMESERIES", recruit_dist_endyr[["Settle#"]]) - 1
					recruit_dist_endyr <- recruit_dist_endyr[seq_len(tmp_brk_line), ]
				}
			}
			recruitment_dist <- list(recruit_dist=recruitment_dist, 
				recruit_dist_Bmark=recruit_dist_Bmark, recruit_dist_endyr=recruit_dist_endyr)
		}
	}
	stats[["maximum_gradient_component"]] <- 
		as.numeric(match_report_table("Convergence_Level", 0, "Convergence_Level", 0, cols=2))
	if ("Gradient" %in% names(parameters)) {
		if (any(!is.na(parameters[["Gradient"]]))) {
			ngrads <- min(5, max(parameters[["Active_Cnt"]], na.rm=TRUE))
			stats[["parameters_with_highest_gradients"]] <- 
				head(parameters[order(abs(parameters[["Gradient"]]), decreasing=TRUE), c("Value", "Gradient")], n=5)
		}
	}
	if (SS_versionNumeric >= 3.3 | substring(SS_version, 1, 9) %in% paste0("SS-V3.24", LETTERS[21:26]) |
		substring(SS_version, 1, 10) %in% paste0("SS-V3.24A", LETTERS)) {
		last_row_index <- 11
	}
	else {
		last_row_index <- 10
	}
	srhead <- match_report_table("SPAWN_RECRUIT", 0, "SPAWN_RECRUIT", last_row_index, cols=1:6)
	if (all(srhead[7, ] == "")) {
		last_row_index <- 12
		srhead <- match_report_table("SPAWN_RECRUIT", 0, "SPAWN_RECRUIT", last_row_index, cols=1:6)
	}
	if (is.null(srhead)) {
		rmse_table <- NULL
		breakpoints_for_bias_adjustment_ramp <- NULL
		sigma_R_in <- parameters["SR_sigmaR", "Value"]
	}
	else {
		rmse_table <- as.data.frame(srhead[-(1:(last_row_index - 1)), 1:5])
		rmse_table <- rmse_table[!grepl("SpawnBio", rmse_table[, 2]), ]
		rmse_table <- type.convert(rmse_table, as.is=TRUE)
		names(rmse_table) <- srhead[last_row_index - 1, 1:5]
		names(rmse_table)[4] <- "RMSE_over_sigmaR"
		row.names(rmse_table) <- NULL
		sigma_R_in <- as.numeric(srhead[grep("sigmaR", srhead[, 2]), 1])
		if (any(srhead[1, ] == "RecDev_method:")) {
			RecDev_method <- as.numeric(srhead[1, which(srhead[1,] == "RecDev_method:") + 1])
		}
		else {
			RecDev_method <- NULL
		}
		biascol <- grep("breakpoints_for_bias", srhead)
		breakpoints_for_bias_adjustment_ramp <- srhead[grep("breakpoints_for_bias", srhead[, biascol]), 1:5]
		colnames(breakpoints_for_bias_adjustment_ramp) <- c("last_yr_early", 
			"first_yr_full", "last_yr_full", "first_yr_recent", "max_bias_adj")
		rownames(breakpoints_for_bias_adjustment_ramp) <- NULL
	}
	raw_recruit <- match_report_table("SPAWN_RECRUIT", last_row_index + 1)
	if (!is.null(raw_recruit) && raw_recruit[1, 1] == "S/Rcurve") {
		raw_recruit <- match_report_table("SPAWN_RECRUIT", last_row_index)
	}
	if (!is.null(raw_recruit) && nrow(raw_recruit) < length(startyr:endyr)) {
		raw_recruit <- match_report_table("SPAWN_RECRUIT", last_row_index + 1, which_blank=2)
		if (raw_recruit[1, 1] == "S/Rcurve") {
			raw_recruit <- match_report_table("SPAWN_RECRUIT", last_row_index, which_blank=2)
		}
	}
	if (is.null(raw_recruit)) {
		recruit <- NULL
	}
	else {
		names(raw_recruit) <- raw_recruit[1, ]
		raw_recruit[raw_recruit == "_"] <- NA
		raw_recruit <- raw_recruit[-(1:2), ]
		recruit <- raw_recruit[-(1:2), ]
		recruit[["dev"]][recruit[["dev"]] == "-nan(ind)"] <- NA
		recruit <- type.convert(recruit, as.is=TRUE)
		recruit <- df.rename(recruit, oldnames=c("year", "spawn_bio", 
			"adjusted", "biasadj"), newnames=c("Yr", "SpawnBio", "bias_adjusted", "biasadjuster"))
	}
	SPAWN_RECR_CURVE <- NULL
	if (!is.na(match_report_line("Full_Spawn_Recr_Curve"))) {
		SPAWN_RECR_CURVE <- match_report_table("Full_Spawn_Recr_Curve", 1, header=TRUE, type.convert=TRUE)
	}
	if (!is.na(match_report_line("SPAWN_RECR_CURVE"))) {
		SPAWN_RECR_CURVE <- match_report_table("SPAWN_RECR_CURVE", 1, header=TRUE, type.convert=TRUE)
	}
	if (SS_versionNumeric >= 3.3) {
		fit_len_comps <- match_report_table("FIT_LEN_COMPS", 1, header=TRUE)
	}
	else {
		fit_len_comps <- NULL
	}
	if (!is.null(dim(fit_len_comps)) && nrow(fit_len_comps) > 0) {
		fit_len_comps[fit_len_comps == "_"] <- NA
		fit_len_comps <- type.convert(fit_len_comps, as.is=TRUE)
	}
	else {
		fit_len_comps <- NULL
	}
	if (SS_versionNumeric < 3.3) {
		lenntune <- match_report_table("FIT_AGE_COMPS", -(nfleets + 
			2), "FIT_AGE_COMPS", -1, cols=1:10, header=TRUE)
		names(lenntune)[10] <- "FleetName"
		lenntune[lenntune == "_"] <- NA
		lenntune <- lenntune[lenntune[["N"]] > 0, c(10, 1, 4:9)]
		lenntune$"MeaneffN/MeaninputN"[lenntune$"MeaneffN/MeaninputN" == "-1.#IND"] <- NA
		lenntune <- type.convert(lenntune, as.is=TRUE)
		lenntune$"HarMean/MeanInputN" <- lenntune$"HarMean(effN)"/lenntune$"mean(inputN*Adj)"
	}
	else {
		lenntune <- match_report_table("Length_Comp_Fit_Summary", 1, header=TRUE)
		if (!is.null(lenntune)) {
			lenntune <- df.rename(lenntune, oldnames=c("FleetName", 
				"Factor", "HarMean_effN"), newnames=c("Fleet_name", "Data_type", "HarMean"))
			if ("Data_type" %in% names(lenntune)) {
				lenntune <- type.convert(lenntune, as.is=TRUE)
			}
			else {
				lenntune <- lenntune[lenntune[["Nsamp_adj"]] > 0, ]
				lenntune <- type.convert(lenntune, as.is=TRUE)
				lenntune$"HarMean(effN)/mean(inputN*Adj)" <- lenntune$HarMean/lenntune$"mean_inputN*Adj"
				lenntune <- df.rename(lenntune, oldnames=c("HarMean", 
					"mean_inputN*Adj"), newnames=c("HarMean(effN)", "mean(inputN*Adj)"))
				lenntune <- lenntune[, names(lenntune) != "mean_effN"]
				end.names <- c("Recommend_Var_Adj", "Fleet_name")
				lenntune <- lenntune[, c(which(!names(lenntune) %in% end.names), which(names(lenntune) %in% end.names))]
			}
		}
	}
	stats[["Length_Comp_Fit_Summary"]] <- lenntune
	fit_age_comps <- match_report_table("FIT_AGE_COMPS", 1, header=TRUE)
	if (!is.null(dim(fit_age_comps)) && nrow(fit_age_comps) > 0) {
		fit_age_comps[fit_age_comps == "_"] <- NA
		fit_age_comps <- type.convert(fit_age_comps, as.is=TRUE)
	}
	else {
		fit_age_comps <- NULL
	}
	if (SS_versionNumeric < 3.3) {
		agentune <- match_report_table("FIT_SIZE_COMPS", -(nfleets + 2), "FIT_SIZE_COMPS", -2, cols=1:10, header=TRUE)
	}
	else {
		start <- match_report_line("Age_Comp_Fit_Summary")
		if (is.na(start)) {
			agentune <- NULL
		}
		else {
			if (rawrep[start + 1, 1] == "") {
				adjust1 <- 2
				which_blank <- 2
			}
			else {
				adjust1 <- 1
				which_blank <- 1
			}
			agentune <- match_report_table("Age_Comp_Fit_Summary", 
				adjust1=adjust1, header=TRUE, which_blank=which_blank)
		}
	}
	agentune <- df.rename(agentune, oldnames=c("FleetName", 
		"N", "Factor", "HarMean_effN"), newnames=c("Fleet_name", "Nsamp_adj", "Data_type", "HarMean"))
	if ("Data_type" %in% names(agentune)) {
		agentune <- type.convert(agentune, as.is=TRUE)
	}
	else {
		if (!is.null(dim(agentune))) {
			names(agentune)[ncol(agentune)] <- "Fleet_name"
			agentune[agentune == "_"] <- NA
			agentune <- agentune[!is.na(agentune[["Nsamp_adj"]]) & agentune[["Nsamp_adj"]] > 0, ]
			agentune$"MeaneffN/MeaninputN"[agentune$"MeaneffN/MeaninputN" == "-1.#IND"] <- NA
			agentune <- type.convert(agentune, as.is=TRUE)
			agentune$"HarMean(effN)/mean(inputN*Adj)" <- agentune$"HarMean(effN)"/agentune$"mean(inputN*Adj)"
			agentune[["Recommend_Var_Adj"]] <- agentune[["Var_Adj"]] * agentune$"HarMean(effN)/mean(inputN*Adj)"
			badnames <- c("mean_effN", "Mean(effN/inputN)", "MeaneffN/MeaninputN")
			agentune <- agentune[, !names(agentune) %in% badnames]
			agentune <- agentune[, c(which(names(agentune) != "Fleet_name"), which(names(agentune) == "Fleet_name"))]
			agentune <- df.rename(agentune, oldnames=c("Var_Adj"), newnames=c("Curr_Var_Adj"))
		}
		else {
			agentune <- NULL
		}
	}
	stats[["Age_Comp_Fit_Summary"]] <- agentune
	fit_size_comps <- NULL
	if (SS_versionNumeric >= 3.3) {
		if (!is.na(match_report_line("FIT_SIZE_COMPS"))) {
			fit_size_comps <- match_report_table("FIT_SIZE_COMPS", 1, header=FALSE, blank_lines=rep_blank_lines)
			if (!is.null(dim(fit_size_comps)) && nrow(fit_size_comps) > 0 && fit_size_comps[1, 1] != "#_none") {
				names(fit_size_comps) <- fit_size_comps[2, ]
				fit_size_comps[["Method"]] <- NA
				fit_size_comps[["Units"]]  <- NA
				fit_size_comps[["Scale"]]  <- NA
				fit_size_comps[["Add_to_comp"]] <- NA
				method_lines <- grep("#Method:", fit_size_comps[, 1])
				method_info <- fit_size_comps[method_lines, ]
				if (any(grepl("Size_Comp_Fit_Summary", fit_size_comps[, 1]))) {
					tune_lines <- grep("Size_Comp_Fit_Summary", fit_size_comps[, 1]) + 1
				}
				else {
					tune_lines <- grep("Factor", fit_size_comps[, 1])
				}
				sizentune <- NULL
				for (imethod in seq_along(method_lines)) {
					start <- method_lines[imethod]
					if (imethod != length(method_lines)) {
						end <- method_lines[imethod + 1] - 1
					}
					else {
						end <- nrow(fit_size_comps)
					}
					fit_size_comps[["Method"]][start:end] <- method_info[imethod, 2]
					fit_size_comps[["Units"]][start:end] <- method_info[imethod, 4]
					fit_size_comps[["Scale"]][start:end] <- method_info[imethod, 6]
					fit_size_comps[["Add_to_comp"]][start:end] <- method_info[imethod, 8]
					sizentune <- rbind(sizentune, fit_size_comps[tune_lines[imethod]:end, ])
				}
				goodcols <- c(1:grep("name", tolower(sizentune[1, ])), grep("Method", names(sizentune)))
				sizentune[1, max(goodcols)] <- "Method"
				sizentune <- sizentune[, goodcols]
				names(sizentune) <- sizentune[1, ]
				sizentune <- df.rename(sizentune, oldnames=c("Factor", "HarMean_effN"), newnames=c("Data_type", "HarMean"))
				sizentune <- sizentune[nchar(sizentune[["Data_type"]]) == 1, ]
				sizentune <- type.convert(sizentune, as.is=TRUE)
				stats[["Size_Comp_Fit_Summary"]] <- sizentune
				fit_size_comps <- dplyr::filter(fit_size_comps, Fleet_Name %in% FleetNames & Fleet %in% 1:nfleets)
			}
		}
		else {
			fit_size_comps <- match_report_table("FIT_SIZE_COMPS", 1, "Size_Comp_Fit_Summary", -(nfleets + 2), header=TRUE)
		}
	}
	if (!is.null(dim(fit_size_comps)) && nrow(fit_size_comps) > 0) {
		fit_size_comps[fit_size_comps == "_"] <- NA
		fit_size_comps <- type.convert(fit_size_comps, as.is=TRUE)
	}
	if (SS_versionNumeric >= 3.3) {
		if (!exists("sizentune")) {
			sizentune <- match_report_table("Size_Comp_Fit_Summary", 1, "OVERALL_COMPS", -1, cols=1:10, header=TRUE)
			if (!is.null(dim(sizentune))) {
				sizentune[, 1] <- sizentune[, 10]
				sizentune <- sizentune[sizentune[["Npos"]] > 0, c(1, 3, 4, 5, 6, 8, 9)]
			}
			else {
				sizentune <- NULL
			}
		}
		stats[["Size_comp_Eff_N_tuning_check"]] <- sizentune
	}
	age_data_info <- NULL
	len_data_info <- NULL
	if (nrow(DM_pars) > 0) {
		if (!is.null(Length_comp_error_controls) | !is.null(Age_comp_error_controls)) {
			if (comp) {
				if (nrow(lendbase) > 0) {
					fit_len_comps_select <- 
						dplyr::select(dplyr::rename(fit_len_comps, Like_sum=Like), Fleet, Time, Sexes, Part, Nsamp_DM)
					lendbase <- dplyr::left_join(lendbase, fit_len_comps_select)
				}
				if (nrow(agedbase) > 0) {
					fit_age_comps_select <- 
						dplyr::select(dplyr::rename(fit_age_comps, Like_sum=Like), Fleet, Time, Sexes, Part, Nsamp_DM)
					agedbase <- dplyr::left_join(agedbase, fit_age_comps_select)
				}
				if (nrow(condbase) > 0) {
					fit_cond_age_select <- 
						dplyr::select(dplyr::rename(fit_age_comps, Like_sum=Like), Fleet, Time, Sexes, Part, Nsamp_DM)
					condbase <- dplyr::left_join(condbase, fit_cond_age_select)
				}
				if (nrow(sizedbase) > 0) {
					fit_size_comps_select <- 
						dplyr::select(dplyr::rename(dplyr::rename(fit_size_comps, Like_sum=Like), method=Method), Fleet, Time, Sexes, Part, Nsamp_DM, method)
					sizedbase <- dplyr::left_join(sizedbase, fit_size_comps_select)
				}
			}
		}
		else {
			if (verbose) {
				message("Reading data.ss_new (or data_echo.ss_new) for info on Dirichlet-Multinomial parameters")
			}
			datname <- get_dat_new_name(dir)
			datfile <- SS_readdat(file=file.path(dir, datname), 
				verbose=verbose, )
			if (is.null(datfile)) {
				starter <- r4ss:::SS_readstarter(file=file.path(dir, "starter.ss"), verbose=verbose)
				datfile <- SS_readdat(file=file.path(dir, starter[["datfile"]]), verbose=verbose, version="3.30")
			}
			age_data_info <- datfile[["age_info"]]
			len_data_info <- datfile[["len_info"]]
			if (!is.null(age_data_info) & !is.null(len_data_info)) {
				age_data_info[["CompError"]] <- as.numeric(age_data_info[["CompError"]])
				age_data_info[["ParmSelect"]] <- as.numeric(age_data_info[["ParmSelect"]])
				len_data_info[["CompError"]] <- as.numeric(len_data_info[["CompError"]])
				len_data_info[["ParmSelect"]] <- as.numeric(len_data_info[["ParmSelect"]])
				if (!any(age_data_info[["CompError"]] > 0) & !any(len_data_info[["CompError"]] > 0)) {
					stop("Problem with Dirichlet-Multinomial parameters: \n", 
						"  Report file indicates parameters exist, but no CompError values\n", 
						"  in data.ss_new are > 0.")
				}
			}
			get_DM_sample_size <- function(CompError, f, sub, data_info, dbase) {
				ipar <- data_info[["ParmSelect"]][f]
				if (ipar %in% 1:nrow(DM_pars)) {
					if (CompError == 1) {
						Theta <- DM_pars[["Theta"]][ipar]
					}
					if (CompError == 2) {
						beta <- DM_pars[["Theta"]][ipar]
					}
				}
				else {
					stop("Issue with Dirichlet-Multinomial parameter:", "Fleet=", f, "and ParmSelect=", ipar)
				}
				if (CompError == 1) {
					Nsamp_DM <- 1/(1 + Theta) + dbase[["Nsamp_adj"]][sub] * Theta/(1 + Theta)
				}
				if (CompError == 2) {
					Nsamp_DM <- dbase[["Nsamp_adj"]][sub] * (1 + beta)/(dbase[["Nsamp_adj"]][sub] + beta)
				}
				Nsamp_DM
			}
			if (comp) {
				if (nrow(agedbase) > 0) {
					agedbase[["Nsamp_DM"]] <- NA
				}
				if (nrow(lendbase) > 0) {
					lendbase[["Nsamp_DM"]] <- NA
				}
				if (nrow(condbase) > 0) {
					condbase[["Nsamp_DM"]] <- NA
				}
				for (f in unique(agedbase[["Fleet"]])) {
					if (age_data_info[["CompError"]][f] > 0) {
						sub <- agedbase[["Fleet"]] == f
						agedbase[["Nsamp_DM"]][sub] <- get_DM_sample_size(CompError=age_data_info[["CompError"]][f], 
							f=f, sub=sub, data_info=age_data_info, dbase=agedbase)
					}
				}
				for (f in unique(lendbase[["Fleet"]])) {
					if (len_data_info[["CompError"]][f] > 0) {
						sub <- lendbase[["Fleet"]] == f
						lendbase[["Nsamp_DM"]][sub] <- get_DM_sample_size(CompError=len_data_info[["CompError"]][f], 
							f=f, sub=sub, data_info=len_data_info, dbase=lendbase)
					}
				}
				for (f in unique(condbase[["Fleet"]])) {
					if (age_data_info[["CompError"]][f] > 0) {
						sub <- condbase[["Fleet"]] == f
						condbase[["Nsamp_DM"]][sub] <- get_DM_sample_size(CompError=age_data_info[["CompError"]][f], 
							f=f, sub=sub, data_info=age_data_info, dbase=condbase)
					}
				}
			}
		}
	}
	jitter_info <- parameters[!is.na(parameters[["Active_Cnt"]]) & 
		!is.na(parameters[["Min"]]), c("Value", "Min", "Max", "Init")]
	jitter_info[["sigma"]] <- (jitter_info[["Max"]] - jitter_info[["Min"]])/(2 * qnorm(0.999))
	jitter_info[["CV"]] <- jitter_info[["sigma"]]/jitter_info[["Init"]]
	jitter_info[["InitLocation"]] <- pnorm(q=jitter_info[["Init"]], 
		mean=(jitter_info[["Max"]] + jitter_info[["Min"]])/2, sd=jitter_info[["sigma"]])
	if (verbose) {
		message("Finished primary run statistics list")
	}
	flush.console()
	if (SS_versionNumeric <= 3.24) {
		returndat[["definitions"]] <- fleetdefs
		returndat[["fleet_ID"]] <- fleet_ID
		returndat[["fleet_area"]] <- fleet_area
		returndat[["catch_units"]] <- catch_units
		returndat[["catch_error"]] <- catch_error
	}
	if (SS_versionNumeric >= 3.3) {
		returndat[["definitions"]] <- fleetdefs
		returndat[["fleet_ID"]] <- fleet_ID
		returndat[["fleet_type"]] <- fleet_type
		returndat[["fleet_timing"]] <- fleet_timing
		returndat[["fleet_area"]] <- fleet_area
		returndat[["catch_units"]] <- catch_units
		if (exists("catch_se")) {
			returndat[["catch_se"]] <- catch_se
			returndat[["equ_catch_se"]] <- equ_catch_se
		}
		else {
			returndat[["catch_se"]] <- NA
			returndat[["equ_catch_se"]] <- NA
		}
	}
	return.def <- function(x) {
		if (exists(x)) {
			get(x)
		}
		else {
			NULL
		}
	}
	returndat[["mcmc"]] <- mcmc
	returndat[["survey_units"]] <- survey_units
	returndat[["survey_error"]] <- survey_error
	returndat[["IsFishFleet"]] <- IsFishFleet
	returndat[["nfishfleets"]] <- nfishfleets
	returndat[["nfleets"]] <- nfleets
	returndat[["nsexes"]] <- nsexes
	returndat[["ngpatterns"]] <- ngpatterns
	returndat[["lbins"]] <- lbins
	returndat[["Lbin_method"]] <- Lbin_method
	returndat[["nlbins"]] <- nlbins
	returndat[["lbinspop"]] <- lbinspop
	returndat[["nlbinspop"]] <- nlbinspop
	returndat[["sizebinlist"]] <- sizebinlist
	returndat[["age_data_info"]] <- age_data_info
	returndat[["len_data_info"]] <- len_data_info
	returndat[["agebins"]] <- agebins
	returndat[["nagebins"]] <- nagebins
	returndat[["accuage"]] <- accuage
	returndat[["nareas"]] <- nareas
	returndat[["startyr"]] <- startyr
	returndat[["endyr"]] <- endyr
	returndat[["nseasons"]] <- nseasons
	returndat[["seasfracs"]] <- seasfracs
	returndat[["seasdurations"]] <- seasdurations
	returndat[["N_sub_seasons"]] <- return.def("N_sub_seasons")
	returndat[["Spawn_month"]] <- return.def("Spawn_month")
	returndat[["Spawn_seas"]] <- return.def("Spawn_seas")
	returndat[["Spawn_timing_in_season"]] <- return.def("Spawn_timing_in_season")
	returndat[["Retro_year"]] <- return.def("Retro_year")
	returndat[["N_forecast_yrs"]] <- return.def("N_forecast_yrs")
	returndat[["Empirical_wt_at_age"]] <- return.def("Empirical_wt_at_age")
	returndat[["N_bio_patterns"]] <- return.def("N_bio_patterns")
	returndat[["N_platoons"]] <- return.def("N_platoons")
	returndat[["NatMort_option"]] <- return.def("NatMort_option")
	returndat[["GrowthModel_option"]] <- return.def("GrowthModel_option")
	returndat[["Maturity_option"]] <- return.def("Maturity_option")
	returndat[["Fecundity_option"]] <- return.def("Fecundity_option")
	returndat[["Start_from_par"]] <- return.def("Start_from_par")
	returndat[["Do_all_priors"]] <- return.def("Do_all_priors")
	returndat[["Use_softbound"]] <- return.def("Use_softbound")
	returndat[["N_nudata"]] <- return.def("N_nudata")
	returndat[["Max_phase"]] <- return.def("Max_phase")
	returndat[["Current_phase"]] <- return.def("Current_phase")
	returndat[["Jitter"]] <- return.def("Jitter")
	returndat[["ALK_tolerance"]] <- return.def("ALK_tolerance")
	returndat[["Length_comp_error_controls"]] <- Length_comp_error_controls
	returndat[["Age_comp_error_controls"]] <- Age_comp_error_controls
	returndat[["Size_comp_error_controls"]] <- Size_comp_error_controls
	returndat[["nforecastyears"]] <- nforecastyears
	returndat[["morph_indexing"]] <- morph_indexing
	returndat[["MGparmAdj"]] <- MGparmAdj
	returndat[["forecast_selectivity"]] <- forecast_selectivity
	returndat[["SelSizeAdj"]] <- SelSizeAdj
	returndat[["SelAgeAdj"]] <- SelAgeAdj
	returndat[["recruitment_dist"]] <- recruitment_dist
	returndat[["recruit"]] <- recruit
	returndat[["SPAWN_RECR_CURVE"]] <- SPAWN_RECR_CURVE
	returndat[["breakpoints_for_bias_adjustment_ramp"]] <- breakpoints_for_bias_adjustment_ramp

	biology <- match_report_table("BIOLOGY", adjust1=ifelse(custom, 2, 1), header=TRUE, type.convert=TRUE)
	biology <- df.rename(biology, oldnames=c("Low", "Mean_Size", "Wt_len", "Wt_len_F", "Mat_len", "Spawn",
		"Wt_len_M", "Fecundity"), newnames=c("Len_lo", "Len_mean", "Wt_F", "Wt_F", "Mat", "Mat*Fec", "Wt_M", "Fec"))
	FecType <- 0
	pl <- parameters[["Label"]]
	FecGrep1 <- grep("Eggs/kg_slope_wt_Fem", pl)
	FecGrep2 <- grep("Eggs_exp_len_Fem", pl)
	FecGrep3 <- grep("Eggs_exp_wt_Fem", pl)
	FecGrep4 <- grep("Eggs_slope_len_Fem", pl)
	FecGrep5 <- grep("Eggs_slope_Wt_Fem", pl)
	if (length(FecGrep1) > 0) {
		FecType <- 1
		FecPar1name <- grep("Eggs/kg_inter_Fem", pl, value=TRUE)[1]
		FecPar2name <- pl[FecGrep1[1]]
	}
	if (length(FecGrep2) > 0) {
		FecType <- 2
		FecPar1name <- grep("Eggs_scalar_Fem", pl, value=TRUE)[1]
		FecPar2name <- pl[FecGrep2[1]]
	}
	if (length(FecGrep3) > 0) {
		FecType <- 3
		FecPar1name <- grep("Eggs_scalar_Fem", pl, value=TRUE)[1]
		FecPar2name <- pl[FecGrep3[1]]
	}
	if (length(FecGrep4) > 0) {
		FecType <- 4
		FecPar1name <- grep("Eggs_intercept_Fem", pl, value=TRUE)[1]
		FecPar2name <- pl[FecGrep4[1]]
	}
	if (length(FecGrep5) > 0) {
		FecType <- 5
		FecPar1name <- grep("Eggs_intercept_Fem", pl, value=TRUE)[1]
		FecPar2name <- pl[FecGrep5[1]]
	}
	if (is.na(lbinspop[1])) {
		lbinspop <- biology[["Len_lo"]][biology[["GP"]] == 1]
	}
	if (length(returndat[["FecPar1"]]) > 1) {
		warning("Plots will only show fecundity and related quantities", "for Growth Pattern 1")
		returndat[["FecPar1"]] <- returndat[["FecPar1"]][1]
		returndat[["FecPar2"]] <- returndat[["FecPar2"]][2]
	}
	if (!is.null(biology)) {
		if (nsexes == 1 && is.na(biology[["Fec"]][1]) && "Wt_M" %in% 
			names(biology)) {
			biology[["Fec"]] <- biology[["Wt_M"]]
			biology <- biology[, !names(biology) %in% "Wt_M"]
		}
		returndat[["SpawnOutputUnits"]] <- ifelse(!is.null(biology[["Fec"]][1]) && 
			!is.na(biology[["Fec"]][1]) && any(biology[["Wt_F"]] != biology[["Fec"]]), "numbers", "biomass")
	}
	returndat[["biology"]] <- biology
	returndat[["FecType"]] <- FecType
	returndat[["FecPar1name"]] <- FecPar1name
	returndat[["FecPar2name"]] <- FecPar2name
	returndat[["FecPar1"]] <- parameters[["Value"]][parameters[["Label"]] == FecPar1name]
	returndat[["FecPar2"]] <- parameters[["Value"]][parameters[["Label"]] == FecPar2name]

	adjust1 <- ifelse(custom, 2, 1)
	M_type <- rawrep[match_report_line("Natural_Mortality") + adjust1 - 1, 2]
	M_type <- as.numeric(gsub(pattern=".*([0-9]+)", replacement="\\1", x=M_type))
	Natural_Mortality <- match_report_table("Natural_Mortality", adjust1=adjust1, header=TRUE, type.convert=TRUE)
	Natural_Mortality_Bmark <- match_report_table("Natural_Mortality_Bmark", adjust1=1, header=TRUE, type.convert=TRUE)
	Natural_Mortality_endyr <- match_report_table("Natural_Mortality_endyr", adjust1=1, header=TRUE, type.convert=TRUE)
	returndat[["M_type"]] <- M_type
	returndat[["Natural_Mortality"]] <- Natural_Mortality
	returndat[["Natural_Mortality_Bmark"]] <- Natural_Mortality_Bmark
	returndat[["Natural_Mortality_endyr"]] <- Natural_Mortality_endyr

	Growth_Parameters <- match_report_table("Growth_Parameters", 1, 
		"Growth_Parameters", 1 + ngpatterns * nsexes, header=TRUE, type.convert=TRUE)
	returndat[["Growth_Parameters"]] <- Growth_Parameters

	Seas_Effects <- match_report_table("Seas_Effects", 1, header=TRUE, type.convert=TRUE)
	returndat[["Seas_Effects"]] <- Seas_Effects

	growthCVtype <- match_report_table("Biology_at_age", 0, "Biology_at_age", 1, header=FALSE)
	growthCVtype <- grep("endyr_with_", unlist(growthCVtype), value=TRUE)
	if (length(growthCVtype) > 0) {
		returndat[["growthCVtype"]] <- strsplit(growthCVtype, split="endyr_with_")[[1]][2]
	}
	else {
		returndat[["growthCVtype"]] <- "unknown"
	}

	growdat <- match_report_table("Biology_at_age", adjust1=ifelse(custom, 2, 1), header=TRUE, type.convert=TRUE)
	if (!is.null(growdat)) {
		growdat <- df.rename(growdat, oldnames=c("Gender"), newnames=c("Sex"))
		nmorphs <- max(growdat[["Morph"]])
		midmorphs <- c(c(0, nmorphs/nsexes) + ceiling(nmorphs/nsexes/2))
	}
	returndat[["endgrowth"]] <- growdat

	test <- match_report_table("MEAN_BODY_WT(", 0, "MEAN_BODY_WT(", 1, header=FALSE)
	wtatage_switch <- length(grep("wtatage.ss", test)) > 0
	returndat[["wtatage_switch"]] <- wtatage_switch

	mean_body_wt <- match_report_table("MEAN_BODY_WT(begin)", 1, header=TRUE, type.convert=TRUE)
	returndat[["mean_body_wt"]] <- mean_body_wt

	mean_size <- match_report_table("MEAN_SIZE_TIMESERIES", 1, "mean_size_Jan_1", -2, cols=1:(4 + accuage + 1), header=TRUE, type.convert=TRUE)
	growthvaries <- FALSE
	if (!is.null(mean_size)) {
		if (SS_versionNumeric < 3.3) {
			mean_size <- mean_size[mean_size[["Beg"]] == 1 & mean_size[["Yr"]] >= startyr & mean_size[["Yr"]] < endyr, ]
		}
		else {
			mean_size <- mean_size[mean_size[["SubSeas"]] == 1 & mean_size[["Yr"]] >= startyr & mean_size[["Yr"]] < endyr, ]
		}
		if (nseasons > 1) {
			mean_size <- mean_size[mean_size[["Seas"]] == 1, ]
		}
		for (morph in unique(mean_size[["Morph"]])) {
			if (sum(!duplicated(mean_size[mean_size[["Morph"]] == morph, paste(0:(accuage - 1))])) > 1) {
				growthvaries <- TRUE
			}
		}
		returndat[["growthseries"]] <- mean_size
		returndat[["growthvaries"]] <- growthvaries
	}

	if (!forecast) {
		sizeselex <- sizeselex[sizeselex[["Yr"]] <= endyr, ]
	}
	returndat[["sizeselex"]] <- sizeselex

	ageselex <- match_report_table("COMBINED_ALK*selL*selA", 1, header=TRUE)
	maximum_ASEL2 <- match_report_table("maximum_ASEL2", adjust1=1, header=TRUE, type.convert=TRUE)
	if (!is.null(ageselex)) {
		if (any(grepl("COMBINED_ALK", names(ageselex)))) {
			ageselex <- match_report_table("AGE_SELEX", 5, header=TRUE)
		}
		ageselex <- df.rename(ageselex, oldnames=c("fleet", "year", "seas", "gender", "morph", "label", "factor"), 
			newnames=c("Fleet", "Yr", "Seas", "Sex", "Morph", "Label", "Factor"))
		if (!forecast) {
			ageselex <- ageselex[ageselex[["Yr"]] <= endyr, ]
		}
		ageselex <- type.convert(ageselex, as.is=TRUE)
	}
	returndat[["ageselex"]] <- ageselex
	returndat[["maximum_ASEL2"]] <- maximum_ASEL2

	exploitation_head <- match_report_table("EXPLOITATION", 1, "EXPLOITATION", 20, header=FALSE)
	#if (exploitation_head[1, 1] == "Info:") { ## buggy
	if (any(exploitation_head[1, 1] == c("Info:","NOTE:"))) {  ## (RH 241217)
		exploitation <- match_report_table("EXPLOITATION", which(exploitation_head[, 1] == "Yr"), header=TRUE, blank_lines=rep_blank_lines)
		exploitation <- exploitation[-grep(":", exploitation[, 1]), ]
		F_method_info <- exploitation_head[grep("F_Method:", exploitation_head[, 2]), 2]
		F_method_info <- gsub(pattern=".", replacement=" ", x=F_method_info, fixed=TRUE)
		F_method_info <- strsplit(F_method_info, split=";", fixed=TRUE)[[1]]
		F_method <- as.numeric(strsplit(F_method_info[[1]], split="=", fixed=TRUE)[[1]][2])
	}
	else {
		exploitation <- match_report_table("EXPLOITATION", 5, header=TRUE)
		F_method <- as.numeric(rawrep[match_report_line("F_Method"), 2])
	}
	returndat[["F_method"]] <- F_method
	if (!is.null(exploitation)) {
		exploitation[exploitation == "_"] <- NA
		exploitation[["Yr"]][exploitation[["Yr"]] %in% c("INIT", "init_yr")] <- startyr - 1
		exploitation <- type.convert(exploitation, as.is=TRUE)
	}
	returndat[["exploitation"]] <- exploitation

	catch <- match_report_table("CATCH", adjust1=ifelse(rawrep[match_report_line("CATCH") + 1, 1] == "#", 2, 1), substr1=FALSE, header=TRUE)
	if (!is.null(catch)) {
		catch <- df.rename(catch, oldnames=c("Name", "Yr.frac"), newnames=c("Fleet_Name", "Time"))
		catch[["Like"]][catch[["Like"]] == "-1.#IND"] <- NA
		catch[["Yr"]][tolower(catch[["Yr"]]) == "init"] <- startyr - 1
		catch <- type.convert(catch, as.is=TRUE)
	}
	returndat[["catch"]] <- catch

	summary_age <- rawrep[match_report_line("TIME_SERIES"), ifelse(custom, 3, 2)]
	summary_age <- as.numeric(substring(summary_age, nchar("BioSmry_age:_") + 1))
	returndat[["summary_age"]] <- summary_age

	timeseries <- match_report_table("TIME_SERIES", 1, header=TRUE)
	timeseries <- timeseries[timeseries[["Seas"]] != "recruits", ]
	timeseries[timeseries == "_"] <- NA
	timeseries <- type.convert(timeseries, as.is=TRUE)
	returndat[["timeseries"]] <- timeseries

	## Get the Summary Biomass age (RH 250102)
	zsb = grep("^TIME_SERIES$",rawrep[,1])
	sba = grep("BioSmry", rawrep[zsb,], value=TRUE)
	returndat[["BioSmry_age"]] = as.numeric(rev(strsplit(sba, split="_")[[1]])[1])

	## New report 61  (RH 241217)
	if (SS_versionFull >= "3.30.23") {
		annual_ts <- match_report_table("ANNUAL_TIME_SERIES", 9, header=TRUE)
		annual_ts[annual_ts == "_"] <- NA
		annual_ts <- type.convert(annual_ts, as.is=TRUE)
		returndat[["annual_ts"]] <- annual_ts
	}

	if (!exists("spawnseas")) {
		spawnseas <- unique(timeseries[["Seas"]][!is.na(timeseries[["SpawnBio"]])])
		if (length(spawnseas) == 0) {
			spawnseas <- NA
		}
	}
	returndat[["spawnseas"]] <- spawnseas

	if (is.null(morph_indexing)) {
		mainmorphs <- NULL
	}
	else {
		if (SS_versionNumeric >= 3.3) {
			temp <- morph_indexing[morph_indexing[["BirthSeas"]] == first_seas_with_recruits & 
				morph_indexing[["Platoon_Dist"]] == max(morph_indexing[["Platoon_Dist"]]), ]
			mainmorphs <- min(temp[["Index"]][temp[["Sex"]] == 1])
			if (nsexes == 2) {
				mainmorphs <- c(mainmorphs, min(temp[["Index"]][temp[["Sex"]] == 2]))
			}
		}
		if (SS_versionNumeric < 3.3) {
			temp <- morph_indexing[morph_indexing[["BirthSeas"]] == first_seas_with_recruits & 
				morph_indexing[["Sub_Morph_Dist"]] == max(morph_indexing[["Sub_Morph_Dist"]]), ]
			mainmorphs <- min(temp[["Index"]][temp[["Sex"]] == 1])
			if (nsexes == 2) {
				mainmorphs <- c(mainmorphs, min(temp[["Index"]][temp[["Sex"]] == 2]))
			}
		}
		if (length(mainmorphs) == 0) {
			warning("Error with morph indexing")
		}
	}
	returndat[["mainmorphs"]] <- mainmorphs

	birthseas <- sort(unique(timeseries[["Seas"]][timeseries[["Recruit_0"]] > 
		0]))
	if (length(birthseas) == 0) {
		birthseas <- sort(unique(morph_indexing[["BirthSeas"]]))
	}
	returndat[["birthseas"]] <- birthseas

	timeseries[["Yr"]] <- timeseries[["Yr"]] + (timeseries[["Seas"]] - 1)/nseasons
	ts <- timeseries[timeseries[["Yr"]] <= endyr + 1, ]
	tsyears <- ts[["Yr"]][ts[["Seas"]] == 1]
	tsspaw_bio <- ts[["SpawnBio"]][ts[["Seas"]] == spawnseas & ts[["Area"]] == 1]
	if (nareas > 1) {
		for (a in 2:nareas) {
			tsspaw_bio <- tsspaw_bio + ts[["SpawnBio"]][ts[["Seas"]] == spawnseas & ts[["Area"]] == a]
		}
	}
	if (nsexes == 1) {
		tsspaw_bio <- tsspaw_bio/2
	}
	depletionseries <- tsspaw_bio/tsspaw_bio[1]
	stats[["SBzero"]] <- tsspaw_bio[1]
	stats[["current_depletion"]] <- depletionseries[length(depletionseries)]
	ls <- nrow(ts) - 1
	totretainedmat <- as.matrix(ts[, substr(names(ts), 1, nchar("retain(B)")) == "retain(B)"])
	ts[["totretained"]] <- 0
	ts[["totretained"]][3:ls] <- rowSums(totretainedmat)[3:ls]
	totcatchmat <- as.matrix(ts[, substr(names(ts), 1, nchar("enc(B)")) == "enc(B)"])
	ts[["totcatch"]] <- 0
	ts[["totcatch"]][3:ls] <- rowSums(totcatchmat)[3:ls]
	if (F_method == 1) {
		stringmatch <- "Hrate:_"
	}
	else {
		stringmatch <- "F:_"
	}
	Hrates <- as.matrix(ts[, substr(names(ts), 1, nchar(stringmatch)) == stringmatch])
	fmax <- max(Hrates)
	depletion_basis <- as.numeric(rawrep[match_report_line("Depletion_basis"), 2])
	if (is.na(depletion_basis)) {
		depletion_basis <- as.numeric(rawrep[grep("Depletion_basis",rawrep[,2]),3])  ## (RH 241217)
	}
	if (is.na(depletion_basis)) {
		depletion_basis <- as.numeric(rawrep[match_report_line("Depletion_method"), 2])
	}
	if (depletion_basis %in% c(1, 3:4)) {
		if (file.exists(file.path(dir, "starter.ss"))) {
			starter <- r4ss:::SS_readstarter(file=file.path(dir, "starter.ss"), verbose=verbose)
			depletion_multiplier <- starter[["depl_denom_frac"]]
		}
		else {
			depletion_multiplier <- NULL
		}
	}
	else {
		depletion_multiplier <- 1
	}
	Bratio_denominator <- rawrep[match_report_line("B_ratio_denominator"), 2]
	if (Bratio_denominator == "no_depletion_basis") {
		Bratio_label <- "no_depletion_basis"
	}
	else {
		if (is.null(depletion_multiplier)) {
			depletion_multiplier <- as.numeric(strsplit(Bratio_denominator, "%")[[1]][1])/100
		}
		if (grepl(pattern="100", x=Bratio_denominator)) {
			Bratio_label <- paste0("B/", substring(Bratio_denominator, 6))
		}
		else {
			Bratio_label <- paste0("B/(", Bratio_denominator, ")")
		}
		if (Bratio_label == "B/Virgin_Biomass") {
			Bratio_label <- "B/B_0"
		}
	}
	returndat[["depletion_basis"]] <- depletion_basis
	returndat[["depletion_multiplier"]] <- depletion_multiplier
	returndat[["Bratio_denominator"]] <- Bratio_denominator
	returndat[["Bratio_label"]] <- Bratio_label

	if (SS_versionNumeric < 3.2) {
		DF_discard <- rawrep[match_report_line("DISCARD_OUTPUT"), 3]
		if (length(grep("T_distribution", DF_discard)) > 0) {
			DF_discard <- as.numeric(strsplit(DF_discard, "=_")[[1]][2])
		}
		if (length(grep("_normal_with_Std_in_as_CV", DF_discard)) > 0) {
			DF_discard <- 0
		}
		if (length(grep("_normal_with_Std_in_as_stddev", DF_discard)) > 0) {
			DF_discard <- -1
		}
		if (length(grep("_lognormal", DF_discard)) > 0) {
			DF_discard <- -2
		}
		shift <- 2
		discard_spec <- NULL
	}
	else {
		DF_discard <- NA
		shift <- 1
		discard_header <- match_report_table("DISCARD_SPECIFICATION", 1, "DISCARD_SPECIFICATION", 20)
		if (!is.null(discard_header)) {
			discard_spec <- match_report_table("DISCARD_SPECIFICATION", which(discard_header[, 3] == "errtype"), header=TRUE, type.convert=TRUE)
			discard_spec <- type.convert(discard_spec, as.is=TRUE)
			names(discard_spec)[1] <- "Fleet"
		}
		else {
			discard_spec <- NULL
		}
	}
	discard <- match_report_table("DISCARD_OUTPUT", shift, header=TRUE)
	if (!is.null(discard) && names(discard)[1] != "Fleet") {
		discard <- match_report_table("DISCARD_OUTPUT", shift, header=FALSE)
		names(discard) <- c("Fleet", "Yr", "Seas", "Obs", "Exp", "Std_in", "Std_use", "Dev")
	}
	discard <- df.rename(discard, oldnames=c("Name", "Yr.frac"), newnames=c("Fleet_Name", "Time"))
	if (!is.null(discard) && nrow(discard) > 1) {
		discard[discard == "_"] <- NA
		if (SS_versionNumeric <= 3.23) {
			discard <- type.convert(discard, as.is=TRUE)
			if (!"Fleet_Name" %in% names(discard)) {
				discard[["Fleet_Name"]] <- discard[["Fleet"]]
			}
			discard[["Fleet"]] <- NA
			for (i in 1:nrow(discard)) {
				discard[["Fleet"]][i] <- strsplit(discard[["Fleet_Name"]][i], "_")[[1]][1]
				discard[["Fleet_Name"]][i] <- substring(discard[["Fleet_Name"]][i], nchar(discard[["Fleet"]][i]) + 2)
			}
			discard_tuning_info <- NULL
		}
		else {
			discard <- type.convert(discard, as.is=TRUE)
			discard_tuning_info <- calc_var_adjust(discard, type="sd")
		}
	}
	else {
		discard <- NA
		discard_tuning_info <- NULL
	}
	returndat[["discard"]] <- discard
	returndat[["discard_spec"]] <- discard_spec
	returndat[["discard_tuning_info"]] <- discard_tuning_info
	returndat[["DF_discard"]] <- DF_discard

	DF_mnwgt <- rawrep[match_report_line("log(L)_based_on_T_distribution"), 1]
	if (!is.na(DF_mnwgt)) {
		DF_mnwgt <- as.numeric(strsplit(DF_mnwgt, "=_")[[1]][2])
		mnwgt <- match_report_table("MEAN_BODY_WT_OUTPUT", 2, header=TRUE)
		mnwgt <- df.rename(mnwgt, oldnames=c("Name"), newnames=c("Fleet_Name"))
		mnwgt[mnwgt == "_"] <- NA
		if (SS_versionNumeric <= 3.23) {
			mnwgt <- type.convert(mnwgt, as.is=TRUE)
			if (!"Fleet_Name" %in% names(mnwgt)) {
				mnwgt[["Fleet_Name"]] <- mnwgt[["Fleet"]]
			}
			mnwgt[["Fleet"]] <- NA
			for (i in 1:nrow(mnwgt)) {
				mnwgt[["Fleet"]][i] <- strsplit(mnwgt[["Fleet_Name"]][i], "_")[[1]][1]
				mnwgt[["Fleet_Name"]][i] <- substring(mnwgt[["Fleet_Name"]][i], nchar(mnwgt[["Fleet_Name"]][i]) + 2)
			}
			mnwgt_tuning_info <- NULL
		}
		else {
			mnwgt <- type.convert(mnwgt, as.is=TRUE)
			mnwgt_tuning_info <- calc_var_adjust(mnwgt, type="CV")
		}
	}
	else {
		DF_mnwgt <- NA
		#mnwgt <- NA
		mnwgt <- match_report_table("MEAN_BODY_WT(Begin)", 2, header=TRUE)  ## (RH 241217) beginning mean weight
		mnwgt_tuning_info <- NULL
	}
	returndat[["mnwgt"]] <- mnwgt
	returndat[["mnwgt_tuning_info"]] <- mnwgt_tuning_info
	returndat[["DF_mnwgt"]] <- DF_mnwgt

	#spr <- match_report_table("SPR_SERIES", 5, header=TRUE)  ## buggy
	spr <- match_report_table("SPR_SERIES", ifelse(SS_versionFull >= "3.30.23", 6, 5), header=TRUE)  ## (RH 241217)
	if (is.null(spr)) {
		spr <- match_report_table("SPR_series", 5, header=TRUE)
	}
	if (!is.null(spr)) {
		names(spr) <- gsub(pattern="SPB", replacement="SSB", names(spr))
		spr <- df.rename(spr, oldnames=c("Year", "spawn_bio", "SPR_std", "Y/R", "F_std"), 
			newnames=c("Yr", "SpawnBio", "SPR_report", "YPR", "F_report"))
		spr[spr == "_"] <- NA
		spr[spr == "&"] <- NA
		spr[spr == "-1.#IND"] <- NA
		spr <- type.convert(spr, as.is=TRUE)
		spr[["spr"]] <- spr[["SPR"]]
		stats[["last_years_SPR"]] <- spr[["spr"]][nrow(spr)]
		stats[["SPRratioLabel"]] <- managementratiolabels[1, 2]
		stats[["last_years_SPRratio"]] <- spr[["SPR_std"]][nrow(spr)]
	}
	returndat[["sprseries"]] <- spr
	returndat[["managementratiolabels"]] <- managementratiolabels
	returndat[["F_report_basis"]] <- managementratiolabels[["Label"]][2]
	returndat[["sprtarg"]] <- sprtarg
	returndat[["btarg"]] <- btarg

	if (!is.na(btarg) & btarg == 0.4 & startyr == 1966 & sprtarg == 0.4 & accuage == 20 & wtatage_switch) {
		if (verbose) {
			message("Setting minimum biomass threshhold to 0.10", 
				" because this looks like the Pacific Hake model.", 
				" You can replace or override in SS_plots via the", " 'minbthresh' input.")
		}
		minbthresh <- 0.1
	}
	returndat[["minbthresh"]] <- minbthresh

	if (length(grep("Kobe_Plot", rawrep[, 1])) != 0) {
		Kobe_head <- match_report_table("Kobe_Plot", 0, "Kobe_Plot", 5, header=TRUE)
		shift <- grep("^Y(ea)?r", Kobe_head[, 1])
		if (length(shift) == 0) {
			shift <- grep("MSY_basis:_Y(ea)?r", Kobe_head[, 1])
			if (length(shift) == 0) {
				stop("Bug: r4ss cannot find the start of table for the Kobe plot.")
			}
		}
		Kobe_warn <- NA
		Kobe_MSY_basis <- NA
		if (length(grep("_basis_is_not", Kobe_head[1, 1])) > 0) {
			Kobe_warn <- Kobe_head[1, 1]
		}
		if (length(grep("MSY_basis", Kobe_head[2, 1])) > 0) {
			Kobe_MSY_basis <- Kobe_head[2, 1]
		}
		Kobe <- match_report_table("Kobe_Plot", shift, header=TRUE)
		Kobe[Kobe == "_"] <- NA
		Kobe[Kobe == "1.#INF"] <- NA
		Kobe[Kobe == "-1.#IND"] <- NA
		names(Kobe) <- gsub("/", ".", names(Kobe), fixed=TRUE)
		Kobe[, 1:3] <- lapply(Kobe[, 1:3], as.numeric)
	}
	else {
		Kobe <- NA
		Kobe_warn <- NA
		Kobe_MSY_basis <- NA
	}
	returndat[["Kobe_warn"]] <- Kobe_warn
	returndat[["Kobe_MSY_basis"]] <- Kobe_MSY_basis
	returndat[["Kobe"]] <- Kobe

	flush.console()
	INDEX_1 <- match_report_table("INDEX_1", 1, "INDEX_1", (nfleets + 1), header=TRUE)
	INDEX_1 <- df.rename(INDEX_1, oldnames=c("NoName", "fleetname"), newnames=c("Name", "Name"))
	if (SS_versionNumeric >= 3.3) {
		ncpue_column <- 11
		INDEX_1 <- match_report_table("INDEX_1", 1, "INDEX_3", -4, header=TRUE)
		INDEX_1 <- INDEX_1[substr(INDEX_1[["Fleet"]], 1, 1) != "#", ]
		ncpue <- sum(as.numeric(INDEX_1[["N"]]), na.rm=TRUE)
	}
	else {
		ncpue_column <- 11
		ncpue <- sum(as.numeric(rawrep[match_report_line("INDEX_1") + 1 + 1:nfleets, ncpue_column]))
	}
	returndat[["index_variance_tuning_check"]] <- INDEX_1

	cpue <- match_report_table("INDEX_2", 1, "INDEX_2", ncpue + 1, header=TRUE)
	cpue[cpue == "_"] <- NA
	if (length(cpue) > 0) {
		cpue <- df.rename(cpue, oldnames=c("Yr.S", "Yr.frac", "Supr_Per", "Name"), 
			newnames=c("Time", "Time", "SuprPer", "Fleet_name"))
		if (SS_versionNumeric < 3.24) {
			cpue[["Name"]] <- NA
			for (i in 1:nrow(cpue)) {
				cpue[["Fleet"]][i] <- strsplit(cpue[["Fleet"]][i], "_")[[1]][1]
				cpue[["Name"]][i] <- substring(cpue[["Fleet"]][i], nchar(cpue[["Fleet"]][i]) + 2)
			}
		}
		if (any(cpue[["Exp"]] == "1.#QNAN")) {
			cpue[["Exp"]][cpue[["Exp"]] == "1.#QNAN"] <- NA
			cpue[["Calc_Q"]][cpue[["Calc_Q"]] == "1.#QNAN"] <- NA
			cpue[["Eff_Q"]][cpue[["Eff_Q"]] == "1.#QNAN"] <- NA
		}
		badrows <- which(cpue[["Use"]] == "")
		if (length(badrows) > 0) {
			columns <- which(names(cpue) == "SE_input"):which(names(cpue) == "Use")
			cpue[badrows, columns] <- cpue[badrows, columns - 1]
			cpue[badrows, "SE_input"] <- NA
		}
		cpue <- type.convert(cpue, as.is=TRUE)
	}
	else {
		cpue <- NULL
	}
	returndat[["cpue"]] <- cpue

	natage <- match_report_table("NUMBERS_AT_AGE", 1, substr1=FALSE, header=TRUE, type.convert=TRUE)
	if (is.null(natage) || nrow(natage) == 0) {
		natage <- NULL
	}
	else {
		natage <- df.rename(natage, oldnames=c("Gender", "SubMorph"), newnames=c("Sex", "Platoon"))
	}
	returndat[["natage"]] <- natage
	natage_annual_1_no_fishery <- match_report_table("NUMBERS_AT_AGE_Annual_1", 1, header=TRUE, type.convert=TRUE)
	natage_annual_2_with_fishery <- match_report_table("NUMBERS_AT_AGE_Annual_2", 1, header=TRUE, type.convert=TRUE)
	returndat[["natage_annual_1_no_fishery"]] <- natage_annual_1_no_fishery
	returndat[["natage_annual_2_with_fishery"]] <- natage_annual_2_with_fishery

	batage <- match_report_table("BIOMASS_AT_AGE", 1, substr1=FALSE, header=TRUE, type.convert=TRUE)
	returndat[["batage"]] <- batage

	col.adjust <- 12
	if (SS_versionNumeric < 3.3) {
		col.adjust <- 11
	}
	natlen <- match_report_table("NUMBERS_AT_LENGTH", 1, substr1=FALSE, header=TRUE, type.convert=TRUE)
	natlen <- df.rename(natlen, oldnames=c("Gender", "SubMorph"), newnames=c("Sex", "Platoon"))
	returndat[["natlen"]] <- natlen

	batlen <- match_report_table("BIOMASS_AT_LENGTH", 1, substr1=FALSE, header=TRUE, type.convert=TRUE)
	returndat[["batlen"]] <- batlen
	fatage <- match_report_table("F_AT_AGE", 1, header=TRUE, type.convert=TRUE)
	returndat[["fatage"]] <- fatage

	discard_at_age <- match_report_table("DISCARD_AT_AGE", 1, header=TRUE, type.convert=TRUE)
	returndat[["discard_at_age"]] <- discard_at_age

	catage <- match_report_table("CATCH_AT_AGE", ifelse(SS_versionFull >= "3.30.23", 2, 1), header=TRUE, type.convert=TRUE)  ## (RH 241217)
	returndat[["catage"]] <- catage

	movement <- match_report_table("MOVEMENT", 1, substr1=FALSE, header=TRUE)
	if (!is.null(movement)) {
		names(movement) <- c(names(movement)[1:6], paste("age", names(movement)[-(1:6)], sep=""))
		movement <- df.rename(movement, oldnames=c("Gpattern"), newnames=c("GP"))
		for (i in 1:ncol(movement)) {
			movement[, i] <- as.numeric(movement[, i])
		}
	}
	returndat[["movement"]] <- movement

	tagreportrates <- match_report_table("Reporting_Rates_by_Fishery", 1, 
		"See_composition_data_output", -1, substr2=TRUE, header=TRUE, type.convert=TRUE)
	returndat[["tagreportrates"]] <- tagreportrates

	tagrelease <- match_report_table("TAG_Recapture", 1, "Tags_Alive", -1, cols=1:10)
	if (!is.null(tagrelease)) {
		tagfirstperiod <- as.numeric(tagrelease[1, 1])
		tagaccumperiod <- as.numeric(tagrelease[2, 1])
		names(tagrelease) <- tagrelease[4, ]
		tagrelease <- tagrelease[-(1:4), ]
		tagrelease <- type.convert(tagrelease, as.is=TRUE)
	}
	else {
		tagrelease <- NULL
		tagfirstperiod <- NULL
		tagaccumperiod <- NULL
	}
	returndat[["tagrelease"]] <- tagrelease
	returndat[["tagfirstperiod"]] <- tagfirstperiod
	returndat[["tagaccumperiod"]] <- tagaccumperiod

	tagsalive <- match_report_table("Tags_Alive", 1, "Total_recaptures", -1)
	if (!is.null(tagsalive)) {
		tagcols <- ncol(tagsalive)
		names(tagsalive) <- c("TG", paste0("period", 0:(tagcols - 2)))
		tagsalive[tagsalive == ""] <- NA
		tagsalive <- type.convert(tagsalive, as.is=TRUE)
	}
	returndat[["tagsalive"]] <- tagsalive

	tagtotrecap <- match_report_table("Total_recaptures", 1)
	if (!is.null(tagtotrecap)) {
		tagcols <- ncol(tagtotrecap)
		names(tagtotrecap) <- c("TG", paste0("period", 0:(tagcols - 2)))
		tagtotrecap[tagtotrecap == ""] <- NA
		tagtotrecap <- type.convert(tagtotrecap, as.is=TRUE)
	}
	returndat[["tagtotrecap"]] <- tagtotrecap
	sdsize_lines <- grep("^sdsize", rawrep[, 1])
	if (length(sdsize_lines) > 0) {
		which_blank <- 1 + length(rep_blank_or_hash_lines[rep_blank_or_hash_lines > 
			match_report_line("AGE_LENGTH_KEY") & rep_blank_or_hash_lines < max(sdsize_lines)])
		rawALK <- match_report_table("AGE_LENGTH_KEY", 4, cols=1:max(6, accuage + 2), header=FALSE, which_blank=which_blank)
		if (length(rawALK) > 1 && length(grep("AGE_AGE_KEY", 
			rawALK[, 1])) == 0) {
			morph_col <- 5
			if (SS_versionNumeric < 3.3 & length(grep("Sub_Seas", rawALK[, 3])) == 0) {
				morph_col <- 3
			}
			starts <- grep("Morph:", rawALK[, morph_col]) + 2
			ends <- grep("mean", rawALK[, 1]) - 1
			N_ALKs <- length(starts)
			ALK <- array(NA, c(nlbinspop, accuage + 1, N_ALKs))
			dimnames(ALK) <- list(Length=rev(lbinspop), TrueAge=0:accuage, Matrix=1:N_ALKs)
			for (i in 1:N_ALKs) {
				ALKtemp <- rawALK[starts[i]:ends[i], 2 + 0:accuage]
				ALKtemp <- type.convert(ALKtemp, as.is=TRUE)
				ALK[, , i] <- as.matrix(ALKtemp)
				Matrix.Info <- rawALK[starts[i] - 2, ]
				Matrix.Info <- Matrix.Info[Matrix.Info != ""]
				dimnames(ALK)$Matrix[i] <- paste(Matrix.Info, collapse=" ")
			}
			returndat[["ALK"]] <- ALK
		}
	}
	rawAAK <- match_report_table("AGE_AGE_KEY", 1)
	if (!is.null(rawAAK)) {
		if (rawAAK[[1]][1] == "no_age_error_key_used" | is.null(dim(rawAAK))) {
			N_ageerror_defs <- 0
		}
		else {
			starts <- grep("KEY:", rawAAK[, 1])
			N_ageerror_defs <- length(starts)
			if (N_ageerror_defs > 0) {
				nrowsAAK <- nrow(rawAAK)/N_ageerror_defs - 3
				AAK <- array(NA, c(N_ageerror_defs, nrowsAAK, accuage + 1))
				age_error_mean <- age_error_sd <- data.frame(age=0:accuage)
				for (i in 1:N_ageerror_defs) {
					AAKtemp <- rawAAK[starts[i] + 2 + 1:nrowsAAK, -1]
					rownames.tmp <- rawAAK[starts[i] + 2 + 1:nrowsAAK, 1]
					AAKtemp <- type.convert(AAKtemp, as.is=TRUE)
					AAK[i, , ] <- as.matrix(AAKtemp)
					age_error_mean[[paste("type", i, sep="")]] <- as.numeric((rawAAK[starts[i] + 1, -1]))
					age_error_sd[[paste("type", i, sep="")]] <- as.numeric((rawAAK[starts[i] + 2, -1]))
				}
				if (!is.null(AAK)) {
					dimnames(AAK) <- list(AgeingErrorType=1:N_ageerror_defs, ObsAgeBin=rownames.tmp, TrueAge=0:accuage)
				}
				returndat[["AAK"]] <- AAK
				returndat[["age_error_mean"]] <- age_error_mean
				returndat[["age_error_sd"]] <- age_error_sd
			}
		}
		returndat[["N_ageerror_defs"]] <- N_ageerror_defs
	}
	if (SS_versionNumeric >= 3.3) {
		yieldraw <- match_report_table("SPR/YPR_Profile", 1, "Finish", -2)
	}
	else {
		yieldraw <- match_report_table("SPR/YPR_Profile", 1)
	}
	if (!is.null(yieldraw)) {
		names <- yieldraw[1, ]
		names[names == "SSB/Bzero"] <- "Depletion"
		yielddat <- yieldraw[c(2:(as.numeric(length(yieldraw[, 1]) - 1))), ]
		yielddat[yielddat == "-nan(ind)"] <- NA
		names(yielddat) <- names
		if ("SPRloop" %in% names) {
			yielddat <- dplyr::filter(yielddat, SPRloop != "ready")
		}
		yielddat <- type.convert(yielddat, as.is=TRUE)
	}
	else {
		yielddat <- NA
	}
	returndat[["equil_yield"]] <- yielddat
	Z_at_age <- match_report_table("Z_AT_AGE_Annual_2", 1, header=TRUE)
	if (!is.null(Z_at_age)) {
		Z_at_age[Z_at_age == "_"] <- NA
		Z_at_age[Z_at_age == "-1.#INF"] <- NA
		Z_at_age <- type.convert(Z_at_age, as.is=TRUE)
	}
	returndat[["Z_at_age"]] <- Z_at_age
	if (!is.na(match_report_line("Report_Z_by_area_morph_platoon"))) {
		M_at_age <- match_report_table("Z_AT_AGE_Annual_1", 1, header=TRUE)
	}
	else {
		M_at_age <- match_report_table("Z_AT_AGE_Annual_1", 1, "-ln(Nt+1", -1, matchcol2=5, header=TRUE)
	}
	if (!is.null(M_at_age)) {
		M_at_age[M_at_age == "_"] <- NA
		M_at_age[M_at_age == "-1.#INF"] <- NA
		M_at_age <- type.convert(M_at_age, as.is=TRUE)
	}
	returndat[["M_at_age"]] <- M_at_age
	if (is.na(match_report_line("Report_Z_by_area_morph_platoon"))) {
		Z_by_area <- NULL
		M_by_area <- NULL
	}
	else {
		if (!is.na(match_report_line("Report_Z_by_area_morph_platoon_2"))) {
			Z_by_area <- match_report_table("Report_Z_by_area_morph_platoon_2", adjust1=1, header=TRUE, type.convert=TRUE)
			M_by_area <- match_report_table("Report_Z_by_area_morph_platoon_1", adjust1=1, adjust2=-3, header=TRUE, type.convert=TRUE)
		}
		else {
			Report_Z_by_area_morph_platoon <- match_report_table("Report_Z_by_area_morph_platoon", adjust1=1, header=FALSE)
			Z_by_area <- match_report_table("With_fishery", adjust1=1, "No_fishery_for_Z=M", adjust2=-1, matchcol1=2, 
				matchcol2=2, obj=Report_Z_by_area_morph_platoon, header=TRUE, type.convert=TRUE)
			M_by_area <- match_report_table("No_fishery_for_Z=M", blank_lines=nrow(Report_Z_by_area_morph_platoon) + 
				  1, adjust1=1, matchcol1=2, obj=Report_Z_by_area_morph_platoon, header=TRUE, type.convert=TRUE)
		}
		returndat["Z_by_area"] <- list(Z_by_area)
		returndat["M_by_area"] <- list(M_by_area)
	}
	Dynamic_Bzero <- match_report_table("Spawning_Biomass_Report_2", 1)
	Dynamic_Bzero2 <- match_report_table("Spawning_Biomass_Report_1", 1)
	if (!is.null(Dynamic_Bzero)) {
		Dynamic_Bzero <- cbind(Dynamic_Bzero, Dynamic_Bzero2[, -(1:2)])
		Dynamic_Bzero <- type.convert(Dynamic_Bzero[-(1:2), ], as.is=TRUE)
		if (ncol(Dynamic_Bzero) == 4) {
			names(Dynamic_Bzero) <- c("Yr", "Era", "SSB", "SSB_nofishing")
		}
		if (nareas > 1 & !is.null(ngpatterns) && ngpatterns == 1) {
			names(Dynamic_Bzero) <- c("Yr", "Era", paste0("SSB_area", 1:nareas), paste0("SSB_nofishing_area", 1:nareas))
			Dynamic_Bzero[["SSB"]] <- apply(Dynamic_Bzero[, 2 + 1:nareas], 1, sum)
			Dynamic_Bzero[["SSB_nofishing"]] <- apply(Dynamic_Bzero[, 2 + nareas + 1:nareas], 1, sum)
		}
	}
	returndat[["Dynamic_Bzero"]] <- Dynamic_Bzero
	if (comp) {
		returndat[["comp_data_exists"]] <- TRUE
		returndat[["lendbase"]] <- lendbase
		returndat[["sizedbase"]] <- sizedbase
		returndat[["agedbase"]] <- agedbase
		returndat[["condbase"]] <- condbase
		returndat[["ghostagedbase"]] <- ghostagedbase
		returndat[["ghostcondbase"]] <- ghostcondbase
		returndat[["ghostlendbase"]] <- ghostlendbase
		returndat[["ladbase"]] <- ladbase
		returndat[["wadbase"]] <- wadbase
		returndat[["tagdbase1"]] <- tagdbase1
		returndat[["tagdbase2"]] <- tagdbase2
		returndat[["morphcompdbase"]] <- morphcompdbase
	}
	else {
		returndat[["comp_data_exists"]] <- FALSE
	}
	returndat[["len_comp_fit_table"]] <- fit_len_comps
	returndat[["age_comp_fit_table"]] <- fit_age_comps
	returndat[["size_comp_fit_table"]] <- fit_size_comps
	returndat[["derived_quants"]] <- der
	returndat[["parameters"]] <- parameters
	returndat[["Dirichlet_Multinomial_pars"]] <- DM_pars
	returndat[["FleetNames"]] <- FleetNames
	returndat[["repfiletime"]] <- repfiletime
	SRRtype <- rawrep[match_report_line("SPAWN_RECRUIT"), 3]
	if (!is.na(SRRtype) && SRRtype == "Function:") {
		SRRtype <- as.numeric(rawrep[match_report_line("SPAWN_RECRUIT"), 4])
	}
	returndat[["SRRtype"]] <- SRRtype
	SSB_final_Label <- paste0("SSB_", endyr + 1)
	if (SSB_final_Label %in% der[["Label"]]) {
		SSB_final_EST <- der[["Value"]][der[["Label"]] == SSB_final_Label]
		SSB_final_SD <- der[["StdDev"]][der[["Label"]] == SSB_final_Label]
		returndat[["Pstar_sigma"]] <- sqrt(log((SSB_final_SD/SSB_final_EST)^2 + 1))
	}
	else {
		returndat[["Pstar_sigma"]] <- NULL
	}
	OFL_final_Label <- paste0("OFLCatch_", endyr + 1)
	if (OFL_final_Label %in% der[["Label"]]) {
		OFL_final_EST <- der[["Value"]][der[["Label"]] == OFL_final_Label]
		OFL_final_SD <- der[["StdDev"]][der[["Label"]] == OFL_final_Label]
		returndat[["OFL_sigma"]] <- sqrt(log((OFL_final_SD/OFL_final_EST)^2 + 1))
	}
	else {
		returndat[["OFL_sigma"]] <- NULL
	}
	if (covar) {
		returndat[["CoVar"]] <- CoVar
		returndat[["stdtable"]] <- stdtable
	}
	recdevEarly <- parameters[substring(parameters[["Label"]], 1, 13) == "Early_RecrDev", ]
	early_initage <- parameters[substring(parameters[["Label"]], 1, 13) == "Early_InitAge", ]
	main_initage <- parameters[substring(parameters[["Label"]], 1, 12) == "Main_InitAge", ]
	recdev <- parameters[substring(parameters[["Label"]], 1, 12) == "Main_RecrDev", ]
	recdevFore <- parameters[substring(parameters[["Label"]], 1, 8) == "ForeRecr", ]
	recdevLate <- parameters[substring(parameters[["Label"]], 1, 12) == "Late_RecrDev", ]
	recruitpars <- NULL
	if (nrow(early_initage) > 0) {
		early_initage[["type"]] <- "Early_InitAge"
		early_initage[["Yr"]] <- startyr - as.numeric(substring(early_initage[["Label"]], 15))
		recruitpars <- rbind(recruitpars, early_initage)
	}
	if (nrow(recdevEarly) > 0) {
		recdevEarly[["type"]] <- "Early_RecrDev"
		recdevEarly[["Yr"]] <- as.numeric(substring(recdevEarly[["Label"]], 15))
		recruitpars <- rbind(recruitpars, recdevEarly)
	}
	if (nrow(main_initage) > 0) {
		main_initage[["type"]] <- "Main_InitAge"
		main_initage[["Yr"]] <- startyr - as.numeric(substring(main_initage[["Label"]], 14))
		recruitpars <- rbind(recruitpars, main_initage)
	}
	if (nrow(recdev) > 0) {
		recdev[["type"]] <- "Main_RecrDev"
		recdev[["Yr"]] <- as.numeric(substring(recdev[["Label"]], 14))
		recruitpars <- rbind(recruitpars, recdev)
	}
	if (nrow(recdevFore) > 0) {
		recdevFore[["type"]] <- "ForeRecr"
		recdevFore[["Yr"]] <- as.numeric(substring(recdevFore[["Label"]], 10))
		recruitpars <- rbind(recruitpars, recdevFore)
	}
	if (nrow(recdevLate) > 0) {
		recdevLate[["type"]] <- "Late_RecrDev"
		recdevLate[["Yr"]] <- as.numeric(substring(recdevLate[["Label"]], 14))
		recruitpars <- rbind(recruitpars, recdevLate)
	}
	if (!is.null(recruitpars)) {
		recruitpars <- recruitpars[order(recruitpars[["Yr"]]), c("Value", "Parm_StDev", "type", "Yr")]
	}
	returndat[["recruitpars"]] <- recruitpars
	if (is.null(recruitpars)) {
		sigma_R_info <- NULL
	}
	else {
		sigma_R_info <- data.frame(period=c("Main", "Early+Main", "Early+Main+Late"), 
			N_devs=0, SD_of_devs=NA, Var_of_devs=NA, mean_SE=NA, mean_SEsquared=NA)
		subset <- recruitpars[["type"]] %in% c("Main_InitAge", "Main_RecrDev")
		within_period <- sigma_R_info[["period"]] == "Main"
		sigma_R_info[["N_devs"]][within_period] <- sum(subset)
		sigma_R_info[["SD_of_devs"]][within_period] <- sd(recruitpars[["Value"]][subset])
		sigma_R_info[["mean_SE"]][within_period] <- mean(recruitpars[["Parm_StDev"]][subset])
		sigma_R_info[["mean_SEsquared"]][within_period] <- mean((recruitpars[["Parm_StDev"]][subset])^2)
		subset <- recruitpars[["type"]] %in% c("Early_RecrDev", "Early_InitAge", "Main_InitAge", "Main_RecrDev")
		within_period <- sigma_R_info[["period"]] == "Early+Main"
		sigma_R_info[["N_devs"]][within_period] <- sum(subset)
		sigma_R_info[["SD_of_devs"]][within_period] <- sd(recruitpars[["Value"]][subset])
		sigma_R_info[["mean_SE"]][within_period] <- mean(recruitpars[["Parm_StDev"]][subset])
		sigma_R_info[["mean_SEsquared"]][within_period] <- mean((recruitpars[["Parm_StDev"]][subset])^2)
		subset <- recruitpars[["type"]] %in% c("Early_RecrDev", "Early_InitAge", "Main_InitAge", "Main_RecrDev", "Late_RecrDev")
		within_period <- sigma_R_info[["period"]] == "Early+Main+Late"
		sigma_R_info[["N_devs"]][within_period] <- sum(subset)
		sigma_R_info[["SD_of_devs"]][within_period] <- sd(recruitpars[["Value"]][subset])
		sigma_R_info[["mean_SE"]][within_period] <- mean(recruitpars[["Parm_StDev"]][subset])
		sigma_R_info[["mean_SEsquared"]][within_period] <- mean((recruitpars[["Parm_StDev"]][subset])^2)
		sigma_R_info[["Var_of_devs"]] <- sigma_R_info[["SD_of_devs"]]^2
		sigma_R_info[["sqrt_sum_of_components"]] <- sqrt(sigma_R_info[["Var_of_devs"]] + sigma_R_info[["mean_SEsquared"]])
		sigma_R_info[["SD_of_devs_over_sigma_R"]] <- sigma_R_info[["SD_of_devs"]]/sigma_R_in
		sigma_R_info[["sqrt_sum_over_sigma_R"]] <- sigma_R_info[["sqrt_sum_of_components"]]/sigma_R_in
		sigma_R_info[["alternative_sigma_R"]] <- sigma_R_in * sigma_R_info[["sqrt_sum_over_sigma_R"]]
		sigma_R_info[["alternative_sigma_R"]][sigma_R_info[["mean_SE"]] == 0] <- "needs_Hessian"
	}
	stats[["sigma_R_in"]] <- sigma_R_in
	stats[["sigma_R_info"]] <- sigma_R_info
	stats[["rmse_table"]] <- rmse_table
	stats[["RecDev_method"]] <- RecDev_method
	RecrDistpars <- parameters[substring(parameters[["Label"]], 1, 8) == "RecrDist", ]
	returndat[["RecrDistpars"]] <- RecrDistpars
	returndat[["wtatage"]] <- wtatage
	returndat[["jitter_info"]] <- jitter_info
	returndat <- c(returndat, stats)
	returndat[["seldev_pars"]] <- seldev_pars
	returndat[["seldev_matrix"]] <- seldev_matrix
	if (printstats) {
		message("\nStatistics shown below (to turn off, change input to printstats=FALSE)")
		stats[["likelihoods_used"]] <- format(stats[["likelihoods_used"]], scientific=20)
		stats[["estimated_non_dev_parameters"]] <- format(stats[["estimated_non_dev_parameters"]], scientific=20)
		print(stats)
	}
	returndat[["logfile"]] <- logfile
	inputs <- list()
	inputs[["dir"]] <- dir
	inputs[["repfile"]] <- repfile
	inputs[["forecast"]] <- forecast
	inputs[["warn"]] <- warn
	inputs[["covar"]] <- covar
	inputs[["verbose"]] <- verbose
	returndat[["inputs"]] <- inputs
	if (verbose) {
		message("completed getSS.output")
	}
	invisible(returndat)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~getSS.output


## make.multifig------------------------2023-09-13
##  Make a multi-panel figure
##  Source: R package 'r4ss' v.1.39.1
## Based on 'r4ss::make_multifig'
## ----------------------------------------r4ss|RH
make.multifig <- function (ptsx, ptsy, yr, linesx=0, linesy=0, ptsSD=0, 
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
   colsex = c("orange","limegreen","slategray3"),
   colvec=c(rgb(1,0,0,0.7), rgb(0,0,1,0.7), rgb(0.1,0.1,0.1,0.7)),
   linescol=c(rgb(0,0.8,0,0.7), rgb(1,0,0,0.7), rgb(0,0,1,0.7)),
   lty=1, lwd=2, pch=1, nlegends=3, 
   legtext=list("yr", "sampsize", "effN"), legx="default", 
   legy="default", legadjx="default", legadjy="default", 
   legsize=c(1.2, 1), legfont=c(2, 1), venusmars=TRUE, 
   sampsizeline=FALSE, effNline=FALSE, sampsizemean=NULL, 
   effNmean=NULL, ipage=0, scalebins=FALSE, sexvec=NULL, 
   multifig_colpolygon=c("grey60","grey80","grey70"), 
   multifig_oma=NULL, lang="e", obsN=0, ...) 
{
#browser();return()
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
			bin.width.table <- data.frame(bin=allbins.obs, 
				width=diffs)
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
	if (length(sampsize) == 1) { sampsize <- 0 }
	if (length(effN) == 1) { effN <- 0 }
	if (length(obsN) == 1) { obsN <- 0 }
	par_old <- par()
	if (is.null(multifig_oma)) {
		if (main == "") {
			multifig_oma <- c(3.5, 5, 1, 1) + 0.1  ## (RH 200825)
		}
		else {
			multifig_oma <- c(3.5, 5, 5, 1) + 0.1
		}
	}
	#par(mfcol=c(nrows, ncols), mar=rep(0, 4), oma=multifig_oma, ...)
	if (npanels==1)
		par(mfrow=c(nrows, ncols), mar=multifig_oma, oma=rep(0, 4),  mgp = c(2,0.5,0), ...)  ## (RH 201208)
	else 
		par(mfrow=c(nrows, ncols), mar=rep(0, 4), oma=multifig_oma,  mgp = c(2,0.5,0), ...)  ## (RH 201207)
	panelrange <- 1:npanels
	if (npages > 1 & ipage != 0) {
		panelrange <- intersect(panelrange, 1:(nrows * ncols) + nrows * ncols * (ipage - 1))
	}
#browser();return()
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
		#abline(h=0, col="grey")
		if (linepos == 2) {
			lines(linesx_i0, linesy_i0, col=linescol[1], lwd=lwd, lty=lty)
			lines(linesx_i1, linesy_i1, col=linescol[2], lwd=lwd, lty=lty)
			lines(linesx_i2, linesy_i2, col=linescol[3], lwd=lwd, lty=lty)
		}
		if (bub) {
			if (length(z_i0) > 0) {
				bubble3(x=ptsx_i0, y=ptsy_i0, z=z_i0, col=rep(colvec[3], length(z_i0)), cexZ1=cexZ1, legend.yadj=1.5, legend=linguaFranca(bublegend,lang), legendloc="topright", maxsize=maxsize, minnbubble=minnbubble, allopen=allopen, add=TRUE)
			}
			if (length(z_i1) > 0) {
				bubble3(x=ptsx_i1, y=ptsy_i1, z=z_i1, col=rep(colvec[1], length(z_i1)), cexZ1=cexZ1, legend.yadj=1.5, legend=linguaFranca(bublegend,lang), legendloc="topright", maxsize=maxsize, minnbubble=minnbubble, allopen=allopen, add=TRUE)
			}
			if (length(z_i2) > 0) {
				bubble3(x=ptsx_i2, y=abs(ptsy_i2), z=z_i2, col=rep(colvec[2], length(z_i2)), cexZ1=cexZ1, legend.yadj=1.5, legend=linguaFranca(bublegend,lang), legendloc="topright", maxsize=maxsize, minnbubble=minnbubble, allopen=allopen, add=TRUE)
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
		} else {
			if (length(ptsx_i0) > 0) {
				if (bars) {
					#drawBars(ptsx_i0, ptsy_i0, col="slategray3", fill="slategray3")  ## (RH 200825)
					drawBars(ptsx_i0, ptsy_i0, col=colvec[3], fill=colvec[3], width=1, lwd=0.5)  ## (RH 210426)
				} else {
					if (polygons)
						polygon(c(ptsx_i0[1], ptsx_i0, tail(ptsx_i0,1)), c(0,ptsy_i0, 0), col=multifig_colpolygon[1])
					points(ptsx_i0, ptsy_i0, type=type, lwd=1,  pch=16, cex=0.7, col=ptscol)
				}
			}
			if (length(ptsx_i1) > 0) {
				if (bars) {
					#drawBars(ptsx_i1, ptsy_i1, col=lucent("pink",0.75), fill=lucent("pink",0.75))  ## (RH 200826) -- females
					drawBars(ptsx_i1, ptsy_i1, col=colvec[1], fill=colvec[1], width=1, lwd=0.5)  ## (RH 210426) -- females
				} else {
					if (polygons)
						polygon(c(ptsx_i1[1], ptsx_i1, tail(ptsx_i1,1)), c(0,ptsy_i1, 0), col=multifig_colpolygon[1])
					points(ptsx_i1, ptsy_i1, type=type, lwd=1,  pch=16, cex=0.7, col=ptscol)
				}
			}
			if (length(ptsx_i2) > 0) {
				if (bars) {
					#drawBars(ptsx_i2, ptsy_i2, col=lucent("skyblue",0.75), fill=lucent("skyblue",0.75))  ## (RH 200826) -- males
					drawBars(ptsx_i2, ptsy_i2, col=colvec[2], fill=colvec[2], width=1, lwd=0.5)  ## (RH 210426) -- males
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
			lines(linesx_i0, linesy_i0, col=linescol[3], lwd=0.5, lty=lty)
			#lines(linesx_i1, linesy_i1, col=linescol[1], lwd=lwd, lty=lty)
			#lines(linesx_i2, linesy_i2, col=linescol[2], lwd=lwd, lty=lty)
			#lines(linesx_i1, linesy_i1, col=darkenRGB(colsex[1],0.5), lwd=lwd, lty=lty)
			#lines(linesx_i2, linesy_i2, col=darkenRGB(colsex[2],0.5), lwd=lwd, lty=lty)
			lines(linesx_i1, linesy_i1, col="red", lwd=lwd, lty=lty)
			lines(linesx_i2, linesy_i2, col="blue", lwd=lwd, lty=lty)
		}
		abline(h=0, col="grey", lwd=0.5) ;box()
		
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
						if (length(obsN)>1) {
							obs <- unique(obsN[sexvec == sex & yr == yr_i])
							text_i <- paste(paste0("N obs = ", round(obs, sampsizeround)), text_i, sep="\n")
						}
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
			if (legadjy[1] == "default") {
				adjy0 = 1.1
				adjy <- ifelse(i < 3, adjy0, adjy0 + adjy0 * (i - 2)) ## (RH 230913)
				#adjy <- ifelse(i < 3, 1.3, 1.3 + 1.3 * (i - 2))
			}
			else {
				adjy <- legadjy[i]
			}
			text(x=textx, y=texty, labels=linguaFranca(text_i,lang), adj=c(adjx, adjy), cex=legsize[i], font=legfont[i])
#if (i==2) {browser();return()}
			if (length(obsN)==1 && text_i2 != text_i && text_i2 != "") {
				text(x=textx, y=texty2, labels=linguaFranca(text_i,lang), adj=c(adjx, -adjy), cex=legsize[i], font=legfont[i])
			}
			if (twosex & !bub & venusmars) {
				pu <- par("usr")
				xval <- pu[2]
				if (length(ptsx_i0) > 0) {
					text(xval, 0.5 * yrange[2], "\\VE+\\MA", vfont=c("serif", "plain"), cex=2, col=darkenRGB(colsex[3],0.5), pos=2)
				}
				if (length(ptsx_i1) > 0) {
					#text(xval, 0.5 * yrange[2], "\\VE", vfont=c("serif", "plain"), cex=2, col=darkenRGB(colsex[1],0.5), pos=2)
					text(xval, 0.5 * yrange[2], "\\VE", vfont=c("serif", "plain"), cex=2, col="red", pos=2)
				}
				if (length(ptsx_i2) > 0) {
					#text(xval, -0.5 * yrange[2], "\\MA", vfont=c("serif", "plain"), cex=2, col=darkenRGB(colsex[2],0.5), pos=2)
					text(xval, -0.5 * yrange[2], "\\MA", vfont=c("serif", "plain"), cex=2, col="blue", pos=2)
				}
			}
		}
		mfg <- par("mfg")
		if( getNpan() %in% rev(panelrange)[1:mfg[4]] ) { ## RH 201207
		#if (mfg[1] == mfg[3] | ipanel == npanels) {
			axis(side=1, at=axis1, labels=linguaFranca(axis1labs,lang), las=xlas)
#browser();return()
		}
		if (mfg[2] == 1) {
			axis(side=2, at=axis2, las=ylas)
			if (twosex) {
				axis(side=2, at=-axis2[axis2 > 0], labels=linguaFranca(format(axis2[axis2 > 0]),lang), las=ylas)
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
				title(main=linguaFranca(main,lang), line=c(2, 0, 3, 3), outer=TRUE, cex.main=cex.main * fixcex)
				#title(xlab=xlab, outer=TRUE, cex.lab=fixcex)
				#title(ylab=ylab, line=ifelse(ylas %in% 1:2, max(3, 2 + 0.4 * maxchar_yaxs), 3.5), outer=TRUE, cex.lab=fixcex)
				mtext(text=linguaFranca(xlab,lang), side=1, line=2, outer=TRUE, cex=1.5)
				mtext(text=linguaFranca(ylab,lang), side=2, line=ifelse(ylas %in% 1:2, max(3, 2 + 0.4 * maxchar_yaxs), 3.5), outer=TRUE, cex=1.5)
			}
			else {
				#title(main=main, xlab=xlab, ylab=ylab, outer=FALSE, cex.main=cex.main)
				title(main=linguaFranca(main,lang), outer=FALSE, cex.main=cex.main)
				mtext(text=linguaFranca(xlab,lang), side=1, line=2, outer=FALSE, cex=1.5)
				mtext(text=linguaFranca(ylab,lang), side=2, line=3.5, outer=FALSE, cex=1.5)
#browser();return()
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


## plotSS.comparisons-------------------2025-10-15
##  Borrowed r4ss code 'SSplotComparisons' to plot retrospectives.
##  Creates a user-chosen set of plots comparing model output from a summary of multiple models,
##  where the collection was created using the SSsummarize function. (r4ss)
## ----------------------------------------r4ss|RH
plotSS.comparisons <- function (summaryoutput, subplots=1:20, plot=TRUE, print=FALSE, 
	png=print, pdf=FALSE, models="all", endyrvec=NULL, 
	indexfleets=NULL, indexUncertainty=TRUE, indexQlabel=TRUE, 
	indexQdigits=4, indexSEvec=NULL, indexPlotEach=FALSE, 
	labels=c("Year", "Spawning biomass (t)", "Fraction of unfished", 
		"Age-0 recruits (1,000s)", "Recruitment deviations", 
		"Index", "Log index", "SPR-related quantity", "Density", 
		"Management target", "Minimum stock size threshold", 
		"Spawning output", "Harvest rate"), col=NULL, shadecol=NULL, 
	pch=NULL, lty=1, lwd=2, spacepoints=10, staggerpoints=1, 
	initpoint=0, tickEndYr=FALSE, shadeForecast=TRUE, xlim=NULL, 
	ylimAdj=1.05, xaxs="i", yaxs="i", type="o", uncertainty=FALSE, 
	shadealpha=0.1, legend=TRUE, legendlabels=NULL, legendloc="topright", 
	legendorder=NULL, legendncol=1, sprtarg=NULL, btarg=NULL, 
	minbthresh=NULL, pwidth=8.5, pheight=6, punits="in", 
	res=400, ptsize=12, plotdir=NULL, filenameprefix="", 
	densitynames=c("SSB_Virgin", "R0"), densityxlabs=NULL, 
	rescale=TRUE, densityscalex=1, densityscaley=1, densityadjust=1, 
	densitysymbols=TRUE, densitytails=TRUE, densitymiddle=FALSE, 
	densitylwd=1, fix0=TRUE, new=TRUE, add=FALSE, 
	par=list(mar=c(5,4,1,1) + 0.1), verbose=TRUE, mcmcVec=FALSE, 
	show_equilibrium=TRUE, lang="e")
{
	meanRecWarning <- TRUE
	ymax_vec <- rep(NA, 17)
	save_png_comparisons <- function(file, lang="e") {  ## needs to do one language at a time
		plotdir = sub("/$","",plotdir)  ## check for and remove ending slash
		#createFdir(lang, dir=plotdir)
		file = paste0(filenameprefix, file)
		plotdir = switch(lang, 'e'=plotdir, 'f'=paste0(plotdir,"/french"))
		if (!dir.exists(plotdir))
			dir.create(plotdir)
		fout = file.path(plotdir, file)
		changeLangOpts(L=lang)
#browser();return()
		clearFiles(fout)
		png(filename=fout, width=pwidth, height=pheight, units=punits, res=res, pointsize=ptsize)
		par(par)
	}
	if (png) {
		print <- TRUE
	}
	if (png & is.null(plotdir)) {
		stop("to print PNG files, you must supply a directory as 'plotdir'")
	}
	if (pdf & png) {
		stop("To use 'pdf', set 'print' or 'png' to FALSE.")
	}
	if (pdf) {
		if (is.null(plotdir)) {
			stop("to write to a PDF, you must supply a directory as 'plotdir'")
		}
		pdffile <- file.path(plotdir, paste0(filenameprefix, 
			"SSplotComparisons_", format(Sys.time(), "%d-%b-%Y_%H.%M"), 
			".pdf"))
		pdf(file=pdffile, width=pwidth, height=pheight)
		if (verbose) {
			message("PDF file with plots will be:", pdffile)
		}
		par(par)
	}
	n <- summaryoutput[["n"]]
	nsexes <- summaryoutput[["nsexes"]]
	startyrs <- summaryoutput[["startyrs"]]
	endyrs <- summaryoutput[["endyrs"]]
	pars <- summaryoutput[["pars"]]
	parsSD <- summaryoutput[["parsSD"]]
	parphases <- summaryoutput[["parphases"]]
	quants <- summaryoutput[["quants"]]
	quantsSD <- summaryoutput[["quantsSD"]]
	SpawnBio <- summaryoutput[["SpawnBio"]]
	SpawnBioLower <- summaryoutput[["SpawnBioLower"]]
	SpawnBioUpper <- summaryoutput[["SpawnBioUpper"]]
	Bratio <- summaryoutput[["Bratio"]]
	BratioLower <- summaryoutput[["BratioLower"]]
	BratioUpper <- summaryoutput[["BratioUpper"]]
	SPRratio <- summaryoutput[["SPRratio"]]
	SPRratioLower <- summaryoutput[["SPRratioLower"]]
	SPRratioUpper <- summaryoutput[["SPRratioUpper"]]
	Fvalue <- summaryoutput[["Fvalue"]]
	FvalueLower <- summaryoutput[["FvalueLower"]]
	FvalueUpper <- summaryoutput[["FvalueUpper"]]
	recruits <- summaryoutput[["recruits"]]
	recruitsLower <- summaryoutput[["recruitsLower"]]
	recruitsUpper <- summaryoutput[["recruitsUpper"]]
	recdevs <- summaryoutput[["recdevs"]]
	recdevsLower <- summaryoutput[["recdevsLower"]]
	recdevsUpper <- summaryoutput[["recdevsUpper"]]
	indices <- summaryoutput[["indices"]]
	mcmc <- summaryoutput[["mcmc"]]
	lowerCI <- summaryoutput[["lowerCI"]]
	upperCI <- summaryoutput[["upperCI"]]
	SpawnOutputUnits <- summaryoutput[["SpawnOutputUnits"]]
	btargs <- summaryoutput[["btargs"]]
	minbthreshs <- summaryoutput[["minbthreshs"]]
	sprtargs <- summaryoutput[["sprtargs"]]
	SPRratioLabels <- summaryoutput[["SPRratioLabels"]]
	FvalueLabels <- summaryoutput[["FvalueLabels"]]

	if (legendloc %in% c("topright")) {
		xleg=0.95; yleg=0.95; xjust=1; yjust=1
	} else if (legendloc %in% c("bottomleft")) {
		xleg=0.05; yleg=0.05; xjust=0; yjust=0
	} else if (legendloc %in% c("topleft")) {
		xleg=0.05; yleg=0.80; xjust=0; yjust=1
		if (strSpp=="418" && subplot==3) { xleg=0.02; yleg=0.88 }
		if (strSpp=="405" && subplot==1) { xleg=0.05; yleg=0.75 }
		if (strSpp=="405" && subplot %in% c(3)) { xleg=0.65; yleg=0.95 }
		if (strSpp=="405" && subplot %in% c(13)) { xleg=0.02; yleg=0.975 }
	} ## add more if neeed
#browser();return()

	if (is.null(btarg)) {
		btarg <- unique(btargs)
		if (length(btarg) > 1) {
			warning("setting btarg=-999 because models don't have matching values")
			btarg <- -999
		}
	}
	if (is.null(minbthresh)) {
		minbthresh <- unique(minbthreshs)
		if (length(minbthresh) > 1) {
			warning("setting minbthresh=-999 because models don't have matching values")
			minbthresh <- -999
		}
	}
	if (is.null(sprtarg)) {
		sprtarg <- unique(sprtargs)
		if (length(sprtarg) > 1) {
			warning("setting sprtarg=-999 because models don't have matching values")
			sprtarg <- -999
		}
	}
	SPRratioLabel <- unique(SPRratioLabels)
	if (length(SPRratioLabel) > 1) {
		warning("setting label for SPR plot to 8th element of input 'labels' ", 
			"because the models don't have matching labels")
		SPRratioLabel <- labels[8]
	}
	FvalueLabel <- unique(FvalueLabels)  ## this makes no sense to me
	if (length(FvalueLabel) > 1) {
		warning("setting label for F plot to 13th element of input 'labels' ", 
			"because the models don't have matching labels")
		FvalueLabel <- labels[13]
	}
	else {
		FvalueLabel <- gsub("_", " ", FvalueLabel)
	}
	if (!is.logical(uncertainty) & is.numeric(uncertainty)) {
		if (any(!uncertainty %in% 1:n)) {
			stop("'uncertainty' should be a subset of the integers\n", 
				" 1-", n, ", where n=", n, " is the number of models.\n", 
				"  Or it can be a single TRUE/FALSE value.\n", 
				"  Or a vector of TRUE/FALSE, of length n=", 
				n)
		}
		else {
			uncertainty <- 1:n %in% uncertainty
		}
	}
	if (is.logical(uncertainty) & length(uncertainty) == 1) {
		uncertainty <- rep(uncertainty, n)
	}
	if (length(uncertainty) != n) {
		stop("'uncertainty' as TRUE/FALSE should have length 1 or n.\n", 
			"  length(uncertainty)=", length(uncertainty))
	}
	if (all(uncertainty)) {
		message("showing uncertainty for all models")
	}
	if (!any(uncertainty)) {
		message("not showing uncertainty for any models")
	}
	if (any(uncertainty) & !all(uncertainty)) {
		message("showing uncertainty for model", ifelse(sum(uncertainty) > 
			1, "s: ", " "), paste(which(uncertainty), collapse=","))
	}
	for (i in 1:n) {
		if (all(is.na(quantsSD[, i]) | quantsSD[, i] == 0)) {
			message("No uncertainty available for model ", i)
			uncertainty[i] <- FALSE
		}
	}
	if (length(unique(nsexes)) > 1) {
		warning("SSplotComparisons no longer divides SpawnBio by 2 for single-sex models\n", 
			"to get female-only spawning biomass output by SS for a single-sex model,\n", 
			"use the new Nsexes=-1 option in the data file.")
	}
	if (models[1] == "all") {
		models <- 1:n
	}
	nlines <- length(models)
	if (any(mcmcVec) & length(mcmc) == 0) {
		mcmcVec <- FALSE
		warning("Setting mcmcVec=FALSE because summaryoutput[['mcmc']] is empty")
	}
	if (nlines > 1 & length(mcmcVec) == 1) {
		mcmcVec <- rep(mcmcVec, nlines)
	}
	if (nlines != length(mcmcVec)) {
		stop("Input 'mcmcVec' must equal 1 or the number of models.\n")
	}
	if (any(subplots %in% 13:14) & !is.null(indices) && nrow(indices) > 0) {
		if (is.null(indexfleets)) {
			indexfleets <- list()
			for (imodel in 1:n) {
				indexfleets[[paste0("model", imodel)]] <- sort(unique(indices[["Fleet"]][indices[["imodel"]] == imodel]))
			}
		}
		else {
			if (!is.null(indexfleets)) {
				if (is.vector(indexfleets) & length(indexfleets) == 1) {
					indexfleets <- rep(indexfleets, n)
				}
				if (length(indexfleets) != n) {
					warning("Skipping index plots: length(indexfleets) should be 1 or n=", n, ".")
					indexfleets <- NULL
				}
			}
		}
		if (!length(unique(lapply(indexfleets, FUN=length))) == 1) {
			message("indexfleets:")
			print(indexfleets)
			warning("Skipping index plots:\n", "Fleets have different numbers of indices listed in 'indexfleets'.")
			indexfleets <- NULL
		}
		index_plot_suffix <- rep("", length(indexfleets))
		if (length(indexfleets[[1]]) > 1) {
			for (iindex in 1:length(indexfleets[[1]])) {
				fleets <- as.numeric(data.frame(indexfleets)[iindex, ])
				if (length(unique(fleets)) == 1) {
					index_plot_suffix[iindex] <- paste0("_flt", fleets[1])
				}
				else {
					index_plot_suffix[iindex] <- paste0("_index", iindex)
				}
			}
		}
	}
	if (is.null(col) & nlines > 3) 
		col <- rich.colors.short(nlines)  ## want year0 to be black
		#col <- rich.colors.short(nlines + 1)[-1]
	if (is.null(col) & nlines < 3) 
		col <- rich.colors.short(nlines)
	if (is.null(col) & nlines == 3) 
		col <- c("blue", "red", "green3")
	if (is.null(shadecol)) {
		shadecol <- adjustcolor(col, alpha.f=shadealpha)
	}
	if (is.null(pch)) {
		pch <- rep(1:25, 10)[1:nlines]
	}
	if (length(col) < nlines) 
		col <- rep(col, nlines)[1:nlines]
	if (length(pch) < nlines) 
		pch <- rep(pch, nlines)[1:nlines]
	if (length(lty) < nlines) 
		lty <- rep(lty, nlines)[1:nlines]
	if (length(lwd) < nlines) 
		lwd <- rep(lwd, nlines)[1:nlines]
	if (!is.expression(legendlabels[1]) && is.null(legendlabels)) {
		legendlabels <- paste("model", 1:nlines)
	}
	if (plot & new & !pdf) {
		dev.new(width=pwidth, height=pheight, pointsize=ptsize, record=TRUE)
		par(par)
	}
	for (iline in (1:nlines)[mcmcVec]) {
		imodel <- models[iline]
		cols <- imodel
		SpawnBioLower[, cols] <- SpawnBioUpper[, cols] <- SpawnBio[, 
			cols] <- NA
		BratioLower[, cols] <- BratioUpper[, cols] <- Bratio[, 
			cols] <- NA
		SPRratioLower[, cols] <- SPRratioUpper[, cols] <- SPRratio[, 
			cols] <- NA
		recruitsLower[, cols] <- recruitsUpper[, cols] <- recruits[, 
			cols] <- NA
		recdevsLower[, cols] <- recdevsUpper[, cols] <- recdevs[, 
			cols] <- NA
		tmp <- grep("SSB", names(mcmc[[imodel]]))
		tmp2 <- c(grep("SSB_unfished", names(mcmc[[imodel]]), 
			ignore.case=TRUE), grep("SSB_Btgt", names(mcmc[[imodel]]), 
			ignore.case=TRUE), grep("SSB_SPRtgt", names(mcmc[[imodel]]), 
			ignore.case=TRUE), grep("SSB_MSY", names(mcmc[[imodel]]), 
			ignore.case=TRUE))
		tmp <- setdiff(tmp, tmp2)
		if (length(tmp) > 0) {
			mcmc.tmp <- mcmc[[imodel]][, tmp]
			mcmclabs <- names(mcmc.tmp)
			lower <- apply(mcmc.tmp, 2, quantile, prob=lowerCI, 
				na.rm=TRUE)
			med <- apply(mcmc.tmp, 2, quantile, prob=0.5, na.rm=TRUE)
			upper <- apply(mcmc.tmp, 2, quantile, prob=upperCI, 
				na.rm=TRUE)
			SpawnBio[, imodel] <- med[match(SpawnBio[["Label"]], 
				mcmclabs)]
			SpawnBioLower[, imodel] <- lower[match(SpawnBioLower[["Label"]], 
				mcmclabs)]
			SpawnBioUpper[, imodel] <- upper[match(SpawnBioUpper[["Label"]], 
				mcmclabs)]
		}
		tmp <- grep("Bratio", names(mcmc[[imodel]]))
		if (length(tmp) > 0) {
			mcmc.tmp <- mcmc[[imodel]][, tmp]
			mcmclabs <- names(mcmc.tmp)
			lower <- apply(mcmc.tmp, 2, quantile, prob=lowerCI, 
				na.rm=TRUE)
			med <- apply(mcmc.tmp, 2, quantile, prob=0.5, na.rm=TRUE)
			upper <- apply(mcmc.tmp, 2, quantile, prob=upperCI, 
				na.rm=TRUE)
			Bratio[, imodel] <- med[match(Bratio[["Label"]], 
				mcmclabs)]
			BratioLower[, imodel] <- lower[match(BratioLower[["Label"]], 
				mcmclabs)]
			BratioUpper[, imodel] <- upper[match(BratioUpper[["Label"]], 
				mcmclabs)]
		}
		tmp <- grep("SPRratio", names(mcmc[[imodel]]))
		if (length(tmp) > 0) {
			mcmc.tmp <- mcmc[[imodel]][, tmp]
			mcmclabs <- names(mcmc.tmp)
			lower <- apply(mcmc.tmp, 2, quantile, prob=lowerCI, 
				na.rm=TRUE)
			med <- apply(mcmc.tmp, 2, quantile, prob=0.5, na.rm=TRUE)
			upper <- apply(mcmc.tmp, 2, quantile, prob=upperCI, 
				na.rm=TRUE)
			SPRratio[, imodel] <- med[match(SPRratio[["Label"]], 
				mcmclabs)]
			SPRratioLower[, imodel] <- lower[match(SPRratioLower[["Label"]], 
				mcmclabs)]
			SPRratioUpper[, imodel] <- upper[match(SPRratioUpper[["Label"]], 
				mcmclabs)]
		}
		tmp <- grep("^Recr_", names(mcmc[[imodel]]))
		tmp2 <- grep("Recr_unfished", names(mcmc[[imodel]]), 
			ignore.case=TRUE)
		tmp <- setdiff(tmp, tmp2)
		if (length(tmp) > 0) {
			mcmc.tmp <- mcmc[[imodel]][, tmp]
			mcmclabs <- names(mcmc.tmp)
			lower <- apply(mcmc.tmp, 2, quantile, prob=lowerCI, 
				na.rm=TRUE)
			med <- apply(mcmc.tmp, 2, quantile, prob=0.5, na.rm=TRUE)
			mean <- apply(mcmc.tmp, 2, mean, na.rm=TRUE)
			upper <- apply(mcmc.tmp, 2, quantile, prob=upperCI, 
				na.rm=TRUE)
			if (!meanRecWarning) {
				message("note: using mean recruitment from MCMC instead of median,\n", 
				  "because it is more comparable to MLE\n")
				meanRecWarning <- TRUE
			}
			recruits[, imodel] <- mean[match(recruits[["Label"]], 
				mcmclabs)]
			recruitsLower[, imodel] <- lower[match(recruitsLower[["Label"]], 
				mcmclabs)]
			recruitsUpper[, imodel] <- upper[match(recruitsUpper[["Label"]], 
				mcmclabs)]
		}
		tmp <- unique(c(grep("_RecrDev_", names(mcmc[[imodel]])), 
			grep("_InitAge_", names(mcmc[[imodel]])), grep("ForeRecr_", 
				names(mcmc[[imodel]]))))
		if (length(tmp) > 0) {
			mcmc.tmp <- mcmc[[imodel]][, tmp]
			mcmclabs <- names(mcmc.tmp)
			lower <- apply(mcmc.tmp, 2, quantile, prob=lowerCI, 
				na.rm=TRUE)
			med <- apply(mcmc.tmp, 2, quantile, prob=0.5, na.rm=TRUE)
			upper <- apply(mcmc.tmp, 2, quantile, prob=upperCI, 
				na.rm=TRUE)
			recdevs[, imodel] <- med[match(recdevs[["Label"]], 
				mcmclabs)]
			recdevsLower[, imodel] <- lower[match(recdevsLower[["Label"]], 
				mcmclabs)]
			recdevsUpper[, imodel] <- upper[match(recdevsUpper[["Label"]], 
				mcmclabs)]
		}
	}
	if (is.null(endyrvec)) {
		endyrvec <- endyrs + 1
	}
	if (length(endyrvec) == 1) {
		endyrvec <- rep(endyrvec, nlines)
	}
	recdevs <- recdevs[!is.na(recdevs[["Yr"]]), ]
	recdevsLower <- recdevsLower[!is.na(recdevsLower[["Yr"]]), 
		]
	recdevsUpper <- recdevsUpper[!is.na(recdevsUpper[["Yr"]]), 
		]
	if (!is.null(endyrvec)) {
		for (iline in 1:nlines) {
			endyr <- endyrvec[iline]
			imodel <- models[iline]
			SpawnBio[SpawnBio[["Yr"]] > endyr, imodel] <- NA
			SpawnBioLower[SpawnBio[["Yr"]] > endyr, imodel] <- NA
			SpawnBioUpper[SpawnBio[["Yr"]] > endyr, imodel] <- NA
			Bratio[Bratio[["Yr"]] > endyr, imodel] <- NA
			BratioLower[Bratio[["Yr"]] > endyr, imodel] <- NA
			BratioUpper[Bratio[["Yr"]] > endyr, imodel] <- NA
			SPRratio[SPRratio[["Yr"]] >= endyr, imodel] <- NA
			SPRratioLower[SPRratio[["Yr"]] >= endyr, imodel] <- NA
			SPRratioUpper[SPRratio[["Yr"]] >= endyr, imodel] <- NA
			Fvalue[Fvalue[["Yr"]] >= endyr, imodel] <- NA
			FvalueLower[Fvalue[["Yr"]] >= endyr, imodel] <- NA
			FvalueUpper[Fvalue[["Yr"]] >= endyr, imodel] <- NA
			recruits[recruits[["Yr"]] > endyr, imodel] <- NA
			recruitsLower[recruits[["Yr"]] > endyr, imodel] <- NA
			recruitsUpper[recruits[["Yr"]] > endyr, imodel] <- NA
			if (!is.null(recdevs)) {
				recdevs[recdevs[["Yr"]] > endyr, imodel] <- NA
				recdevsLower[recdevs[["Yr"]] > endyr, imodel] <- NA
				recdevsUpper[recdevs[["Yr"]] > endyr, imodel] <- NA
			}
		}
	}
	addpoly <- function(yrvec, lower, upper) {
		lower[lower < 0] <- 0
		for (iline in (1:nlines)[uncertainty]) {
			imodel <- models[iline]
			good <- !is.na(lower[, imodel]) & !is.na(upper[, 
				imodel])
			polygon(x=c(yrvec[good], rev(yrvec[good])), y=c(lower[good, 
				imodel], rev(upper[good, imodel])), border=NA, 
				col=shadecol[iline])
		}
	}
	plotSpawnBio <- function(show_uncertainty=TRUE, lang="e") {
		if (!any(uncertainty)) {
			show_uncertainty <- FALSE
		}
		if (is.null(xlim)) {
			if (show_equilibrium) {
				xlim <- range(SpawnBio[["Yr"]])
			}
			else {
				xlim <- range(SpawnBio[["Yr"]][-c(1, 2)])
			}
			if (!is.null(endyrvec) & all(endyrvec < max(xlim))) {
				xlim[2] <- max(endyrvec)
			}
		}
		ylim <- ylimAdj * range(0, SpawnBio[SpawnBio[["Yr"]] >= xlim[1] & SpawnBio[["Yr"]] <= xlim[2], models], na.rm=TRUE)
		if (show_uncertainty) {
			ylim <- range(ylim, ylimAdj * SpawnBioUpper[SpawnBio[["Yr"]] >= xlim[1] & SpawnBio[["Yr"]] <= xlim[2], models[uncertainty]], na.rm=TRUE)
		}
		xlim = xlim + c(-2,2)
		if (length(unique(SpawnOutputUnits)) != 1) {
			warning("Some models may have different units", " for spawning output than others")
		}
		if (any(SpawnOutputUnits == "numbers")) {  ## this seems incorrect
			ylab <- labels[12]
		}
		else {
			ylab <- labels[2]
		}
		ylab = labels[2]  ## ad hoc fix (RH 230626)
		yunits <- 1
		if (rescale & ylim[2] > 1000 & ylim[2] < 1e+06) {
			yunits <- 1000
			ylab <- gsub("(t)", "(x1000 t)", ylab, fixed=TRUE)
			ylab <- gsub("eggs", "x1000 eggs", ylab, fixed=TRUE)
		}
		if (rescale & ylim[2] > 1e+06) {
			yunits <- 1e+06
			ylab <- gsub("(t)", "(million t)", ylab, fixed=TRUE)
			ylab <- gsub("eggs", "millions of eggs", ylab, fixed=TRUE)
		}
		if (rescale & ylim[2] > 1e+09) {
			yunits <- 1e+09
			ylab <- gsub("million", "billion", ylab, fixed=TRUE)
		}
		if (!add) {
			plot(0, type="n", xlim=xlim, ylim=ylim, xlab=linguaFranca(labels[1],lang), ylab=linguaFranca(ylab,lang), xaxs=xaxs, yaxs=yaxs, axes=FALSE)
			#abline(v=endyrvec,col="gainsboro")
			#text(endyrvec,0,endyrvec,srt=90,cex=0.8,adj=c(0,0))
		}
		if (legend) {
			x0 = rep(SpawnBio[["Yr"]][2], nlines)
			y0 = as.numeric(SpawnBio[3, models])
			points(x0, y0, col=col, pch=pch, cex=1, lwd=2)
#browser();return()
			ylast = as.vector(sapply(SpawnBio[-(1:2), models],function(x){rev(x[!is.na(x)])[1]}))
			segments(x0=endyrvec, y0=0, x1=endyrvec, y1=ylast, lty=1, col=lucent("black",0.25), lwd=0.5)
			ylow =  par()$usr[3] + 0.01 * diff(par()$usr[3:4])
			show = seq(1,length(endyrvec),2)
			text(endyrvec[show], rep(ylow,length(show)), endyrvec[show], col=col[show], font=2, srt=90, cex=0.6, adj=c(0,0))
			addLegend(xleg, yleg, pch=pch, col=col, legend=endyrvec, cex=0.8, lwd=2, xjust=xjust, yjust=yjust, bty="n")
			#add_legend(legendlabels, legendloc=legendloc, legendorder=legendorder, legendncol=legendncol, col=col, pch=pch, lwd=lwd, lty=lty)
		}
		if (show_uncertainty) {
			addpoly(yrvec=SpawnBio[["Yr"]][-(1:2)], lower=SpawnBioLower[-(1:2),], upper=SpawnBioUpper[-(1:2),])
			xEqu <- SpawnBio[["Yr"]][2] - (1:nlines)/nlines
		}
		else {
			xEqu <- rep(SpawnBio[["Yr"]][2], nlines)
		}
		if (spacepoints %in% c(0, 1, FALSE)) {
			matplot(SpawnBio[["Yr"]][-(1:2)], SpawnBio[-(1:2), models], col=col, pch=pch, lty=lty, lwd=lwd, type=type, ylim=ylim, add=TRUE)
		}
		else {
			matplot(SpawnBio[["Yr"]][-(1:2)], SpawnBio[-(1:2), models], col=col, lty=lty, lwd=lwd, type="l", ylim=ylim, add=TRUE)
			#SpawnBio2 <- SpawnBio
			#for (iline in 1:nlines) {
			#	imodel <- models[iline]
			#	SpawnBio2[(SpawnBio2[["Yr"]] - initpoint)%%spacepoints != (staggerpoints * iline)%%spacepoints, imodel] <- NA
			#}
			#matplot(SpawnBio2[["Yr"]][-(1:2)], SpawnBio2[-(1:2), models], col=col, pch=pch, lwd=lwd, type="p", ylim=ylim, add=TRUE)
		}
		if (show_equilibrium) {
			old_warn <- options()$warn
			options(warn=-1)
			if (show_uncertainty) {
				arrows(x0=xEqu[models[uncertainty]], y0=as.numeric(SpawnBioLower[1, 
					models[uncertainty]]), x1=xEqu[models[uncertainty]], 
					y1=as.numeric(SpawnBioUpper[1, models[uncertainty]]), 
					length=0.01, angle=90, code=3, col=col[uncertainty], lwd=2)
			}
			options(warn=old_warn)
			points(x=xEqu, SpawnBio[1, models], col=col, pch=pch, cex=1, lwd=2)
		}
		if (!add) {
			abline(h=0, col="grey")
			if (tickEndYr) {
				ticks <- graphics::axTicks(1)
				axis(1, at=c(ticks[ticks < max(endyrvec)], max(endyrvec)))
			}
			else {
				allyrs = SpawnBio$Yr[1]:endyrvec[1]
				axis(1, at=intersect(allyrs, seq(1900,2050,5)), tcl=-0.2, labels=F)
				axis(1, at=intersect(allyrs, seq(1900,2050,10)), tcl=-0.4, labels=T)
			}
			if (!is.null(endyrvec) & max(endyrvec) > 1 + max(endyrs) & shadeForecast) {
				rect(xleft=max(endyrs) + 1, ybottom=par()$usr[3], xright=par()$usr[2], ytop=par()$usr[4], col=gray(0, alpha=0.1), border=NA)
			}
			yticks <- pretty(ylim,n=16)
			ytklab = yticks[seq(1,length(yticks),2)]
			axis(2, at=yticks, labels=FALSE, tcl=-0.2)
			axis(2, at=ytklab, labels=format(ytklab/yunits), tcl=-0.4, las=1)
		}
		box()
#browser();return()
		return(ylim[2])
	}
	plotBratio <- function(show_uncertainty=TRUE, lang="e") {
		## RH redefine Bratio as Bt/B0
		BtB0 = SpawnBio[-(1:2),]
		BtB0$Label = sub("SSB","BtB0",BtB0$Label)
		B0 = as.numeric(BtB0[1,models])
		Bratio = BtB0
		Bratio[,models] = sweep(BtB0[,models],2,B0,"/")

		if (!any(uncertainty)) {
			show_uncertainty <- FALSE
		}
		if (is.null(xlim)) {
			xlim <- range(Bratio[["Yr"]])
			if (!is.null(endyrvec) & all(endyrvec < max(xlim))) {
				xlim[2] <- max(endyrvec)
			}
		}
		xlim = xlim + c(-2,2)

		ylim <- ylimAdj * range(0, Bratio[Bratio[["Yr"]] >= xlim[1] & Bratio[["Yr"]] <= xlim[2], models], na.rm=TRUE)
		if (show_uncertainty) {
			ylim <- ylimAdj * range(ylim/ylimAdj, BratioUpper[Bratio[["Yr"]] >= xlim[1] & Bratio[["Yr"]] <= xlim[2], models[uncertainty]], na.rm=TRUE)
		}
		if (!add) {
			#plot(0, type="n", xlim=xlim, ylim=ylim, xlab=labels[1], ylab=labels[3], xaxs=xaxs, yaxs=yaxs, axes=FALSE)
			plot(0, type="n", xlim=xlim, ylim=ylim, xlab=linguaFranca(labels[1],lang), ylab=expression(italic(B)[italic(t)]/italic(B)[0]), xaxs=xaxs, yaxs=yaxs, axes=FALSE)
		}
		if (legend) {
#browser();return()
			ylast = as.vector(sapply(Bratio[, models],function(x){rev(x[!is.na(x)])[1]}))
			points(endyrvec, ylast, col=col, pch=pch, cex=1.2, lwd=2)
			addLegend(xleg, yleg, pch=pch, col=col, legend=endyrvec, cex=0.8, lwd=2, xjust=xjust, yjust=yjust, bty="n")
		}
		if (show_uncertainty) {
			addpoly(Bratio[["Yr"]], lower=BratioLower, upper=BratioUpper)
		}
		if (spacepoints %in% c(0, 1, FALSE)) {
			matplot(Bratio[["Yr"]], Bratio[, models], col=col, pch=pch, lty=lty, lwd=lwd, type=type, ylim=ylim, add=TRUE)
		}
		else {
			matplot(Bratio[["Yr"]], Bratio[, models], col=col, pch=pch, lty=lty, lwd=lwd, type="l", ylim=ylim, add=TRUE)
			if (type != "l") {
				Bratio2 <- Bratio
				for (iline in 1:nlines) {
					imodel <- models[iline]
					Bratio2[(Bratio2[["Yr"]] - initpoint)%%spacepoints != (staggerpoints * iline)%%spacepoints, imodel] <- NA
				}
				#matplot(Bratio2[["Yr"]], Bratio2[, models], col=col, pch=pch, lty=lty, lwd=lwd, type="p", ylim=ylim, add=TRUE)
			}
		}
		yticks <- pretty(par()$yaxp[1:2])
		btarg = c(0.4, 0.32, 0.16)
		if (any(btarg > 0)) {
			abline(h=btarg, col=c("green3","blue","red"), lty=2)
			#text(min(Bratio[["Yr"]]) + 4, btarg + 0.03, labels[10], adj=0)
			text(min(Bratio[["Yr"]]) + 4, btarg + 0.02, paste0(btarg,"B0 ~ ",linguaFranca(c("TRP","USR","LRP"),lang)), adj=0)
			yticks <- sort(c(btarg, yticks))
		}
		if (minbthresh > 0) {
			abline(h=minbthresh, col="red", lty=2)
			text(min(Bratio[["Yr"]]) + 4, minbthresh + 0.03, 
				labels[11], adj=0)
			yticks <- sort(c(minbthresh, yticks))
		}
		if (!add) {
			#abline(h=0, col="grey")
			abline(h=1, col="grey", lty=2)
			if (tickEndYr) {
				ticks <- graphics::axTicks(1)
				axis(1, at=c(ticks[ticks < max(endyrvec)],  max(endyrvec)))
			}
			else {
				allyrs = Bratio$Yr[1]:endyrvec[1]
				axis(1, at=intersect(allyrs, seq(1900,2050,5)), tcl=-0.2, labels=F)
				axis(1, at=intersect(allyrs, seq(1900,2050,10)), tcl=-0.4, labels=T)
			}
#browser();return()
			if (!is.null(endyrvec) & max(endyrvec) > 1 + max(endyrs) & 
				shadeForecast) {
				rect(xleft=max(endyrs) + 1, ybottom=par()$usr[3], 
				  xright=par()$usr[2], ytop=par()$usr[4], 
				  col=gray(0, alpha=0.1), border=NA)
			}
			axis(2, at=seq(0,1,0.1), tcl=-0.2, labels=F)
			axis(2, at=yticks, tcl=-0.4, las=1)
			box()
		}
		#if (legend) {
		#	add_legend(legendlabels, legendloc=legendloc, legendorder=legendorder, legendncol=legendncol, col=col, pch=pch, lwd=lwd, lty=lty)
		#}
#browser();return()
		return(ylim[2])
	}
	plotSPRratio <- function(show_uncertainty=TRUE, lang="e") {
		if (!any(uncertainty)) {
			show_uncertainty <- FALSE
		}
		if (is.null(xlim)) {
			xlim <- range(SPRratio[["Yr"]])
			if (!is.null(endyrvec) & all(endyrvec < max(xlim))) {
				xlim[2] <- max(endyrvec)
			}
		}
		ylim <- ylimAdj * range(0, SPRratio[SPRratio[["Yr"]] >= 
			xlim[1] & SPRratio[["Yr"]] <= xlim[2], models], na.rm=TRUE)
		if (show_uncertainty) {
			ylim <- ylimAdj * range(ylim/ylimAdj, SPRratioUpper[SPRratio[["Yr"]] >= 
				xlim[1] & SPRratio[["Yr"]] <= xlim[2], models[uncertainty]], 
				na.rm=TRUE)
		}
		par(par)
		if (!add) {
			if (isTRUE(!is.na(SPRratioLabel) && SPRratioLabel == 
				paste0("(1-SPR)/(1-SPR_", floor(100 * sprtarg), 
				  "%)"))) {
				newmar <- oldmar <- par()$mar
				newmar[4] <- newmar[2]
				par(mar=newmar)
			}
			plot(0, type="n", xlim=xlim, ylim=ylim, xlab=labels[1], 
				ylab="", xaxs=xaxs, yaxs=yaxs, las=1, 
				axes=FALSE)
			axis(2)
		}
		if (show_uncertainty) {
			addpoly(SPRratio[["Yr"]], lower=SPRratioLower, 
				upper=SPRratioUpper)
		}
		if (spacepoints %in% c(0, 1, FALSE)) {
			matplot(SPRratio[["Yr"]], SPRratio[, models], col=col, 
				pch=pch, lty=lty, lwd=lwd, type=type, 
				ylim=ylim, add=TRUE)
		}
		else {
			matplot(SPRratio[["Yr"]], SPRratio[, models], col=col, 
				pch=pch, lty=lty, lwd=lwd, type="l", 
				ylim=ylim, add=TRUE)
			if (type != "l") {
				SPRratio2 <- SPRratio
				for (iline in 1:nlines) {
				  imodel <- models[iline]
				  SPRratio2[(SPRratio2[["Yr"]] - initpoint)%%spacepoints != 
					(staggerpoints * iline)%%spacepoints, imodel] <- NA
				}
				matplot(SPRratio2[["Yr"]], SPRratio2[, models], 
				  col=col, pch=pch, lty=lty, lwd=lwd, 
				  type="p", ylim=ylim, add=TRUE)
			}
		}
		abline(h=0, col="grey")
		if (sprtarg > 0) {
			if (isTRUE(SPRratioLabel == "1-SPR")) {
				abline(h=sprtarg, col="red", lty=2)
				text(SPRratio[["Yr"]][1] + 4, (sprtarg + 0.03), 
				  labels[10], adj=0)
				mtext(side=2, text=SPRratioLabel, line=par()$mgp[1], 
				  col=par()$col.lab, cex=par()$cex.lab)
			}
			else {
				yticks <- pretty(ylim)
				if (isTRUE(!is.na(SPRratioLabel) && SPRratioLabel == 
				  paste0("(1-SPR)/(1-SPR_", floor(100 * sprtarg), 
					"%)"))) {
				  abline(h=1, col="red", lty=2)
				  text(SPRratio[["Yr"]][1] + 4, 1 + 0.03, labels[10], 
					adj=0)
				  axis(4, at=yticks, labels=yticks * (1 - 
					sprtarg), las=1)
				  mtext(side=4, text="1 - SPR", line=par()$mgp[1], 
					col=par()$col.lab, cex=par()$cex.lab)
				  mtext(side=2, text=paste("(1-SPR)/(1-SPR_", 
					100 * sprtarg, "%)", sep=""), line=par()$mgp[1], 
					col=par()$col.lab, cex=par()$cex.lab)
				}
				else {
				  message("No line added to SPR ratio plot, ", 
					"as the settings used in this model ", "have not yet been configured in SSplotComparisons.")
				  mtext(side=2, text=SPRratioLabel, line=par()$mgp[1], 
					col=par()$col.lab, cex=par()$cex.lab)
				}
			}
		}
		else {
			mtext(side=2, text=SPRratioLabel, line=par()$mgp[1], 
				col=par()$col.lab, cex=par()$cex.lab)
		}
		if (!add) {
			if (tickEndYr) {
				ticks <- graphics::axTicks(1)
				axis(1, at=c(ticks[ticks < max(endyrvec)],  max(endyrvec)))
			}
			else {
				axis(1)
			}
			if (!is.null(endyrvec) & max(endyrvec) > 1 + max(endyrs) & 
				shadeForecast) {
				rect(xleft=max(endyrs) + 1, ybottom=par()$usr[3], 
				  xright=par()$usr[2], ytop=par()$usr[4], 
				  col=gray(0, alpha=0.1), border=NA)
			}
		}
		if (legend) {
			add_legend(legendlabels, legendloc=legendloc, legendorder=legendorder, 
				legendncol=legendncol, col=col, pch=pch, 
				lwd=lwd, lty=lty)
		}
		box()
		if (exists("oldmar")) {
			par(mar=oldmar)
		}
		return(ylim[2])
	}
	plotF <- function(show_uncertainty=TRUE, F2u=TRUE, lang="e") {
		if (!any(uncertainty)) {
			show_uncertainty <- FALSE
		}
		if (F2u && all(substring(Fvalue$Label,1,1)=="F") ) { ## convert fishing mortality to harvest rates
			Fvalue[,models] = 1 - exp(-Fvalue[,models])
			FvalueLabel = "Harvest rate"
		} else {
			FvalueLabel = "Fishing mortality"
		}
		if (is.null(xlim)) {
			xlim <- range(Fvalue[["Yr"]])
			if (!is.null(endyrvec) & all(endyrvec < max(xlim))) {
				xlim[2] <- max(endyrvec)
			}
		}
		ylim <- ylimAdj * range(0, Fvalue[Fvalue[["Yr"]] >= xlim[1] & Fvalue[["Yr"]] <= xlim[2], models], na.rm=TRUE)
		if (show_uncertainty) {
			ylim <- ylimAdj * range(ylim/ylimAdj, FvalueUpper[Fvalue[["Yr"]] >= xlim[1] & Fvalue[["Yr"]] <= xlim[2], models[uncertainty]], na.rm=TRUE)
		}
		par(par)
		if (!add) {
			plot(0, type="n", xlim=xlim, ylim=ylim, xlab=linguaFranca(labels[1],lang), ylab="", xaxs=xaxs, yaxs=yaxs, las=1, axes=FALSE)
			if (tickEndYr) {
				ticks <- graphics::axTicks(1)
				axis(1, at=c(ticks[ticks < max(endyrvec)],  max(endyrvec)))
			}
			else {
				axis(1)
			}
			axis(2)
		}
		if (show_uncertainty) {
			addpoly(Fvalue[["Yr"]], lower=FvalueLower, upper=FvalueUpper)
		}
		if (spacepoints %in% c(0, 1, FALSE)) {
			matplot(Fvalue[["Yr"]], Fvalue[, models], col=col, pch=pch, lty=lty, lwd=lwd, type=type, ylim=ylim, add=TRUE)
		}
		else {
			matplot(Fvalue[["Yr"]], Fvalue[, models], col=col, pch=pch, lty=lty, lwd=lwd, type="l", ylim=ylim, add=TRUE)
#browser();return()
			if (type != "l") {
				Fvalue2 <- Fvalue
				for (iline in 1:nlines) {
					imodel <- models[iline]
					Fvalue2[Fvalue2[["Yr"]]%%spacepoints != (staggerpoints * iline)%%spacepoints, imodel] <- NA
				}
				matplot(Fvalue2[["Yr"]], Fvalue2[, models], col=col, pch=pch, lty=lty, lwd=lwd, type="p", ylim=ylim, add=TRUE)
			}
		}
		abline(h=0, col="grey")
		mtext(side=2, text=linguaFranca(FvalueLabel,lang), line=par()$mgp[1], col=par()$col.lab, cex=par()$cex.lab)
		box()
		if (legend) {
			addLegend(xleg, yleg, col=col, legend=linguaFranca(legendlabels,lang), cex=0.8, lwd=3, xjust=0, yjust=1, bty="n")
			#add_legend(legendlabels, legendloc=legendloc, legendorder=legendorder, legendncol=legendncol, col=col, pch=pch, lwd=lwd, lty=lty)
		}
		return(ylim[2])
	} ## end plotF
	plotRecruits <- function(show_uncertainty=TRUE, recruit_lines=TRUE, lang="e") {
		if (!any(uncertainty)) {
			show_uncertainty <- FALSE
		}
		if (is.null(xlim)) {
			if (show_equilibrium) {
				xlim <- range(recruits[["Yr"]])
			}
			else {
				xlim <- range(recruits[["Yr"]][-c(1, 2)])
			}
			if (!is.null(endyrvec) & all(endyrvec < max(xlim))) {
				xlim[2] <- max(endyrvec)
			}
		}
		xlim = xlim + c(-2,2)

		ylim <- ylimAdj * range(0, recruits[recruits[["Yr"]] >= xlim[1] & recruits[["Yr"]] <= xlim[2], models], na.rm=TRUE)
		if (show_uncertainty) {
			ylim <- ylimAdj * range(ylim/ylimAdj, recruitsUpper[recruits[["Yr"]] >= 
				xlim[1] & recruits[["Yr"]] <= xlim[2], models[uncertainty]], 
				na.rm=TRUE)
		}
		ylab <- labels[4]
		yunits <- 1
		if (ylim[2] > 1000 & ylim[2] < 1e+06) {
			yunits <- 1000
			ylab <- gsub("1,000s", "millions", ylab)
		}
		if (ylim[2] > 1e+06) {
			yunits <- 1e+06
			ylab <- gsub("1,000s", "billions", ylab)
		}
		if (spacepoints %in% c(0, 1, FALSE)) {
			matplot(recruits[["Yr"]][-(1:2)], recruits[-(1:2), models], col=col, pch=pch, lty=lty, lwd=lwd, type=type, xlim=xlim, ylim=ylim, xlab=labels[1], ylab=ylab, xaxs=xaxs, yaxs=yaxs, axes=FALSE, add=add)
		}
		else {
			matplot(recruits[["Yr"]][-(1:2)], recruits[-(1:2), models], col=col, pch=pch, lty=lty, lwd=lwd, type="l", xlim=xlim, ylim=ylim, xlab=linguaFranca(labels[1],lang), ylab=linguaFranca(ylab,lang), xaxs=xaxs, yaxs=yaxs, axes=FALSE, add=add)
			if (type != "l") {
				recruits2 <- recruits
				for (iline in 1:nlines) {
					imodel <- models[iline]
					recruits2[(recruits2[["Yr"]]%%spacepoints - initpoint) != (staggerpoints * iline)%%spacepoints, imodel] <- NA
				}
				matplot(recruits2[["Yr"]][-(1:2)], recruits2[-(1:2), models], col=col, pch=pch, lty=lty, lwd=lwd, type="p", xlab=labels[1], ylab=ylab, xaxs=xaxs, yaxs=yaxs, axes=FALSE, ylim=ylim, add=TRUE)
			}
		}
		if (legend) {
			#ylast = as.vector(sapply(recruits[, models],function(x){rev(x[!is.na(x)])[1]}))
			#points(endyrvec, ymax, col=col, pch=pch, cex=1.2, lwd=2)
			zmax = as.vector(sapply(recruits[, models],function(x){ findPV(max(x,na.rm=T),x[!is.na(x)]) }))
			xmax = recruits$Yr[zmax]
			ymax = diag(as.matrix(recruits[zmax,models]))
			segments(xmax-2, ymax, xmax+2, ymax, col="grey10", lty=3, lwd=1)
			points(xmax, ymax, col=col, pch=pch, cex=1.2, lwd=2)
			addLegend(xleg, yleg, pch=pch, col=col, legend=endyrvec, cex=0.8, lwd=2, xjust=0, yjust=1, bty="n")
		}
		if (show_uncertainty) {
			xEqu <- recruits[["Yr"]][2] - (1:nlines)/nlines
		}
		else {
			xEqu <- rep(recruits[["Yr"]][1], nlines)
		}
		if (show_equilibrium) {
			points(x=xEqu, y=recruits[1, models], col=col, 
				pch=pch, cex=1.2, lwd=lwd)
		}
		if (show_uncertainty) {
			for (iline in 1:nlines) {
				imodel <- models[iline]
				if (uncertainty[imodel]) {
				  xvec <- recruits[["Yr"]]
				  if (nlines > 1) 
					xvec <- xvec + 0.4 * iline/nlines - 0.2
				  old_warn <- options()$warn
				  options(warn=-1)
				  arrows(x0=xvec[-c(1, 2)], y0=pmax(as.numeric(recruitsLower[-c(1, 
					2), imodel]), 0), x1=xvec[-c(1, 2)], y1=as.numeric(recruitsUpper[-c(1, 
					2), imodel]), length=0.01, angle=90, 
					code=3, col=col[imodel])
				  options(warn=old_warn)
				  if (show_equilibrium) {
					arrows(x0=xEqu[imodel], y0=pmax(as.numeric(recruitsLower[1, 
					  imodel]), 0), x1=xEqu[imodel], y1=as.numeric(recruitsUpper[1, 
					  imodel]), length=0.01, angle=90, code=3, 
					  col=col[imodel])
				  }
				}
			}
		}
		#abline(h=0, col="grey")
		if (legend) {
			#add_legend(legendlabels, legendloc=legendloc, legendorder=legendorder, legendncol=legendncol, col=col, pch=pch, lwd=lwd, lty=lty)
		}
		if (!add) {
			if (tickEndYr) {
				ticks <- graphics::axTicks(1)
				axis(1, at=c(ticks[ticks < max(endyrvec)], max(endyrvec)))
			}
			else {
				allyrs = recruits$Yr[1]:endyrvec[1]
				axis(1, at=intersect(allyrs, seq(1900,2050,5)), tcl=-0.2, labels=F)
				axis(1, at=intersect(allyrs, seq(1900,2050,10)), tcl=-0.4, labels=T)
			}
			if (!is.null(endyrvec) & max(endyrvec) > 1 + max(endyrs) & shadeForecast) {
				rect(xleft=max(endyrs) + 1, ybottom=par()$usr[3], xright=par()$usr[2], ytop=par()$usr[4], col=gray(0, alpha=0.1), border=NA)
			}
			yticks <- pretty(ylim,n=10)
			yticks = yticks[yticks<=ylim[2]]
			#axis(2, at=sort(unique(c(yticks,rev(rev((yticks+diff(yticks)[1]/2))[-1])))), tcl=-0.2, labels=F)
			axis(2, at=yticks, labels=format(yticks/yunits), las=1)
			box()
		}
#browser();return()
		return(ylim[2])
	}
	plotRecDevs <- function(show_uncertainty=TRUE, lang="e") {
		if (any(is.na(recdevs[["Yr"]]))) {
			warning("Recdevs associated with initial age structure may not be shown")
		}
		if (!any(uncertainty)) {
			show_uncertainty <- FALSE
		}
		if (is.null(xlim)) {
			xlim <- range(recdevs[["Yr"]], na.rm=TRUE)
			if (!is.null(endyrvec) & all(endyrvec < max(xlim))) {
				xlim[2] <- max(endyrvec)
			}
		}
		ylim <- ylimAdj * range(recdevs[recdevs[["Yr"]] >= xlim[1] & recdevs[["Yr"]] <= xlim[2], models], na.rm=TRUE)
		if (any(is.infinite(ylim))) {
			warning("Skipping recdev plots. Infinite ylim may indicate ", "all values are NA in summaryoutput[[\"recdevs\"]]")
			return(ylim[2])
		}
		xlim = xlim + c(-2,2)

		if (show_uncertainty) {
			if (all(is.na(recdevsLower[, models]))) {
				return(invisible(NA))
			}
			ylim <- ylimAdj * range(recdevsLower[recdevs[["Yr"]] >= 
				xlim[1] & recdevs[["Yr"]] <= xlim[2], models], 
				recdevsUpper[recdevs[["Yr"]] >= xlim[1] & recdevs[["Yr"]] <= 
				  xlim[2], models], na.rm=TRUE)
		}
		#ylim <- range(-ylim, ylim)
		ylim = extendrange(ylim)
		if (!add) {
			plot(0, xlim=xlim, ylim=ylim, axes=FALSE, type="n", xlab=linguaFranca(labels[1],lang), ylab=linguaFranca(labels[5],lang), xaxs=xaxs, yaxs=yaxs, las=1)
			allvals =  seq(-50,50,0.1)
			axis(2, at=allvals[allvals>=ylim[1] & allvals<=ylim[2]], tcl=-0.2, labels=F)
			axis(2, las=1, tcl=-0.4)
			abline(h=0, col="grey")
			rdper = list()
			rdper[["rd1"]] = grep("Early_InitAge", recdevs$Label)
			rdper[["rd2"]] = grep("Early_RecrDev", recdevs$Label)
			rdper[["rd3"]] = grep("Main_RecrDev", recdevs$Label)
			rdper[["rd4"]] = grep("Late_RecrDev", recdevs$Label)
			rdyrs = lapply(rdper,function(x){recdevs[,"Yr"][range(x)]})
			rdbrk = sapply(rdyrs, function(x){x[1]})
			abline(v=rdbrk[-1], lty=2)
			text(sapply(rdyrs,mean), rep(ylim[2]-0.02*diff(ylim),length(rdyrs)), switch(lang, 'e'=c("Init","Early","Main","Late"), 'f'=c("Initiale", "Tot", "Principale", "Tard")), font=2, adj=c(0.5,1), cex=1.5)
		}
		if (legend) {
			zmax = as.vector(sapply(recdevs[, models],function(x){ findPV(max(x,na.rm=T),x[!is.na(x)]) }))
			xmax = recdevs$Yr[zmax]
			ymax = diag(as.matrix(recdevs[zmax,models]))
			segments(xmax-2, ymax, xmax+2, ymax, col="grey10", lty=3, lwd=1)
			points(xmax, ymax, col=col, pch=pch, cex=1.2, lwd=2)
			addLegend(xleg, yleg, pch=pch, col=col, legend=endyrvec, cex=0.8, lwd=2, xjust=0, yjust=1, bty="n")
		}
		if (show_uncertainty) {
			for (iline in 1:nlines) {
				imodel <- models[iline]
				if (uncertainty[imodel]) {
					xvec <- recdevs[["Yr"]]
					if (nlines > 1) 
					xvec <- xvec + 0.4 * iline/nlines - 0.2
					arrows(x0=xvec, y0=as.numeric(recdevsLower[, imodel]), x1=xvec, y1=as.numeric(recdevsUpper[, imodel]), length=0.01, angle=90, code=3, col=col[iline])
				}
			}
		}
		for (iline in 1:nlines) {
			imodel <- models[iline]
			yvec <- recdevs[, imodel]
			xvec <- recdevs[["Yr"]]
			#points(xvec, yvec, pch=pch[iline], lwd=lwd[iline], col=col[iline])
			lines(xvec, yvec, lwd=lwd[iline], col=col[iline])
		}
		if (!add) {
			if (tickEndYr) {
				ticks <- graphics::axTicks(1)
				axis(1, at=c(ticks[ticks < max(endyrvec)], max(endyrvec)))
			}
			else {
				allyrs = recdevs$Yr[1]:endyrvec[1]
				axis(1, at=intersect(allyrs, seq(1900,2050,5)), tcl=-0.2, labels=F)
				axis(1, at=intersect(allyrs, seq(1900,2050,10)), tcl=-0.4, labels=T)
			}
			if (!is.null(endyrvec) & max(endyrvec) > 1 + max(endyrs) & shadeForecast) {
				rect(xleft=max(endyrs) + 1, ybottom=par()$usr[3],  xright=par()$usr[2], ytop=par()$usr[4], col=gray(0, alpha=0.1), border=NA)
			}
			box()
		}
		if (legend) {
			#add_legend(legendlabels, legendloc=legendloc, legendorder=legendorder, legendncol=legendncol, col=col, pch=pch, lwd=lwd, lty=lty)
		}
#browser();return()
		return(ylim[2])
	}
	plotPhase <- function(show_uncertainty=TRUE, lang="e") {
		if (!any(uncertainty)) {
			show_uncertainty <- FALSE
		}
		xlim <- range(0, ylimAdj * Bratio[, models], na.rm=TRUE)
		ylim <- range(0, ylimAdj * SPRratio[, models], na.rm=TRUE)
		if (!add) {
			plot(0, type="n", xlim=xlim, ylim=ylim, xlab=labels[3], 
				ylab=SPRratioLabel, xaxs=xaxs, yaxs=yaxs, 
				las=1)
		}
		goodyrs <- intersect(Bratio[["Yr"]], SPRratio[["Yr"]])
		lastyr <- max(goodyrs)
		for (iline in 1:nlines) {
			imodel <- models[iline]
			xvals <- Bratio[Bratio[["Yr"]] %in% goodyrs, imodel]
			yvals <- SPRratio[SPRratio[["Yr"]] %in% goodyrs, imodel]
			lines(xvals, yvals, col=col[iline], lty=lty[iline], lwd=lwd[iline], type="l")
			points(tail(xvals, 1), tail(yvals, 1), col=col[iline], 
				pch=pch[iline], lwd=lwd[iline])
		}
		abline(h=1, v=1, col="grey", lty=2)
		if (btarg > 0) 
			abline(v=btarg, col="red", lty=2)
		if (sprtarg > 0) 
			abline(h=sprtarg, col="red", lty=2)
		if (legend) {
			add_legend(legendlabels, legendloc=legendloc, legendorder=legendorder, 
				legendncol=legendncol, col=col, pch=pch, 
				lwd=lwd, lty=lty)
		}
		return(ylim[2])
	}
	plotIndices <- function(log=FALSE, iindex, lang="e") {
		indices2 <- NULL
		for (iline in 1:nlines) {
			imodel <- models[iline]
			subset2 <- indices[["imodel"]] == imodel & indices[["Yr"]] <= endyrvec[iline] & indices[["Fleet"]] == indexfleets[[imodel]][iindex]
			indices2 <- rbind(indices2, indices[subset2, ]) ## fleet-specific data
		}
		yr <- indices2[["Yr"]]
		obs <- indices2[["Obs"]]
		exp <- indices2[["Exp"]]
		imodel <- indices2[["imodel"]]
		Q <- indices2[["Calc_Q"]]
		if (log) {
			obs <- log(obs)
			exp <- log(exp)
			ylab <- labels[7]
		}
		else {
			ylab <- labels[6]
		}
		if (indexUncertainty) {
			if (indexPlotEach) {
				if (is.null(indexSEvec)) {
					indexSEvec <- indices2[["SE"]]
				}
				y <- obs
				if (log) {
					upper <- qnorm(0.975, mean=y, sd=indexSEvec)
					lower <- qnorm(0.025, mean=y, sd=indexSEvec)
				}
				else {
					upper <- qlnorm(0.975, meanlog=log(y), sdlog=indexSEvec)
					lower <- qlnorm(0.025, meanlog=log(y), sdlog=indexSEvec)
				}
			}
			else {
				subset <- indices2[["imodel"]] == models[1]
				if (is.null(indexSEvec)) {
					indexSEvec <- indices2[["SE"]][subset]
				}
				y <- obs[subset]  ## added subset (RH 220603)
				if (log) {
					upper <- qnorm(0.975, mean=y, sd=indexSEvec)
					lower <- qnorm(0.025, mean=y, sd=indexSEvec)
				}
				else {
					upper <- qlnorm(0.975, meanlog=log(y), sdlog=indexSEvec)
					lower <- qlnorm(0.025, meanlog=log(y), sdlog=indexSEvec)
				}
			}
		}
		else {
			upper <- NULL
			lower <- NULL
		}
		sub <- !is.na(indices2[["Like"]])
		ylim <- range(exp, obs[sub], lower[sub], upper[sub], na.rm=TRUE)
		if (!any(sub)) {
			ylim <- range(exp, obs, lower, upper, na.rm=TRUE)
		}
		if (!log) {
			ylim <- c(0, ylimAdj * ylim[2])
		}
		else {
			ylim <- ylim + c(-1, 1) * (ylimAdj - 1) * diff(ylim)
		}
		meanQ <- rep(NA, nlines)
		if (!add) {
			if (!is.null(endyrvec)) {
				xlim <- c(min(yr), max(yr))
			}
			else {
				xlim <- range(yr)
			}
			plot(0, type="n", xlim=xlim, yaxs=yaxs, ylim=ylim, xlab=linguaFranca("Year",lang), ylab="", axes=FALSE)
		}
		if (!log & yaxs != "i") {
			abline(h=0, col="grey")
		}
		if (indexPlotEach) {
			for (iline in (1:nlines)[!mcmcVec]) {
				adj <- 0.2 * iline/nlines - 0.1
				imodel <- models[iline]
				if (any(is.na(indices2[["like"]]))) {
					warning("NA's found in likelihood, may cause issues with index plots")
				}
				subset <- indices2[["imodel"]] == imodel & !is.na(indices2[["Like"]])
				if (indexUncertainty) {
					arrows(x0=yr[subset] + adj, y0=lower[subset], x1=yr[subset] + adj, y1=upper[subset], length=0.01, angle=90, code=3, col=adjustcolor(col, alpha.f=0.7)[iline])
				}
				points(yr[subset] + adj, obs[subset], pch=21, cex=1.5, col=1, bg=adjustcolor(col, alpha.f=0.7)[iline])
			}
		}
		else {
			imodel <- models[which(endyrvec == max(endyrvec))[1]]
			subset <- indices2[["imodel"]] == imodel & !is.na(indices2[["Like"]])
			if (indexUncertainty) {
				arrows(x0=yr[subset], y0=lower[subset], x1=yr[subset], y1=upper[subset], angle=90, code=3, length=0.025, col=1, lwd=2)
			}
			points(yr[subset], obs[subset], pch=15, col="red", cex=1)
		}
		Qtext <- rep("(Q =", nlines)
		for (iline in (1:nlines)[!mcmcVec]) {
			imodel <- models[iline]
			subset <- indices2[["imodel"]] == imodel
			meanQ[iline] <- mean(Q[subset])
			if (indexQlabel && any(Q[subset] != mean(Q[subset]))) {
				Qtext[iline] <- "(mean Q ="
			}
			x <- yr[subset]
			y <- exp[subset]
			lines(x, y, pch=pch[iline], lwd=lwd[iline], lty=lty[iline], col=col[iline], type=type)
		}
		legendlabels2 <- legendlabels
		if (indexQlabel) {
			legendlabels2 <- paste(legendlabels, Qtext, format(meanQ, digits=indexQdigits), ")")
		}
		if (legend) {
			#add_legend(legendlabels, legendloc=legendloc, legendorder=legendorder, legendncol=legendncol, col=col, pch=pch, lwd=lwd, lty=lty)
			xlast  = sapply(split(yr,indices2$imodel),function(x){rev(x)[1]})
			ylast  = sapply(split(exp,indices2$imodel),function(x){rev(x)[1]})
			points(xlast, ylast, col=col, pch=pch, cex=1.2, lwd=2)
			xpos = 0.95
			#if (iindex %in% c(1,3,4,8:50)) xpos = 0.2
			#if (iindex %in% c(2,7,79)) xpos = 0.3
			addLegend(xleg, yleg, pch=pch, col=col, legend=linguaFranca(endyrvec,lang), cex=0.8, lwd=2, xjust=0, yjust=1, bty="n", title=linguaFranca(gsub("_", " ", indices2$Fleet_name[1]), lang) )
		}
		if (!add) {
			#xticks <- pretty(xlim)
			#axis(1, at=xticks, labels=format(xticks))
			allyrs = min(indices2$Yr):max(indices$Yr)
			axis(1, at=intersect(allyrs, seq(1900,2050,1)), tcl=-0.2, labels=F)
			axis(1, at=intersect(allyrs, seq(1900,2050,5)), tcl=-0.4, labels=T)
			if (tickEndYr) {
				axis(1, at=max(endyrvec))
			}
			yticks <- pretty(ylim,n=10)
			ytkbig = yticks[seq(1,length(yticks),2)]
			if (min(ytkbig[-1])>1e3){
				ytklab = ytkbig/1000; ylab = "Index (1000s)"
			} else {
				ytklab = ytkbig ; yalb ="Index"
			}
			axis(2, at=yticks, tcl=-0.2, labels=F)
			axis(2, at=ytkbig, tcl=-0.4, labels=ytklab, las=1)
			mtext(linguaFranca(ylab,lang), side=2, line=2, cex=1.5)
			box()
		}
#browser();return()
		return(ylim[2])
	}  ## end plotIndices

	plotDensities <- function(parname, xlab, denslwd, limit0=TRUE, cumulative=FALSE, lang="e") {
		if (any(!mcmcVec)) {
			vals <- rbind(pars[pars[["Label"]] == parname, names(pars) != "recdev"], quants[quants[["Label"]] == parname, ])
			if (nrow(vals) != 1) {
				warn <- paste("problem getting values for parameter:", parname, "")
				if (nrow(vals) == 0) {
					warn <- paste(warn, "no Labels match in either parameters or derived quantities")
				}
				if (nrow(vals) > 0) {
				  warn <- paste(warn, "Too many matching Labels:", 
					pars[["Label"]][pars[["Label"]] == parname], 
					quants[["Label"]][quants[["Label"]] == parname])
				}
				warning(warn)
				return(NULL)
			}
			valSDs <- rbind(parsSD[pars[["Label"]] == parname, ], quantsSD[quants[["Label"]] == parname, ])
		}
		xmax <- xmin <- ymax <- NULL
		mcmcDens <- vector(mode="list", length=nlines)
		good <- rep(TRUE, nlines)
		for (iline in 1:nlines) {
			imodel <- models[iline]
			if (mcmcVec[iline]) {
				mcmcColumn <- grep(parname, colnames(mcmc[[imodel]]), fixed=TRUE)
				if (length(mcmcColumn) == 0) {
					message("No columns selected from MCMC for '", parname, "' in model ", imodel)
					good[iline] <- FALSE
				}
				if (length(mcmcColumn) > 1) {
					message("Too many columns selected from MCMC for model ", imodel, ":")
					print(names(mcmc[[imodel]])[mcmcColumn])
					warning("Please specify a unique label in the mcmc dataframe", "or specify mcmcVec=FALSE for model ", imodel, " (or mcmcVec=FALSE applying to all models).")
					good[iline] <- FALSE
				}
				if (good[iline]) {
					mcmcVals <- mcmc[[imodel]][, mcmcColumn]
					xmin <- min(xmin, quantile(mcmcVals, 0.005, na.rm=TRUE))
					if (limit0) {
						xmin <- max(0, xmin)
					}
					if (fix0 & !grepl("R0", parname)) {
						xmin <- 0
					}
					xmax <- max(xmax, quantile(mcmcVals, 0.995, na.rm=TRUE))
					z <- density(mcmcVals, cut=0, adjust=densityadjust)
					z[["x"]] <- z[["x"]][c(1, 1:length(z[["x"]]), length(z[["x"]]))]
					z[["y"]] <- c(0, z[["y"]], 0)
					ymax <- max(ymax, max(z[["y"]]))
					mcmcDens[[iline]] <- z
				}
			}
			else {
				parval <- vals[1, imodel]
				parSD <- valSDs[1, imodel]
				if (!is.numeric(parval)) 
					parval <- -1
				if (!is.na(parSD) && parSD > 0) {
					xmin <- min(xmin, qnorm(0.005, parval, parSD))
					if (limit0) 
						xmin <- max(0, xmin)
					if (fix0 & !grepl("R0", parname)) 
						xmin <- 0
					xmax <- max(xmax, qnorm(0.995, parval, parSD))
					x <- seq(xmin, xmax, length=500)
					mle <- dnorm(x, parval, parSD)
					mlescale <- 1/(sum(mle) * mean(diff(x)))
					mle <- mle * mlescale
					ymax <- max(ymax, max(mle))
				}
				else {
					xmin <- min(xmin, parval)
					xmax <- max(xmax, parval)
				}
			}
		}
		if (grepl("Bratio", parname)) {
			xmin <- 0
		}
		if (limit0) {
			xmin <- max(0, xmin)
		}
		if (fix0 & !grepl("R0", parname)) {
			xmin <- 0
		}
		xlim <- c(xmin, xmin + (xmax - xmin) * densityscalex)
		x <- seq(xmin, xmax, length=500)
		xunits <- 1
		if (rescale & xmax > 1000 & xmax < 3e+06) {
			xunits <- 1000
			xlab2 <- "'1000 t"
		}
		if (rescale & xmax > 3e+06) {
			xunits <- 1e+06
			xlab2 <- "million t"
		}
		if (is.null(ymax)) {
			message("  skipping plot of ", parname, " because it seems to not be estimated in any model")
		}
		else {
			par(par)
			if (!add) {
				if (cumulative) {
					plot(0, type="n", xlim=xlim, axes=FALSE, xaxs="i", yaxs=yaxs, ylim=c(0, 1), xlab=xlab, ylab="")
				}
				else {
					plot(0, type="n", xlim=xlim, axes=FALSE, xaxs="i", yaxs=yaxs, ylim=c(0, 1.1 * ymax * densityscaley), xlab=xlab, ylab="")
				}
			}
			if (grepl("Bratio", parname)) {
				if (btarg > 0) {
					abline(v=btarg, col="red", lty=2)
					text(btarg + 0.03, par()$usr[4], labels[10], adj=1.05, srt=90)
				}
				if (minbthresh > 0) {
					abline(v=minbthresh, col="red", lty=2)
					text(minbthresh + 0.03, par()$usr[4], labels[11], adj=1.05, srt=90)
				}
			}
			symbolsQuants <- c(0.025, 0.125, 0.25, 0.5, 0.75, 0.875, 0.975)
			for (iline in (1:nlines)[good]) {
				imodel <- models[iline]
				if (mcmcVec[iline]) {
					mcmcColumn <- grep(parname, colnames(mcmc[[imodel]]), fixed=TRUE)
					mcmcVals <- mcmc[[imodel]][, mcmcColumn]
					x2 <- quantile(mcmcVals, symbolsQuants, na.rm=TRUE)
					x <- mcmcDens[[iline]][["x"]]
					if (!cumulative) {
						y <- mcmcDens[[iline]][["y"]]
						yscale <- 1/(sum(y) * mean(diff(x)))
						y <- y * yscale
					}
					else {
						y <- cumsum(mcmcDens[[iline]][["y"]])/sum(mcmcDens[[iline]][["y"]])
					}
					y2 <- NULL
					for (ii in x2) {
						y2 <- c(y2, min(y[abs(x - ii) == min(abs(x - ii))]))
					}
					if (!cumulative) {
						polygon(c(x[1], x, rev(x)[1]), c(0, y, 0), col=shadecol[iline], border=NA)
					}
					else {
						polygon(c(x[1], x, rev(x)[c(1, 1)]), c(0, y, 1, 0), col=shadecol[iline], border=NA)
					}
					lines(x, y, col=col[iline], lwd=2)
					if (!cumulative) {
						if (densitysymbols) {
							points(x2, y2, col=col[iline], pch=pch[iline])
						}
						lines(rep(x2[median(1:length(x2))], 2), c(0, y2[median(1:length(x2))]), col=col[iline])
					}
					else {
						if (densitysymbols) {
							points(x2, symbolsQuants, col=col[iline], pch=pch[iline])
						}
						lines(rep(median(mcmcVals), 2), c(0, 0.5), col=col[iline])
					}
				}
				else {
					parval <- vals[1, imodel]
					parSD <- valSDs[1, imodel]
					if (!is.na(parSD) && parSD > 0) {
						xmin <- min(xmin, qnorm(0.005, parval, parSD))
						if (limit0) {
							xmin <- max(0, xmin)
						}
						if (fix0 & !grepl("R0", parname)) {
							xmin <- 0
						}
						x <- seq(xmin, max(xmax, xlim), length=500)
						x2 <- qnorm(symbolsQuants, parval, parSD)
						if (cumulative) {
							y <- mle <- pnorm(x, parval, parSD)
							y2 <- mle2 <- pnorm(x2, parval, parSD)
						}
						else {
							mle <- dnorm(x, parval, parSD)
							mle2 <- dnorm(x2, parval, parSD)
							mlescale <- 1/(sum(mle) * mean(diff(x)))
							y <- mle <- mle * mlescale
							y2 <- mle2 <- mle2 * mlescale
						}
						polygon(c(x[1], x, rev(x)[1]), c(0, mle, 0), col=shadecol[iline], border=NA)
						lines(x, mle, col=col[iline], lwd=2)
						if (!cumulative) {
							if (densitysymbols) {
								points(x2, mle2, col=col[iline], pch=pch[iline], cex=1.2)
							}
							lines(rep(parval, 2), c(0, dnorm(parval, parval, parSD) * mlescale), col=col[iline], lwd=ifelse(iline==1,3,2), lty=ifelse(iline==1,1,3) )
#if(iline==4) {browser();return()}
						}
						else {  ## cumulative
							if (densitysymbols) {
								points(x2, symbolsQuants, col=col[iline], pch=pch[iline])
							}
							lines(rep(parval, 2), c(0, 0.5), col=col[iline], lwd=denslwd)
						}
					}
					else {
						abline(v=parval, col=col[iline], lwd=denslwd)
					}
				}
				if (densitytails & densitymiddle) {
					warning("You are shading both tails and central 95% of density plots", "which is illogical")
				}
				doShade <- FALSE
				if (mcmcVec[iline]) {
					doShade <- TRUE
				}
				else {
					if (!is.na(parSD) && parSD > 0) {
						doShade <- TRUE
					}
				}
				if (densitytails & doShade) {
					x.lower <- x[x <= x2[1]]
					y.lower <- y[x <= x2[1]]
					x.upper <- x[x >= rev(x2)[1]]
					y.upper <- y[x >= rev(x2)[1]]
					polygon(c(x.lower[1], x.lower, rev(x.lower)[1]), c(0, y.lower, 0), col=shadecol[iline], border=NA)
					polygon(c(x.upper[1], x.upper, rev(x.upper)[1]), c(0, y.upper, 0), col=shadecol[iline], border=NA)
				}
				if (densitymiddle & doShade) {
					x.middle <- x[x >= x2[1] & x <= rev(x2)[1]]
					y.middle <- y[x >= x2[1] & x <= rev(x2)[1]]
					polygon(c(x.middle[1], x.middle, rev(x.middle)[1]), c(0, y.middle, 0), col=shadecol[iline], border=NA)
				}
			}
			if (!add) {
				abline(h=0, col="grey")
				xticks <- pretty(xlim)
				axis(1, at=xticks, labels=format(xticks/xunits))
				theLine <- par()$mgp[1]
				if (cumulative) {
					axis(2, at=symbolsQuants, labels=format(symbolsQuants), cex.axis=0.9)
					mtext(side=2, line=theLine, text="Cumulative Probability", col=par()$col.lab, cex=par()$cex.lab)
				}
				else {
					mtext(side=2, line=theLine, text=labels[9], col=par()$col.lab, cex=par()$cex.lab)
				}
				box()
			}
			if (xunits != 1) {
				message("x-axis for ", parname, " in density plot has been divided by ", xunits, " (so may be in units of ", xlab2, ")")
			}
			if (legend) {
				addLegend(xleg, yleg, pch=pch, col=col, legend=legendlabels, cex=0.8, lwd=2, xjust=0, yjust=1, bty="n")
				#add_legend(legendlabels, legendloc=ifelse(cumulative, "topleft", legendloc), legendorder=legendorder, legendncol=legendncol, col=col, pch=pch,  lwd=lwd, lty=lty)
			}
		}
		return(NA)
	} ## end plotDensities

	uncertaintyplots <- intersect(c(2, 4, 6, 8, 10, 12), subplots)
	if (!any(uncertainty) & length(uncertaintyplots) > 0) {
		message("skipping plots with uncertainty:", paste(uncertaintyplots, collapse=","))
	}
	par(par)
	if (1 %in% subplots) {
		if (verbose) {
			message("subplot 1: spawning biomass")
		}
		if (plot) {
			for (l in lang) {
				changeLangOpts(L=lang)
				ymax_vec[1] <- plotSpawnBio(show_uncertainty=FALSE, lang=l)
			}; eop()
		}
		if (print) {
			for (l in lang) {
				save_png_comparisons("compare1_spawnbio.png", lang=l)
				ymax_vec[1] <- plotSpawnBio(show_uncertainty=FALSE, lang=l)
				dev.off()
			}; eop()
		}
	}
	if (2 %in% subplots) {
		if (any(uncertainty)) {
			if (verbose) {
				message("subplot 2: spawning biomass with uncertainty intervals")
			}
			if (plot) {
				for (l in lang) {
					changeLangOpts(L=lang)
					ymax_vec[2] <- plotSpawnBio(show_uncertainty=TRUE, lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang) {
					save_png_comparisons("compare2_spawnbio_uncertainty.png", lang=l)
					ymax_vec[2] <- plotSpawnBio(show_uncertainty=TRUE, lang=l)
					dev.off()
				}; eop()
			}
		}
	}

	if (3 %in% subplots) {
		if (verbose) {
			message("subplot 3: biomass ratio (hopefully equal to fraction of unfished)")
		}
		if (plot) {
			for (l in lang) {
				changeLangOpts(L=lang)
				ymax_vec[3] <- plotBratio(show_uncertainty=FALSE, lang=l)
			}; eop()
		}
		if (print) {
			for (l in lang) {
				save_png_comparisons("compare3_Bratio.png", lang=l)
				ymax_vec[3] <- plotBratio(show_uncertainty=FALSE, lang=l)
				dev.off()
			}; eop()
		}
	}
	if (4 %in% subplots) {
		if (any(uncertainty)) {
			if (verbose) {
				message("subplot 4: biomass ratio with uncertainty")
			}
			if (plot) {
				for (l in lang){
					changeLangOpts(L=lang)
					ymax_vec[4] <- plotBratio(show_uncertainty=TRUE, lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang){
					save_png_comparisons("compare4_Bratio_uncertainty.png", lang=l)
					ymax_vec[4] <- plotBratio(show_uncertainty=TRUE, lang=l)
					dev.off()
				}; eop()
			}
		}
	}
	if (5 %in% subplots) {
		if (verbose) {
			message("subplot 5: SPR ratio")
		}
		if (plot) {
			for (l in lang) {
				changeLangOpts(L=lang)
				ymax_vec[5] <- plotSPRratio(show_uncertainty=FALSE, lang=l)
			}; eop()
		}
		if (print) {
			for (l in lang) {
				save_png_comparisons("compare5_SPRratio.png", lang=l)
				ymax_vec[5] <- plotSPRratio(show_uncertainty=FALSE, lang=l)
				dev.off()
			}; eop()
		}
	}
	if (6 %in% subplots) {
		if (any(uncertainty)) {
			if (verbose) {
				message("subplot 6: SPR ratio with uncertainty")
			}
			if (plot) {
				for (l in lang) {
					changeLangOpts(L=lang)
					ymax_vec[6] <- plotSPRratio(show_uncertainty=TRUE, lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang) {
					save_png_comparisons("compare6_SPRratio_uncertainty.png", lang=l)
					ymax_vec[6] <- plotSPRratio(show_uncertainty=TRUE, lang=l)
					dev.off()
				}; eop()
			}
		}
	}
	if (7 %in% subplots) {
		if (verbose) {
			message("subplot 7: F value")
		}
		if (plot) {
			for (l in lang) {
				changeLangOpts(L=lang)
				ymax_vec[7] <- plotF(show_uncertainty=FALSE, lang=l)
			}; eop()
		}
		if (print) {
			for (l in lang) {
				save_png_comparisons("compare7_Fvalue.png", lang=l)
				ymax_vec[7] <- plotF(show_uncertainty=FALSE, lang=l)
				dev.off()
			}; eop()
		}
	}
	if (8 %in% subplots) {
		if (any(uncertainty)) {
			if (verbose) {
				message("subplot 8: F value with uncertainty")
			}
			if (plot) {
				for (l in lang) {
					changeLangOpts(L=lang)
					ymax_vec[8] <- plotF(show_uncertainty=TRUE, lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang) {
					save_png_comparisons("compare8_Fvalue_uncertainty.png", lang=l)
					ymax_vec[8] <- plotF(show_uncertainty=TRUE, lang=l)
					dev.off()
				}; eop()
			}
		}
	}
	if (9 %in% subplots) {
		if (verbose) {
			message("subplot 9: recruits")
		}
		if (plot) {
			for (l in lang) {
				ymax_vec[9] <- plotRecruits(show_uncertainty=FALSE, lang=l)
			}; eop()
		}
		if (print) {
			for (l in lang) {
				save_png_comparisons("compare9_recruits.png", lang=l)
				ymax_vec[9] <- plotRecruits(show_uncertainty=FALSE, lang=l)
				dev.off()
			}; eop()
		}
	}
	if (10 %in% subplots) {
		if (any(uncertainty)) {
			if (verbose) {
				message("subplot 10: recruits with uncertainty")
			}
			if (plot) {
				for (l in lang) {
					changeLangOpts(L=lang)
					ymax_vec[10] <- plotRecruits(lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang) {
					save_png_comparisons("compare10_recruits_uncertainty.png", lang=l)
					ymax_vec[10] <- plotRecruits(lang=l)
					dev.off()
				}; eop()
			}
		}
	}
	if (11 %in% subplots) {
		if (verbose) 
			message("subplot 11: recruit devs")
		if (is.null(recdevs)) {
			message("No recdevs present in the model summary, skipping plot.")
		}
		else {
			if (plot) {
				for (l in lang) {
					ymax_vec[11] <- plotRecDevs(show_uncertainty=FALSE, lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang) {
					save_png_comparisons("compare11_recdevs.png", lang=l)
					ymax_vec[11] <- plotRecDevs(show_uncertainty=FALSE, lang=l)
					dev.off()
				}; eop()
			}
		}
	}
	if (12 %in% subplots) {
		if (any(uncertainty)) {
			if (verbose) {
				message("subplot 12: recruit devs with uncertainty")
			}
			if (plot) {
				for (l in lang) {
					changeLangOpts(L=lang)
					ymax_vec[12] <- plotRecDevs(lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang) {
					save_png_comparisons("compare12_recdevs_uncertainty.png", lang=l)
					ymax_vec[12] <- plotRecDevs(lang=l)
					dev.off()
				}; eop()
			}
		}
	}
	if (13 %in% subplots & !is.null(indices) && nrow(indices) > 0) {
		if (verbose) {
			message("subplot 13: index fits")
		}
		for (iindex in 1:length(indexfleets[[1]])) {
			#if (iindex!=1) next
			if (plot) {
				for (l in lang) {
					changeLangOpts(L=lang)
					ymax_vec[13] <- plotIndices(log=FALSE, iindex=iindex, lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang) {
					save_png_comparisons(paste0("compare13_indices", index_plot_suffix[iindex], ".png"), lang=l)
					ymax_vec[13] <- plotIndices(log=FALSE, iindex=iindex, lang=l)
					dev.off()
				}; eop()
			}
		}
	}
	if (14 %in% subplots & !is.null(indices) && nrow(indices) > 
		0) {
		if (verbose) {
			message("subplot 14: index fits on a log scale")
		}
		for (iindex in 1:length(indexfleets[[1]])) {
			if (plot) {
				for (l in lang) {
					changeLangOpts(L=lang)
					ymax_vec[14] <- plotIndices(log=TRUE, iindex=iindex, lang=l)
				}; eop()
			}
			if (print) {
				for (l in lang) {
					save_png_comparisons(paste0("compare14_indices_log", index_plot_suffix[iindex], ".png"), lang=l)
					ymax_vec[14] <- plotIndices(log=TRUE, iindex=iindex, lang=l)
					dev.off()
				}; eop()
			}
		}
	}
	if (15 %in% subplots) {
		if (verbose) {
			message("subplot 15: phase plot")
		}
		if (plot) {
			for (l in lang) {
				changeLangOpts(L=lang)
				ymax_vec[15] <- plotPhase(lang=l)
			}; eop()
		}
		if (print) {
			for (l in lang) {
				save_png_comparisons("compare15_phase_plot.png", lang=l)
				ymax_vec[15] <- plotPhase(lang=l)
				dev.off()
			}; eop()
		}
	}
	if (16 %in% subplots | 17 %in% subplots) {
		if (any(uncertainty)) {
			if (verbose) {
				message("subplots 16 and 17: densities")
			}
			expandednames <- NULL
			for (i in 1:length(densitynames)) {
				matchingnames <- c(pars[["Label"]], quants[["Label"]])[grep(densitynames[i], c(pars[["Label"]], quants[["Label"]]), fixed=TRUE)]
				expandednames <- c(expandednames, matchingnames)
			}
			if (length(expandednames) == 0) {
				warning("  No parameter/quantity names matching 'densitynames' input.")
			}
			else {
				message("  parameter/quantity names matching 'densitynames' input:")
				print(expandednames)
				ndensities <- length(expandednames)
				densitytable <- data.frame(name=expandednames, 
				  label=expandednames, stringsAsFactors=FALSE)
				if (!is.null(densityxlabs) && length(densityxlabs) == 
				  ndensities) {
				  densitytable[["label"]] <- densityxlabs
				  message("  table of parameter/quantity labels with associated", 
					" x-axis label:")
				  print(densitytable)
				}
				else {
				  if (!is.null(densityxlabs)) {
					warning("length of 'densityxlabs' doesn't match the number of values ", 
					  "matching 'densitynames' so parameter labels will be used instead")
				  }
				}
				if (16 %in% subplots) {
					for (iplot in 1:ndensities) {
						name <- densitytable[iplot, 1]
						xlab <- densitytable[iplot, 2]
						if (plot) {
							for (l in lang) {
								changeLangOpts(L=lang)
								ymax_vec[16] <- plotDensities(parname=name, xlab=xlab, denslwd=densitylwd, lang=l)
							}; eop()
						}
						if (print) {
							for (l in lang) {
								save_png_comparisons(paste("compare16_densities_", name, ".png", sep=""), lang=l)
								ymax_vec[16] <- plotDensities(parname=name, xlab=xlab, denslwd=densitylwd, lang=l)
								dev.off()
							}; eop()
						}
					}
				}
				if (17 %in% subplots) {
					for (iplot in 1:ndensities) {
						name <- densitytable[iplot, 1]
						xlab <- densitytable[iplot, 2]
						if (plot) {
							for (l in lang) {
								changeLangOpts(L=lang)
								ymax_vec[17] <- plotDensities(parname=name, xlab=xlab, denslwd=densitylwd, cumulative=TRUE, lang=l)
							}; eop()
						}
						if (print) {
							for (l in lang) {
								save_png_comparisons(paste("compare17_densities_", name, ".png", sep=""), lang=l)
								ymax_vec[17] <- plotDensities(parname=name, xlab=xlab, denslwd=densitylwd, cumulative=TRUE, lang=l)
								dev.off()
							}; eop()
						}
					}
				}
			}
		}
	}
	if (pdf) 
		dev.off()
	return(invisible(ymax_vec))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.comparisons


## plotSS.comps-------------------------2023-09-13
##  Plot age proportions and model fit
##  Source: R package 'r4ss' v.1.39.1
##  Based on 'r4ss::SSplotComps'
## ----------------------------------------r4ss|RH
plotSS.comps <- function (replist, subplots=c(1:21, 24), kind="LEN", sizemethod=1, 
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
	fart <- function(opar) { if (any("windows"%in%names(dev.list()))) par(opar); eop() }
	on.exit(fart(oldpar))
	changeLangOpts(L=lang)
	if (missing(outnam))
		outnam = NULL
	ttput(outnam)

	if (!is.element(subplots, c(1:21, 24))) {
		mess = c(" \n","Choose another subplot from:",
			"  1  = multipanel (by year) fits to age proportions for multiple fleets",
			"  2  = single panel (by fleet) bubble plots of age proportions by year",
			"  3  = multi-panel bubble plots for conditional age-at-length",
			"  4  = multi-panel plot of fit to conditional age-at-length for specific years",
			"  5  = Pearson residuals for A-L key",
			"  6  = multi-panel plot of point and line fit to conditional age-at-length for specific length bins",
			"  7  = sample size plot",
			"  8  = TA1.8 Francis weighting plot",
			"  9  = TA1.8 Francis weighting plot for conditional data",
			"  10 = Andre's mean age and std. dev. in conditional AAL",
			"  11-20 = Qu'est-ce qu'il mange en hiver?",
			"  21 = composition by fleet aggregating across years",
			"  22 = composition by fleet aggregating across years within each season",
			"  23 = composition by fleet aggregating across seasons within a year",
			"  24 = bubble plot comparison of length or age residuals"
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
									sslab = "N obs = "
									obsN  = 0
								} else {
									sslab = "N wtd = "
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
				#if (7 %in% subplots & samplesizeplots & !datonly & !("DM_effN" %in% names(dbase) && any(!is.na(dbase$DM_effN))) & !(kind %in% c("GSTAGE", "GSTLEN", "L@A", "W@A"))) {
				if (7 %in% subplots & samplesizeplots & !datonly & !(kind %in% c("GSTAGE", "GSTLEN", "L@A", "W@A"))) {
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
						addLabel(0.98,0.95,txt=fleets.lab[unique(dbasegood$Fleet)], adj=c(1,0), cex=1, col="red")
#browser();return()

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


## plotSS.francis-----------------------2025-09-11
##  Plot mean age fits using Francis (2011) methodology
##  Modified function 'r4ss::SSMethod.TA1.8' :
##  Apply Francis composition weighting method TA1.8
## ----------------------------------------r4ss|RH
plotSS.francis <- function(
   fit, type, fleet, part=0:2, sexes=0:3, seas=NULL, method=NULL, 
   plotit=TRUE, printit=TRUE, datonly=FALSE, plotadj=!datonly,
   maxpanel=1000, fleetnames=NULL, label.part=TRUE, label.sex=TRUE, 
   set.pars=TRUE, col.obs=c("green3","green"), col.fit=lucent(c("blue2","cyan"),0.5),
   png=FALSE, pngres=400, PIN=c(8,9), outnam, lang=c("e","f") )
{
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

	oldpar = par(no.readonly=TRUE)
	fart <- function(opar) { if (any("windows"%in%names(dev.list()))) par(opar) }
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
	sex.flag <- type!='con' & max(tapply(dbase$'Sexes', dbase$Fleet, function(x)length(unique(x)))) > 1
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
			#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )

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
			if (Npanel < 6) {
				Nr = Npanel; Nc = 1
			} else {
				Nc = 2; Nr = ceiling(Npanel/2)
			}
			if(set.pars){
				# save current graphical parameters
				par_current <- par()
				# set new parameters
				if (sex.flag) {
					par(mfrow=c(Nr,Nc), mar=c(1.5,2,1,1)+0.1, mgp=c(0,0.5,0), oma=c(1.2,1.2,0,0), las=1)
				} else {
					par(mfcol=c(Nr,Nc), mar=c(1.5,2,1,1)+0.1, mgp=c(0,0.5,0), oma=c(1.2,1.2,0,0), las=1)
				}
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
				## Lines
				fl <- fleetnames[subpldat[1,'Fleet']]
				yr <- paste(subpldat[1,'Yr'])
				lab <- if(type=='con') ifelse(Nfleet>1,paste(yr,fl),yr) else fl
#browser();return()
				if(sex.flag & label.sex){
					lab <- paste(lab,ifelse(subpldat[1,'sexes']==0,'comb','sex'))
				}
				if(method.flag){
					lab <- paste(lab,'meth',subpldat[1,'method'])
				}
				if(label.part){
					lab <- paste(lab,partition.labels)
				}
				if(sex.flag & label.sex){
					lab = sub("sex \\(whole", paste0(switch(subpldat[1,'sexes'], "(female","(male","(	other"), " whole"), lab)
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
				par(mfrow=par_current$mfrow, mfcol=par_current$mfcol, mar=par_current$mar, mgp=par_current$mgp, oma=par_current$oma, las=par_current$las)
			}
		}; eop()
	} ## end if plotit
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
		Output <- list(agedat=pldat, w=Nmult, lo=confint[1,], hi=confint[2,], w.francis=wj)
		Outs <- paste0("Francis Weights - ", type, ": ", fleetnames[fleet], ": ", round(Output[["w"]],4), " (", round(Output[["lo"]],4), "-", round(Output[["hi"]],4), ")")
#browser();return()
		
		
		#tmp <- matrix(sample(pldat[,'Std.res'],1000*nrow(pldat),replace=TRUE),nrow(pldat))
		#confint <- as.vector(quantile(apply(tmp,2,function(x)1/var(x,na.rm=TRUE)), c(0.025,0.975), na.rm=TRUE))
		#Output <- c(w=Nmult,lo=confint[1],hi=confint[2])
		#Outs <- paste("Francis Weights - ", type, ": ", fleetnames[fleet],": ", round(Nmult,4), " (",round(confint[1],4),"-",round(confint[2],4),")", sep="")
		if(printit){
			print(Outs)
		}
		w.francis = Output; ttput(w.francis)
		return(Output)
	} ## end if not data only
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.francis


## plotSS.index-------------------------2025-07-24
##  Plot SS model fit to abundance index series
##  Code based on r4ss function 'SSplotIndices'.
## ----------------------------------------r4ss|RH
plotSS.index <- function (replist, subplots=c(1:10, 12), plot=TRUE, print=FALSE, 
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
		fout = switch(lang, 'e' = file.path(plotdir,"english", file), 'f' = file.path(plotdir,"french", file) )
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
	names(labels.default)=1:length(labels.default)
	if (is.null(names(labels)))
		names(labels)=1:length(labels)
	labels.new=labels.default
	labels.new[names(labels)]=labels
	labels=labels.new
	nlabs =sapply(labels,countVec) ## how many labels are available
#browser();return()
	get.lab <- function(i, j, k=nlabs)
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
				lower_total <- qlnorm(0.025, meanlog=log(y[include]), sdlog=cpueuse$SE[include])
				upper_total <- qlnorm(0.975, meanlog=log(y[include]), sdlog=cpueuse$SE[include])
			}
			else {
				lower_total <- qnorm(0.025, mean=log(y[include]), sd=cpueuse$SE[include])
				upper_total <- qnorm(0.975, mean=log(y[include]), sd=cpueuse$SE[include])
			}
		}
		if (error == -1) {
			lower_total <- qnorm(0.025, mean=y[include], sd=cpueuse$SE[include])
			upper_total <- qnorm(0.975, mean=y[include], sd=cpueuse$SE[include])
		}
		if (error > 0) {
			lower_total <- log(y[include]) + qt(0.025, df=error) * cpueuse$SE[include]
			upper_total <- log(y[include]) + qt(0.975, df=error) * cpueuse$SE[include]
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
				#ylim <- c(0, 1.05 * min(max(upper_total, zmax, na.rm=TRUE), max(maximum_ymax_ratio * y)))
				ylim <- c(0, 1.05 * min(max(upper_total, max(z), na.rm=TRUE), max(maximum_ymax_ratio * y)))
			}
			if (log) {
				ylim <- range(c(lower_total, upper_total), na.rm=TRUE)
			}
			plot(x=x[include], y=y[include], type="n", xlim=xlim, ylim=ylim, yaxs=ifelse(log, "r", "i"), 
				xlab=linguaFranca(get.lab(1,ifleet),l), ylab=linguaFranca("Abundance Index",l), main=linguaFranca(main,l), cex.main=cex.main, cex.axis=1.2, cex.lab=1.5, ...)
			panlab = ifelse(!log, get.lab(2,ifleet), get.lab(5,ifleet))
			addLabel(0.925, 0.975, linguaFranca(panlab,l), adj=c(1,1), cex=1.5, col="blue")
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
			segments(x[include], lower_input, x[include], upper_input, col=colvec1[s], lwd=3, lend=1)
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
#browser();return()
		}
		else {
			points(x=x[include], y=log(y[include]), pch=pch1, cex=cex, bg=bg, col=colvec1[s])
			if (addexpected) {
				lines(x, log(z), lwd=2, col=col3)
			}
		}
		if (legend & length(colvec1) > 1) {
			legend(x=legendloc, legend=seasnames, pch=pch1, col=colvec1, cex=cex)
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
		#ylim <- c(-1.05, 1.05) * max(abs(y[include]))
		ylim <- c(1.1,1.4) * c(pmin(-1, min(y[include])), pmax(1, max(y[include])))
		if (!add) {
			plot(x=x[include], y=y[include], type="n", xlim=xlim, ylim=ylim, yaxs="i", xlab=linguaFranca(get.lab(1,ifleet),l), ylab=linguaFranca(ylab,l), main=linguaFranca(main,l), cex.main=cex.main, cex.axis=1.25, cex.lab=1.5, ...)
		}
		ii = match(ifleet,fleetvec)  ##grep(ifleet,fleetvec) ## (RH 250912) grep grabs 1, 10, 11, etc.
		#points(x=x[include], y=y[include], pch=pch1, cex=cex, bg=adjustcolor(colvec1[s], alpha.f=0.7), col=adjustcolor(colvec1[ii], alpha.f=0.7))
		nx = length(x[include])
		abline(h=-10:10, lty=3, col="slategray")  ## (RH 250724)
		segments (x0=x[include], y0=rep(0,nx), x1=x[include], y1=y[include], lty=1, col=colvec1[ii], lwd=1)  ## (RH 250724)
		abline(h=0, lty=1, col="grey60")
		points(x=x[include], y=y[include], pch=21, cex=1.5, col=colvec1[ii], bg=colvec2[ii])
		addLabel(0.975, 0.95, linguaFranca(gsub("_+","  ",fleetnames[ifleet]),l), col="grey30", cex=1.25, adj=c(1,0.25))
#browser();return()
#if(ifleet==4) {browser();return()}
		if (legend & length(colvec1) > 1) {
			legend(x=legendloc, legend=linguaFranca(seasnames,l), pch=pch1, pt.bg=colvec1, col=colvec1, cex=cex)
		}
	} ## end index_resids.fn

	obs_vs_exp.fn <- function(log=FALSE, l="e", ...) {
		main <- paste(get.lab(2,ifleet), Fleet, sep=" ")
		if (!mainTitle) {
			main <- ""
		}
		if (!add) {
			if (!log) {
				plot(y[include], z[include], type="n", xlab=linguaFranca(get.lab(3,ifleet),l), 
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
				fout = switch(lang, 'e' = file.path(plotdir, "english", onefile), 'f' = file.path(plotdir,"french", onefile) )
				clearFiles(fout)
				png(filename=fout, units="in", res=res, width=PIN[1], height=PIN[2])
			}
			rc=.findSquare(length(fleetvec))
			expandGraph(mfcol=rc, mar=c(3,3.5,0.5,1), oma=c(1,0,0,0), mgp = c(2,0.5,0))
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
				colvec1 <- rep(col1,length(fleetvec))[1:length(fleetvec)]
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
				colvec2 <- rep(col2,length(fleetvec))[1:length(fleetvec)]
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
#browser();return()
			if (exists("Q_extraSD_info") && ifleet %in% Q_extraSD_info$Fleet) {
				cpueuse$SE_input <- cpueuse$SE - Q_extraSD_info$Value[Q_extraSD_info$Fleet == ifleet]
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
				} ## end plot
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
						file <- paste0("index8_q_vs_vuln_bio_", Fleet, ".png")
						caption <- paste0("Catchability vs. vulnerable biomass for fleet ", 
						Fleet, "<br> \n", "This plot should illustrate curvature of nonlinear catchability relationship<br> \n", 
						"or reveal patterns associated with random-walk catchability.")
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						q_vs_vuln_bio.fn(l=lang)
						dev.off(); eop()
					}
				} ## end print
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
				} ## edn plot
				if (print) {
					if (10 %in% subplots && !onepage) {
						file <- paste0("index10_resids_SE_total_", Fleet, ".png")
						caption <- paste0("Residuals of fit to index for ", Fleet, ".")
						if (error == 0) {
							caption <- paste0(caption, "<br>Values are (log(Obs) - log(Exp))/SE ", 
							"where SE is the total standard error including any ", "estimated additional uncertainty.")
						}
						else {
							caption <- paste0(caption, "<br>Values are based on the total standard error ", 
							"including any estimated additional uncertainty.")
						}
						plotinfo <- pngfun(file=file, caption=caption, lang=lang)
						index_resids.fn(option=1, l=lang)
						dev.off(); eop()
					} else if (10 %in% subplots && onepage) {
						index_resids.fn(option=1, l=lang)
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
						file <- paste0("index12_resids_SE_total_", Fleet, ".png")
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
				} ## end print
			} ## end include
		} ## end ifleet
		if (onepage && print) {dev.off(); eop()}
	} ## end j in  subplots
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
				plot(0, type="n", xlab=get.lab(1,ifleet), main=main, cex.main=cex.main, col=usecols[1], ylab=get.lab(8,ifleet), xlim=xlim, ylim=ylim)
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


## plotSS.pars--------------------------2025-08-20
##  Plot parameter fits and priors.
##  Based on 'r4ss::SSplotPars'
## ----------------------------------------r4ss|RH
plotSS.pars <- function (replist, plotdir=NULL, xlab="Parameter value", 
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
		fout = switch(lang, 'e' = file.path(plotdir,"english", file), 'f' = file.path(plotdir,"french", file) )
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
			abline(v=c(Pmin,Pmax),lty=3,col="darkgrey")
			addLabel(0.975,0.975, linguaFranca(gsub("_"," ",header),l), adj=c(1,1), cex=1, col="darkgreen", font=2)
			#addLabel(0.975,0.9, paste0("Init=",initval), adj=c(1,1), cex=1, col="darkgreen", font=2)
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
			par(xpd=TRUE)  #par(xpd=NA)
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
				if (parname=="BH_h") {
					#B0 = diag.mpd$B[1]; R0 = exp(diag.mpd$P[grep("R0",names(diag.mpd$P))]); h=finalval
					#shape1 = ((1-h)*B0) / (4*h*R0) ## alpha
					#shape2 = (5*h-1) / (4*h*R0)    ## beta
					#xmin = qbeta(0.21, shape1=shape2, shape2=shape1)
					#xmax = qbeta(0.99, shape1=shape2, shape2=shape1)
					xmin = 0.2; xmax=1
#if (parname=="BH_h"){browser();return()}
				} else {
					xmin <- qnorm(p=0.01, mean=finalval, sd=ifelse(parsd<=abs(finalval), parsd, abs(finalval)))
					xmax <- qnorm(p=0.99, mean=finalval, sd=ifelse(parsd<=abs(finalval), parsd, abs(finalval)))
				}
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
		postparname <- parname
		if (showpost) {
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
#.flush.cat(postparname, "   ", parname, "\n")
		if (is.null(xlim)) {
			if (fitrange & ((!is.na(parsd) && parsd != 0) | showpost)) {
				if (fitnudge<=0) {
					xmin <- max(Pmin, xmin, na.rm=TRUE)
					xmax <- min(Pmax, xmax, na.rm=TRUE)
				} else {
					xspan = abs(diff(c(xmin,xmax)))
					xmin  = xmin - (fitnudge * xspan)
					xmax  = xmax + (fitnudge * xspan)
#if (postparname=="BH_h"){browser();return()}
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
			ylim2 <- c(0, 1.2 * ymax)
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
			par(mfrow=c(nrows, ncols), mar=c(1,0.5,1,0), oma=c(2,2,0,1), mgp=c(2,0.5,0))
		}
		if (print | plot) {
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


## plotSS.profile-----------------------2024-04-03
## Commandeer r4ss function SSplotProfile
## Run 'repeatMPDs.r' before running this function.
## ----------------------------------------r4ss|RH
plotSS.profile <- function (summaryoutput, plot=TRUE, print=FALSE, models="all", 
	profile.string="steep", profile.label="Spawner-recruit steepness (h)", 
	exact=FALSE, ylab="Change in -log-likelihood", components=c("TOTAL", 
		"Catch", "Equil_catch", "Survey", "Discard", "Mean_body_wt", 
		"Length_comp", "Age_comp", "Size_at_age", "SizeFreq", 
		"Morphcomp", "Tag_comp", "Tag_negbin", "Recruitment", 
		"InitEQ_Regime", "Forecast_Recruitment", "Parm_priors", 
		"Parm_softbounds", "Parm_devs", "F_Ballpark", "Crash_Pen"), 
	component.labels=c("Total", "Catch", "Equilibrium catch", 
		"Index data", "Discard", "Mean body weight", "Length data", 
		"Age data", "Size-at-age data", "Generalized size data", 
		"Morph composition data", "Tag recapture distribution", 
		"Tag recapture total", "Recruitment", "Initital equilibrium recruitment", 
		"Forecast recruitment", "Priors", "Soft bounds", "Parameter deviations", 
		"F Ballpark", "Crash penalty"), minfraction=0.01, sort.by.max.change=TRUE, 
	col="default", pch="default", lty=1, lty.total=1, 
	lwd=2, lwd.total=3, cex=1, cex.total=1.5, xlim="default", 
	ymax="default", xaxs="r", yaxs="r", type="o", legend=TRUE, 
	legendloc="topright", pwidth=8, pheight=6, punits="in", 
	res=400, ptsize=10, cex.main=1, plotdir=NULL, add_cutoff=FALSE, 
	cutoff_prob=0.95, verbose=TRUE, lang=c("f","e"), ...) 
{
	if (print) {
		if (is.null(plotdir)) {
			stop("to print PNG files, you must supply a directory as 'plotdir'")
		}
		if (!file.exists(plotdir)) {
			if (verbose) 
				cat("creating directory:", plotdir, "\n")
			dir.create(plotdir, recursive=TRUE)
		}
	}
	if (length(components) != length(component.labels)) {
		stop("Inputs 'components' and 'component.labels' should have equal length")
	}
	n <- summaryoutput$n
	likelihoods <- summaryoutput$likelihoods
	if (is.null(likelihoods)) {
		stop("Input 'summaryoutput' needs to be a list output from SSsummarize\n", "and have an element named 'likelihoods'.")
	}
	pars <- summaryoutput$pars
	if (models[1] == "all") {
		models <- 1:n
	}
	else {
		if (!all(models %in% 1:n)) 
			stop("Input 'models' should be a vector of values from 1 to n=", n, " (for your inputs).\n")
	}
	if (exact) {
		parnumber <- match(profile.string, pars$Label)
	}
	else {
		parnumber <- grep(profile.string, pars$Label)
	}
	if (length(parnumber) <= 0) {
		stop("No parameters matching profile.string='", profile.string, "'", sep="")
	}
	parlabel <- pars$Label[parnumber]
#browser();return()
	if (length(parlabel) > 1) {
		stop("Multiple parameters matching profile.string='", profile.string, "':\n", paste(parlabel, collapse=", "), "\nYou may need to use 'exact=TRUE'.", sep="")
	}
	parvec <- as.numeric(pars[pars$Label == parlabel, models])
	cat("Parameter matching profile.string='", profile.string, "': '", parlabel, "'\n", sep="")
	cat("Parameter values (after subsetting based on input 'models'):\n")
	print(parvec)
	if (xlim[1] == "default") 
		xlim <- range(parvec)
	prof.table <- data.frame(t(likelihoods[likelihoods$Label %in% components, models]))
	names(prof.table) <- likelihoods[likelihoods$Label %in% components, ncol(likelihoods)]
	component.labels.good <- rep("", ncol(prof.table))
	for (icol in 1:ncol(prof.table)) {
		ilabel <- which(components == names(prof.table)[icol])
		component.labels.good[icol] <- component.labels[ilabel]
	}
	subset <- parvec >= xlim[1] & parvec <= xlim[2]
	for (icol in 1:ncol(prof.table)) {
		prof.table[, icol] <- prof.table[, icol] - min(prof.table[subset, icol])
	}
	if (ymax == "default") 
		ymax <- 1.05 * max(prof.table[subset, ])
	ylim <- c(0, ymax)
	column.max <- apply(prof.table[subset, ], 2, max)
	change.fraction <- column.max/column.max[1]
	include <- change.fraction >= minfraction  ## determines which columns to use
	nlines <- sum(include)
	message("\nLikelihood components showing max change as fraction of total change.\n", "To change which components are included, change input 'minfraction'.\n")
	print(data.frame(frac_change=round(change.fraction, 4), include=include, label=component.labels.good))
#browser();return()

	if (nlines == 0) {
		stop("No components included, 'minfraction' should be smaller.")
	}
	component.labels.used <- component.labels.good[include]
	prof.table            <- prof.table[order(parvec), include]	
	parvec                <- parvec[order(parvec)]
	change.fraction       <- change.fraction[include]
	if (nlines > 1) {
		if (sort.by.max.change) {
			neworder <- c(1, 1 + order(change.fraction[-1], decreasing=TRUE))
#browser();return()
			prof.table <- prof.table[, neworder]
			component.labels.used <- component.labels.used[neworder]
			col =  col[names(prof.table)]
		}
	}
	if (col[1] == "default") {
		col <- rich.colors.short(nlines)
	} else {
		col = rep(col,nlines)[1:nlines]
	}
	if (pch[1] == "default") {
		pch <- 1:nlines
	}
	lwd <- c(lwd.total, rep(lwd, nlines - 1))
	cex <- c(cex.total, rep(cex, nlines - 1))
	lty <- c(lty.total, rep(lty, nlines - 1))
	dots = list(...)
	if (is.null(dots$log))
		log = FALSE
	else
		log = dots$log

	## Main plotting routine
	plotprofile <- function(log, lang="e", print) {
		fout.e = paste0("profile_plot_likelihood-[",profile.string,"]",ifelse(log,"-ylog2",""),".png")
		for (l in lang) {
			changeLangOpts(L=l)
			if (print) {
				createFdir(lang=l, dir=plotdir)
				fout = switch(l, 'e' = file.path(plotdir, fout.e), 'f' = file.path(plotdir,"french", fout.e) )
				clearFiles(fout)
				png(filename=fout, width=pwidth, height=pheight, units=punits, res=res, pointsize=ptsize)
				expandGraph()
			}
#browser();return()
			if (log) {
				ylim = log2(GT0(ylim,eps=0.05))
				ylab = sub("^Change", "Log2 change", ylab)
			}
			plot(0, type="n", xlim=xlim, ylim=ylim, xlab=linguaFranca(profile.label,l), ylab=linguaFranca(ylab,l), yaxs=yaxs, xaxs=xaxs)
			abline(h=ifelse(log,log2(GT0(0,eps=0.05)),0), col="grey")
			if (add_cutoff) {
				abline(h=0.5 * qchisq(p=cutoff_prob, df=1), lty=2)
			}
			#matplot(parvec, prof.table, type=type, pch=pch, col=lucent(col,0.75), cex=cex, lty=lty, lwd=lwd, add=T)
			ytab = prof.table[,ncol(prof.table):1]
			if (log)
				ytab  = log2(GT0(ytab,eps=0.05))
			#matplot(parvec, ytab, type=type, pch=rev(pch), col=rev(lucent(col,0.75)), cex=rev(cex), lty=rev(lty), lwd=rev(lwd), add=T)
			matplot(parvec, ytab, type=type, pch=rev(pch), col=rev(col), cex=rev(cex), lty=rev(lty), lwd=rev(lwd), add=T)
			if (legend) 
				legend(legendloc, bty="n", legend=linguaFranca(component.labels.used,l), lwd=lwd, pt.cex=cex, lty=lty, pch=pch, col=col)
			box()
			if (print) dev.off(); eop()
		}
#browser();return()
	}
	if (plot) {
		expandGraph()
		plotprofile(log=log, lang=lang, print=print)
	}
	out <- data.frame(parvec=parvec, prof.table)
	names(out)[1] <- parlabel
	return(invisible(out))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.profile


## plotSS.rdevs-------------------------2023-08-15
##  Plot reruitment deviations
##  Based on 'r4ss::SSplotPars'
## ----------------------------------------r4ss|RH
plotSS.rdevs <- function (replist, subplots=1:3, plot=TRUE, print=FALSE, 
   add=FALSE, uncertainty=TRUE, minyr=-Inf, maxyr=Inf, 
   forecastplot=FALSE, col1="black", col2="blue", col3="green3", 
   col4="red", legendloc="topleft", labels=c("Year", "Asymptotic standard error estimate", 
   "Log recruitment deviation", "Bias adjustment fraction, 1 - stddev^2 / sigmaR^2"), 
   pwidth=8, pheight=6, punits="in", res=400, ptsize=10, 
   cex.main=1, plotdir="default", verbose=TRUE, outnam, lang="e") 
{
	oldpar = par(no.readonly=TRUE)
	fart <- function(opar) { if (any("windows"%in%names(dev.list()))) par(opar); eop() }
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
	recdevFore <- recdevForeTrue <- parameters[substring(parameters$Label, 1, 8) == "ForeRecr", ]
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
#browser();return()
			ylim <- range(c(recdevEarly$Value, recdev$Value, 
				recdevFore$Value)[goodyrs], na.rm=TRUE)

			## Recruitment deviation function
			recdevfunc <- function(uncertainty, l="e") {
				alldevs <- rbind(recdevEarly, recdev, recdevFore)[goodyrs,]
				#colvec <- c(rep(col3, nrow(recdevEarly)), rep(col1,nrow(recdev)), rep(col4, nrow(recdevFore)))[goodyrs]
				colvec <- c(rep(col3, nrow(recdevEarly)), rep(col1,nrow(recdev)), rep(col2,nrow(recdevLate)), rep(col4, nrow(recdevForeTrue)))[goodyrs]
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
				plot(Yr, Yr, type="n", xlab=linguaFranca(labels[1],l), ylab=linguaFranca(labels[3],l), xlim=xlim, ylim=ylim, cex.axis=1.2, cex.lab=1.5)
				abline(h=0, col="grey")
				if (uncertainty)
					arrows(Yr, recdev_lo, Yr, recdev_hi, length=0.03, code=3, angle=90, lwd=1.2, col=colvec)
				lines(Yr, val, lty=3)
				points(Yr, val, pch=16, col=colvec)
#browser();return()
				R0  = replist$parameters[grep("R0",rownames(replist$parameters)),"Value"]
				R0a = show0(round(R0,3),3)
				R0b = formatC(exp(R0),format="d",big.mark=switch(l,'e'=",",'f'=" "))
				addLabel(0.98,0.96,bquote(LN~italic(R)[0]== .(R0a)), cex=1, adj=c(1,0), col="blue")
				addLabel(0.98,0.93,bquote(italic(R)[0]== .(R0b)), cex=1, adj=c(1,0), col="blue")
				recdevs =  alldevs[,c("Yr","Value","Parm_StDev")]
				ttput(recdevs)
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


## plotSS.selex-------------------------2025-07-21
##  Plot selectivity curves with maturity ogive.
##  Source: R package 'r4ss' v.1.39.1
##  Based on 'r4ss::SSplotSelex'
## ----------------------------------------r4ss|RH
plotSS.selex <- function (replist, infotable=NULL, fleets="all", fleetnames="default", 
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
	fart <- function(opar) { if (any("windows"%in%names(dev.list()))) par(opar) }
	on.exit(fart(oldpar))

	## drag race herstory for browsing wayward data
	sumtingwong <- function(mess, debug, elvira) {
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
		clearFiles(file.path(plotdir, fnam))
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
		plt.selectivity(adbase, mainTitle=species.name, ptypes=ptypes, pngres=pngres, lang=lang, sobj=sobj, maxage=maxage, PIN=c(pwidth,pheight))
#browser();return()
	}
	if (!is.null(plotinfo)) 
		plotinfo$category <- "Sel"
	return(invisible(list(infotable=infotable2, plotinfo=plotinfo)))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.selex


## plotSS.stdres------------------------2026-02-20
## Plot standardised residuals -- three plots on one page.
## Modified from PBSawatea code in 'PBSscape.r'.
## Resolved sex combinations (for F & M only)
## Based on 'r4ss::SSplotComps'
## -----------------------------------------------
plotSS.stdres <- function(replist, kind="AGE", fleets="all",
   fleetnames="default", sexes="both", datonly=FALSE, aggregates_by_mkt = FALSE,
   labels=c("Length (cm)", "Age (yr)", "Year", "Observed sample size",
   "Effective sample size", "Proportion", "cm", "Frequency", "Weight", "Length",
   "(mt)", "(numbers x1000)", "Stdev (Age)", "Conditional AAL plot, ", "Size bin"),
   plot=TRUE, print=FALSE, type="Multinomial", useOSA=TRUE, usePearson=FALSE,
   ptypes="png", pngres=400, PIN=c(7,9), outnam, lang="e", ...)
{
	oldpar = par(no.readonly=TRUE)
	fart <- function(opar) { if (any("windows"%in%names(dev.list()))) par(opar); eop() }
	on.exit(fart(oldpar))

	if (sum(useOSA,usePearson)==0)
		stop("Choose either OSA or Pearson residuals (or both)")

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
	## Figure out sex combinations
	if (sexes[1] %in% c("all","comb")) {
		nsexes   = 1
		sexes    = list('comb'=1:2)
		titlesex = "(F+M)"
	} else if (sexes[1] %in% c("both","MF","FM")) {
		nsexes   = 2
		sexes    = list('female'=1, 'male'=2)
		titlesex = "(female left, male right)"
	} else if (nsexes==1 || sexes[1] %in% c("one","single")) {
		#sexes <- 0:nsexes
		nsexes   = 1
		sexes    = list('asex' = 1)
		titlesex = "(single sex)"
	} else if (nsexes > 1 && length(sexes) == 1) {
		nsexes = 1
		if (as.character(sexes) %in% c("0","comb","all")) {
			sexes  = list('comb'=1:2)
			titlesex <- "(sexes combined)"
			filesex <- "sex0"
		} else if (as.character(sexes) %in% c("1","F","female")) {
			sexes  = list('female'=1)
			titlesex <- "(female)"
			filesex <- "sex1"
		} else if (as.character(sexes) %in% c("2","M","male")) {
			sexes  = list('male'=2)
			titlesex <- "(male)"
			filesex <- "sex2"
		}
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
	dbase_kind <- dbase_kind[dbase_kind$Fleet %in% fleets & dbase_kind$sex %in% .su(unlist(sexes)), ]
	if (nrow(dbase_kind)==0) {
		message("Sum Ting Wong"); browser(); return()
	}

	## Replace std.res. generated by r4ss using calculation in function 'calcStdRes'
	## Note: Best to use r4ss calculation for Multinomial and Dirichlet-Multinomial
	#if (type %in% c("Multinomial","Fournier","Coleraine")) {
	if (type %in% c("Fournier","Coleraine")) {
		dbase_kind$Pearson_orig = dbase_kind$Pearson
		dbase_temp = calcStdRes(dbase_kind,type=type)
		dbase_kind$Pearson = dbase_temp$stdRes
	}
	## Try out 'one-step-ahead' (OSA) aka forecast quantile residuals [Trijoulet et al. 2023] (RH 240405)
	if (useOSA) {
		require(compResidual) ## only installed in R develop
#browser();return()

		## Process OSA residuals by fleet
		fbase = split(dbase_kind, dbase_kind$Fleet)
		fappy = lapply(fbase, function(xbase) {
			## Following the example on github for Gulf of Maine haddock,
			##  data need to be in year by age residual format
			OSAlist = list()
			ybase = split(xbase, xbase$Yr)
			yappy = lapply(ybase, function(sbase) {
				#sbase = data.frame()
				#for (ss in .su(abase$Sex)) {
 					#sbase = abase[is.element(abase$Sex, ss),]
				sexy = c("F","M"); names(sexy)=1:2
				slab = paste0(sexy[as.character(sbase$Sex)], pad0(sbase$Bin,2))
				sobs = sbase$Obs; names(sobs) = paste0("obs",slab)
				sexp = sbase$Exp; names(sexp) = paste0("pred",slab)
				sOSA = data.frame(Fleet=.su(sbase$Fleet), Year=.su(sbase$Yr), ESS=round(.su(sbase$effN)), t(sobs), t(sexp))
				#}
			})
			yOSA = do.call("rbind", lapply(yappy, data.frame, stringsAsFactors=FALSE))
			ifleet = paste0("fleet_", .su(yOSA$Fleet))
			fleet.name = fleets.all[.su(yOSA$Fleet)]
			OSAlist[[ifleet]][["osa_dat"]] = yOSA
			sexy = unlist(sexes)
			sexy = c(sexes, list(all=c(1,2)))  ## YTR 2024 : Nick Fisch requested that OSA be done on both sexes
			for (s in 1:length(sexy)) {
				ss  = sexy[[s]]
				sss = toupper(substring(names(sexy)[s],1,1))
				if (sss=="A"){ 
					sss = c("F", "M"); ii = "A"; iii = "sex_A"
				} else {
					ii = sss; iii = paste0("sex_", sss)
				}
				isex = paste0("sex_",sss)
				## extract observations
				obs <- yOSA[, grep(paste0("^obs",sss,collapse="|"), colnames(yOSA))]
				## Need to standardise? This changes OSA substantially.
				## Decided not to standardise obs because it's multiplied by ESS (one value for both sexes)
				#obs <- sweep(obs, 1, apply(obs,1,sum), "/")

				if (length(sss)>1) {
					mess = paste0("obs <- ",paste0(paste0("obs[,grep(\"", sss, "\",colnames(obs))]"), collapse=" + "))
					eval(parse(text=mess))
					colnames(obs) = sub(sss[1], "A", colnames(obs))
				}
#if (s==3) {browser();return()}
				## multiply by effective sample size and round:
				obs <- round(obs * yOSA[, "ESS"])
				#obs <- pmax(obs,1)
				## extract predictions
				pred <- yOSA[, grep(paste0("^pred",sss,collapse="|"), colnames(yOSA))] 
#browser(); return()
				if (length(sss)>1) {
					mess = paste0("pred <- ",paste0(paste0("pred[,grep(\"", sss, "\",colnames(pred))]"), collapse=" + "))
					eval(parse(text=mess))
					colnames(pred) = sub(sss[1], "A", colnames(pred))
				}
				## Need to standardise? (not sure yet)  PJS and I decided not to standardise (241007)
				#pred <- sweep(pred, 1, apply(pred,1,sum), "/")
				
				set.seed(123)
#browser();return()
				## calculate residuals:
				res <- resMulti(t(obs), t(pred))
				## Add names to sample number for plotting:
				colnames(res) <- yOSA[,"Year"]
				OSAlist[[ifleet]][["oas_obs"]][[iii]]  = obs
				OSAlist[[ifleet]][["oas_pred"]][[iii]] = pred
				OSAlist[[ifleet]][["oas_res"]][[iii]]  = res
				## Plot results if figures are requested
				#fout.e = paste0("osa.residuals.",sub("_","",ifleet),".",sub("_","",isex))
				fout.e = paste0("osa.residuals.",sub("_","",ifleet),".",sub("_","",iii))
				#fout.e = paste0("osa.residuals.",sub("_","",ifleet),".sex",paste0(sss,collapse=""))
				for (l in lang) {
					changeLangOpts(L=l)
					#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
					fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
					for (p in ptypes) {
						if (print && p=="eps") {
							clearFiles(paste0(fout,".eps"))
							postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
						} else if (print && p=="png"){
							clearFiles(paste0(fout,".png"))
							png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
						}
						lab.main = paste0(fleet.name, " - ", switch(ii, 'A'="Females + Males", 'B'="Both sexes", 'F'="Females", 'M'="Males", "Unknown"))
						plot_cres(res, maxLag=20, oma=c(0,0,1,0), mar=c(3,3,1,1), main=lab.main, bscale=0.4, lang=l)
#browser();return()
						#mtext(switch(ii, 'A'="Females + Males", 'B'="Both sexes", 'F'="Females", 'M'="Males", "Unknown"), outer=TRUE, side=3, line=0, cex=1.5)
						if (print && p %in% c("eps","png")) dev.off()
					} ## end p (ptypes) loop
				}; eop()
			} ## end s (sex) loop
			return(OSAlist)
		}) ## end fappy
		OSA = unlist(fappy, recursive=FALSE) ## remove redundant 1st level
		names(OSA) = sub("^[0-9]+\\.","",names(OSA))
		if (exists("redoData") && redoData)  ## (RH 250912) added to interact with 'sweaveMPD.Tnw'
			save("OSA", file=paste0(sub("MPD","OSA",basename(getwd())),".rda"))
	} ## end if useOSA (calculations)
#browser(); return()

	for (f in fleets) {
		fdbase = dbase_kind[is.element(dbase_kind$Fleet,f),]
		## Already has Pearson residuals computed
		if (is.null(ttcall(outnam)))
			outnam = "ageresFleet"
#browser();return()
		fnam = paste0(outnam,f) ##,sub("mf$","",gsub("\\+","",tolower(extract.between(titlesex,"(",")"))))) wtf?
		fout.e = fnam
		for (l in lang) {
			changeLangOpts(L=l)
			#fout = switch(l, 'e' = fout.e, 'f' = paste0("./french/",fout.e) )
			fout = switch(l, 'e' = paste0("./english/",fout.e), 'f' = paste0("./french/",fout.e) )
			for (p in ptypes) {
				if (print && p=="eps") {
					clearFiles(paste0(fout,".eps"))
					postscript(paste0(fout,".eps"), width=PIN[1], height=PIN[2], horizontal=FALSE,  paper="special")
				} else if (print && p=="png"){
					clearFiles(paste0(fout,".png"))
					png(paste0(fout,".png"), units="in", res=pngres, width=PIN[1], height=PIN[2])
				}
#browser();return()
				par(mfcol=c(3, sum(useOSA,usePearson)*nsexes), mai=c(0.45,0.3,0.1,0.1), omi=c(0,0.25,0.4,0), mgp=c(2,0.75,0))
				#par(mfcol=c(3,ifelse(useOSA,2,1)*nsexes), mai=c(0.45,0.3,0.1,0.1), omi=c(0,0.25,0.4,0), mgp=c(2,0.75,0))
				#par(mfcol=c(3,2), mai=c(0.45,0.3,0.1,0.1), omi=c(0,0.25,0.4,0), mgp=c(2,0.75,0))
				for ( s in 1:nsexes) {
					ss = sexes[[s]]
					sfdbase = fdbase[is.element(fdbase$sex,ss),]
					if (usePearson) {
						plt.ageResids(sfdbase, resfld="Pearson", main="", lang=l, ...)   ## by age class  (RH 240408 eclipse)
						addLabel(0.05,0.95,"Pearson",col="red",cex=1.2,adj=0)
						plt.yearResids(sfdbase, resfld="Pearson", lang=l, ...)           ## by year
						plt.cohortResids(sfdbase, resfld="Pearson", lang=l, ...)         ## by cohort (year of birth)
					}
					if (useOSA) {
						mess = paste0("osares = OSA$fleet_", f, "$oas_res$sex_", toupper(substring(names(sexes)[s],1,1)) )
						eval(parse(text=mess))
						sfdbase$OSA = as.vector(rbind(osares, rep(NA, ncol(osares))))
						plt.ageResids(sfdbase, resfld="OSA", main="", lang=l, ...)   ## by age class  (RH 240408 eclipse)
						addLabel(0.05,0.95,"OSA",col="blue",cex=1.2,adj=0)
						plt.yearResids(sfdbase, resfld="OSA", lang=l, ...)           ## by year
						plt.cohortResids(sfdbase, resfld="OSA", lang=l, ...)         ## by cohort (year of birth)
					}
				}
#browser();return()
				#title = paste(toUpper(tolower(gsub("[[:punct:]]+", " ", FleetNames[f]))),titlesex)
				#title = paste(gsub("URVEY", "urvey", gsub("ISHERY", "ishery", gsub("YNOPTIC", "ynoptic", gsub("[[:punct:]]+", " ", FleetNames[f])))),titlesex)
				title = paste(gsub("_", " ", FleetNames[f]),titlesex)
				mtext(linguaFranca(title,l), side=3, outer=TRUE, line=0.25, cex=1.5)
				mtext(linguaFranca("Standardised Residuals",l), side=2, outer=TRUE, line=0, cex=1.5)
#browser();return()
				if (print && p %in% c("eps","png")) dev.off()
			} ## end p (ptypes) loop
		}; eop()
	}
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.stdres


## plotSS.stock.recruit-----------------2021-08-26
## Plot stock-recruitment function (based on MPDs)
## Based on 'r4ss::SSplotSpawnrecruit'
## ----------------------------------------r4ss|RH
## xLimSR and yLimSR fixed here for YMR to have Run 26 and 27 figs
##  on same scales. Use these first two values to scale to data:
# xLimSR =c(0, max(obj$B$SB))
# yLimSR=c(0, max(obj$B$R, na.rm=TRUE))
#xLimSR=c(0, max(c(max(obj$B$SB),45000)))   # so it draw bigger if necessary
#yLimSR=c(0, max(c(max(obj$B$R, na.rm=TRUE),55000)))
## -----------------------------------------------
plotSS.stock.recruit <- function (replist, subplot = 1:3, add = FALSE, plot = TRUE, print = FALSE, 
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

	## Additional columns (RH 230227)
	B$Rexp = srFun(B$SB)
#browser();return()

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


## plotSS.ts----------------------------2025-08-22
##  Plot SS time series (Bt, BtB0, Age-0 recruits)
##  Modified r4ss function 'SSplotTimeseries'.
## ----------------------------------------r4ss|RH
plotSS.ts <- function (replist, subplot, add=FALSE, areas="all", areacols="default", 
   areanames="default", forecastplot=TRUE, uncertainty=TRUE, 
   bioscale=1, minyr=-Inf, maxyr=Inf, plot=TRUE, print=FALSE, 
   plotdir="default", verbose=TRUE, btarg="default", minbthresh="default", 
   xlab="Year", labels=NULL, punits="in", ptsize=10, cex.main=1, sobj=NULL,
   res=400, outnam, PIN=c(9,9), lang="e") 
{
	oldpar = par(no.readonly=TRUE)
	fart <- function(opar) { if (any("windows"%in%names(dev.list()))) par(opar) }
	on.exit(fart(oldpar))

	if (missing(subplot)) 
		stop("'subplot' input required")
	if (length(subplot) > 1) 
		stop("function can only do 1 subplot at a time")
	options(big.mark=",")

	pngfun <- function(file, caption=NA, lang="e") {
		metas = c("\\", "/", ":", "*", "?", "\"", "<", ">", "|")
		meta  = paste0(c("[",metas,"]"),collapse="")
		fnam  = gsub("[_ ]+", "_", gsub(meta,"",file))
		createFdir(lang, dir=plotdir)
		changeLangOpts(L=lang)
		fout = switch(lang, 'e' = file.path(plotdir,"english", file), 'f' = file.path(plotdir,"french", file) )
		clearFiles(fout)
		##.flush.cat("Figure file:", fout, "\nsaved to:", plotdir, "\n")
		png(filename=fout, width=PIN[1], height=PIN[2], units=punits, res=res, pointsize=ptsize)
		plotinfo <- rbind(plotinfo, data.frame(file=fout, caption=caption))
		return(plotinfo)
	}
	plotinfo <- NULL
	if (is.null(labels)) {
		## 1 Total biomass (mt) with forecast
		## 2 Total biomass by area (spatial models only)
		## 3 Total biomass (mt) at beginning of spawning season with forecast
		## 4 Summary biomass (mt) with forecast
		## 5 Summary biomass (mt) by area (spatial models only)
		## 6 Summary biomass (mt) at beginning of season 1 with forecast
		## 7 Spawning output with forecast with ~95% asymptotic intervals
		## 8 Spawning output by area (spatial models only)
		## 9 Relative spawning output with forecast with ~95% asymptotic intervals
		## 10 Relative spawning output by area (spatial models only)
		## 11 Age-0 recruits (1,000s) with forecast with ~95% asymptotic intervals
		## 12 Age-0 recruits by area (spatial models only)
		## 13 Fraction of recruits by area (spatial models only)
		## 14 Age-0 recruits (1,000s) by birth season with forecast
		## 15 Fraction of total Age-0 recruits by birth season with forecast
		labels <- c("Total biomass (t)", "Total biomass (t) by area", "Total biomass (t) at beginning of season",
			"Summary biomass (t)", "Summary biomass (t) by area", "Summary biomass (t) at beginning of season",
			"Spawning biomass (t)", "Spawning biomass (t) by area", "Depletion (Bt / B0)",
			"Depletion (Bt / B0) by area", "Age-0 recruits (1000s)", "Age-0 recruits (1000s) by area",
			"Fraction recruits by area", "Age-0 recruits (1000s) by birth season", "Fraction Age-0 recruits by birth season")
	}
	SS_versionshort <- replist$SS_versionshort
	timeseries      <- replist$timeseries
	timeseries2     <- replist$catch          ## (RH 250325) vulnerable biomass is reported here (250821: but it looks sketchy by Area and Fleet)
	nseasons        <- replist$nseasons
	spawnseas       <- replist$spawnseas
	birthseas       <- replist$birthseas
	startyr         <- replist$startyr
	endyr           <- replist$endyr
	projyr          <- endyr + replist$N_forecast_yrs  ## use if forecastplot is TRUE (RH 250102)
	nsexes          <- replist$nsexes
	nareas          <- replist$nareas
	derived_quants  <- replist$derived_quants
	seasfracs       <- replist$seasfracs
	B_ratio_denominator <- replist$B_ratio_denominator
	recruitment_dist    <- replist$recruitment_dist
	if (btarg == "default") 
		btarg <- replist$btarg
	if (length(minbthresh)==1 && minbthresh == "default") 
		minbthresh <- replist$minbthresh
	if (areacols[1] == "default") {
		#if (exists("narea")) nareas = narea  ## hokey fix for now
		areacols <- areabgs <- rich.colors.short(nareas)
		if (nareas <= 5) {
			defcols <- c("blue", "green4", "red", "darkgreen", "purple")  ## 5ABC 3CD 5DE (POP 2023)
			defbgs  <- c("cyan", "green", "pink", "yellow", "thistle")    ## 5ABC 3CD 5DE (POP 2023)
			defcols <- c("blue", "red", "green4", "orange", "purple")     ## 5ABC 5DE 3CD (SGR 2025)
			defbgs  <- c("cyan", "pink", "green", "yellow", "thistle")    ## 5ABC 5DE 3CD (SGR 2025)
			areacols = defcols[1:nareas]
			areabgs  = defbgs[1:nareas]
			gearcols = defcols[1:ngear]
			gearbgs  = defbgs[1:ngear]
		}
		if (nareas > 5) {  ## narea = proxy for # fisheries (for now)
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
	ts  <- timeseries
	ts2 <- timeseries2 ## (RH 250326)
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
			ylab <- labels[7]
#browser();return()
			#ylab <- labels[subplot]
		}
		if (subplot %in% 9:10) {
			B0 = ts$SpawnBio[grep("TIME",ts$Era)[1]]  ## do not want 'VIRG' or 'INIT'
			#B0 = ts$SpawnBio[!is.na(ts$SpawnBio)][1]
			yvals <- ts[,"SpawnBio",drop=FALSE] / B0 # $SpawnBio/ts$SpawnBio[!is.na(ts$SpawnBio)][1]
			ylab <- labels[9]
		}
		if (subplot %in% 11:15) {
			yvals <- ts[,"Recruit_0",drop=FALSE] ##$Recruit_0
			ylab <- labels[11]
			if (all(yvals[ts$Era == "VIRG",] == 0 & max(ts$Seas == 1))) {
				yvals[ts$Era == "VIRG",] <- derived_quants["Recr_Virgin", "Value"]
			}
			if (all(yvals[ts$Era == "INIT",] == 0 & max(ts$Seas == 1))) {
				yvals[ts$Era == "INIT",] <- derived_quants["Recr_Unfished", "Value"]
			}
		}
		if (subplot %in% c(13, 15)) 
			ylab <- labels[13]
		if (subplot %in% c(101)){
			ylab = "Biomass comparisons (tonnes)"
			ysub  = ts[,c("Area","Bio_all","SmryBio_SX:1_GP:1","SmryBio_SX:2_GP:1","SpawnBio")] * bioscale
			ylist = split(ysub[,-1,drop=FALSE], ysub[,"Area"])
			yvals = sapply(1:ncol(ylist[[1]]), function(i){Reduce("+", lapply(ylist, "[[", i))})
			yvals = as.data.frame(yvals)
			dimnames(yvals) = list(year=.su(ts$Yr),vals=colnames(ylist[[1]]))
		}

		if (any(grepl("^Hrate:",colnames(ts)))) { Fmeth = 1 } else { Fmeth = 3 }  ## temporary fix

		if (subplot %in% c(102,103)){
			cfleets  = unique(grep("^Hrate:|^F:",colnames(ts),value=T))  ## commercial  fleets
			ncfleet  = gsub("[Fu]:_", "", cfleets)
			ncfleets = length(cfleets)
			if (narea>1)
				ts$Fleet = as.numeric(rep(ncfleet, each=length(.su(ts$Yr))))
			for (j in 1:ncfleets) {
				jj = ncfleet[j]
				if (Fmeth==1)
					ts[,paste0("u",jj)]  = ts[,paste0("Hrate:_",jj)]
				else 
					ts[,paste0("u",jj)]  = 1 - exp(-ts[,paste0("F:_",jj)])
				uzero = ts[,paste0("u",jj)] <= 0 | is.na(ts[,paste0("u",jj)]) ## (RH 250321) added is.na condition
				## (RH 250821) disable VB from catch (looks sketchy); use narea>100 to disable
				if (narea>100 && exists("ts2") && !is.null(ts2) && "vuln_bio"%in%colnames(ts2)) {  ## (RH 250326) get vuln_bio from the catch series
					zf1 = is.element(ts$Fleet,jj)
					zy1 = is.element(ts$Yr,.su(ts2$Yr))
					zf2 = is.element(ts2$Fleet,jj)
					zy2 = is.element(ts2$Yr,.su(ts$Yr))
					ts[zf1 & zy1, paste0("V",jj)] = ts2[zf2 & zy2, "vuln_bio"]
					zna = is.na(ts[,paste0("V",jj)])
					ts[zna,paste0("V",jj)] = 0
				} else {
					ts[uzero,paste0("V",jj)]  = 0
					## Note: 'sel(B)' appears to be the catch (='obs_cat') in 'ts' object; V=C/u
					ts[!uzero,paste0("V",jj)] = ts[!uzero,paste0("sel(B):_",jj)] / ts[!uzero,paste0("u",jj)]
					ts[1,paste0("V",jj)] = ts[2,paste0("V",jj)] = 0
				}
			}
#browser();return()
			if (subplot %in% c(102)){
				ysub  = ts[,c("Area","Bio_all",paste0("V",ncfleet))] * bioscale
				ylab = "Biomass (tonnes)"
			}
			if (subplot %in% c(103)){
				ysub  = ts[,c("Area", paste0("u",ncfleet))]
				ylab = "Harvest Rate"
			}
			ylist = split(ysub[,-1,drop=FALSE],ysub[,"Area"])
			yvals = sapply(1:ncol(ylist[[1]]), function(i){Reduce("+", lapply(ylist, "[[", i))})
			yvals = as.data.frame(yvals)
			dimnames(yvals) = list(year=.su(ts$Yr),vals=colnames(ylist[[1]]))
			yvals[yvals<=0 | is.na(yvals)] = NA
#browser();return()
		}
		if (!is.element(subplot, c(1:15,101:103))) {
			stop("subplot should be a value from 1 to 15 (r4ss) or 101 to 103 (PBSsynth)")
		}
		main = ylab
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
#browser();return()
			}
			if (subplot %in% c(1, 4, 7, 11, 13)) {
				#yvals2 <- rep(NA, length(ts$YrSeas))
				yvals2 <- as.data.frame(array(NA, dim=c(length(ts$YrSeas),1), dimnames=list(rownames(ts),"sumting")))
				for (iyr in 1:nrow(yvals)) {
					y <- ts$YrSeas[iyr]
					yvals2[iyr,] <- sum(yvals[ts$YrSeas == y,])
				}
				if (subplot == 13) {
					## yvals = "Recruit_0"
#browser();return()
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
			ymax <- max(yvals, na.rm=TRUE)
			if (subplot == 10) {
				ymax = 1
				for (iarea in 1:nareas) {
					yvals <- ts[ts$Area == iarea,"SpawnBio",drop=FALSE]/(ts$SpawnBio[ts$Area == iarea & ts$Seas == spawnseas][1])
					ymax <- max(ymax, max(yvals, na.rm=TRUE))
				}
#browser();return()
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
#browser();return()
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
#browser();return()
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
			expandGraph(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75,0.5,0))
			#plot(yrvals, yvals[plot1 | plot2 | plot3], type="n", xlab=xlab, xlim=xlim, ylim=c(0, 1.05 * ymax), yaxs="i", ylab=ylab, main=main, cex.main=cex.main, font.main=1)
#browser();return()
			if (subplot %in% c(103)) ylim = c(-ymax*0.075, 1.02*ymax)
			else                     ylim = c(0, 1.05 * ymax)
			plot(0, 0, type="n", xlim=xlim, ylim=ylim, yaxs=ifelse(subplot%in%c(103),"i","i"), xlab=linguaFranca(xlab,l), ylab=linguaFranca(ylab,l), main=NULL, cex.main=cex.main, cex.lab=1.4, font.main=1, yaxt="n")
			axis(1, at=intersect(seq(1900,3000,5),xlim[1]:xlim[2]), tcl=-0.2, labels=FALSE)
			axis(2, at=pretty(ylim,n=10), labels=format(pretty(ylim,n=10), big.mark=options()$big.mark,trim=T), cex=1.4)  ## nMark nZero -- make sure Options()$big.mark has been set
		}
		if (subplot %in% c(9, 10)) {
			addtarg <- function() {
				if (btarg > 0 & btarg < 1) {
					abline(h=btarg, col="green4",lty=2)
					text(max(startyr, minyr) + 1, btarg + 0.02, linguaFranca(ifelse(subplot%in%c(9,10),paste0(btarg,"B0"),labels[10]),l), adj=0, col="green4")
				}
				if (all(minbthresh > 0) & all(minbthresh < 1)) {
					for (i in 1:length(minbthresh)) {
						ii = sort(minbthresh)[i]
#browser();return()
						icol = switch(i, "red", "blue", "green4")
						abline(h=ii, col=icol, lty=2)
						text(max(startyr, minyr) + 1, ii + 0.02, linguaFranca(ifelse(subplot%in%c(9,10),paste0(ii,"B0"),labels[10]),l), adj=0, col=icol)
					}
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
				points(ts$Yr[plot1], yvals[plot1,], pch=17, col=mycol)
				lines(ts$Yr[plot2], yvals[plot2,], type=mytype, col=mycol)
				points(ts$Yr[plot3], yvals[plot3,], pch=19, col=mycol)
			}
			#legend("topright", legend=linguaFranca(paste("Season", birthseas),l), lty=1, pch=1, col=seascols, bty="n")
			addLegend(0.975, 0.975, legend=linguaFranca(paste("Season", birthseas),l), lty=1, pch=1, col=seascols, bty="n", xjust=1, yjust=1)
		}
		else {
			if (subplot %in% c(1,4,7,9,11,14,15, 101:103)) 
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
				}
				if (subplot %in% 9:10) {
					plot1 <- NULL
					#plot2[3] <- FALSE ## WTF?
				}
				mycol <- areacols[iarea]
				mybg  <- areabgs[iarea]
#print(c(iarea, mycol, mybg))
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
					for (yy in 1:ncol(yvals)) {
						#points(ts$YrSeas[plot1], yvals[plot1,yy], pch=19, col=mycol)
						#lines(ts$YrSeas[plot2], yvals[plot2,yy], type=mytype, col=mycol)
						points(ts$YrSeas[plot1], yvals[plot1,yy], pch=24, col=mycol, bg=mybg, cex=1.2)
						lines(ts$YrSeas[plot2], yvals[plot2,yy], type="l", col=mycol)
						points(ts$YrSeas[plot2], yvals[plot2,yy], pch=21, col=mycol, bg=mybg)
						points(ts$YrSeas[plot3], yvals[plot3,yy], pch=17, col=mycol)
					}
					legtxt = c(legtxt, paste0("SS3 ", switch(sp, '7'="Bt", '8'="Bt by area", '9'="Bt/B0", "age-0 recruits")))
					leglty = c(leglty, 1)
					legcol = c(legcol, mycol)
#if(iarea==3) {browser();return()}
				}
				else {
					for (j in 1:ncol(yvals)) {
						jj = colnames(yvals)[j]
						if (grepl("SX:",jj)) {
							col.sex = ifelse(grepl("SX:1",jj),"orange",.colBlind["bluegreen"])
							points(ts$YrSeas[plot1], yvals[plot1,j], pch=24, col=col.sex, bg=colorspace::lighten(col.sex, amount=0.25)) #lucent(mycol,0.2))
							lines(ts$YrSeas[plot2|plot3], yvals[plot2|plot3,j], col=col.sex, lwd=2)
						} else {
							#col.pch = ifelse(grepl("Bio_all",jj),"black",mycol)
							#col.fleet = c("blue","green3","red","purple") ## may need more colours
							#col.pch = switch(jj, 'Bio_all'="black", 'V1'=col.fleet[1], 'V2'=col.fleet[2], 'V3'=col.fleet[3], 'V4'=col.fleet[4], 'u1'=col.fleet[1], 'u2'=col.fleet[2], 'u3'=col.fleet[3], 'u4'=col.fleet[4], mycol)
							col.fleet = gearcols
							bg.fleet  = gearbgs
							## This is so problematic:
							col.pch = switch(jj, 'Bio_all'="black", 'V1'=col.fleet[1], 'V2'=col.fleet[2], 'V3'=col.fleet[3], 'V10'=col.fleet[4], 'V11'=col.fleet[5], 'u1'=col.fleet[1], 'u2'=col.fleet[2], 'u3'=col.fleet[3], 'u10'=col.fleet[4], 'u11'=col.fleet[5], mycol)
							col.bg = switch(jj, 'Bio_all'="gainsboro", 'V1'=bg.fleet[1], 'V2'=bg.fleet[2], 'V3'=bg.fleet[3], 'V10'=bg.fleet[4], 'V11'=bg.fleet[5], 'u1'=bg.fleet[1], 'u2'=bg.fleet[2], 'u3'=bg.fleet[3], 'u10'=bg.fleet[4], 'u11'=bg.fleet[5], mybg)
							#col.pch = col.fleet[j]
							#col.bg  = bg.fleet[j]
#print(c(col.pch, mycol))
							#points(ts$YrSeas[plot1], yvals[plot1,j], pch=22, col=col.pch, bg=colorspace::lighten(col.pch, amount=0.25)) #lucent(mycol,0.2))
							#lines(ts$YrSeas[plot2], yvals[plot2,j], col="gainsboro", lwd=2)
							#points(ts$YrSeas[plot2], yvals[plot2,j], pch=21, col=col.pch, bg=colorspace::lighten(col.pch, amount=0.25)) #lucent(mycol,0.2))
							#points(ts$YrSeas[plot3], yvals[plot3,j], pch=24, cex=0.8, col=col.pch, bg=colorspace::lighten(col.pch, amount=0.25)) #lucent(mycol,0.2))
							points(ts$YrSeas[plot1], yvals[plot1,j], pch=24, col=col.pch, bg=col.bg) ## equilbrium start
							lines(ts$YrSeas[plot2], yvals[plot2,j], col=col.pch, lwd=2)
							points(ts$YrSeas[plot2], yvals[plot2,j], pch=21, col=col.pch, bg=col.bg) ## main reconstruction
							points(ts$YrSeas[plot3], yvals[plot3,j], pch=24, cex=0.8, col=col.pch, bg=col.bg) ## forecast
#if(j==4){browser();return()}
#browser();return()
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
						if (subplot %in% c(7,9) || colnames(yvals)[j] %in% c("SpawnBio")) {
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
					} ## end j loop
				}
			} ## end iarea loop
			if (nareas > 1 & subplot %in% c(2, 3, 5, 6, 8, 10, 12, 13)) {
				if (exists("area.names") && narea==length(areas))  ## need to source 'initialise.r' first
					areanames = area.names
				addLegend(0.975, 0.975, legend=linguaFranca(areanames[areas],l), lty=1, seg.len=3, pch=21, col=areacols[areas], pt.bg=areabgs[areas], bty="n", xjust=1, yjust=1)
#browser();return()
				#addLegend (0.975, 0.975, legend=linguaFranca(legtxt,l), lty=leglty, col=legcol, lwd=2, seg.len=3, bty="n", xjust=1, yjust=1)
			}
			if (subplot %in% c(101:103)) {
				years = if(forecastplot) startyr:projyr else startyr:endyr  ## (RH 250102)
				cdat  = ts[is.element(ts$Yr,years),c("Area",grep("retain\\(B)",colnames(ts),value=TRUE))]
				clist = split(cdat[,-1,drop=FALSE],cdat[,"Area"])
				catch = sapply(1:ncol(clist[[1]]), function(i){Reduce("+", lapply(clist, "[[", i))})
				catch =  as.data.frame(catch[,,drop=FALSE])
#browser();return()
				if (ncol(catch)>1)
					cumcat = t(apply(catch[,,drop=FALSE],1,cumsum))
				else
					cumcat = catch
				dimnames(catch) = dimnames(cumcat) = list(year=years,gear=gear.names)
				catch$total = apply(catch,1,sum)
				#ts$catch = apply(ts[,grep("retain\\(B)",colnames(ts))],1,sum)
				base = 0
				if (subplot %in% c(103)){
					cumcat.orig = cumcat
					#catch = sapply(catch, scaleVec, Tmin=ylim[1], Tmax=0)
					#catch = t(apply(catch, 1, scaleVec, Tmin=ylim[1], Tmax=0))
					#dimnames(catch) = dimnames(catch.orig)
					cumcat = scaleVec(cumcat.orig, Tmin=ylim[1], Tmax=0)
					base = ylim[1]
				}
#browser();return()
				#drawBars(ts$Yr[plot2|plot3], ts$catch[plot2|plot3], col="gainsboro", fill="hotpink", width=1, base=base)
				## Quick fix for now -- plot total catch then first biggest fleet catch
				#drawBars(ts$Yr[plot2|plot3], catch[plot2|plot3,"total"], col="red", fill="red", width=1, base=base)
				#drawBars(ts$Yr[plot2|plot3], catch[plot2|plot3,1], col="red", fill="pink", width=1, base=base)
				if (subplot %in% c(101)) {
					drawBars(years, catch[,"total"], col="grey30", fill="yellow", width=1, base=base)
				} else {
					for (a in 1:ngear) {
						abase = if (a==1) base else cumcat[,gear.names[a-1]]
						ycat  = cumcat[,gear.names[a]]
						drawBars(x=years, y=ycat, col=gearcols[a], fill=gearbgs[a], width=1, base=abase)
					}
				}
				legtxt = NULL
				if (subplot %in% 101) {
					sumage  = replist$BioSmry_age
					legtxt = c("Total biomass","Female summary biomass", "Male summary biomass", "Spawning biomass", "Coastwide catch")
					if (!is.null(sumage)) {
						legtxt[2:3] = sub("Female", convUTF("\\u{2640}"), legtxt[2:3])
						legtxt[2:3] = sub("Male", convUTF("\\u{2642}"), legtxt[2:3])
						legtxt[2:3] = paste0(legtxt[2:3], " (a", convUTF("\\u{2265}"), sumage, "y)")
					}
					legcol = c("black", "orange", .colBlind["bluegreen"], mycol, "black")
					legbg  = c("gainsboro", NA, NA, mybg, "yellow")
					#legpch = c(21,95,95,21,22)
					legpch = c(21,NA,NA,21,22)
					leglty = c(NA,1,1,NA,NA)
					#legpch = c(21,151,151,21,22)  ## newer R versions don't support as many symbols
				}
				else if (subplot %in% 102) {
					gear.area.names = paste0(gear.names, " (", area.names, ")")  ## too wordy (confusing)
					legmain = ifelse(nareas==1, "Single-area model", "Multi-area model")
					legtxt = c("Total biomass", paste0("Vulnerable biomass - ", gear.names), paste0("Catch - ", gear.names))
					legcol = c("black", col.fleet[1:ncfleets], col.fleet[1:ncfleets])
					legbg  = c("gainsboro", bg.fleet[1:ncfleets], bg.fleet[1:ncfleets])
					legpch = c(21,rep(21,ngear),rep(22,ngear))
				}
				else if (subplot %in% 103) {
					gear.area.names = paste0(gear.names, " (", area.names, ")")  ## too wordy (confusing)
					legmain = ifelse(nareas==1, "Single-area model", "Multi-area model")
					legtxt = c(paste0("Harvest rate - ", gear.names), paste0("Catch (t) - ", gear.names))
					legcol = c(col=col.fleet[1:ncfleets], col.fleet[1:ncfleets])
					legbg  = c(bg.fleet[1:ncfleets], bg.fleet[1:ncfleets])
					legpch = c(rep(21,ngear),rep(22,ngear))
				}
				if (!is.null(legtxt)){
					if (subplot %in% c(101)) {
						if (forecastplot) {
							zfor = grep("Total|Spawning", legtxt); zoff = grep("Total|Spawning", legtxt, invert=T)
							forpch = legpch; forpch[zfor] = 21; forpch[zoff] = NA
							fortxt = rep(paste0(rep(" ",2*max(sapply(legtxt,nchar))),collapse=""),length(legtxt))
							addLegend(ifelse(print,0.965,0.95), 0.975, legend=fortxt, pch=forpch, col=legcol, pt.bg=legbg, pt.cex=1.5, pt.lwd=1,  bty="n", xjust=1, yjust=1, cex=1.2)
							legpch[zfor] = 24
							addLegend(0.975, 0.975, legend=linguaFranca(legtxt,l), pch=legpch, lty=leglty, lwd=2, col=legcol, pt.bg=legbg, pt.cex=1.5, pt.lwd=1, bty="n", xjust=1, yjust=1, cex=1.2)
#browser();return()
						} else {
							addLegend(0.97, 0.975, legend=linguaFranca(legtxt,l), pch=legpch, col=legcol, pt.bg=legbg, bty="n", xjust=1, yjust=1, cex=1.2)
						}
					}
					if (subplot %in% c(102, 103)) {
						if (forecastplot) {
							if (subplot %in% c(102)) {
								zfor = grep("Total|Vulnerable", legtxt); zoff = grep("Total|Vulnerable", legtxt, invert=T)
							}
							if (subplot %in% c(103)) {
								zfor = grep("Harvest", legtxt); zoff = grep("Harvest", legtxt, invert=T)
							}
							forpch = legpch; forpch[zfor] = 21; forpch[zoff] = NA
							#fortxt = rep(paste0(rep(" ",1.75*max(sapply(legtxt,nchar))),collapse=""),length(legtxt))
							#addLegend(0.925, 0.995, legend=fortxt, pch=forpch, col=legcol, pt.bg=legbg, bty="n", xjust=1, yjust=1, cex=1, title=legmain, title.font=2, title.adj=0)
							#legpch[zfor] = 24
							#addLegend(0.975, 0.995, legend=linguaFranca(legtxt,l), pch=legpch, col=legcol, pt.bg=legbg, bty="n", xjust=1, yjust=1, cex=1, title="")
							## Trick is to plot legend twice, using white text the first time
							addLegend(0.960, 0.995, legend=linguaFranca(legtxt,l), pch=forpch, col=legcol, pt.bg=legbg, bty="n", xjust=1, yjust=1, cex=1, title=legmain, title.font=2, title.adj=0, title.col="black", text.col="white")
#browser();return()
							legpch[zfor] = 24
							addLegend(0.975, 0.995, legend=linguaFranca(legtxt,l), pch=legpch, col=legcol, pt.bg=legbg, bty="n", xjust=1, yjust=1, cex=1, title="")
						} else {
							addLegend(0.975, 0.995, legend=linguaFranca(legtxt,l), pch=legpch, col=legcol, pt.bg=legbg, bty="n", xjust=1, yjust=1, cex=1, title=legmain, title.font=2, title.adj=0, title.col="black", text.col="black")
						}
					}
				}
				box()
			}
		}
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


## plotSS.yield ------------------------2026-02-17
## Plot the equilibrium yield curve for an MPD run
## r4ss Function: SSplotYield
## ----------------------------------------r4ss|RH
plotSS.yield=function (replist, subplots=1:5, 
   refpoints=c("MSY", "Btgt", "SPR", "Current", "LRP", "USR"), 
   add=FALSE, plot=TRUE, print=FALSE, labels=c("Fraction unfished",
   "Equilibrium yield (t)", "Total biomass (t)", "Surplus production (t)",
   "Yield per recruit (kg)", "Spawning output"), col="black", col2="black",
   lty=1, lwd=2, cex.main=1, pwidth=6.5, pheight=5, punits="in", 
   res=400, ptsize=10, plotdir="./english", verbose=TRUE, lang=c("f","e")) 
{
	if (print)
		createFdir(lang)
	plotinfo <- NULL
	equil_yield <- replist[["equil_yield"]]
	equil_yield <- equil_yield[equil_yield[["SPRloop"]] != 3, ]
	equil_yield <- equil_yield[order(equil_yield[["Depletion"]], decreasing=FALSE), ]
	if ("Tot_Catch" %in% names(equil_yield)) {
		equil_yield[["Catch"]] <- equil_yield[["Tot_Catch"]]
	}
	nareas <- replist[["nareas"]]
	nseasons <- replist[["nseasons"]]
	timeseries <- replist[["timeseries"]]
	SSB0 <- replist[["derived_quants"]]["SSB_Virgin", "Value"]
	yieldfunc <- function(refpoints=NULL) {
		if (!add) {
			expandGraph(mar=c(3.5,3.5,1,1), mgp=c(2,0.5,0))
			plot(0, type="n", xlim=c(0, max(equil_yield[["Depletion"]], 1, na.rm=TRUE)), 
				ylim=c(0, max(equil_yield[["Catch"]], na.rm=TRUE)), xlab=labels[1], ylab=labels[2],
				cex.axis=1.2, cex.lab=1.5)
			abline(h=0, col="grey")
			abline(v=0, col="grey")
		}
		lines(equil_yield[["Depletion"]], equil_yield[["Catch"]], lwd=lwd, col=col, lty=lty)
		#colvec <- c(4, 2, 3, 1)  ## default
		colvec <- c("red", "blue", "purple", "slategray", .colBlind[c("redpurple", "bluegreen")] )  ## MSY, Btgt, SPR, Current, LRP, USR
		ltyvec <- c(1, 5, 3, 4, 3, 2)
		if ("MSY" %in% refpoints) {
			MSY = replist[["derived_quants"]]["SSB_MSY","Value"]/SSB0
			lines(x=rep(MSY,2), y=c(0, replist[["derived_quants"]]["Dead_Catch_MSY","Value"]),
				col=colvec[1], lwd=2, lty=ltyvec[1])
		}
		if ("Btgt" %in% refpoints) {
			TRP  = round(replist[["derived_quants"]]["SSB_Btgt","Value"]/SSB0, 5)
			YTRP = replist[["derived_quants"]]["Dead_Catch_Btgt","Value"]
			LRP  = 0.4 * TRP
			USR  = 0.8 * TRP
			YLRP = approx(x=equil_yield[["Depletion"]], y=equil_yield[["Catch"]], xout = LRP)$y
			YUSR = approx(x=equil_yield[["Depletion"]], y=equil_yield[["Catch"]], xout = USR)$y
			lines(x=rep(TRP, 2), y=c(0, YTRP), col=colvec[2], lwd=2, lty=ltyvec[2])
			if ("LRP" %in% refpoints)
				lines(x=rep(LRP, 2), y=c(0, YLRP), col=colvec[5], lwd=2, lty=ltyvec[5])
			if ("USR" %in% refpoints)
				lines(x=rep(USR, 2), y=c(0, YUSR), col=colvec[6], lwd=2, lty=ltyvec[6])
		}
		if ("SPR" %in% refpoints) {
		    lines(x=rep(replist[["derived_quants"]]["SSB_SPR", 
		        "Value"]/SSB0, 2), y=c(0, replist[["derived_quants"]]["Dead_Catch_SPR", 
		        "Value"]), col=colvec[3], lwd=2, lty=ltyvec[3])
		}
		if ("Current" %in% refpoints) {
			Current = replist[["current_depletion"]]
			which_val <- which(abs(equil_yield[["Depletion"]] - 
				replist[["current_depletion"]]) == min(abs(equil_yield[["Depletion"]] - 
				replist[["current_depletion"]])))[1]
			lines(x=rep(Current,2), y=c(0, equil_yield[["Catch"]][which_val]), 
				col=colvec[4], lwd=2, lty=ltyvec[4])
		}
		which_lines <- c("MSY" %in% refpoints, "Btgt" %in% refpoints, "SPR" %in% refpoints,
			"Current" %in% refpoints, "LRP" %in% refpoints, "USR" %in% refpoints)
		if (any(which_lines)) {
			z = match(c("MSY","SPR","Current","Btgt","LRP","USR"),refpoints)
			z = na.omit(z)
			legtxt = sub("Current", paste0("Current = ",round(Current,2),"B0"),
				sub("MSY", paste0("MSY = ",round(MSY,2),"B0"),
				sub("USR", paste0("USR = ",USR,"B0"), 
				sub("LRP", paste0("LRP = ",LRP,"B0"), 
				sub("Btgt", paste0("TRP = ",show0(TRP,2),"B0"), refpoints)))))
			legend("topright", bty="n", lwd=2, lty=ltyvec[which_lines][z], col=colvec[which_lines][z], legend=legtxt[z], seg.len=4, inset=0.01)
		}
	}
	if (1 %in% subplots | 2 %in% subplots) {
		if (!is.null(equil_yield[[1]][1]) && any(!is.na(equil_yield[[1]]))) {
			if (any(!is.na(equil_yield[["Depletion"]])) & any(!is.na(equil_yield[["Catch"]])) & 
				any(!is.infinite(equil_yield[["Depletion"]]))) {
				if (1 %in% subplots) {
					if (plot) {
						yieldfunc()
					}
					if (print) {
						file <- "yield1_yield_curve.png"
						caption <- "Yield curve"
						plotinfo <- save_png(plotinfo=plotinfo, file=file, plotdir=plotdir, pwidth=pwidth, 
							pheight=pheight, punits=punits, res=res, ptsize=ptsize, caption=caption)
						yieldfunc()
						dev.off()
					}
				}
				if (2 %in% subplots & !is.null(refpoints)) {
					if (plot) {
						yieldfunc(refpoints=refpoints)
					}
					if (print) {
						file <- "yield2_yield_curve_with_refpoints.png"
						caption <- "Yield curve with reference points"
						plotinfo <- save_png(plotinfo=plotinfo, file=file, plotdir=plotdir, pwidth=pwidth, 
							pheight=pheight, punits=punits, res=res, ptsize=ptsize, caption=caption)
						yieldfunc(refpoints=refpoints)
						dev.off()
					}
				}
			}
			else {
				message("Skipped equilibrium yield plots: equil_yield has all NA values")
			}
		}
		else {
			message("Skipped equilibrium yield plots: no equil_yield results in this model")
		}
	}
	df <- dplyr::summarise(dplyr::group_by(dplyr::ungroup(dplyr::summarise(dplyr::group_by(dplyr::mutate(dplyr::filter(timeseries, 
		!Era %in% c("VIRG", "FORE")), catch_tot=rowSums(pick(starts_with("dead(B)")), 
		na.rm=TRUE)), Yr, Seas), sum_Bio_all=sum(Bio_all), 
		sum_SpawnBio=sum(SpawnBio), sum_catch_tot=sum(catch_tot))), 
		Yr), mean_Bio_all=mean(sum_Bio_all), mean_SpawnBio=mean(sum_SpawnBio, 
		na.rm=TRUE), catch_tot=sum(sum_catch_tot))
	Nyrs <- nrow(df)
	df[["sprod"]] <- NA
	df[["sprod"]][1:(Nyrs - 1)] <- df[["mean_Bio_all"]][2:Nyrs] - 
		df[["mean_Bio_all"]][1:(Nyrs - 1)] + df[["catch_tot"]][1:(Nyrs - 1)]
	df <- dplyr::filter(df, !is.na(sprod))
	sprodfunc <- function(bio_col, xlab) {
		x <- df[[bio_col]]
		y <- df[["sprod"]]
		xlim <- c(0, max(x, na.rm=TRUE))
		ylim <- c(min(0, y, na.rm=TRUE), max(y, na.rm=TRUE))
		if (!add) {
			plot(0, ylim=ylim, xlim=xlim, xlab=xlab, ylab=labels[4], type="n")
		}
		lines(x, y, col=col2)
		old_warn <- options()[["warn"]]
		options(warn=-1)
		s <- seq(length(y) - 1)
		arrows(x[s], y[s], x[s + 1], y[s + 1], length=0.06, angle=20, col=col2, lwd=1.2)
		options(warn=old_warn)
		abline(h=0, col="grey")
		abline(v=0, col="grey")
		points(x[1], y[1], col=col2, bg="white", pch=21)
	}
	YPR_timeseries <- function() {
		if ("Era" %in% names(sprseries)) {
			sub <- sprseries[["Era"]] != "FORE"
		} else {
			sub <- sprseries[["Yr"]] <= replist[["endyr"]]
		}
		plot(x=sprseries[["Yr"]][sub], y=sprseries[["YPR"]][sub], 
			ylim=c(0, 1.1 * max(sprseries[["YPR"]][sub], na.rm=TRUE)), 
			xlab="Year", ylab=labels[5], type="l", lwd=2, 
			col="blue", yaxs="i")
	}
	if (3 %in% subplots) {
		if (plot) {
			sprodfunc(bio_col="mean_Bio_all", xlab=labels[3])
		}
		if (print) {
			file <- "yield3_surplus_production.png"
			caption <- paste("Surplus production vs. total biomass plot. For interpretation, see<br>", 
				"<blockquote>Walters, Hilborn, and  Christensen, 2008,", 
				"Surplus production dynamics in declining and", 
				"recovering fish populations. <i>Can. J. Fish. Aquat. Sci.</i>", 
				"65: 2536-2551. <a href='https://doi.org/10.1139/F08-170'>https://doi.org/10.1139/F08-170</a>.</blockquote>")
			plotinfo <- save_png(plotinfo=plotinfo, file=file, plotdir=plotdir, pwidth=pwidth, pheight=pheight, 
				punits=punits, res=res, ptsize=ptsize, caption=caption)
			sprodfunc(bio_col="mean_Bio_all", xlab=labels[3])
			dev.off()
		}
	}
	if (4 %in% subplots) {
		if (plot) {
			sprodfunc(bio_col="mean_SpawnBio", xlab=labels[6])
		}
		if (print) {
			file <- "yield4_surplus_production.png"
			caption <- paste("Surplus production vs. spawning biomass plot. For interpretation, see<br>", 
				"<blockquote>Forrest, Kronlund, Cleary, and Grinnell. 2023. An", 
				"evidence-based approach for selecting a limit reference point for Pacific", 
				"herring (<i>Clupea pallasii</i>) stocks in British Columbia, Canada. <i>Can. J. Fish.", 
				"Aquat. Sci.</i> 80: 1071-1083.", "<a href='https://doi.org/10.1139/cjfas-2022-0168'>https://doi.org/10.1139/cjfas-2022-0168</a>.</blockquote>")
			plotinfo <- save_png(plotinfo=plotinfo, file=file, plotdir=plotdir, pwidth=pwidth, pheight=pheight, 
				punits=punits, res=res, ptsize=ptsize, caption=caption)
			sprodfunc(bio_col="mean_SpawnBio", xlab=labels[6])
			dev.off()
		}
	}
	if (5 %in% subplots) {
		sprseries <- replist[["sprseries"]]
		if (is.null(sprseries)) {
			if (verbose) {
				message("Skipping yield per recruit plot because SPR_SERIES not in output")
			}
		} else {
			if (plot) {
				YPR_timeseries()
			}
			if (print) {
				file <- "yield5_YPR_timeseries.png"
				caption <- "Time series of yield per recruit (kg)"
				plotinfo <- save_png(plotinfo=plotinfo, file=file, plotdir=plotdir, pwidth=pwidth, pheight=pheight, 
					punits=punits, res=res, ptsize=ptsize, caption=caption)
				YPR_timeseries()
				dev.off()
			}
		}
	}
	if (!is.null(plotinfo)) 
		plotinfo[["category"]] <- "Yield"
	return(invisible(plotinfo))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plotSS.yield
