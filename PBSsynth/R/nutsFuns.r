##==========================================================
## PBS Stock Synthesis (nutsFuns):
##  Functions created by Chris Grandin to run on a Linux server.
## ---------------------------------------------------------
## run_adnuts............Run adnuts MCMC for the model found in the directory given
## fix.posteriors........Remove multiple header/datarows found in posteriors files to leave just one
## load_ss_files.........Load all the SS files for output and input, and return the model object
## calc.mcmc.............Return a list of mcmc calculations, e.g. quantiles for various values
## create.key.nuisance_posteriors_files......Probably does what function name suggests
## load_models...........Load models from files created using 'create_rds_file()'
## modify_starter_mcmc_type...Change the MCMC setting in an SS starter file so the extra MCMC files are not created
## write_psv2............Write an ADMB PSV file
## write_admb_cov2.......Write a covariance matrix to 'admodel.cov'
## stop_quietly..........Stop the program without issuing an error message
## ---------------------------------------------CG

##==========================================================

## run_adnuts---------------------------2023-09-18
## Run adnuts MCMC for the model found in the directory given
## ---------------------------------------------CG
##' Run adnuts MCMC for the model found in the directory given
##'
##' @details
##' `path` is the directory in which the MLE will be run, a subdirectory of
##' this called `mcmc` is where the MCMC will be run using the `NUTS`
##' algorithm. Inside the `mcmc` directory, several temporary subdirectories
##' will be created, one for each MCMC chain labeled `chain_*`, AkA CPU
##' number used in the parallel execution. These will disappear once the
##' run has completed and the output has been merged.
##'
##' @param path Directory where the model files reside
##' @param num_chains The number of chains to run in parallel. If `NULL`,
##' 1 less than the number of cores on the machine will be used
##' 1 less than the number of cores on the machine will be used
##' @param seed The random seed used to draw the random seeds for each chain
##' @param num_samples The number of samples to output
##' @param num_warmup_samples The warmup samples (equivalent of burnin)
##' @param adapt_delta The target acceptance rate. See [adnuts::sample_admb()]
##' @param run_extra_mcmc If `TRUE`, run SS extra mcmc option which outputs
##' files into the `sso` subdirectory. If `FALSE`, those files will not be
##' created and the `posteriors.sso` and `dervied_posteriors.sso` files
##' will be in the running directory
##' @param fn_exe The name of the executable which was built using ADMB
##' @param overwrite Logical. If `TRUE`, don't ask user if they want to
##' overwrite if the directory already exists, just do it
##' @param input_files The input files for SS
##' @param hess_step Logical. If `TRUE`, use the `hess_step` algorithm`
##' @param fn_logfile The filename of the logfile
##'
##' @return Nothing
##' @export
#
##ss_executable="/bin/ss3/build/ss3"  ## does not work because 'file(fn, "wb")'
###   evaluates to '/home/haigh/POP2023/models/Run12/mcmc//bin/ss3/build/ss3.psv' in 'write-psv.R:24'
#
##source ~/.bashrc  ## explicitly puts ss3 on PATH (see last line)
#                  ## export PATH=/home/haigh:/bin/ss3/build:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin
#system(". ~/.bashrc", intern = TRUE, wait = TRUE)  ## https://phoenixnap.com/kb/linux-source-command

run_adnuts <- function(path, num_chains=NULL, seed=42, num_samples=2000,
  num_warmup_samples=250, adapt_delta=0.9, run_extra_mcmc=FALSE,
  hess_step=TRUE, fn_exe=ifelse(exists("ss_executable"), ss_executable, "ss3"),
  overwrite=FALSE, fn_logfile="model_output.log",
  input_files=c(
  ifelse(exists("starter_file_name"), starter_file_name, "starter.ss"), 
  ifelse(exists("starter_file_name"), forecast_file_name, "forecast.ss"), 
  ifelse(exists("weight_at_age_file_name"), weight_at_age_file_name, "wtatage.ss"), 
  ifelse(exists("control_file_name"), control_file_name, "control.ss"),
  ifelse(exists("data_file_name"), data_file_name, "data.ss")) )
{
  # Determine if the caller is calling from an Rstudio session
  is_rstudio <- Sys.getenv("RSTUDIO") == "1"

  # Get the actual full path of the executable on whatever machine you're on
  fn_exe <- get_model_executable(fn_exe)

  # Chains to run in parallel
  num_machine_cores  <- detectCores()
  if(is.null(num_chains)){
    num_chains <- num_machine_cores - 1
  }
  if(num_machine_cores <= num_chains){
    msg <- paste0("The number of available cores (", num_machine_cores, ") ",
                  "is less than or equal to the number of chains ",
                  "requested (", num_chains, ")")
    if(is_rstudio){
      message(red(symbol$cross, msg))
      stop_quietly(call. = FALSE)
    }else{
      stop(msg, call. = FALSE)
    }
  }

  set.seed(seed)
  seeds      <- sample(1:1e4, size = num_chains)
  mcmc_path  <- file.path(path, "mcmc")
  rdata_file <- file.path(mcmc_path, paste0(basename(path), ".Rdata"))
#browser();return()

  if(dir.exists(mcmc_path)){
    ovrw <- 1
    if(!overwrite && interactive()){
      ovrw <- menu(c("Yes", "No"),
                   title = paste0(mcmc_path, " directory exists. Overwrite?"))
    }
    if(ovrw == 1){
      unlink(mcmc_path, recursive = TRUE, force = TRUE)
    }else{
      msg <- paste0("The `mcmc` directory `", mcmc_path, "` exists and was ",
                    "not modified. Delete it or set `ovrwrite` to `TRUE` if ",
                    "you want to run the procedure.")
      if(is_rstudio){
        message(red(symbol$cross, msg))
        stop_quietly(, call. = FALSE)
      }else{
        stop(msg, call. = FALSE)
      }
    }
  }
  dir.create(mcmc_path, showWarnings = FALSE)

  if(hess_step){
    input_files <- c(input_files, "admodel.cov", "admodel.hes", "ss.bar")
  }

  message("\n")
  if(!is.null(fn_logfile)){
    msg <- paste0("The logfile name `", fn_logfile, "` was supplied so no ",
                  "ADMB model output will appear in the console while the ",
                  "adnuts procedure is running\n")
    if(is_rstudio){
      message(green(symbol$circle_double, msg))
    }else{
      message(msg)
    }
  }

  # Run MLE and optimization MCMC (path) ----
  msg <- paste0("Running optimization MCMC (chain length 15) ",
                "to ensure hessian is good ",
                "and optimize without bias adjustment turned on\n")
  if(is_rstudio){
    message(green(symbol$info, msg))
  }else{
    message(msg)
  }
  # Run MLE to create .ss_new files
  cmd <- paste0("cd ", path, " && ", fn_exe)
  if(!is.null(fn_logfile)){
    cmd <- paste0(cmd, " > ", fn_logfile, " 2>&1")
  }
  ss_new_run <- system_(cmd, intern = TRUE, wait = TRUE)
##  ss_new_run = character()
  if(!is.null(attributes(ss_new_run)) && attr(ss_new_run, "status")){
    msg <- paste0("System command returned an error (status 1):\n", cmd)
    if(is_rstudio){
      message(red(symbol$cross, msg))
      stop_quietly(call. = FALSE)
    }else{
      stop(msg, call. = FALSE)
    }
  }
  cmd <- paste0("cd ", path, " && ", fn_exe,  " -nox -iprint 200 -mcmc 15")
  if(!is.null(fn_logfile)){
    cmd <- paste0(cmd, " > ", fn_logfile, " 2>&1")
  }
  initial_run <- system_(cmd, intern = TRUE, wait = TRUE)
  if(!is.null(attributes(initial_run)) && attr(initial_run, "status")){
    msg <- paste0("System command returned an error (status 1):\n", cmd)
    if(is_rstudio){
      message(red(symbol$cross, msg))
      stop_quietly(call. = FALSE)
    }else{
      stop(msg, call. = FALSE)
    }
  }

  # Run initial MCMC (mcmc_path) ----
  input_files <- file.path(path, input_files)
  input_files <- input_files[file.exists(input_files)]

  file.copy(input_files, mcmc_path, overwrite = TRUE)
  if(run_extra_mcmc){
    dir.create(file.path(mcmc_path, "sso"), showWarnings = FALSE)
    modify_starter_mcmc_type(mcmc_path, 2)
  }else{
    modify_starter_mcmc_type(mcmc_path, 1)
  }
  # rdata_file <- file.path(mcmc_path, "pop.Rdata")
  # Run MLE to create .ss_new files
  cmd <- paste0("cd ", mcmc_path, " && ", fn_exe)
  if(!is.null(fn_logfile)){
    cmd <- paste0(cmd, " > ", fn_logfile, " 2>&1")
  }
  ss_new_run <- system_(cmd, intern = TRUE, wait = TRUE)
  if(!is.null(attributes(ss_new_run)) && attr(ss_new_run, "status")){
    msg <- paste0("System command returned an error (status 1):\n", cmd)
    if(is_rstudio){
      message(red(symbol$cross, msg))
      stop_quietly(call. = FALSE)
    }else{
      stop(msg, call. = FALSE)
    }
  }
  # Copy .ss_new files to .ss files in the mcmc directory
  r4ss::copy_SS_inputs(dir.old = mcmc_path,
                       dir.new = mcmc_path,
                       use_ss_new = TRUE,
                       overwrite = TRUE,
                       verbose = FALSE)

  # The -hbf 1 argument is a technical requirement because NUTS uses a
  # different set of bounding functions and thus the mass matrix will be
  # different
  if(hess_step){
    msg <- paste0("Running re-optimization MCMC (chain length 15) ",
                  "to get the correct mass ",
                  "matrix for the NUTS run (use `-hbf 1`) and the ",
                  "`hess_step` algorithm\n")
    if(is_rstudio){
      message(green(symbol$info, msg))
    }else{
      message(msg)
    }
    cmd <- paste0("cd ", mcmc_path, " && ", fn_exe,
                  " -hbf 1 -nox -iprint 200 -mcmc 15 -hess_step 10 ",
                  "-binp ss.bar")
    if(!is.null(fn_logfile)){
      cmd <- paste0(cmd, " > ", fn_logfile, " 2>&1")
    }
    system_(cmd, intern = TRUE, wait = TRUE)
  }else{
    msg <- paste0("Running re-optimization MCMC (chain length 15) ",
                  "to get the correct mass matrix ",
                  "for the NUTS run (use `-hbf 1`)\n")
    if(is_rstudio){
      message(green(symbol$info, msg))
    }else{
      message(msg)
    }
    cmd <- paste0("cd ", mcmc_path, " && ",
                  fn_exe, " -hbf 1 -nox -iprint 200 -mcmc 15")
    if(!is.null(fn_logfile)){
      cmd <- paste0(cmd, " > ", fn_logfile, " 2>&1")
    }
    system_(cmd, intern = TRUE, wait = TRUE)
  }

  ## Hake uses this formula
  num_iters <- ceiling(((num_chains * num_warmup_samples) + num_samples) / num_chains)

  # Run ADNUTS MCMC
  msg <- paste0("Running Initial NUTS MCMC (chain length ", num_iters, ", warmup ", num_warmup_samples, ") ",
                "with MLE mass matrix\n")
  if(is_rstudio){
    message(green(symbol$info, msg))
  }else{
    message(msg)
  }


  nuts_initial <- sample_admb(model = fn_exe,
                              path = mcmc_path,
                              algorithm = "nuts",
                              num_samples = num_iters,
                              seeds = seeds,
                              num_chains = num_chains,
                              warmup = num_warmup_samples,
                              control = list(metric = "mle",
                                             adapt_delta = adapt_delta),
                              fn_logfile = fn_logfile)
#browser();return()

  # Run again for inference using updated mass matrix.
  # Increase adapt_delta toward 1 if you have divergences (runs will take
  # longer). Note this is in unbounded parameter space
  mass <- nuts_initial$covar_est
  inits <- sample_inits(nuts_initial, num_chains)

  msg <- paste0("Running updated NUTS MCMC for inference, acceptance ",
                "ratio (adapt_delta) = ", adapt_delta, "\n")
  if(is_rstudio){
    message(green(symbol$info, msg))
  }else{
    message(msg)
  }
  nuts_updated <- sample_admb(model = fn_exe,
                              path = mcmc_path,
                              algorithm = "nuts",
                              num_samples = num_iters,
                              init = inits,
                              seeds = seeds,
                              num_chains = num_chains,
                              warmup = num_warmup_samples,
                              mceval = TRUE,
                              control = list(metric = mass,
                                             adapt_delta = adapt_delta),
                              fn_logfile = fn_logfile)

  save(list = ls(all.names = TRUE), file = rdata_file, envir = environment())

  cmd <- paste0("cd ", mcmc_path, " && ", fn_exe, " -mceval")
  if(!is.null(fn_logfile)){
    cmd <- paste0(cmd, " > ", fn_logfile, " 2>&1")
  }
  #system_(cmd, intern = TRUE, wait = TRUE)  ## already run in sample_admb

  msg <- "Finished `run_adnuts()`\n"
  if(is_rstudio){
    message(green(symbol$info, msg))
  }else{
    message(msg)
  }

  invisible(nuts_updated)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~run_adnuts

## Code to activate if runing R from R gui console?
#run_from_R  = FALSE
#if (run_from_R) {
#  root_path   = "/home/haigh/POP2023"
#  models_path = paste0(root_path, "/models")
#  r_path      = paste0(root_path, "/R")
#  num_chains  = 8
#  num_samples = 2000  ## final number of samples
#  prop_warmup = 0.5 
#  ## Offhore rockfish assumes burnin=0 in starter.ss and takes % burn-in during sampling
#  ## a=num_samples, b=burnin explicit, c=num_chains, w=prop_warmup (implicit burnin), W=num_warmup_samples
#  burnin      = 0 ## preceeding the final MCMC samples (specified in starter.ss)
#  icalc       = function(a,b=0,c,w){ ceiling( (a+b)/(c*(1-w)) ) }         ## w = proportion burn-in
#  num_iters   = icalc(a=num_samples, b=burnin, c=num_chains, w=prop_warmup)
#  num_warmup_samples = num_iters * (1 - ((num_samples + burnin) / (num_iters*num_chains)))
#
#  run_extra_mcmc = FALSE
#  adapt_delta    = 0.9
#  verbose        = FALSE
#  models         = "Run04" #"CAR24" #  ## use Canary for debugging CG's code (shorter convergence times)
#  input_files    = c("starter.ss","forecast.ss","control.ss","data.ss")  ## seems incapable of having data and control files with different names
#
#  source(paste0(r_path, "/load-functions.R"))  ## includes call to 'fns2fix.R'
#
#  for (i in 1:length(models)){
#    ii = models[i]
#    model_path = paste0(models_path, "/", ii)
#    run_adnuts(path      = model_path,
#      run_extra_mcmc     = run_extra_mcmc,
#      num_chains         = num_chains,
#      adapt_delta        = adapt_delta,
#      num_samples        = num_samples,
#      num_warmup_samples = num_warmup_samples,
#      overwrite          = TRUE)
#  }
#}

## fix.posteriors-----------------------2023-09-18
##  Remove multiple header/datarows found in 
##  posteriors files to leave just one
## ---------------------------------------------CG
##' Remove multiple header/datarows found in posteriors files to leave just one
##'
##' @param dir Directory where the files reside
##'
##' @return If the file contains only one header/datarow (i.e is of correct format)
##' then return a 1-row dataframe of the contents
##' @export
fix.posteriors <- function(dir)
{
  do.it <- function(file){
    posts <- read.table(file.path(dir, file),
                        header = TRUE,
                        fill = TRUE,
                        stringsAsFactors = FALSE)
    if(all(grepl("^[[:digit:]]", posts[,1]))){
      return(posts)
    }
    write.table(posts[1:(grep("\\D+", posts[,1])[1] - 1),],
                file.path(dir, file),
                quote = FALSE,
                row.names = FALSE)
  }
  do.it(posts_file_name)
  do.it(derposts_file_name)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~fix.posteriors


## load_ss_files------------------------2023-09-18
##  Load all the SS files for output and input, 
##  and return the model object
## ---------------------------------------------CG
##' Load all the SS files for output and input, and return the model object
##'
##' @details If MCMC directory is present, load that and perform calculations for mcmc parameters.
##'
##' @param model.dir Directory the model reesides in
##' @param key_posts Vector of key posteriors used to create key posteriors file
##' @param key_posts_fn Name of the key posteriors file
##' @param nuisance_posts_fn Name of the nuisance posteriors file
##' @param printstats Print info on each model loaded via [r4ss::SS_output()]
##'
##' @return A model object representing the output from the SS model
##' @export
load_ss_files <- function(model_path = NA,
                          key_posts = c("NatM",
                                        "SR_LN",
                                        "SR_BH_steep",
                                        "Q_extraSD",
                                        "ln.EffN_mult._1",
                                        "ln.EffN_mult._2"),
                          key_posts_fn = "keyposteriors.csv",
                          nuisance_posts_fn = "nuisanceposteriors.csv",
                          printstats = FALSE,
                          ...)
{
  stopifnot(!is.na(model_path))

  # Load MPD results
  model <- tryCatch({
    SS_output(dir = model_path,
              verbose = FALSE,
              printstats = printstats,
              covar = FALSE,
              wtfile = "wtatage.ss")
  }, error = function(e){
    SS_output(dir = model_path,
              verbose = FALSE,
              printstats = printstats,
              covar = FALSE,
              forecast = FALSE,
              wtfile = "wtatage.ss")
  })

  # Load the data file and control file for the model
  # Get the file whose name contains "_data.ss" and "_control.ss"
  # If there is not exactly one of each, stop with error.
  model_path_listing <- tolower(dir(model_path))
  dat_fn_ind <- grep("_data.ss", model_path_listing)
  ctl_fn_ind <- grep("_control.ss", model_path_listing)
  par_fn_ind <- grep("ss.par", model_path_listing)
  if(!length(dat_fn_ind)){
    stop("Error in model ", model_path,
         ", there is no data file. A data file is any file whose name end with _data.ss.\n",
         call. = FALSE)
  }
  if(length(dat_fn_ind) > 1){
    stop("Error in model ", model_path,
         ", there is more than one data file. A data file is any file whose name ends with _data.ss.\n\n",
         call. = FALSE)

  }
  if(!length(ctl_fn_ind)){
    stop("Error in model ", model_path,
         ", there is no control file. A control file is any file whose name ends with _control.ss.\n\n",
         call. = FALSE)

  }
  if(length(ctl_fn_ind) > 1){
    stop("Error in model ", model_path,
         ", there is more than one control file. A control file is any file whose name ends with _control.ss.\n\n",
         call. = FALSE)
  }
  model$path <- model_path
  model$dat_file <- file.path(model_path, model_path_listing[dat_fn_ind])
  model$ctl_file <- file.path(model_path, model_path_listing[ctl_fn_ind])
  model$par_file <- file.path(model_path, model_path_listing[par_fn_ind])
  model$dat <- SS_readdat(model$dat_file, verbose = FALSE)
  model$ctl <- readLines(model$ctl_file)
  model$ctl <- gsub("\t", " ", model$ctl)

  # model$par <- readLines(par_fn)
  # Set default mcmc members to NA. Later code depends on this.
  model$mcmc <- NA
  # Set the mcmc and extra mcmc paths and record their existence
  model$mcmc_path <- file.path(model_path, "mcmc")
  model$mcmc_exists <- dir.exists(model$mcmc_path)
  model$extra_mcmc_path <- file.path(model$mcmc_path, "sso")
  model$extra_mcmc_exists <- dir.exists(model$extra_mcmc_path)
  # Save the posterior names from the mcmc output. This is necessary for the function `plot_mcmc_param_stats()`
  posteriors_dir <- ifelse(model$extra_mcmc_exists, model$extra_mcmc_path, model$mcmc_path)

  # If it has an mcmc sub-directory, load that as well
  if(dir.exists(posteriors_dir)){
    tmp <- readLines(file.path(posteriors_dir, "posteriors.sso"), n = 1)
    tmp <- stringr::str_split(tmp, "[:space:]")[[1]]
    # Remove Empty string, Iter and Objective_function as they are not parameters
    model$post_names <- tmp[!tmp %in% c("", "Iter", "Objective_function")]
    fix.posteriors(posteriors_dir)
    model$mcmc <- SSgetMCMC(dir = posteriors_dir,
                            writecsv = FALSE,
                            verbose = FALSE)
    # replace any SPB with SSB
    names(model$mcmc) <- gsub(pattern="SPB", replacement="SSB", names(model$mcmc))
    create.key.nuisance_posteriors_files(model,
                                         key_posts,
                                         key_posts_fn,
                                         nuisance_posts_fn)
    # Do the mcmc calculations, e.g. quantiles for SB, SSB, DEPL, RECR, RECRDEVS
    model$mcmccalcs <- calc.mcmc(model$mcmc)

  }
  model
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~load_ss_files


## calc.mcmc----------------------------2023-09-18
##  Return a list of mcmc calculations, 
##  e.g. quantiles for various values
## ---------------------------------------------CG
##' Return a list of mcmc calculations, e.g. quantiles for various values
##'
##' @param mcmc The output of the [r4ss::SS_getMCMC()] function as a data.frame
##' @param lower Lower quantile value
##' @param upper Upper quantile value
##' @param biomass.scale Scale the biomass by this amount. The default for
##'  2021 and before was 2e6, but is not 1e6 due to changes in SS3 (hake Issue
##'   #866). Biomass will be shown in the millions of tonnes (and is females only).
##' @param recruitment.scale Scale the recruitment by this amount. The default is 1e6
##' because recruitment will be shown in millions of tonnes
##'
##' @return
##' @export
calc.mcmc <- function(mcmc,
                      lower = 0.025,
                      upper = 0.975,
                      biomass.scale = 1e6,
                      recruitment.scale = 1e6)
{
  ssb <- mcmc[,grep("SSB",names(mcmc))] / biomass.scale
  svirg <- quantile(ssb[,names(ssb) == "SSB_Virgin"],
                    c(lower, 0.5, upper))
  sinit <- quantile(ssb[,names(ssb) == "SSB_Initial"],
                    c(lower, 0.5, upper))

  # sinit.post exists so that depletion calculations can be done for each posterior
  sinit.post <- ssb[,names(ssb) == "SSB_Initial"]

  names(ssb) <- gsub("SSB_", "", names(ssb))
  cols.to.strip <- c("Virgin", "Initial")
  ssb <- strip.columns(ssb, cols.to.strip)

  slower <- apply(ssb, 2, quantile, prob = lower, na.rm = TRUE)
  smed   <- apply(ssb, 2, quantile, prob = 0.5, na.rm = TRUE)
  supper <- apply(ssb, 2, quantile, prob = upper, na.rm = TRUE)

  depl   <- apply(ssb, 2, function(x){x / sinit.post})
  dlower <- apply(depl, 2, quantile, prob = lower, na.rm = TRUE)
  dmed   <- apply(depl, 2, quantile, prob = 0.5, na.rm = TRUE)
  dupper <- apply(depl, 2, quantile, prob = upper, na.rm = TRUE)

  ## 1e6 used here because recruitment will be shown in the millions of tonnes
  recr <- mcmc[,grep("Recr_", names(mcmc))] / recruitment.scale
  recr <- recr[,-grep("Fore", names(recr))]
  names(recr) <- gsub("Recr_", "", names(recr))
  rvirg <- quantile(recr[,names(recr) == "Virgin"],
                    c(lower, 0.5, upper))
  rinit <- quantile(recr[,names(recr) == "Initial"],
                    c(lower, 0.5, upper))
  runfished <- quantile(recr[,grepl("unfished", names(recr), ignore.case = TRUE)],
                        c(lower, 0.5, upper))

  cols.to.strip <- c("Virgin", "Initial",
    grep("unfished", names(recr), ignore.case = TRUE, value = TRUE))
  recr <- strip.columns(recr, cols.to.strip)

  rmed <- apply(recr, 2, quantile, prob = 0.5, na.rm = TRUE)
  rmean <- apply(recr, 2, mean, na.rm = TRUE)
  rlower <- apply(recr, 2, quantile,prob = lower, na.rm = TRUE)
  rupper <- apply(recr, 2, quantile,prob = upper, na.rm = TRUE)

  dev <- mcmc[,c(grep("Early_InitAge_", names(mcmc)),
                 grep("Early_RecrDev_", names(mcmc)),
                 grep("Main_RecrDev_", names(mcmc)),
                 grep("Late_RecrDev_", names(mcmc)),
                 grep("ForeRecr_", names(mcmc)))]

  names(dev) <- gsub("Early_RecrDev_", "", names(dev))
  names(dev) <- gsub("Main_RecrDev_", "", names(dev))
  names(dev) <- gsub("Late_RecrDev_", "", names(dev))
  names(dev) <- gsub("ForeRecr_", "", names(dev))

  # Change the Early_Init names to be the correct preceeding years
  start_yr <- as.numeric(min(names(dev)))
  early <- grep("Early_InitAge_", names(dev))
  num.early.yrs <- length(early)
  early.yrs <- seq(start_yr - num.early.yrs, start_yr - 1, 1)
  late.yrs <- names(dev[-early])
  names(dev) <- c(as.character(early.yrs), late.yrs)

  devlower <- apply(dev, 2, quantile, prob = lower, na.rm = TRUE)
  devmed <- apply(dev, 2, quantile, prob = 0.5, na.rm = TRUE)
  devupper <- apply(dev, 2, quantile, prob = upper, na.rm = TRUE)

  spr <- mcmc[,grep("SPRratio_", names(mcmc))]
  names(spr) <- gsub("SPRratio_", "", names(spr))

  plower <- apply(spr, 2, quantile, prob = lower, na.rm = TRUE)
  pmed <- apply(spr, 2, quantile, prob = 0.5, na.rm = TRUE)
  pupper <- apply(spr, 2, quantile, prob = upper, na.rm = TRUE)

  f <- mcmc[,grep("F_", names(mcmc))]
  names(f) <- gsub("F_", "", names(f))
  flower <- apply(f, 2, quantile, prob = lower, na.rm = TRUE)
  fmed   <- apply(f, 2, quantile, prob = 0.5, na.rm = TRUE)
  fupper <- apply(f, 2, quantile, prob = upper, na.rm = TRUE)

  # Calculations for the reference points table
  probs <- c(lower, 0.5, upper)

  unfish.fem.bio <-
    f(round(quantile(mcmc$SSB_Virgin,
                     prob = probs) / biomass.scale, 3) * 1000,
      0)
  unfish.recr <-
    f(round(quantile(mcmc$Recr_Virgin,
                     prob = probs) / recruitment.scale, 3) * 1000,
      0)
  f.spawn.bio.bf40 <-
    f(round(quantile(mcmc$SSB_SPR,
                     prob = probs) / biomass.scale, 3) * 1000,
      0)
  spr.msy.proxy <- c(latex.bold("--"),
                     "40\\%",
                     latex.bold("--"))
  exp.frac.spr <-
    paste0(f(100 * quantile(mcmc$annF_SPR,
                            prob = probs),
             1),
           "\\%")
  yield.bf40 <-
    f(round(quantile(mcmc$Dead_Catch_SPR,
                     prob = probs) / recruitment.scale, 3) * 1000,
      0)
  fem.spawn.bio.b40 <-
    f(round(quantile(mcmc$SSB_Btgt,
                     prob = probs) / biomass.scale, 3) * 1000,
      0)
  spr.b40 <-
    paste0(f(100 * quantile(mcmc$SPR_Btgt,
                            prob = probs),
             1),
           "\\%")
  exp.frac.b40 <-
    paste0(f(100 * quantile(mcmc$annF_Btgt,
                            prob = probs),
             1),
           "\\%")
  yield.b40 <-
    f(round(quantile(mcmc$Dead_Catch_Btgt,
                     prob = probs) / recruitment.scale, 3) * 1000,
      0)
  fem.spawn.bio.bmsy <-
    f(round(quantile(mcmc$SSB_MSY,
                     prob = probs) / biomass.scale, 3) * 1000,
      0)
  spr.msy <- paste0(f(100 * quantile(mcmc$SPR_MSY,
                                     prob = probs),
                      1),
                    "\\%")
  exp.frac.sprmsy <-
    paste0(f(100 * quantile(mcmc$annF_MSY,
                            prob = probs),
             1),
           "\\%")
  msy <-
    f(round(quantile(mcmc$Dead_Catch_MSY,
                     prob = probs) / recruitment.scale, 3) * 1000,
      0)

  # Return a list of the calculated values
  sapply(c("svirg",
           "sinit",
           "slower",
           "smed",
           "supper",
           "dlower",
           "dmed",
           "dupper",
           "rvirg",
           "rinit",
           "runfished",
           "rlower",
           "rmed",
           "rupper",
           "rmean",
           "devlower",
           "devmed",
           "devupper",
           "plower",
           "pmed",
           "pupper",
           "flower",
           "fmed",
           "fupper",
           ## Reference points
           "unfish.fem.bio",
           "unfish.recr",
           "f.spawn.bio.bf40",
           "spr.msy.proxy",
           "exp.frac.spr",
           "yield.bf40",
           "fem.spawn.bio.b40",
           "spr.b40",
           "exp.frac.b40",
           "yield.b40",
           "fem.spawn.bio.bmsy",
           "spr.msy",
           "exp.frac.sprmsy",
           "msy"),
         function(x){get(x)})
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calc.mcmc


## -------------------------------------2023-09-18
##  Probably does what function name suggests.
## ---------------------------------------------CG
create.key.nuisance_posteriors_files <- function(model,
                                                 posterior_regex,
                                                 key_post_file,
                                                 nuisance_post_file)
{
  ## Creates the two files for key and nuisance posteriors
  if(model$extra_mcmc_exists){
    key_file <- here::here(model$extra_mcmc_path, key_post_file)
    nuisance_file <- here::here(model$extra_mcmc_path, nuisance_post_file)
  }else{
    key_file <- here::here(model$mcmc_path, key_post_file)
    nuisance_file <- here::here(model$mcmc_path, nuisance_post_file)
  }

  mc <- model$mcmc
  mc_names <- names(mc)
  mcmc_grep <- unique(grep(paste(posterior_regex, collapse="|"), mc_names))
  mcmc_names <- mc_names[mcmc_grep]
  keys <- mc[, mcmc_grep]
  nuisances <- mc[, -mcmc_grep]
  write.csv(keys, key_file, row.names = FALSE)
  write.csv(nuisances, nuisance_file, row.names = FALSE)
}
##~~~~~~~~~~~~create.key.nuisance_posteriors_files


## load_models--------------------------2023-09-18
##  Load models from files created using [create_rds_file()]
## ---------------------------------------------CG
##' Load models from files created using [create_rds_file()]
##'
##' @details Load model(s) and return as a list if more than one. If only one,
##'  return that object or if ret.single.list is TRUE, return a 1-element list.
##'
##' @param model_dirs A vector of model directory names
##' @param ret_single_list See details
##'
##' @return A list of model objects
##' @export
##'
##' @examples
##' base <- load_models("base")
load_models <- function(model_dirs, ret_single_list=FALSE)
{
  ret_list <- NULL
  model_rds_files <- file.path(rootd_models, model_dirs, paste0(model_dirs, ".rds"))
  if(!all(file.exists(model_rds_files))){
    stop("The following files do not exist, run build_rds() on the associated directories:\n",
         paste(model_rds_files[!file.exists(model_rds_files)], collapse = "\n"),
         call. = FALSE)
  }
  for(i in 1:length(model_rds_files)){
    small_file <- file.path(rootd_models, model_dirs[i], paste0("small_", model_dirs[i], ".rds"))
    if(file.exists(small_file)){
      #message("Trying ", small_file)
      #gzfile(small_file, "rb")
      ret_list[[i]] <- readRDS(small_file)
      message("Loaded small RDS file: ", small_file)
    }else{
      ret_list[[i]] <- readRDS(model_rds_files[i])
      message("Loaded large RDS file: ", model_rds_files[i])
    }
  }
  if(length(model_dirs) == 1){
    if(ret_single_list){
      ret_list
    }else{
      ret_list[[1]]
    }
  }else{
    ret_list
  }
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~load_models


## modify_starter_mcmc_type-------------2023-09-18
##  Change the MCMC setting in an SS starter file
##  so the extra MCMC files are not created.
## ---------------------------------------------CG
#' Change the MCMC setting in an SS starter file so the extra MCMC files are
#' not created
#'
#' @param path The path where `starter.ss` resides
#' @param value Value to change the MCMC setting to. 1 = Do not create extra
#' MCMC files, 2 or 3 - create extra MCMC files
#' @param starter_fn The starter filename
#'
#' @return Nothing
#' @export
modify_starter_mcmc_type <- function(path,
                                     value,
                                     starter_fn = "starter.ss"){

  # Make a modification to the starter file
  if(!dir.exists(path)){
    stop("The directory ", path, " does not exist",
         call. = FALSE)
  }
  if(!file.exists(file.path(path, starter_fn))){
    stop("The file ", file.path(path, starter_fn), " does not exist",
         call. = FALSE)
  }
  starter_contents <- readLines(file.path(path,
                                          starter_fn))
  mcmc_output_ind <- grep("MCMC output detail|MCMC_output_detail",
                          starter_contents)
  mcmc_output_val <- starter_contents[mcmc_output_ind]
  mcmc_output_val <- gsub("^.*(#.*)",
                          "\\1",
                          mcmc_output_val)
  mcmc_output_val <- paste0(value,
                            " ",
                            mcmc_output_val,
                            " - *Modified by modify_starter_mcmc_type()*")
  starter_contents[mcmc_output_ind] <- mcmc_output_val
  writeLines(starter_contents, file.path(path, starter_fn))
}
##~~~~~~~~~~~~~~~~~~~~~~~~modify_starter_mcmc_type


## write_psv2---------------------------2023-09-18
##  Write an ADMB PSV file
## ---------------------------------------------CG
##' Write an ADMB PSV file
##'
##' @details Useful to combine multiple MCMC runs together into a single
##' `.psv file` which can then be executed with '-mceval'.
##'
##' @param path Directory to write the PSV file in
##' @param fn_psv Model PSV file name with or without the PSV extension
##' @param samples A matrix or data.frame of samples, each column is a
##' parameter, each row a sample.
write_psv2 <- function(path, fn_psv, samples)
{
  if(!dir.exists(path)){
    stop("Directory `path` = `", path, "` does not exist",
         call. = FALSE)
  }
  fn_psv = basename(fn_psv)  ## get rid of pathways (RH 230314)
  if(!length(grep("\\.(psv)|(PSV)$", fn_psv))){
    fn_psv <- paste0(fn_psv, ".psv")
  }
  fn <- file.path(path, fn_psv)

  samples <- as.matrix(samples)
  con <- file(fn, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(object = ncol(samples), con)
  writeBin(object = as.vector(t(samples)), con)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~write_psv2


## write_admb_cov2----------------------2018-09-18
##  Write a covariance matrix to 'admodel.cov'
## ---------------------------------------------CG
##' Write a covariance matrix to `admodel.cov`.
##'
##' @param path Path to model.
##' @param cov_unbounded The covariance matrix in unbounded space.
##' @param hbf The hybrid_bounded_flag value. Use `hbf = 1` for HMC.
##' @param fn_cov Name of the ADMB covariance file
write_admb_cov2 <- function(path = NULL,
                           cov_unbounded,
                           hbf = NULL,
                           fn_cov = "admodel.cov")
{
  if(!dir.exists(path)){
    stop("`path` does not exist",
         call. = FALSE)
  }
  fn_cov = basename(fn_cov)  ## get rid of pathways (RH 230314)
  fn <- file.path(path, fn_cov)
  if(!file.exists(fn)){
    stop("The file `", fn, "` does not exist",
         call. = FALSE)
  }

  dest_fn <- file.path(path, "admodel_original.cov")
  tmp <- file.copy(from = fn,
                   to = dest_fn)

  # Read in the output files
  results <- read_admb_cov(path)
  if(is.null(hbf)){
    hbf = results$hybrid_bounded_flag
  }
  scale <- results$scale
  num_pars <- results$num_pars
  if(nrow(cov_unbounded) != num_pars)
    stop("Invalid size of covariance matrix, has ", nrow(cov_unbounded),
         " rows but should have ", num_pars,
         call. = FALSE)

  # Write it to file using original scales, although these are ignored.
  file_new <- file(fn, "wb")
  on.exit(close(file_new))
  writeBin(as.integer(num_pars), file_new)
  writeBin(as.vector(as.numeric(cov_unbounded)), file_new)
  writeBin(as.integer(hbf), file_new)
  writeBin(as.vector(scale), file_new)
  invisible()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~write_admb_cov2


## stop_quietly-------------------------2023-09-18
##  Stop the program without issuing an error message
## ---------------------------------------------CG
##' Stop the program without issuing an error message
##'
##' @param ... Arguments passed in to be written out as a message
##'
##' @return Nothing
##' @export
stop_quietly <- function(...) {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  message(unlist(list(...)))
  stop()
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~stop_quietly

