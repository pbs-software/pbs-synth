\name{runADnuts}
\alias{run_adnuts}
\alias{fix.posteriors}
\alias{load_ss_files}
\alias{calc.mcmc}
\alias{create.key.nuisance_posteriors_files}
\alias{load_models}
\alias{modify_starter_mcmc_type}
\alias{write_psv2}
\alias{write_admb_cov2}
\alias{stop_quietly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Run ADMB NUTS
}
\description{
  Run Automatic Differentiation Model Builder (ADMB) using no-U-turn (NUTS) sampling.
}
\usage{
run_adnuts(path, num_chains=NULL, seed=19, num_samples=2000, 
  num_warmup_samples=250, adapt_delta=0.9, run_extra_mcmc=FALSE, 
  hess_step=TRUE, fn_exe=ifelse(exists("ss_executable"), ss_executable, "ss3"), 
  overwrite=FALSE, fn_logfile="model_output.log", 
  input_files=c(ifelse(exists("starter_file_name"), starter_file_name, 
  "starter.ss"), ifelse(exists("starter_file_name"), forecast_file_name, 
  "forecast.ss"), ifelse(exists("weight_at_age_file_name"), 
  weight_at_age_file_name, "wtatage.ss"), ifelse(exists("control_file_name"), 
  control_file_name, "control.ss"), ifelse(exists("data_file_name"), 
  data_file_name, "data.ss")))

fix.posteriors(dir)

load_ss_files(model_path=NA,
  key_posts=c("NatM", "SR_LN", "SR_BH_steep", "Q_extraSD", 
  "ln.EffN_mult._1", "ln.EffN_mult._2"), key_posts_fn="keyposteriors.csv", 
  nuisance_posts_fn="nuisanceposteriors.csv", printstats=FALSE, ...)

calc.mcmc(mcmc, lower=0.025, upper=0.975, biomass.scale=1e6,
  recruitment.scale=1e6)

create.key.nuisance_posteriors_files(
  model, posterior_regex,  key_post_file, nuisance_post_file)

load_models(model_dirs, ret_single_list=FALSE)

modify_starter_mcmc_type(path, value, starter_fn="starter.ss")

write_psv2(path, fn_psv, samples)

write_admb_cov2(path=NULL, cov_unbounded, hbf=NULL, fn_cov="admodel.cov")

stop_quietly(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{\code{character} -- directory where model input files reside}
  \item{num_chains}{\code{numeric} -- number of chains to use in parallel}
  \item{seed}{\code{numeric} -- random seed start value}
  \item{num_samples}{\code{numeric} -- number of MCMC samples wanted}
  \item{num_warmup_samples}{\code{numeric} -- number of warm-up samples per chain}
  \item{adapt_delta}{\code{numeric} -- quality of convergence}
  \item{run_extra_mcmc}{\code{logical} -- if \code{TRUE}, return Report files
    for each MCMC sample retained}
  \item{hess_step}{\code{logical} -- if \code{TRUE}, use the \code{'hess_step'} algorithm
    (aside: Hessian = square matrix of second-order partial derivatives of a scalar-valued function)}
  \item{fn_exe}{\code{character} -- name (with optional location) of the compiled SS3 executable}
  \item{overwrite}{\code{logical} -- if \code{TRUE}, overwrite previous model fits}
  \item{fn_logfile}{\code{character} -- name (with optional location) of the model's log file}
  \item{input_files}{\code{character} -- vector of input file names; the basic four comprise
    \code{'starter.ss'}, \code{'data.ss'}, \code{'control.ss'}, and \code{'forecast.ss'}}
  \cr\bold{Args for extra alias funs:}\cr--------------------\cr 
  \bold{fix.posteriors}
  \item{dir}{\code{character} -- directory where the files reside}
  \bold{load_ss_files}
  \item{model_path}{\code{character} -- directory where model is being evaluated}
  \item{key_posts}{\code{character} -- names of the key parameters}
  \item{key_posts_fn}{\code{character} -- file name for posteriors of key parameters}
  \item{nuisance_posts_fn}{\code{character} -- file name for posteriors of nuisance parameters}
  \item{printstats}{\code{logical} -- if \code{TRUE}, print the statistics}
  \item{...}{\code{dots} -- additional arguments to function \code{'load_ss_files'}}
  \bold{calc.mcmc}
  \item{mcmc}{\code{list} -- output of the \code{'r4ss::SS_getMCMC()'} function as a data.frame}
  \item{lower}{\code{numeric} -- lower quantile value}
  \item{upper}{\code{numeric} -- upper quantile value}
  \item{biomass.scale}{\code{numeric} -- scale the biomass by this amount;
    default for 2021 and before was 2e6, but is now 1e6 due to changes in SS3
    (hake Issue #866); biomass will be shown in the millions of tonnes (and is females only)}
  \item{recruitment.scale}{\code{numeric} -- scale the recruitment by this amount;
    default is 1e6 because recruitment will be shown in millions of tonnes}
  \bold{create.key.nuisance_posteriors_files}
  \item{model}{\code{list} -- probably data frame from \code{'r4ss::SS_getMCMC()'}}
  \item{posterior_regex}{\code{character} -- a Jamaican posterior?}
  \item{key_post_file}{\code{character} -- file name for posteriors of key parameters}
  \item{nuisance_post_file}{\code{character} -- file name for posteriors of nuisance parameters}
  \bold{load_models}
  \item{model_dirs}{\code{character} -- vector of model directory names}
  \item{ret_single_list}{\code{logical} -- if \code{TRUE}, return a single-element list}
  \bold{modify_starter_mcmc_type}
  \item{path}{\code{character} -- path where \code{'starter.ss'} resides}
  \item{value}{\code{numeric} -- value to change the MCMC setting to;
    \code{1} = do not create extra MCMC files, \code{2} or \code{3} = create extra MCMC files}
  \item{starter_fn}{\code{character} -- starter file name}
  \bold{write_psv2}
  \item{path}{\code{character} -- directory to write the PSV file in}
  \item{fn_psv}{\code{character} -- model PSV file name with or without the PSV extension}
  \item{samples}{\code{numeric} --  matrix or data frame of samples;
    each column is a parameter, each row a sample}
  \bold{write_admb_cov2}
  \item{path}{\code{character} -- path to model}
  \item{cov_unbounded}{\code{numeric} -- covariance matrix in unbounded space}
  \item{hbf}{\code{numeric} -- \code{'hybrid_bounded_flag'} value; use \code{'hbf=1'} for HMC}
  \item{fn_cov}{\code{character} -- name of the ADMB covariance file}
  \bold{stop_quietly}
  \item{...}{\code{dots} -- arguments passed to \code{'stop_quietly'} to be written out as a message}
}
\details{
  \code{'run_adnuts()'}\cr \code{'path'} is the directory in which the MLE will be run;
  a subdirectory of this, called \code{'mcmc'}, is where the MCMC will
  be run using the \code{'NUTS'} algorithm.
  Inside the \code{'mcmc'} directory, several temporary subdirectories
  will be created, one for each MCMC chain labeled \code{'chain_*'},
  aka CPU number used in the parallel execution. 
  These will disappear once the run has completed and the output has been merged.

  \code{'fix.posteriors()'}\cr Remove multiple header/datarows found in posteriors
  files to leave just one.

  \code{'load_ss_files()'}\cr Load all the SS files for output and input,
  and return the model object.

  \code{'calc.mcmc()'}\cr Return a list of mcmc calculations, e.g., quantiles
  for various values.

  \code{'create.key.nuisance_posteriors_files()'}\cr Undocumented mystery.

  \code{'load_models()'}\cr Load models from files created using \code{'create_rds_file()'}.

  \code{'modify_starter_mcmc_type()'}\cr Change the MCMC setting in an 
  SS starter file so the extra MCMC files are not created.

  \code{'write_psv2()'}\cr Write an ADMB PSV file.

  \code{'write_admb_cov2()'}\cr Write a covariance matrix to \code{'admodel.cov'}.

  \code{'stop_quietly()'}\cr Stop the program without issuing an error message.
}
\value{
  Nothing in particular.
}
\references{
  Monnahan, C.C. and Kristensen, K. 2018. 
  \href{https://doi.org/10.1371/journal.pone.0197954}{%
  No-U-turn sampling for fast Bayesian inference in ADMB and TMB: Introducing the adnuts and tmbstan R packages}.
  \emph{PLoS ONE} \bold{13}(5). e0197,954.
}
\author{
  \href{mailto:chrisgrandin@gmail.com}{Chris Grandin}, Program Head -- Pacific Hake (PBS), Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Remote office, Vancouver BC\cr
  Last modified \code{Rd: 2024-08-26}
}
\note{
  This suite of functions, created by Chris Grandin, works on a Linux
  server maintained by Chris Grandin.
  Heaven forbid that he should vanish.

  The R package \code{'rstan'} imports \code{'rstan'}; however,
  installation of the latter on Windows using \code{'R-devel'} remains elusive.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSsynth}:\cr
  \code{\link[PBSsynth]{calcRhat}},
  \code{\link[PBSsynth]{gatherMCMC}},
  \code{\link[PBSsynth]{load_extra_mcmc}},
  \code{\link[PBSsynth]{plotSS.compo}},
  \code{\link[PBSsynth]{plotSS.pmcmc}},
  \code{\link[PBSsynth]{plotSS.senso}},
  \code{\link[PBSsynth]{prepCP}},
  \code{\link[PBSsynth]{tabSS.decision}}\cr
In package \pkg{r4ss}:\cr
  \code{\link[r4ss]{SSgetMCMC}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{nutsFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.
