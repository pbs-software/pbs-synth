# Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States of America.
# Foreign copyrights may apply. See copyright.txt for more information.
# user support available at: NMFS.Stock.Synthesis@noaa.gov
# user info available at: https://vlab.ncep.noaa.gov/group/stock-synthesis
# V3.30.18; 2021-xx-xx safe; Stock Synthesis by Richard Methot (NOAA) using ADMB 12.3
# C input names of data and control files
data.30.00.ss #_dat_file -- data file name
control.30.00.ss #_ctl_file -- control file name
0 #_init_par_vals -- 0=use init values in control file; 1=use ss.par
0 #_run_disp_det (0=none other than ADMB outputs; 1=one brief line per iteration; 2=fuller display per iteration)
1 #_det_age_str_rep -- detailed output (0=minimal for data-limited; 1=high (w/ wtatage.ss_new); 2=brief)
0 #_write_iter1_det -- write 1st iteration details to echoinput.sso file (0;1)
0 #_write_par_trace -- write parm values to ParmTrace.sso (0=no; 1=good;active; 2=good;all; 3=every_iter;all_parms; 4=every;active)
0 #_write_cum_rep -- write to cumreport.sso (0=no; 1=like&timeseries; 2=add survey fits)
1 #_calc_full_priors -- Include prior_like for non-estimated parameters (0;1)
1 #_use_soft_bounds -- Use Soft Boundaries to aid convergence (0;1) (recommended)
1 #_data_file_out -- Number of data files to produce: 1st is input; 2nd is estimates; 3rd and higher are bootstrap
8 #_turn_off_est -- Turn off estimation for parameters entering after this phase
0 #_mcmc_burn -- MCeval burn interval
1 #_mcmc_thin -- MCeval thin interval
0 #_jit_init_par -- jitter initial parm value by this fraction
-1 #_sd_rep_syr -- min yr for sdreport outputs (-1 for styr)
-2 #_sd_rep_eyr -- max yr for sdreport outputs (-1 for endyr+1; -2 for endyr+Nforecastyrs
0 #_sd_rep_xyr -- N individual STD years
# COND_xyr -- vector of year values
0.0001 #_fin_conv_crit -- final convergence criteria (e.g. 1.0e-04)
0 #_retro_yr -- retrospective year relative to end year (e.g. -4)
1 #_min_age_b_sum -- for calc of summary biomass
1 #_dp_bias -- Depletion basis:  denom is: 0=skip; 1=rel X*SPB0; 2=rel SPBmsy; 3=rel X*SPB_styr; 4=rel X*SPB_endyr
1 #_dp_denom_frac -- Fraction (X) for Depletion denominator (e.g. 0.4)
2 #_spr_rep_basis -- SPR_report_basis:  0=skip; 1=(1-SPR)/(1-SPR_tgt); 2=(1-SPR)/(1-SPR_MSY); 3=(1-SPR)/(1-SPR_Btarget); 4=rawSPR
1 #_ann_f_units -- Annual_F_units: 0=skip; 1=exploitation(Bio); 2=exploitation(Num); 3=sum(Apical_F's); 4=true F for range of ages; 5=unweighted avg. F for range of ages
# 10 15 #_COND_f_units -- min and max age over which average F will be calculated with F_reporting=4 or 5
0 #_f_rep_basis: 0=raw_annual_F; 1=F/Fspr; 2=F/Fmsy ; 3=F/Fbtgt; where F means annual_F
0 #_mcmc_out_det -- MCMC output detail: integer part (0=default; 1=adds obj func components); and decimal part (added to SR_LN(R0) on first call to mcmc)
0 #_alk_tol -- ALK tolerance (example 0.0001)
# -1 #_COND_seed_boot -- random number seed for bootstrap data (-1 to use long(time) as seed): 1585083947
3.3 #_eof -- check value for end of starter file and for version control
