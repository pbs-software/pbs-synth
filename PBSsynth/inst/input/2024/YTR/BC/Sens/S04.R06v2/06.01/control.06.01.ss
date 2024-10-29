# Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States of America.
# Copyright (c) 2008-2021 ADMB Foundation and Regents of the University of California
# user support available at: NMFS.Stock.Synthesis@noaa.gov
# user info available at: https://vlab.ncep.noaa.gov/group/stock-synthesis
# V3.30.22.01; 2022-09-30; Stock Synthesis 3 by Richard Methot (NOAA) using ADMB-13.1 safe libraries was compiled using GNU C++ 12.2.0 64bit on Apr 27 2023
# C growth parameters are estimated
# C spawner-recruitment bias adjustment Not tuned For optimality
# data_and_control_files: rebsn_data.ss // rebsn_control.ss
0 #_use_watage -- 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1 #_npatt_growth -- N Growth Patterns (Growth Patterns; Morphs; Bio Patterns; GP are terms used interchangeably in SS)
1 #_nplat_growth -- N platoons within each growth pattern
# 1 #_COND_ctl_plat_sd_ratio -- platoon_within/between_stdev_ratio (no read if N_platoons=1)
# 1 #_COND_ctl_plat_vector_dist -- vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
4 #_r_dist_meth -- Recruitment distribution method for parameters: 2=main effects for GP; settle timing; and area; 3=each settle entity; 4=none; no parameters (only if growth pattern x settlement x area = 1)
1 #_spr_opts -- Spawner-Recruitment options (not implement yet; but required); Future usage: 1=global; 2=by area
1 #_r_nsettle -- Number of recruitment settlement assignments
0 #_future_useless -- Future feature; not implement yet but required
#_GP Month Area Age at settlement #_(for each settlement assignment)
1 1 1 0 #_patt_mon_area_age
# 1 1 2 0 #_patt_mon_area_age
# 1 1 3 0 #_patt_mon_area_age
#
# 0 #_COND_nmove_defn -- N_movement_definitions goes here if Nareas > 1
# 1 #_COND_first_age_move -- first age that moves (real age at begin of season; not integer) also cond on do_migration>0
# 1 1 2 4 10 #_COND_move_defn -- example move definition for seas=1; morph=1; source=1 dest=2; age1=4; age2=10
#
#_Time_Blocks
0 #_nblock_patts #_Number of block patterns.
# 1 #_blocks_per_pattern
# #_begin_end_yrs_blocks
# 1970 1970
#
# Auto-generation: controls for all timevary parameters
1 #_env_block_dev -- Environmental/Block/Deviation: adjust method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
1 1 1 1 1 #_auto_gen_vals -- Auto-generation: 5 values control AG for par block sections: 1=biology; 2=SR; 3=Q; 4=tag (future); 5=selectivity
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line from control file; 2 = read then autogen if parm min==-12345
#
# Available timevary codes
# Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
# Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
# EnvLinks: 1: P(y)=P_base*exp(TVP*env(y)); 2: P(y)=P_base+TVP*env(y); 3: null; 4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
# DevLinks: 1: P(y)*=exp(dev(y)*dev_se; 2: P(y)+=dev(y)*dev_se; 3: random walk; 4: zero-reverting random walk with rho; 21-24 keep last dev for rest of years
#
# Prior_types: 0=none; 1=symmetric beta; 2=full beta (CASAL); 3=lognormal w/out bias adjustment; 4=lognormal w/ bias adjustment; 5=gamma; 6=normal
#
# #_BIOLOGY----------------------------------------------p.80
#
0 #_bio_m_opts -- natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
# 2 #_bio_n_brkpts -- COND 1: number of breakpoints; then read a vector of ages for these breakpoints. Later, per sex x GP, read N parameters for the natural mortality at each breakpoint.
# 9 10 #_bio_vec_brkpts -- COND 1: vector of age breakpoints
1 #_bio_growth_mod -- GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1 #_bio_gmod_amin -- Reference age for first size-at-age (post-settlement) parameter
45 #_bio_gmod_amax -- Reference age for second size-at-age parameter (999 to use as L infinity)
-998 #_bio_gmod_exp_decay -- exponential decay for growth above maxage (plus group: fixed at 0.20 in SS v.3.24; should approx initial Z); -998 disable growth above maxage (plus group); -999 replicates calculation in SS v.3.24
0 #_bio_gmod_holder -- placeholder for future growth feature
0 #_bio_sd_add_laa -- SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_bio_cv_patt -- CV_Growth_Pattern:0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
3 #_bio_mat_opt -- maturity_option (starts at age 0): 1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
0 0 0 0 0 0.148148154 0.215348908 0.302217729 0.405996853 0.518911222 0.629926425 0.728716702 0.809126071 0.869955096 0.913471789 0.943374407 0.963357732 0.976464824 0.984956668 0.990414604 0.993904578 0.996128844 0.997543464 0.998441954 0.999012143 0.999373794 0.999603098 0.999748457 0.999840589 0.999898979 0.999935983 0.999959433 0.999974294 0.99998371 0.999989678 0.999993459 0.999995855 0.999997373 0.999998336 0.999998945 0.999999332 0.999999577 0.999999732 0.99999983 0.999999892 0.999999932 #_bio_mat_ogive
# Age_Fecundity by growth pattern from wt-at-age.ss now invoked by read bodywt flag
9.06336065 #_bio_age_mature -- first mature age (value is overridden if maturity option is 3 or 4)
1 #_bio_fec_opt -- fecundity option: (1) eggs=Wt*(a+b*Wt); (2) eggs=a*L^b; (3) eggs=a*Wt^b; (4) eggs=a+b*L; (5) eggs=a+b*W
0 #_bio_herm_opt -- hermaphroditism option:0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_bio_par_offset_meth -- parameter offset approach (1=none; 2=M; G; CV_G as offset from female-GP1; 3=like SS2 V1.x)
#
# Mortality_Growth_Parameters (mgp) [p.88] note: gp=Growth Pattern
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var devlink devminyr devmaxyr dev_PH Block Block_Fxn label
# GP1=5ABC, GP2=3CD, GP3=5DE; M: prior means based on estimated values in previous assessments
# BC coast Sex=1 (female) GP=1 Natural Mortality
0.02 0.2 0.08 0.08 0.024 6 4 0 0 0 0 0 0 0 #_mgpF_gp1_M
# 0.02 0.2 0.12 0.12 0.03 0 5 0 0 0 0 0 0 0 #_mgpF_gp1_M2
# Sex=1 (female) GP=1 Growth 0.3
2 20 11.87231054 11.87231054 3.561693161 0 -2 0 0 0 0 0 0 0 #_mgpF_gp1_L_Amin
30 100 55.68063404 55.68063404 16.70419021 0 -2 0 0 0 0 0 0 0 #_mgpF_gp1_L_Amax
0.01 0.3 0.141854641 0.141854641 0.042556392 0 -2 0 0 0 0 0 0 0 #_mgpF_gp1_VB_K
0.03 2 0.034444192 0.034444192 0.010333258 0 -2 0 0 0 0 0 0 0 #_mgpF_gp1_CV_young
0.03 2 0.052008299 0.052008299 0.01560249 0 -2 0 0 0 0 0 0 0 #_mgpF_gp1_CV_old
# Sex=1 (female) GP=1 Weight-Length
-3 3 1.18065E-05 1.18065E-05 3.54196E-06 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_WL_alpha
-3 4 3.071384121 3.071384121 0.921415236 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_WL_beta
# Sex=1 (female) GP=1 Maturity & Fecundity
-3 3 9.06336065 9.06336065 2.719008195 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_Mat_50pc -- Value ignored for maturity option 3, 4, and 6.
-3 3 -0.25 -0.25 -0.075 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_Mat_slope -- Logistic slope (must have negative value). Value ignored for maturity option 3, 4, and 6.
-3 3 1 1 0.3 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_Fec_alpha -- (eggs/kg)
-3 3 0 0 0 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_Fec_beta -- (eggs/kg)
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var devlink devminyr devmaxyr dev_PH Block Block_Fxn label
# GP1=5ABC, GP2=3CD, GP3=5DE; prior means based on estimated values in previous assessments
# BC coast Sex=2 (male) GP=1 Natural Mortality
0.02 0.2 0.08 0.08 0.024 6 4 0 0 0 0 0 0 0 #_mgpM_gp1_M
# 0.02 0.2 0.08 0.08 0.02 0 5 0 0 0 0 0 0 0 #_mgpM_gp1_M2
# Sex=2 (male) GP=1 Growth
2 20 11.02203001 11.02203001 3.306609003 0 -2 0 0 0 0 0 0 0 #_mgpM_gp1_L_Amin
30 100 47.96144006 47.96144006 14.38843202 0 -2 0 0 0 0 0 0 0 #_mgpM_gp1_L_Amax
0.01 0.3 0.199279322 0.199279322 0.059783796 0 -2 0 0 0 0 0 0 0 #_mgpM_gp1_VB_K
0.03 2 0.063668527 0.063668527 0.019100558 0 -2 0 0 0 0 0 0 0 #_mgpM_gp1_CV_young
0.03 2 0.030148524 0.030148524 0.009044557 0 -2 0 0 0 0 0 0 0 #_mgpM_gp1_CV_old
# Sex=2 (male) GP=1 Weight-Length
-3 3 1.09801E-05 1.09801E-05 3.29404E-06 0 -50 0 0 0 0 0 0 0 #_mgpM_gp1_WL_alpha
-3 4 3.093536902 3.093536902 0.928061071 0 -50 0 0 0 0 0 0 0 #_mgpM_gp1_WL_beta
# Recruitment Distribution (only for methods 2 or 3) -- parameter values are the natural log of apportionment weight
# -5 5 0 0 0 0 -4 0 0 0 0 0 0 0 #_bio_r_dist_grp1 # if recr option=3, only need explicit settlement parameters
# -5 5 0 0 1 6 3 0 2 1935 2010 5 0 0 #_bio_r_dist_area1
# -5 5 0 0 1 6 3 0 2 1975 2010 5 0 0 #_bio_r_dist_area2
# -5 5 0 0 1 0 -3 0 0 0 0 0 0 0 #_bio_r_dist_area3
# -5 5 0 0 0 0 -3 0 0 0 0 0 0 0 #_bio_r_dist_month1 # if recr option=3, only need explicit settlement parameters
# Cohort growth dev base
0.2 5 1 1 1 0 -1 0 0 0 0 0 0 0 #_bio_cohort_grow_dev
# Fraction female; by GP
0.001 0.999 0.5 0.5 0.15 0 -50 0 0 0 0 0 0 0 #_bio_gp1_fracF
#
## time varying recruitment distribution parameters (short lines)
#_LO HI INIT PRIOR PR_SD PR_type PHASE
# 0.0001 12 1 1 0.5 0 -5 #_bio_rdist_area1_se
# -0.99 0.99 1 0 0.5 0 -6 #_bio_rdist_area1_rho
# 0.0001 12 1 1 0.5 0 -5 #_bio_rdist_area2_se
# -0.99 0.99 1 0 0.5 0 -6 #_bio_rdist_area2_rho
#
# seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_bio_seas_WL1F bio_seas_WL2F bio_seas_mat1 bio_seas_mat2 bio_seas_fec1 bio_seas_fec2 bio_seas_WL1M bio_seas_WL2M malewtlen2 bio_seas_L1 bio_seas_K
#
# #_SPAWNER-RECRUITMENT----------------------------------p.97
#
3 #_spr_opts -- Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0 #_spr_use_h -- 0/1 to use steepness (h) in initial equilibrium recruitment calculation
0 #_spr_future_sigr -- 0/1 future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var use_dev dev_mnyr dev_mxyr dev_PH Block Blk_Fxn name
1 16 10 10 10 6 1 0 0 0 0 0 0 0 #_spr_ln(r0)
# 0.02 1.601 0.7 0.7 0.418310321 0 -5 0 0 0 0 0 0 0 #_spr_bh_steep #SS mu= 0.430107527 SS tau= 5.143892697 SS alpha= 4.574000224 SS beta= 2.212426966
0.2 1 0.67 0.67 0.17 2 5 0 0 0 0 0 0 0 #_spr_bh_steep #_format used by Hake assessment but substitute Forrest et al. (2010) posterior predictive prior for BH
-15 15 1.2 1 0 0 -1 0 0 0 0 0 0 0 #_spr_sigma_r
-5 5 0 0 99 0 -50 0 0 0 0 0 0 0 #_spr_regime
0 2 0 1 99 0 -50 0 0 0 0 0 0 0 #_spr_autocorr
#
# Recruitment deviation setup
2 #_r_rdev -- 0=none; 1=devvector (sum to zero); 2=simple deviations (sum close to zero); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1935 #_r_rdev_yr1 -- first year of main recr_devs; early devs can preceed this era
2015 #_r_rdev_yrN -- last year of main recr_devs; forecast devs start in following year # POP=change to 1935:2010
2 #_r_rdev_phase -- recdev phase
1 #_r_advance_opts -- (0/1) to read 13 advanced options
1920 #_r_rdev_early_start -- (0=none; neg value makes relative to recdev_start)
2 #_r_rdev_early_phase
0 #_r_fcast_phase -- forecast recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_r_lambda -- for Fcast_recr_like occurring before endyr+1
1942 #_r_yrN_nobias_adj -- last year with no bias adj in MPD; begin of ramp
1978 #_r_yr1_fullbias_adj -- first year with full bias adj in MPD; begin of plateau
2009 #_r_yrN_fullbias_adj -- last year with full bias adj in MPD
2015 #_r_yrN_ramp -- end year for ramp in MPD (can be in forecast to shape ramp; but SS sets bias_adj to 0.0 for fcast yrs)
0.65 #_r_maxbias_adj -- maximum bias adj in MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
0 #_r_period_cycles -- period of cycles in recruitment (N parms read below)
-5 #_r_rdev_min -- minimum recruitment deviation
5 #_r_rdev_max -- maximum recruitment deviation
0 #_r_read_rdev -- read in the recruitment deviations
# end of advanced SR options
#
# placeholder for full parameter lines for recruitment cycles
# read specified recr devs
# Yr Input_value
#
# explicit recruitment deviations (example):
# 1935 1936 1937 1938 ... 2021 2022
# -0.5 0.5 -0.5 0.5 0 -0.5 0.5
# implementation error by year in forecast: 0 0 0 0 0 0 0 0 0 0
#
# #_FISHING MORTALITY METHOD----------------------------p.110
# #
0.1 #_f_init_guess -- F ballpark value in units of annual_F
-1966 #_f_yr_guess -- F ballpark year (neg value to disable)
3 #_f_method -- F Method: 1=Pope's (discrete); 2=Baranov (continuous); 3=Hybrid; 4=Fleet-specific hybrid (preferred)
2.9 #_f_max -- max F or harvest rate; depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
4 #_f_niter_tune -- COND 3: N iterations for tuning F in hybrid method (recommend 3 to 7)
#
# initial Fs for the two fisheries unless starting from unexploited (when catch in year -999 = 0)
# LO HI INIT PRIOR PR_SD PR_type PHASE label
# 0 3 0.01 0 99 0 1 #_f_init_trawl
# 0 3 0.01 0 99 0 1 #_f_init_other
# F rates by fleet (example):
# Yr: 1971 1972 1973 ... 2011 2012
# seas: 1 1 1 1 1 1
# FISHERY 0 0.01 0.02 0 0.1 0.1
#
# #_CATCHABILITY----------------------------------------p.112
#
# Q_setup for fleets with cpue or survey data
# 1: fleet number
# 2: link type: (1=simple q; 1 parm; 2=mirror simple q; 1 mirrored parm; 3=q and power; 2 parm; 4=mirror with offset; 2 parm)
# 3: extra input for link; i.e. mirror fleet#_or dev index number
# 4: 0/1 to select extra sd parameter
# 5: 0/1 for biasadj or not when using informative prior
# 6: 0/1 to float
#_fleet_num link_type link_info extra_sd bias_adjust Q_float
# 1 1 0 0 0 1 #_q_setup_idx_bc
2 1 0 0 0 1 #_q_setup_idx_qcs
3 1 0 0 0 1 #_q_setup_idx_wcvi
4 1 0 0 0 1 #_q_setup_idx_wchg
5 1 0 0 0 1 #_q_setup_idx_hs
6 1 0 0 0 1 #_q_setup_idx_gig
7 1 0 0 0 1 #_q_setup_idx_nmfs
# 8 1 0 0 0 1 #_q_setup_idx_wvh
-9999 0 0 0 0 0 #_end_fleet_q -- terminator for fleet/survey catchability
#
# Q_parms(if_any); Qunits_are_ln(q)
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var use_dev dev_mnyr dev_mxyr dev_PH Block Block_Fxn label
# -15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_bc
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_qcs
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_wcvi
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_wchg
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_hs
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_gig
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_nmfs
# -15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_wvh
# no timevary Q parameters
#
# ##_SELECTIVITY-AND-DISCARD-----------------------------p.116
#
# size_selex_patterns
# Pattern:_0; parm=0; selex=1.0 for all sizes
# Pattern:_1; parm=2; logistic; with 95% width specification
# Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
# Pattern:_15; parm=0; mirror another age or length selex
# Pattern:_6; parm=2+special; non-parm len selex
# Pattern:_43; parm=2+special+2; like 6; with 2 additional param for scaling (average over bin range)
# Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
# Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
# Pattern:_21; parm=2+special; non-parm len selex; read as pairs of size; then selex
# Pattern:_22; parm=4; double_normal as in CASAL
# Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
# Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL); using joiners
# Pattern:_25; parm=3; exponential-logistic in size
# Pattern:_27; parm=3+special; cubic spline
# Pattern:_42; parm=2+special+3; // like 27; with 2 additional param for scaling (average over bin range)
# discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
0 0 0 0 #_s_pat_len_bc
0 0 0 0 #_s_pat_len_qcs
0 0 0 0 #_s_pat_len_wcvi
0 0 0 0 #_s_pat_len_wchg
0 0 0 0 #_s_pat_len_hs
0 0 0 0 #_s_pat_len_gig
0 0 0 0 #_s_pat_len_nmfs
# 0 0 0 0 #_s_pat_len_wvh
#
# age_selex_patterns
# Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
# Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
# Pattern:_11; parm=2; selex=1.0 for specified min-max age
# Pattern:_12; parm=2; age logistic
# Pattern:_13; parm=8; age double logistic
# Pattern:_14; parm=nages+1; age empirical
# Pattern:_15; parm=0; mirror another age or length selex
# Pattern:_16; parm=2; Coleraine - Gaussian
# Pattern:_17; parm=nages+1; empirical as random walk N parameters to read can be overridden by setting special to non-zero
# Pattern:_41; parm=2+nages+1; // like 17; with 2 additional param for scaling (average over bin range)
# Pattern:_18; parm=8; double logistic - smooth transition
# Pattern:_19; parm=6; simple 4-parm double logistic with starting age
# Pattern:_20; parm=6; double_normal;using joiners
# Pattern:_26; parm=3; exponential-logistic in age
# Pattern:_27; parm=3+special; cubic spline in age
# Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
# Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
20 0 3 0 #_s_pat_age_bc TRAWL_FISHERY_BC
20 0 3 0 #_s_pat_age_qcs QCS_SYNOPTIC
20 0 3 0 #_s_pat_age_wcvi WCVI_SYNOPTIC
20 0 3 0 #_s_pat_age_wchg WCHG_SYNOPTIC
20 0 3 0 #_s_pat_age_hs HS_SYNOPTIC
15 0 0 2 #_s_pat_age_gig GIG_HISTORICAL
20 0 3 0 #_s_pat_age_nmfs NMFS_TRIENNIAL
# 15 0 0 4 #_s_pat_age_wvh WCVI_HISTORICAL
#
# AGE_SELECTIVITY_PATTERN_20
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var use_dev dev_mnyr dev_mxyr dev_PH Block Blk_Fxn
# 1_TRAWL_FISHERY_BC AgeSelex
# female_selectivity_pattern_20
5 40 11 11 11 6 3 0 0 0 0 0 0 0 #_s_p1_trawl_bc -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_trawl_bc -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 2 2 2 6 4 0 0 0 0 0 0 0 #_s_p3_trawl_bc -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_trawl_bc -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_trawl_bc -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_trawl_bc -- final: selectivity at last bin; as logistic between 0 and 1
# male_selectivity_option_3
# p1size (age) at which a dogleg occurs
-10 10 0 0 10 6 4 0 0 0 0 0 0 0 #_s_p1_trawl_bc_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_trawl_bc_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_trawl_bc_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_trawl_bc_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_trawl_bc_male -- male parameter 5 is the apical selectivity for males
# 2_QCS_SYNOPTIC AgeSelex
# female_selectivity_pattern_20
5 40 11 11 11 6 3 0 0 0 0 0 0 0 #_s_p1_qcs -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_qcs -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 2 2 2 6 4 0 0 0 0 0 0 0 #_s_p3_qcs -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_qcs -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_qcs -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_qcs -- final: selectivity at last bin; as logistic between 0 and 1
# male_selectivity_option_3
# p1size (age) at which a dogleg occurs
-10 10 0 0 10 6 4 0 0 0 0 0 0 0 #_s_p1_qcs_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_qcs_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_qcs_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_qcs_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_qcs_male -- male parameter 5 is the apical selectivity for males
# 3_WCVI_SYNOPTIC AgeSelex
# female_selectivity_pattern_20
5 40 11 11 11 6 3 0 0 0 0 0 0 0 #_s_p1_wcvi -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_wcvi -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 2 2 2 6 4 0 0 0 0 0 0 0 #_s_p3_wcvi -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_wcvi -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_wcvi -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_wcvi -- final: selectivity at last bin; as logistic between 0 and 1
# male_selectivity_option_3
# p1size (age) at which a dogleg occurs
-10 10 0 0 10 6 4 0 0 0 0 0 0 0 #_s_p1_wcvi_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_wcvi_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_wcvi_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_wcvi_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_wcvi_male -- male parameter 5 is the apical selectivity for males
# 4_WCHG_SYNOPTIC AgeSelex
# female_selectivity_pattern_20
5 40 15 15 15 0 -3 0 0 0 0 0 0 0 #_s_p1_wchg -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_wchg -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 3 3 3 0 -4 0 0 0 0 0 0 0 #_s_p3_wchg -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_wchg -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_wchg -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_wchg -- final: selectivity at last bin; as logistic between 0 and 1
# male_selectivity_option_3
# p1size (age) at which a dogleg occurs
-10 10 0 0 10 0 -4 0 0 0 0 0 0 0 #_s_p1_wchg_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_wchg_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_wchg_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_wchg_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_wchg_male -- male parameter 5 is the apical selectivity for males
# 5_HS_SYNOPTIC AgeSelex
# female_selectivity_pattern_20
5 40 10 10 10 0 -3 0 0 0 0 0 0 0 #_s_p1_hs -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_hs -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 2 2 2 0 -4 0 0 0 0 0 0 0 #_s_p3_hs -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_hs -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_hs -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_hs -- final: selectivity at last bin; as logistic between 0 and 1
# male_selectivity_option_3
# p1size (age) at which a dogleg occurs
-10 10 0 0 10 0 -4 0 0 0 0 0 0 0 #_s_p1_hs_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_hs_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_hs_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_hs_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_hs_male -- male parameter 5 is the apical selectivity for males
# 6_GIG_HISTORICAL AgeSelex # mirrors fleet 2 (QCS synoptic)
# female_selectivity_pattern_20
# 5 40 11.5 11.5 11.5 6 3 0 0 0 0 0 0 0 #_s_p1_gig -- peak: beginning size (or age) for the plateau (=mu in Awatea)
# 0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_gig -- top: width of plateau; as logistic between peak and maximum length (or age)
# -15 15 2.2 2.2 2.2 6 4 0 0 0 0 0 0 0 #_s_p3_gig -- ascending width: parameter value is ln(width) (=log varL in Awatea)
# -15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_gig -- descending width: parameter value is ln(width) (=log varR in Awatea)
# -10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_gig -- initial: selectivity at first bin; as logistic between 0 and 1
# -10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_gig -- final: selectivity at last bin; as logistic between 0 and 1
# male_selectivity_option_3
# p1size (age) at which a dogleg occurs
# -10 10 0 0 10 6 4 0 0 0 0 0 0 0 #_s_p1_gig_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
# -15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_gig_male -- added to the third selectivity parameter (width of ascending side)
# -15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_gig_male -- added to the fourth selectivity parameter (width of descending side)
# -10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_gig_male -- added to the sixth selectivity parameter (selectivity at final size bin)
# 0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_gig_male -- male parameter 5 is the apical selectivity for males
# 7_NMFS_TRIENNIAL AgeSelex
# female_selectivity_pattern_20
1 40 11 11 11 6 3 0 0 0 0 0 0 0 #_s_p1_nmfs -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_nmfs -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 2 2 2 6 4 0 0 0 0 0 0 0 #_s_p3_nmfs -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_nmfs -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_nmfs -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_nmfs -- final: selectivity at last bin; as logistic between 0 and 1
# male_selectivity_option_3
# p1size (age) at which a dogleg occurs
-8 10 0 0 10 6 4 0 0 0 0 0 0 0 #_s_p1_nmfs_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_nmfs_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_nmfs_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_nmfs_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_nmfs_male -- male parameter 5 is the apical selectivity for males
# 8_WCVI_HISTORICAL AgeSelex # mirrors fleet 3 (WCVI synoptic)
#
# Dirichlet-Multinomial parameters controlling age-comp weights (also choose Dirichlet error structure in data file)
# -5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w1_ln(dm_parm)_bc
# -5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w2_ln(dm_parm)_qcs
# -5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w3_ln(dm_parm)_wcvi
# -5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w4_ln(dm_parm)_wchg #_Apparently do not specify fleets w/out AF data in this section
# -5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w5_ln(dm_parm)_hs
# -5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w6_ln(dm_parm)_gig
# -5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w7_ln(dm_parm)_nmfs
# -5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w8_ln(dm_parm)_gbr #_Apparently do not specify fleets w/out AF data in this section
# #_LO HI INIT PRIOR PR_SD PR_type PHASE env-var use_dev dev_mnyr dev_mxyr dev_PH Block Blk_Fxn
# no timevary selex parameters
#
0 #_s_use_2d_ar1 selectivity(0/1)
# no 2D_AR1 selex offset used
#
# #_TAG_RECAPTURE_PARAMETERS----------------------------p.138
#
# Tag loss and Tag reporting parameters go next
0 #_tag_custom -- TG_custom: 0=no read and autogen if tag data exist; 1=read
#_COND -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0 #_placeholder if no parameters
# no timevary parameters
#
# #_VARIANCE_ADJUSTMENT_FACTORS-------------------------p.140
#
# Input variance adjustments factors:
# 1=add_to_survey_CV
# 2=add_to_discard_stddev
# 3=add_to_bodywt_CV
# 4=mult_by_lencomp_N
# 5=mult_by_agecomp_N
# 6=mult_by_size-at-age_N
# 7=mult_by_generalized_sizecomp
# Note: adjustment factor 5 allegedly not needed if Dirichlet-Multinomial likelihood used (see hake example; Chris Grandin)
#_factor fleet value desc
# 1 1 0 #_vadj_cv_bc_cpue ## no CPUE for YTR
1 2 0 #_vadj_cv_qcs_idx
1 3 0 #_vadj_cv_wcvi_idx
1 4 0 #_vadj_cv_wchg_idx
1 5 0 #_vadj_cv_hs_idx
1 6 0 #_vadj_cv_gig_idx
1 7 0 #_vadj_cv_nmfs_idx
# 1 8 0 #_vadj_cv_gbr_idx
# SS_tune_comps(rlist; option="Francis") on Multinomial
5 1 2.50529132374868 #_vadj_af_trawl_bc
5 2 1.90800534408521 #_vadj_af_qcs
5 3 2.50220276676217 #_vadj_af_wcvi
# 5 5 1 #_vadj_af_hs
# 5 6 1 #_vadj_af_gig
5 7 1.54355740876591 #_vadj_af_nmfs
-9999 0 0 #_end_vadj -- terminator_variance_adjustment_factors
#
# #_LAMBDAS_(EMPHASIS_FACTORS)--------------------------p.142
#
1 #_lambda_max_phase -- maximum ambda phase
1 #_lambda_sd_off -- SD offset; must be 1 if any growthCV; sigmaR; or survey extraSD is an estimated parameter
# #_read 3 changes to default Lambdas (default value is 1.0)
# #_Like_comp codes: 1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch;
# #_10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
# #_like_comp fleet phase value sizefreq_method
-9999 1 1 1 1 #_lambda_max_end -- terminator for max lambda phase
#
# lambdas (for info only; columns are phases)
# 0 0 0 0 #_CPUE/survey:_1
# 1 1 1 1 #_CPUE/survey:_2
# 1 1 1 1 #_CPUE/survey:_3
# 1 1 1 1 #_lencomp:_1
# 1 1 1 1 #_lencomp:_2
# 0 0 0 0 #_lencomp:_3
# 1 1 1 1 #_agecomp:_1
# 1 1 1 1 #_agecomp:_2
# 0 0 0 0 #_agecomp:_3
# 1 1 1 1 #_size-age:_1
# 1 1 1 1 #_size-age:_2
# 0 0 0 0 #_size-age:_3
# 1 1 1 1 #_init_equ_catch
# 1 1 1 1 #_recruitments
# 1 1 1 1 #_parameter-priors
# 1 1 1 1 #_parameter-dev-vectors
# 1 1 1 1 #_crashPenLambda
# 0 0 0 0 #_f_ballpark_lambda
#
# #_CONTROLS_FOR VARIANCE_OF_DERIVED_QUANTITIES---------p.144
#
0 #_read_sd_rep -- (0/1) read specs for more stddev reporting
# 0 0 0 0 0 0 0 0 0 #_placeholder for selex_fleet; 1=len/2=age/3=both; year; N selex bins; 0 or Growth pattern; N growth ages; 0 or NatAge_area(-1 for all); NatAge_yr; N Natages
# placeholder for vector of selex bins to be reported
# placeholder for vector of growth ages to be reported
# placeholder for vector of NatAges ages to be reported
999 #_eof -- EOF for control
