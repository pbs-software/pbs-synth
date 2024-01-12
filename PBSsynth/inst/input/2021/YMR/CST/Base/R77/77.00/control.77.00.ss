#_Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States of America.
#_Foreign copyrights may apply. See copyright.txt for more information.
#_user support available at: NMFS.Stock.Synthesis@noaa.gov
#_user info available at: https://vlab.ncep.noaa.gov/group/stock-synthesis
#_V3.30.16.00; 2020-09-03 safe; Stock Synthesis by Richard Methot (NOAA) using ADMB 12.2
#_C growth parameters are estimated
#_C spawner-recruitment bias adjustment Not tuned For optimality
#_data_and_control_files: rebsn_data.ss // rebsn_control.ss
0 #_use_watage -- 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1 #_npatt_growth -- N Growth Patterns (Growth Patterns; Morphs; Bio Patterns; GP are terms used interchangeably in SS)
1 #_nplat_growth -- N platoons within each growth pattern
#_1 #_COND_ctl_plat_sd_ratio -- platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_1 #_COND_ctl_plat_vector_dist -- vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
4 #_r_dist_meth -- Recruitment distribution method for parameters: 2=main effects for GP; settle timing; and area; 3=each settle entity; 4=none; no parameters (only if growth pattern x settlement x area = 1)
1 #_spr_opts -- Spawner-Recruitment options (not implement yet; but required); Future usage: 1=global; 2=by area
1 #_r_nsettle -- Number of recruitment settlement assignments
0 #_future_useless -- Future feature; not implement yet but required
#_Growth Pattern Month Area Age at settlement#_(for each settlement assignment)
1 1 1 0 #_patt_mon_area_age
#
#_0 #_COND_nmove_defn -- N_movement_definitions goes here if Nareas > 1
#_1.0 #_COND_first_age_move -- first age that moves (real age at begin of season; not integer) also cond on do_migration>0
#_1 1 1 2 4 10 #_COND_move_defn -- example move definition for seas=1; morph=1; source=1 dest=2; age1=4; age2=10
#
0 #_nblock_patts
#_1 #_blocks_per_pattern
#_begin_end_yrs_blocks
#_1970 1970
#
#_Auto-generation: controls for all timevary parameters
1 #_env_block_dev -- Environmental/Block/Deviation: adjust method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
1 1 1 1 1 #_auto_gen_vals -- Auto-generation: 5 values control AG for par block sections: 1=biology; 2=SR; 3=Q; 4=tag (future); 5=selectivity
#_where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line from control file; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks: 1: P(y)=P_base*exp(TVP*env(y)); 2: P(y)=P_base+TVP*env(y); 3: null; 4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks: 1: P(y)*=exp(dev(y)*dev_se; 2: P(y)+=dev(y)*dev_se; 3: random walk; 4: zero-reverting random walk with rho; 21-24 keep last dev for rest of years
#
#_Prior_types: 0=none; 1=symmetric beta; 2=full beta (CASAL); 3=lognormal w/out bias adjustment; 4=lognormal w/ bias adjustment; 5=gamma; 6=normal
#
##_BIOLOGY----------------------------------------------p.80
#
0 #_bio_m_opts -- natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
1 #_bio_growth_mod -- GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1 #_bio_gmod_amin -- Reference age for first size-at-age (post-settlement) parameter
60 #_bio_gmod_amax -- Reference age for second size-at-age parameter (999 to use as L infinity)
-998 #_bio_gmod_exp_decay -- exponential decay for growth above maxage (plus group: fixed at 0.20 in SS v.3.24; should approx initial Z); -998 disable growth above maxage (plus group); -999 replicates calculation in SS v.3.24
0 #_bio_gmod_holder -- placeholder for future growth feature
0 #_bio_sd_add_laa -- SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_bio_cv_patt -- CV_Growth_Pattern:0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
3 #_bio_mat_opt -- maturity_option (starts at age 0): 1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
0 0 0 0 0 0.054889304 0.087974358 0.135228401 0.199353167 0.281852434 0.382176397 0.496991952 0.619838106 0.741396635 0.850484437 0.935676137 0.987252309 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 #_bio_mat_ogive
#_Age_Fecundity by growth pattern from wt-at-age.ss now invoked by read bodywt flag
1 #_bio_age_mature -- first mature age (value is overridden if maturity option is 3 or 4)
1 #_bio_fec_opt -- fecundity option: (1) eggs=Wt*(a+b*Wt); (2) eggs=a*L^b; (3) eggs=a*Wt^b; (4) eggs=a+b*L; (5) eggs=a+b*W
0 #_bio_herm_opt -- hermaphroditism option:0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_bio_par_offset_meth -- parameter offset approach (1=none; 2=M; G; CV_G as offset from female-GP1; 3=like SS2 V1.x)
#
#_Mortality_Growth_Parameters (mgp) [p.88] note: gp=Growth Pattern
#_LO HI INIT PRIOR PR_SD PR_type PHASE env_var devlink devminyr devmaxyr dev_PH Block Block_Fxn label
#_Sex:1 (female) Growth Pattern:1 Natural Mortality (YMR Gertseva age 90)
0.02 0.15 0.040 0.040 0.011 0 -3 0 0 0 0 0 0 0 #_mgpF_gp1_M
#_Sex:1 (female) Growth Pattern:1 Growth
2 20 15.91308821 5 1 6 -2 0 0 0 0 0 0 0 #_mgpF_gp1_L_Amin
30 100 48.17227503 45 9 6 -2 0 0 0 0 0 0 0 #_mgpF_gp1_L_Amax
0.01 0.3 0.115689018 0.1 0.025 6 -2 0 0 0 0 0 0 0 #_mgpF_gp1_VB_K
0.03 2 0.084294477 0.65 0.3 0 -2 0 0 0 0 0 0 0 #_mgpF_gp1_CV_young
0.03 2 0.046779459 0.65 0.3 0 -2 0 0 0 0 0 0 0 #_mgpF_gp1_CV_old
#_Sex:1 (female) Growth Pattern:1 Weight-Length (YMR Surv CST 210203)
-3 3 7.83513E-06 7.83513E-06 99 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_WL_alpha
-3 4 3.182609251 3.182609251 99 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_WL_beta
#_Sex:1 (female) Growth Pattern:1 Maturity & Fecundity
-3 3 11.0244863 24 99 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_Mat_50pc -- Value ignored for maturity option 3, 4, and 6.
-3 3 -0.25 -0.25 99 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_Mat_slope -- Logistic slope (must have negative value). Value ignored for maturity option 3, 4, and 6.
-3 3 1 1 99 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_Fec_alpha -- (eggs/kg)
-3 3 0 0 99 0 -50 0 0 0 0 0 0 0 #_mgpF_gp1_Fec_beta -- (eggs/kg)
#_LO HI INIT PRIOR PR_SD PR_type PHASE env_var devlink devminyr devmaxyr dev_PH Block Block_Fxn label
#_Sex:2 (male) Growth Pattern:1 Natural Mortality (YMR Gertseva age 90)
0.02 0.15 0.040 0.040 0.011 0 -3 0 0 0 0 0 0 0 #_mgpM_gp1_M
#_Sex:2 (male) Growth Pattern:1 Growth
2 20 14.94831083 5 1 6 -2 0 0 0 0 0 0 0 #_mgpM_gp1_L_Amin
30 100 46.64549546 45 9 6 -2 0 0 0 0 0 0 0 #_mgpM_gp1_L_Amax
0.01 0.3 0.128827999 0.1 0.025 6 -2 0 0 0 0 0 0 0 #_mgpM_gp1_VB_K
0.03 2 0.080770028 0.65 0.4 0 -2 0 0 0 0 0 0 0 #_mgpM_gp1_CV_young
0.03 2 0.03901777 0.65 0.4 0 -2 0 0 0 0 0 0 0 #_mgpM_gp1_CV_old
#_Sex:2 (male) Growth Pattern:1 Weight-Length (YMR Surv CST 210203)
-3 3 6.45932E-06 6.45932E-06 99 0 -50 0 0 0 0 0 0 0 #_mgpM_gp1_WL_alpha
-3 4 3.241045666 3.241045666 99 0 -50 0 0 0 0 0 0 0 #_mgpM_gp1_WL_beta
#_Sex:2 (male) Growth Pattern:1 Maturity & Fecundity
#_-3 3 24 24 99 0 -50 0 0 0 0 0 0 0 #_mgpM_gp1_Mat_50pc
#_-3 3 -0.25 -0.25 99 0 -50 0 0 0 0 0 0 0 #_mgpM_gp1_Mat_slope
#_-3 3 1 1 99 0 -50 0 0 0 0 0 0 0 #_mgpM_gp1_Fec_alpha
#_-3 3 0 0 99 0 -50 0 0 0 0 0 0 0 #_mgpM_gp1_Fec_beta
#_Recruitment Distribution(only for methods 2 or 3)
#_0 2 1 1 99 0 -50 0 0 0 0 0 0 0 #_bio_r_dist_grp1
#_0 2 1 1 99 0 -50 0 0 0 0 0 0 0 #_bio_r_dist_area1
#_0 2 1 1 99 0 -50 0 0 0 0 0 0 0 #_bio_r_dist_time1
#_Cohort growth dev base
0.2 5 1 1 1 0 -1 0 0 0 0 0 0 0 #_bio_cohort_grow_dev
#_Fraction female; by GP
0.001 0.999 0.5 0.5 0.5 0 -50 0 0 0 0 0 0 0 #_bio_gp1_fracF
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_bio_seas_WL1F bio_seas_WL2F bio_seas_mat1 bio_seas_mat2 bio_seas_fec1 bio_seas_fec2 bio_seas_WL1M bio_seas_WL2M malewtlen2 bio_seas_L1 bio_seas_K
#
##_SPAWNER-RECRUITMENT----------------------------------p.97
#
3 #_spr_opts -- Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0 #_spr_use_h -- 0/1 to use steepness (h) in initial equilibrium recruitment calculation
0 #_spr_future_sigr -- 0/1 future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var use_dev dev_mnyr dev_mxyr dev_PH Block Blk_Fxn name
1 16 8 8 8 6 1 0 0 0 0 0 0 0 #_spr_ln(r0)
0.02 1.601 0.7 0.7 0.418310321 0 -1 0 0 0 0 0 0 0 #_spr_bh_steep   #SS mu= 0.430107527 SS tau= 5.143892697 SS alpha= 4.574000224 SS beta= 2.212426966
-15 15 0.9 1 0 0 -1 0 0 0 0 0 0 0 #_spr_sigma_r
-5 5 0 0 99 0 -50 0 0 0 0 0 0 0 #_spr_regime
0 2 0 1 99 0 -50 0 0 0 0 0 0 0 #_spr_autocorr
#
#_Recruitment deviation setup
1 #_r_rdev -- 0=none; 1=devvector (sum to zero); 2=simple deviations (sum close to zero); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1950 #_r_rdev_yr1 -- first year of main recr_devs; early devs can preceed this era
2012 #_r_rdev_yrN -- last year of main recr_devs; forecast devs start in following year
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
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
#_read specified recr devs
#_Yr Input_value
#
#_explicit recruitment deviations
#1935 1936 1937 1938 1939 1940 1941 1942 1943 1944 1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020
#-0.499095 -0.482613 -0.434728 -0.377828 -0.352243 -0.359967 -0.325885 -0.215498 -0.136413 -0.232993 -0.357066 -0.338344 -0.185062 -0.000775487 0.0987642 -0.0206452 -0.133815 -0.0355136 0.277344 0.806261 0.671974 0.489025 0.73955 0.590117 0.240556 0.205738 0.671658 0.679474 0.0682313 -0.014923 0.143063 -0.000800232 -0.207249 -0.148113 -0.0192354 -0.00582447 0.0983492 -0.00722273 -0.300078 -0.348326 -0.0259654 0.219585 0.121448 0.259239 0.167265 0.138322 0.634521 0.439072 0.365514 1.32991 0.697917 0.382806 0.151452 -0.154913 -0.325421 -0.169825 0.464605 0.707631 -0.0144109 -0.23701 -0.211188 -0.0914514 0.0179154 0.0326441 -0.0444984 -0.108368 -0.14137 -0.167997 -0.193027 -0.212822 -0.229387 -0.241817 -0.250301 -0.256754 -0.261921 -0.266045 -0.269971 -0.272281 -0.27559 -0.276484 -0.277171 -0.277685 -0.279001 -0.279006 -0.279006 -0.279006
#_implementation error by year in forecast: 0 0 0 0 0 0 0 0 0 0
#
##_FISHING MORTALITY METHOD----------------------------p.110
#
0.1 #_f_init_guess -- F ballpark value in units of annual_F
-1966 #_f_yr_guess -- F ballpark year (neg value to disable)
1 #_f_method -- F Method: 1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
0.8 #_f_max -- max F or harvest rate; depends on F_Method (u_max for Pope's method -- use in Awatea)
#_no additional F input needed for Fmethod 1
#_if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
#_if Fmethod=3; read N iterations for tuning for Fmethod 3
#4 #_f_niter_tune -- N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial Fs for the two fisheries unless starting from unexploited (when catch in year -999 = 0)
#_LO HI INIT PRIOR PR_SD PR_type PHASE label
#_0 3 0.01 0 99 0 1 #_f_init_trawl
#_0 3 0.01 0 99 0 1 #_f_init_other
#_F rates by fleet
#_Yr: 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011
#_seas: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#_FISHERY 0 0.00211081 0.010609 0.0107037 0.0217063 0.0333334 0.0459509 0.0599453 0.0757167 0.107737 0.146876 0.162531 0.180868 0.202893 0.230365 0.266192 0.314644 0.338215 0.354481 0.356016 0.338877 0.238035 0.242891 0.250688 0.26355 0.283377 0.227156 0.238194 0.247552 0.252337 0.253174 0.0129829 0.0279253 0.038022 0.0447387 0.0493313 0.0527092 0.0554663 0.0579282 0.0602317 0.0624094
#
##_CATCHABILITY----------------------------------------p.112
#
#_Q_setup for fleets with cpue or survey data
#_1: fleet number
#_2: link type: (1=simple q; 1 parm; 2=mirror simple q; 1 mirrored parm; 3=q and power; 2 parm; 4=mirror with offset; 2 parm)
#_3: extra input for link; i.e. mirror fleet#_or dev index number
#_4: 0/1 to select extra sd parameter
#_5: 0/1 for biasadj or not when using informative prior
#_6: 0/1 to float
#_fleet_num link_type link_info extra_sd bias_adjust Q_float
1 1 0 0 0 1 #_q_setup_cpue_bt
2 1 0 0 0 1 #_q_setup_idx_qcs
3 1 0 0 0 1 #_q_setup_idx_wcvi
4 1 0 0 0 1 #_q_setup_idx_wchg
5 1 0 0 0 1 #_q_setup_idx_gig
-9999 0 0 0 0 0 #_end_fleet_q -- terminator for fleet/survey catchability
#
#_Q_parms(if_any); Qunits_are_ln(q)
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var use_dev dev_mnyr dev_mxyr dev_PH Block Block_Fxn label
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_cpue_bt
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_qcs
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_wcvi
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_wchg
-15 15 -3 -3 6 0 -1 0 0 0 0 0 0 0 #_q_parms_idx_gig
#_no timevary Q parameters
#
##_SELECTIVITY-AND-DISCARD-----------------------------p.116
#
#_size_selex_patterns
#_Pattern:_0; parm=0; selex=1.0 for all sizes
#_Pattern:_1; parm=2; logistic; with 95% width specification
#_Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#_Pattern:_15; parm=0; mirror another age or length selex
#_Pattern:_6; parm=2+special; non-parm len selex
#_Pattern:_43; parm=2+special+2; like 6; with 2 additional param for scaling (average over bin range)
#_Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
#_Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#_Pattern:_21; parm=2+special; non-parm len selex; read as pairs of size; then selex
#_Pattern:_22; parm=4; double_normal as in CASAL
#_Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#_Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL); using joiners
#_Pattern:_25; parm=3; exponential-logistic in size
#_Pattern:_27; parm=3+special; cubic spline
#_Pattern:_42; parm=2+special+3; // like 27; with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
0 0 0 0 #_s_pat_len_fishery
0 0 0 0 #_s_pat_len_qci
0 0 0 0 #_s_pat_len_wcvi
0 0 0 0 #_s_pat_len_wchg
0 0 0 0 #_s_pat_len_gig
#
#_age_selex_patterns
#_Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#_Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#_Pattern:_11; parm=2; selex=1.0 for specified min-max age
#_Pattern:_12; parm=2; age logistic
#_Pattern:_13; parm=8; age double logistic
#_Pattern:_14; parm=nages+1; age empirical
#_Pattern:_15; parm=0; mirror another age or length selex
#_Pattern:_16; parm=2; Coleraine - Gaussian
#_Pattern:_17; parm=nages+1; empirical as random walk N parameters to read can be overridden by setting special to non-zero
#_Pattern:_41; parm=2+nages+1; // like 17; with 2 additional param for scaling (average over bin range)
#_Pattern:_18; parm=8; double logistic - smooth transition
#_Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#_Pattern:_20; parm=6; double_normal;using joiners
#_Pattern:_26; parm=3; exponential-logistic in age
#_Pattern:_27; parm=3+special; cubic spline in age
#_Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#_Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
20 0 3 0 #_s_pat_len_fishery
20 0 3 0 #_s_pat_len_qci
20 0 3 0 #_s_pat_len_wcvi
20 0 3 0 #_s_pat_len_wchg
20 0 3 0 #_s_pat_len_gig
#
#_AGE_SELECTIVITY_PATTERN_20
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var use_dev dev_mnyr dev_mxyr dev_PH Block Blk_Fxn
#_1_TRAWL+_FISHERY AgeSelex
#_female_selectivity_pattern_20
1 40 10.7 10.7 2.14 6 3 0 0 0 0 0 0 0 #_s_p1_fishery -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_fishery -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 1.6 1.6 0.32 6 4 0 0 0 0 0 0 0 #_s_p3_fishery -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_fishery -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_fishery -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_fishery -- final: selectivity at last bin; as logistic between 0 and 1
#_male_selectivity_option_3
#_p1size (age) at which a dogleg occurs
-20 20 0 0 1 0 -4 0 0 0 0 0 0 0 #_s_p1_fishery_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_fishery_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_fishery_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_fishery_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_fishery_male -- male parameter 5 is the apical selectivity for males
#_2_QCS_SYNOPTIC AgeSelex
#_female_selectivity_pattern_20
1 40 15.6 15.6 3.12 6 3 0 0 0 0 0 0 0 #_s_p1_qcs -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_qcs -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 3.72 3.72 0.744 6 4 0 0 0 0 0 0 0 #_s_p3_qcs -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_qcs -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_qcs -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_qcs -- final: selectivity at last bin; as logistic between 0 and 1
#_male_selectivity_option_3
#_p1size (age) at which a dogleg occurs
-20 20 0 0 1 0 -4 0 0 0 0 0 0 0 #_s_p1_fishery_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_qcs_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_qcs_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_qcs_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_qcs_male -- male parameter 5 is the apical selectivity for males
#_3_WCVI_SYNOPTIC AgeSelex
#_female_selectivity_pattern_20
1 40 15.4 15.4 3.08 6 3 0 0 0 0 0 0 0 #_s_p1_wcvi -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_wcvi -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 3.44 3.44 0.688 6 4 0 0 0 0 0 0 0 #_s_p3_wcvi -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_wcvi -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_wcvi -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_wcvi -- final: selectivity at last bin; as logistic between 0 and 1
#_male_selectivity_option_3
#_p1size (age) at which a dogleg occurs
-20 20 0 0 1 0 -4 0 0 0 0 0 0 0 #_s_p1_fishery_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_wcvi_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_wcvi_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_wcvi_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_wcvi_male -- male parameter 5 is the apical selectivity for males
#_4_WCHG_SYNOPTIC AgeSelex
#_female_selectivity_pattern_20
1 40 10.8 10.8 2.16 6 3 0 0 0 0 0 0 0 #_s_p1_wchg -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_wchg -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 2.08 2.08 0.416 6 4 0 0 0 0 0 0 0 #_s_p3_wchg -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_wchg -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_wchg -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_wchg -- final: selectivity at last bin; as logistic between 0 and 1
#_male_selectivity_option_3
#_p1size (age) at which a dogleg occurs
-20 20 0 0 1 0 -4 0 0 0 0 0 0 0 #_s_p1_fishery_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_wchg_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_wchg_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_wchg_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_wchg_male -- male parameter 5 is the apical selectivity for males
#_5_GIG_HISTORICAL AgeSelex
#_female_selectivity_pattern_20
1 40 17.4 17.4 3.48 6 3 0 0 0 0 0 0 0 #_s_p1_gig -- peak: beginning size (or age) for the plateau (=mu in Awatea)
0 60 0 0 99 0 -3 0 0 0 0 0 0 0 #_s_p2_gig -- top: width of plateau; as logistic between peak and maximum length (or age)
-15 15 4.6 4.6 0.92 6 4 0 0 0 0 0 0 0 #_s_p3_gig -- ascending width: parameter value is ln(width) (=log varL in Awatea)
-15 500 200 0 1 0 -4 0 0 0 0 0 0 0 #_s_p4_gig -- descending width: parameter value is ln(width) (=log varR in Awatea)
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p5_gig -- initial: selectivity at first bin; as logistic between 0 and 1
-10 10 -999 0 99 0 -50 0 0 0 0 0 0 0 #_s_p6_gig -- final: selectivity at last bin; as logistic between 0 and 1
#_male_selectivity_option_3
#_p1size (age) at which a dogleg occurs
-20 20 0 0 1 0 -4 0 0 0 0 0 0 0 #_s_p1_fishery_male -- size (age) at which a dogleg (inflection point?) occurs (=Delta in Awatea)
-15 15 0 0 2.5 0 -4 0 0 0 0 0 0 0 #_s_p2_gig_male -- added to the third selectivity parameter (width of ascending side)
-15 500 0 0 0.6 0 -4 0 0 0 0 0 0 0 #_s_p3_gig_male -- added to the fourth selectivity parameter (width of descending side)
-10 10 0 0 99 0 -50 0 0 0 0 0 0 0 #_s_p4_gig_male -- added to the sixth selectivity parameter (selectivity at final size bin)
0 1 1 1 0 0 -50 0 0 0 0 0 0 0 #_s_p5_gig_male -- male parameter 5 is the apical selectivity for males
#_Dirichlet-Multinomial parameters controlling age-comp weights (also choose Dirichlet error structure in data file)
#_-5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w1_ln(dm_parm)_fishery
#_-5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w2_ln(dm_parm)_qcs
#_-5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w3_ln(dm_parm)_wcvi
#_-5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w4_ln(dm_parm)_wchg
#_-5 10 0 0 1.813 6 2 0 0 0 0 0 0 0 #_s_w5_ln(dm_parm)_gig
#_LO HI INIT PRIOR PR_SD PR_type PHASE env-var use_dev dev_mnyr dev_mxyr dev_PH Block Blk_Fxn
#_no timevary selex parameters
#
0 #_s_use_2d_ar1 selectivity(0/1)
#_no 2D_AR1 selex offset used
#
##_TAG_RECAPTURE_PARAMETERS----------------------------p.138
#
#_Tag loss and Tag reporting parameters go next
0 #_tag_custom -- TG_custom: 0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0 #_placeholder if no parameters
#_no timevary parameters
#
##_VARIANCE_ADJUSTMENT_FACTORS-------------------------p.140
#
#_Input variance adjustments factors:
#_1=add_to_survey_CV
#_2=add_to_discard_stddev
#_3=add_to_bodywt_CV
#_4=mult_by_lencomp_N
#_5=mult_by_agecomp_N
#_6=mult_by_size-at-age_N
#_7=mult_by_generalized_sizecomp
#_Note: adjustment factor 5 reportedly not needed if Dirichlet-Multinomial likelihood used (see hake example; Chris Grandin)
#_factor fleet value desc
1 1 0 #_vadj_cv_bt_cpue
1 2 0 #_vadj_cv_qcs_idx
1 3 0 #_vadj_cv_wcvi_idx
1 4 0 #_vadj_cv_wchg_idx
1 5 0 #_vadj_cv_gig_idx
#_SS_tune_comps(rlist; option="Francis") on Multinomial
5 1 1 #_vadj_af_mn_fishery
5 2 1 #_vadj_af_mn_qcs
5 3 1 #_vadj_af_mn_wcvi
5 4 1 #_vadj_af_mn_wchg
5 5 1 #_vadj_af_mn_gig
-9999 0 0 #_end_vadj -- terminator_variance_adjustment_factors
#
##_LAMBDAS_(EMPHASIS_FACTORS)--------------------------p.142
#
1 #_lambda_max_phase -- maximum ambda phase
1 #_lambda_sd_off -- SD offset; must be 1 if any growthCV; sigmaR; or survey extraSD is an estimated parameter
#_read 3 changes to default Lambdas (default value is 1.0)
#_Like_comp codes: 1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch;
#_10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#_like_comp fleet phase value sizefreq_method
-9999 1 1 1 1 #_lambda_max_end -- terminator for max lambda phase
#
#_lambdas (for info only; columns are phases)
#_0 0 0 0 #_CPUE/survey:_1
#_1 1 1 1 #_CPUE/survey:_2
#_1 1 1 1 #_CPUE/survey:_3
#_1 1 1 1 #_lencomp:_1
#_1 1 1 1 #_lencomp:_2
#_0 0 0 0 #_lencomp:_3
#_1 1 1 1 #_agecomp:_1
#_1 1 1 1 #_agecomp:_2
#_0 0 0 0 #_agecomp:_3
#_1 1 1 1 #_size-age:_1
#_1 1 1 1 #_size-age:_2
#_0 0 0 0 #_size-age:_3
#_1 1 1 1 #_init_equ_catch
#_1 1 1 1 #_recruitments
#_1 1 1 1 #_parameter-priors
#_1 1 1 1 #_parameter-dev-vectors
#_1 1 1 1 #_crashPenLambda
#_0 0 0 0 #_f_ballpark_lambda
#
##_CONTROLS_FOR VARIANCE_OF_DERIVED_QUANTITIES---------p.144
#
0 #_read_sd_rep -- (0/1) read specs for more stddev reporting
#_0 0 0 0 0 0 0 0 0 #_placeholder for #_selex_fleet; 1=len/2=age/3=both; year; N selex bins; 0 or Growth pattern; N growth ages; 0 or NatAge_area(-1 for all); NatAge_yr; N Natages
#_placeholder for vector of selex bins to be reported
#_placeholder for vector of growth ages to be reported
#_placeholder for vector of NatAges ages to be reported
999 #_eof -- EOF for control
