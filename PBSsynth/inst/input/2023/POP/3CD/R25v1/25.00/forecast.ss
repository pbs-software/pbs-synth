# Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States of America.
# Foreign copyrights may apply. See copyright.txt for more information.
# user support available at: NMFS.Stock.Synthesis@noaa.gov
# user info available at: https://vlab.ncep.noaa.gov/group/stock-synthesis
# V3.30.20; 2022-09-30; Stock Synthesis 3 by Richard Methot (NOAA) using ADMB-13.0 safe libraries, compiled using GNU C++ 10.3.0 64bit on Aug 10 2022.
# C generic forecast file
# for all year entries except rebuilder; enter either: actual year; -999 for styr; 0 for endyr; neg number for rel. endyr
2 #_bmark_rfpt -- Benchmarks/reference_points: 0=skip; 1=calc F_spr; F_btgt; F_msy; 2=calc F_spr; F0.1; F_msy
2 #_msy_methods -- 1= set to F(SPR); 2=calc F(MSY); 3=set to F(Btgt) or F0.1; 4=set to F(endyr)
0.45 #_spr_target -- spawner per recruit target(e.g. 0.40)
0.4 #_ssb_target -- spawning biomass target (e.g. 0.40)
-999 0 -999 0 -999 0 -999 0 -999 0 #_bmark_yrs -- benchmark years: start and end for B; S; F; R; SPR (enter actual year; or values of 0 or -integer to be rel. endyr)
1 #_bmark_f_basis -- 1 = use year range; 2 = set relF same as forecast below
2 #_fcast -- Forecast (required but ignored if no benchmarks): -1=none; 0=simple_1yr; 1=F(SPR); 2=F(MSY) 3=F(Btgt) or F0.1; 4=Ave F (uses first-last relF yrs); 5=input annual F scalar
# where none and simple require no input after this line; simple sets forecast F same as end year F
11 #_fcast_nyrs -- Number of forecast years (must be >=1)
1 #_scalar_f -- (only used for Do_Forecast==5) such that apical_F(f)=Fmult*relF(f)
-999 0 -999 0 1935 2014 #_fcast_yrs -- forecast years: start and end for S; F; R (enter actual year; or values of 0 or -integer to be rel. endyr)
0 #_fcast_sel_opt -- Forecast selectivity options (0=Fcast selex is mean from year range; 1=Fcast selectivity from annual time-vary parms)
0 #_ctl_rule -- method: 0=none (additional control rule inputs will be ignored), 1=catch as function of SSB, buffer on F; 2=F as function of SSB, buffer on F; 3=catch as function of SSB, buffer on catch (US W Coast GF approach); 4=F is a function of SSB, buffer on catch.
0.002 #_ctl_rule_inflect -- Relative biomass level to unfished biomass above which F is constant at control rule F. If set to -1 the ratio of BMSY to the unfished spawning biomass will automatically be used.
0.001 #_ctl_rule_cutoff -- Relative biomass level to unfished biomass below which F is set to 0 (management threshold).
0.75 #_ctl_rule_buff -- (buffer) enter Control rule target as fraction of Flimit (e.g. 0.75); negative value invokes list of [year; scalar] with filling from year to YrMax
# COND_ctl_buff -1 -- Conditional input for annual control rule buffer (see manual)
3 #_fcast_nloops -- number of forecast loops (1=OFL only; 2=ABC; 3=get F from forecast ABC catch with allocations applied)
3 #_fcast_loop_rec1 -- First forecast loop with stochastic recruitment
1 #_fcast_rec -- Forecast recruitment: 0= spawn_recr; 1=value*spawn_recr_fxn; 2=value*VirginRecr; 3=recent mean from yr range above (need to set phase to -1 in control file to get constant recruitment in MCMC)
1 #_scalar_r -- Scalar (value multiplier for SSR, options 1 or 2)  or number of years of recent main recruitments to average (option 3)
0 #_fcast_loop_ctl5 -- Forecast loop control #5 (reserved for future bells&whistles)
2200 #_fcast_yr1 -- First year for caps and allocations (should be after years with fixed inputs)
0 #_impl_error -- stddev of log(realized catch/target catch) in forecast (set value>0.0 to cause active impl_error)
0 #_rebuild_wcgf -- perform West Coast gfish rebuilder output (0/1)
1996 #_rebuild_catch_yr -- Rebuilder: first year catch could have been set to zero (Ydecl)(-1 to set to 1999)
-1 #_rebuild_start_yr -- Rebuilder: year for current age structure (Yinit) (-1 to set to endyear+1)
1 #_fleet_rel_f -- relative F: 1=use first-last alloc year; 2=read seas; fleet; alloc list below
# Note that fleet allocation is used directly as average F if Do_Forecast=4
2 #_fcast_catch_basis -- basis for fcast catch tuning and for fcast catch caps and allocation (2=deadbio; 3=retainbio; 5=deadnum; 6=retainnum); NOTE: same units for all fleets
# 2 #_COND_rel_f -- Conditional input if relative F choice = 2
# #_enter list of: season; fleet; relF; if used; terminate with season=-9999
# 1 1 1
# -9999 0 0 #_end_rel_f -- terminator for list of relF
-9999 -1 #_end_max_tot_fleet_cat -- enter list of: fleet number; max annual catch for fleets with a max; terminate with fleet=-9999
-9999 -1 #_end_max_tot_area_cat -- enter list of area ID and max annual catch; terminate with area=-9999
-9999 -1 #_end_fleet_assign_alloc -- enter list of fleet number and allocation group assignment; if any; terminate with fleet=-9999
# COND_alloc_grps -- if N allocation groups >0; list year; allocation fraction for each group
# list sequentially because read values fill to end of N forecast
# terminate with -9999 in year field
# no allocation groups
2 #_fcast_cat_basis -- basis for input forecast catch: -1=read basis with each obs; 2=dead catch; 3=retained catch; 99=input apical_F; NOTE: bio vs num based on fleet's catchunits
# COND_==-1 -- Forecasted catches: enter one line per number of fixed forecast year catch
# COND_>0 -- Forecasted catches: enter one line per number of fixed forecast year catch
# enter list of Fcast catches or Fa; terminate with line having year=-9999
# RH: First forecast is current year so need N+1 years of forecast to get N forecast years
# Before 5C adjustment: 5-yr average (2018-2022) for POP: 5ABC=1442, 3CD=840, 5DE=1024
# After 5C adjustment: 5-yr average (2018-2022) for POP: 5ABC=1618, 3CD=840, 5DE=848
#_Yr Seas Fleet Catch
#2024 1 1 1618
#2025 1 1 1618
#2026 1 1 1618
#2027 1 1 1618
#2028 1 1 1618
#2029 1 1 1618
#2030 1 1 1618
#2031 1 1 1618
#2032 1 1 1618
#2033 1 1 1618
#2034 1 1 1618
2024 1 1 840
2025 1 1 840
2026 1 1 840
2027 1 1 840
2028 1 1 840
2029 1 1 840
2030 1 1 840
2031 1 1 840
2032 1 1 840
2033 1 1 840
2034 1 1 840
#2024 1 3 848
#2025 1 3 848
#2026 1 3 848
#2027 1 3 848
#2028 1 3 848
#2029 1 3 848
#2030 1 3 848
#2031 1 3 848
#2032 1 3 848
#2033 1 3 848
#2034 1 3 848
-9999 0 0 0 #_end_catpols -- year, season, fleet, catch (or_F)
999 #_eof -- verify end of input for forecast file
