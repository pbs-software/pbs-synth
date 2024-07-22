Yellowtail Rockfish (YTR) 2024 Stock Assessment
===============================================
Model platform : Stock Synthesis 3
Version        : V3.30.22.01; 2022-09-30; Stock Synthesis 3 by Richard Methot (NOAA) using ADMB-13.1 safe libraries was compiled using GNU C++ 12.2.0 64bit on Apr 27 2023

Base run (use Francis mean-age to reweight AF data, use Multinomial):
---------------------------------------------------------------------
	* Run02v1 -- BC coastwide (240320);
		+ 7 fleets (1=TRAWL_FISHERY_BC, 2=QCS_SYNOPTIC, 3=WCVI_SYNOPTIC, 4=WCHG_SYNOPTIC, 5=HS_SYNOPTIC, 6=GIG_HISTORICAL, 7=NMFS_TRIENNIAL)
		+ apply Francis weights: 2.4746, 1.9031, 2.4378, 1.5631;
		+ MCMCs run on Linux server (240521): v1c, nsamp=20000, nchain=8, burnin=2500/chain, thin=10, extra=FALSE.

Sensitivity runs (based on R02v1)
---------------------------------MCMC
	* S01 Run10v1 (240412) -- split M between ages 9 and 10 (~age 50% maturity) 
		+ used prior N(0.08,0.024) for males and young females;
		+ used prior N(0.12,0.036) for mature females;
		+ widened bounds on M to (0.001, 0.2) for males and young females;
		+ widened bounds on M to (0.001, 0.3) for mature females.
		+ applied Francis mean-age weights: w = 2.4173, 1.9256, 2.4312, 1.6455
		+ ran MCMCs on Linux server (240515): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S02 Run11v1 (240415) -- estimate dome-shaped selectivity for females;
		+ require estimating male selectivity parameters with females as offsets
		+ estimate $\beta_{1,3,4,6}$ and $\Delta_{1,2,3,4}$
		+ applied Francis mean-age weights: w = 2.3506, 2.0577, 3.8136, 1.8349
		+ ran MCMCs on Linux server (240516): v1b, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* s03 Run05v1 (240328)} -- change sigmaR from 0.9 to 0.6
		+ applied Francis mean-age weights: w = 2.3975, 1.9056, 2.3011, 1.6077
		run MCMCs on Linux server (240522): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S04 Run06v1 (240402) -- change sigmaR from 0.9 to 1.2
		+ applied Francis mean-age weights: w = 2.5031, 1.9054, 2.5014, 1.5416
		+ ran MCMCs on Linux server (240522): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S05 Run07v1 (240402) -- estimate sigmaR using N(0.9,0.9) in phase 1
		+ applied Francis mean-age weights: w = 1.0961, 3.6777, 7.0888, 1.1594
		+ ran MCMCs on Linux server (240523): v1b, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
		+ (initial NUTS took >24h, updated NUTS $\sim$5h)
	* S06 Run08v1 (240402) -- use Dirichlet-Mulitinomial for AFs (no reweight)
		+ ran MCMCs on Linux server (240523): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S07 Run12v1 (240521) -- replace ageing error vector AE3 with AE1 (no ageing error)
		+ applied Francis mean-age weights: w = 2.5527, 1.8759, 2.3043, 1.4690
		+ ran MCMCs on Linux server (240524): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S08 Run13v1 (240521) -- replace ageing error vector AE3 with AE5 (smoothed age-precision error)
		+ applied Francis mean-age weights: w = 2.4513, 1.8822, 2.4134, 1.5304
		+ ran MCMCs on Linux server (240524): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S09 Run14v1 (240521) -- replace ageing error vector AE3 with AE6 (CASAL-style 10% CV)
		+ applied Francis mean-age weights: w = 2.4154, 1.9343, 2.4874, 1.5277
		+ ran MCMCs on Linux server (240524): \textbf{v1a} nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S10 Run15v1 (240522) -- reduce catch in years 1965-1995 by 30%
		+ applied Francis mean-age weights: w = 2.4214, 1.9291, 2.4286, 1.6190
		+ ran MCMCs on Linux server (240527): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S11 Run16v1 (240522) -- increase catch in years 1965-1995 by 50%
		+ applied Francis mean-age weights: w = 2.5533, 1.9003, 2.5046, 1.5538
		+ ran MCMCs on Linux server (240527): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S12 Run09v3 (240511) -- use geostatistical abundance indices for the four synoptic surveys
		+ applied Francis mean-age weights: w = 2.4980, 2.1765, 2.6342, 1.4475
		+ ran MCMCs on Linux server (240527): v3a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S13 Run04v2 (240521) -- add the HBLL abundance index series (North: SSID=22, South: SSID=36).
		+ applied Francis mean-age weights: w = 2.5616, 1.9902, 2.6143, 1.5814 
		+ ran MCMCs on Linux server (240527): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
	* S14 Run17v1 (240620) -- split trawl fleet into BT and MW trawl fleets
		+ applied Francis mean-age weights: w = 2.1869, 2.8211, 1.8149, 2.4053, 1.4768
		+ ran MCMCs on Linux server (240620): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
