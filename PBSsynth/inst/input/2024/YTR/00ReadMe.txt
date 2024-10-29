Yellowtail Rockfish (YTR) 2024 Stock Assessment
===============================================
Model platform : Stock Synthesis 3
Version        : V3.30.22.01; 2022-09-30; Stock Synthesis 3 by Richard Methot (NOAA) using ADMB-13.1 safe libraries was compiled using GNU C++ 12.2.0 64bit on Apr 27 2023

Base run (use Francis mean-age to reweight AF data, use Multinomial):
---------------------------------------------------------------------
	* B1 R02v2 -- BC coastwide (241009);
		+ 7 fleets (1=TRAWL_FISHERY_BC, 2=QCS_SYNOPTIC, 3=WCVI_SYNOPTIC, 4=WCHG_SYNOPTIC, 5=HS_SYNOPTIC, 6=GIG_HISTORICAL, 7=NMFS_TRIENNIAL)
		+ apply Francis weights: 2.4770, 1.9064, 2.4379, 1.5656
		+ MCMCs run on Linux server (241009): v2a, nsamp=20000, nchain=8, burnin=2500/chain, thin=10, extra=FALSE.

Sensitivity runs (based on R02v1)
---------------------------------MCMC
	* S01 R10v2 (241010) -- split M between ages 9 and 10 (~age 50% maturity) 
		+ used prior N(0.08,0.024) for males and young females;
		+ used prior N(0.12,0.036) for mature females;
		+ widened bounds on M to (0.001, 0.2) for males and young females;
		+ widened bounds on M to (0.001, 0.3) for mature females.
		+ applied Francis mean-age weights: w = 2.4185, 1.9283, 2.4302, 1.6476
		+ ran MCMCs on Linux server (241011): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S02 R11v2 (241010) -- estimate dome-shaped selectivity for females;
		+ require estimating male selectivity parameters with females as offsets
		+ estimate $\beta_{1,3,4,6}$ and $\Delta_{1,2,3,4}$
		+ applied Francis mean-age weights: w = 2.3518, 2.0599, 3.8144, 1.8358
		+ ran MCMCs on Linux server (241011): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* s03 R05v2 (241010) -- change sigmaR from 0.9 to 0.6
		+ applied Francis mean-age weights: w = 2.4005, 1.9100, 2.3001, 1.6107
		run MCMCs on Linux server (241011): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S04 R06v2 (241010) -- change sigmaR from 0.9 to 1.2
		+ applied Francis mean-age weights: w = 2.5053, 1.9080, 2.5022, 1.5436
		+ ran MCMCs on Linux server (241011): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S05 R07v3 (241010) -- estimate sigmaR using N(0.9,0.9) in phase 1
		+ applied Francis mean-age weights: w = 1.0896, 3.6942, 7.1840, 1.1531
		+ ran MCMCs on Linux server (241012): v3a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE
		+ (initial NUTS took >24h, updated NUTS $\sim$5h)

	* S06 R08v2 (241010) -- use Dirichlet-Mulitinomial for AFs (no reweight)
		+ ran MCMCs on Linux server (241011): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S07 R12v2 (241010) -- replace ageing error vector AE3 with AE1 (no ageing error)
		+ applied Francis mean-age weights: w = 2.5546, 1.8792, 2.3041, 1.4707
		+ ran MCMCs on Linux server (241012): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S08 R13v2 (241010) -- replace ageing error vector AE3 with AE5 (smoothed age-precision error)
		+ applied Francis mean-age weights: w = 2.4559, 1.8872, 2.4153, 1.5423
		+ ran MCMCs on Linux server (241012): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S09 R14v2 (241010) -- replace ageing error vector AE3 with AE6 (CASAL-style 10% CV)
		+ applied Francis mean-age weights: w = 2.4187, 1.9377, 2.4873, 1.5301
		+ ran MCMCs on Linux server (241012): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S10 R15v2 (241010) -- reduce catch in years 1965-1995 by 30%
		+ applied Francis mean-age weights: w = 2.4237, 1.9327, 2.4296, 1.6221
		+ ran MCMCs on Linux server (241012): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S11 R16v2 (241010) -- increase catch in years 1965-1995 by 50%
		+ applied Francis mean-age weights: w = 2.5564, 1.9040, 2.5025, 1.5555
		+ ran MCMCs on Linux server (241012): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S12 R09v4 (241010) -- use geostatistical abundance indices for the four synoptic surveys
		+ applied Francis mean-age weights: w = 2.5007, 2.1805, 2.6335, 1.4501
		+ ran MCMCs on Linux server (241012): v4a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S13 R04v3 (241010) -- add the HBLL abundance index series (North: SSID=22, South: SSID=36).
		+ applied Francis mean-age weights: w = 2.5643, 1.9941, 2.6142, 1.5839
		+ ran MCMCs on Linux server (241012): v3a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S14 R17v2 (241010) -- split trawl fleet into BT and MW trawl fleets
		+ applied Francis mean-age weights: w = 2.1902, 2.8195, 1.8184, 2.4056, 1.4792
		+ ran MCMCs on Linux server (241012): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

Sensitivity runs (based on R02v1)
---------------------------------MPD
	* S15 R18v2 (241010) -- use different priors on M (with higher means)
		+ applied Francis mean-age weights: w = 2.3551, 1.9028, 2.4194, 1.5926
		+ ran MCMCs on Linux server (241013): v2a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

	* S16 R19v1 (240906) -- use flavidus fecundity parameters from Dick et al. (2017) : exp(alpha)=2.875e-07, beta=4.590
		+ applied Francis mean-age weights: w = 2.4783, 1.9081, 2.4380, 1.5669
		+ ran MCMCs on Linux server (241012): v1a, nsamp=10,000, nchain=8, burnin=1250/chain, thin=5, extra=FALSE

