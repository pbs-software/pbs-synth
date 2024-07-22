Pacific Ocean Perch (POP) 2023 Stock Assessment
===============================================
Model platform : Stock Synthesis 3
Version        : V3.30.20; 2022-09-30; Stock Synthesis 3 by Richard Methot (NOAA) using ADMB-13.0 safe libraries, compiled using GNU C++ 10.3.0 64bit on Aug 10 2022.

Base runs (use Francis mean-age to reweight AF data, use Multinomial):
----------------------------------------------------------------------
	* Run 21v3 -- BC coastwide (5ABC, 3CD, 5DE), multi-area model (Sep 8, 2023);
		+ 9 fleets (1=TRAWL_FISHERY_5ABC, 2=TRAWL_FISHERY_3CD, 3=TRAWL_FISHERY_5DE, 4=QCS_SYNOPTIC, 5=WCVI_SYNOPTIC, 6=WCHG_SYNOPTIC, 7=GIG_HISTORICAL, 8=NMFS_TRIENNIAL, 9=WCVI_HISTORICAL)
		+ apply Francis weights: 2.7923, 3.1274, 2.8821, 0.5744, 1.0766, 1.2071, 0.7348, 0.4874;
		+ MCMCs run on Linux server (230909): v.3a, nsamp=20000, nchain=8, burnin=2500/chain, thin=10, extra=TRUE.
	* Run 24v1 -- 5ABC (area 1) single-area model (Sep 14, 2023);
		+ 3 fleets (1=TRAWL_FISHERY_5ABC, 2=QCS_SYNOPTIC, 3=GIG_HISTORICAL)
		+ apply Francis weights: 2.7682, 0.5926, 0.7486;
		+ MCMCs run on Linux server (230915): v.1a, nsamp=20000, nchain=8, burnin=2500/chain, thin=10, extra=TRUE.
	* Run 25v1 -- 3CD  (area 2) single-area model (Sep 14, 2023);
		+ 4 fleets (1=TRAWL_FISHERY_3CD, 2=WCVI_SYNOPTIC, 3==NMFS_TRIENNIAL, 4=WCVI_HISTORICAL)
		+ apply Francis weights: 3.3348, 1.0405, 0.4856;
		+ MCMCs run on Linux server (230915): v.1a, nsamp=20000, nchain=8, burnin=2500/chain, thin=10, extra=TRUE.
	* Run 26v1 -- 5DE  (area 3) single-area model (Sep 14, 2023);
		+ 2 fleets (1=TRAWL_FISHERY_5DE, 2=WCHG_SYNOPTIC)
		+ apply Francis weights: 3.0688, 1.3401;
		+ MCMCs run on Linux server (230915): v.1a, nsamp=20000, nchain=8, burnin=2500/chain, thin=10, extra=TRUE.

Sensitivity runs (based on R21v3)
---------------------------------MCMC
	* Run17v18 -- use Dirichlet-multinomial parameterisation:
		+ no Francis reweights;
		+ MCMCs run on Linux server (230910): v.18a, nsamp=20000, nchain=8, burnin=2500/chain, thin=10, extra=TRUE.
	* Run27v1 -- fix parameter Rdist for 5ABC to 0:
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.7823, 2.8671, 3.0780, 0.5757, 0.9667, 1.2445, 0.7059, 0.4634 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230916): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
	* Run28v1 -- fix parameter Rdist for 3CD to 0:
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.7693, 2.9989, 3.1270, 0.5861, 1.0491, 1.3352, 0.7262, 0.4762 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230916): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
	* Run29v1 -- apply no ageing error:
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.6171, 3.0860, 2.7799, 0.5820, 1.0435, 1.2434, 0.6843, 0.4868 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230916): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
	* Run30v1 -- use smoothed ageing error from age-reader CVs:
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.7593, 3.0975, 2.8381, 0.5741, 1.0673, 1.1788, 0.6984, 0.4863 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230916): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
	* Run31v1 -- use constant-CV ageing error:
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.9101, 3.0522, 2.7957, 0.5767, 1.0773, 1.2077, 0.7569, 0.4860 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230916): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
	* Run32v1 -- reduce commercial catch (1965-95) by 30%:
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.8833, 3.0629, 2.8249, 0.5660, 1.0609, 1.0910, 0.7311, 0.4902 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230916): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
	* Run33v1 -- increase commercial catch (1965-95) by 50%:
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.5216, 3.1830, 2.8400, 0.5822, 1.1006, 1.3723, 0.7304, 0.4862 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230916): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
	* Run34v1 -- reduce sigmaR to 0.6 (from 0.9):
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.7267, 3.2302, 2.9101, 0.5741, 1.0934, 1.2708, 0.7125, 0.4873 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230916): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
	* Run35v1 -- increase sigmaR to 1.2 (from 0.9):
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.8098, 3.0808, 2.8794, 0.5733, 1.0676, 1.1977, 0.7414, 0.4876 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS);
		+ MCMCs run on Linux server (230920): v.1a, nsamp=10000, nchain=8, burnin=1250/chain, thin=5, extra=TRUE.
---------------------------------MPD
	* Run22v2 -- use separate midwater and bottom trawl fleets for 3CD and 5ABC:
		+ separated midwater and bottom trawl catches in 3CD and 5ABC
		+ estimated a shared selectivity based on combined AF data (3CD+5ABC);
		+ 5ABC trawl, 3CD trawl, 5DE trawl, QCS synoptic, WCVI synoptic, WCHG synoptic, GIG historical, NMFS triennial, 3CD midwater
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.5089, 4.7872, 2.8770, 0.5705, 1.0907, 1.2134, 0.7399, 0.5203, 2.7465 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS, MW).
	* Run36v2 -- add Hecate Strait synoptic survey to 5DE data:
		+ link HS selectivity to neighbouring WCHG Synoptic (fleet 6);
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.7889, 3.1266, 2.9207, 0.5744, 1.0726, 1.2861, 0.7342, 0.4862 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS).
	* Run37v1 -- use empirical proportions mature rather than MLE fitted ogive:
		+ requested by Noel Cadigan at RPR;
		+ use Multinomial distribution for estimating AFs;
		+ apply Francis weights: 2.7922, 3.1272, 2.882, 0.5745, 1.0765, 1.207, 0.7348, 0.4875 (5ABC, 3CD, 5DE, QCS, WCVI, WCHG, GIG, NMFS).
