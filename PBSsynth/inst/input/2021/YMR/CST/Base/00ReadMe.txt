Yellowmouth Rockfish (YMR) 2021 Stock Assessment
================================================
Model platform : Stock Synthesis 3
Version        : V3.30.16.00; 2020-09-03 safe; Stock Synthesis by Richard Methot (NOAA) using ADMB 12.2
Base run (reweighted once):
	* BC coastwide with one area (3CD5ABCDE) and
		5 fleets (1=TRAWL+_FISHERY, 2=QCS_SYNOPTIC, 3=WCVI_SYNOPTIC, 4=WCHG_SYNOPTIC, 5=GIG_HISTORICAL).
	* fixed natural mortality M to five levels: 0.04, 0.045, 0.05, 0.055, and 0.06 for a total of five reference models using one axis of uncertainty (M):
		+ B1: Run77 (M=0.04)
		+ B2: Run71 (M=0.045)
		+ B3: Run75 (M=0.05) [central run]
		+ B4: Run72 (M=0.055)
		+ B5: Run76 (M=0.06)
Assumptions:
	* assumed two sexes (females, males);
	* set plus age class A to 60 years;
	* assumed one commercial fishery dominated by trawl (bottom + midwater), 
		with minor removals by halibut longline, sablefish trap, lingcod longline,
		inshore longline, and salmon troll, pooled into a single catch series
		with associated age frequency (AF) data drawn from the trawl fishery;
	* used one commercial bottom trawl fishery abundance index series (bottom trawl CPUE index, 1996-2020);
	* used four survey abundance index series (QCS Synoptic, WCVI Synoptic, WCHG Synoptic, and GIG Historical), with age frequency (AF) data;
	* assumed a wide (weak) normal prior N (8, 8) on log R0 to help stabilise the model;
	* used informed normal priors for the two selectivity parameters (mu_g, v_gL) for
		all fleets (fishery and surveys), and set the male selectivity offset (Delta_g) to 0;
	* estimated recruitment deviations from 1950 to 2012;
	* applied abundance reweighting: added CV process error to index CVs, cp=0.3296
		for the commercial CPUE series and c_p=0 for the surveys;
	* applied composition reweighting: adjusted AF effective sample sizes using a
		harmonic mean ratio method based on McAllister and Ianelli (1997);
	* fixed the standard deviation of recruitment residuals (sigmaR) to 0.9;
	* used an ageing error vector based on the CV of observed lengths at age,
		described in Appendix D, Section D.2.3.
