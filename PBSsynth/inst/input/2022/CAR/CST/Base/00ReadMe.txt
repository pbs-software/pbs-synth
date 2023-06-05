Canary Rockfish (CAR) 2022 Stock Assessment
===========================================
Model platform : Stock Synthesis 3
Version        : # V3.30.18; 2021-xx-xx safe; Stock Synthesis by Richard Methot (NOAA) using ADMB 12.3
Base run:
	* Run 24 -- BC coastwide with one area (3CD5ABCDE) and
		8 fleets (1=TRAWL_FISHERY, 2=OTHER_FISHERY, 3=QCS_SYNOPTIC, 4=WCVI_SYNOPTIC, 5=NMFS_TRIENNIAL, 6=HS_SYNOPTIC, 7=WCHG_SYNOPTIC, 8=GIG_HISTORICAL).
Assumptions:
	* assumed two sexes (females, males);
	* estimated a single mortality M per sex to represent all ages;
	* set plus-age class A to 60 years;
	* assumed two commercial fisheries: 'Trawl' (predominant with ~97% of catch) and 'Other';
		+ Trawl fishery comprised bottom and midwater trawl gear;
		+ Other fishery included non-trawl gear (halibut longline, sablefish trap/longline, dogfish/lingcod troll, and hook & line rockfish);
		+ age frequency (AF) data were only available from the Trawl fishery;
	* used one commercial bottom trawl fishery abundance index series (bottom trawl CPUE index, 1996-2021);
	* used six survey abundance index series (QCS Synoptic, WCVI Synoptic, NMFS Triennial, HS Synoptic,
		WCHG Synoptic, and GIG Historical), with age frequency (AF) data for the first three surveys;
	* assumed a wide (weak) normal prior N (7, 7) on log R0 to help stabilise the model;
	* used informed normal priors for the three primary selectivity parameters (mu_g, v_gL, Delta_g, see Appendix E)
		for all fleets (fishery and surveys) derived from Table J.7 in Stanley et al. (2009);
	* estimated recruitment deviations from 1950 to 2012;
	* applied abundance reweighting: added CV process error to index CVs, c_p=0.178
		for the commercial CPUE series and c_p=0 for the surveys (see Appendix E);
	* used SS3's Dirichlet-Multinomial error distribution to fit AF data instead of applying composition reweighting;
	* fixed the standard deviation of recruitment residuals (sigmaR) to 0.9;
	* used an ageing error vector based on the CV of observed lengths at age, described in Appendix D, Section D.2.3.
