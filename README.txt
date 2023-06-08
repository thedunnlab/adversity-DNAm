README
Documentation for the paper " The timing of childhood adversity associates with epigenetic patterns across childhood and adolescence: Results from a prospective, longitudinal study" published in Lancet Child and Adolescent Health (2023)
Author: 	Alexandre Lussier, PhD
Last updated:	June 8, 2023

In this folder, you can find the following information regarding our manuscript. 

1. Manuscripts
	- Preprint of the manuscript by Lussier and colleagues described here, and supplemental materials (In press at Lancet Child and Adolescent Health).
	- Manuscript by Zhu and colleagues, which describes the application of the SLCMA to epigenome-wide settings. <https://doi.org/10.1093/aje/kwaa246>
	- Manuscript by Dunn and colleagues, which describes the application of the SLCMA to analyze the time-dependent relationship between childhood adversity and DNA methylation. <https://doi.org/10.1016/j.biopsych.2018.12.023>
	- Manuscript by Lussier and colleagues, which describes the findings between childhood adversity and DNA methylation at age 7. <https://doi.org/10.1080/15592294.2022.2028072> and <https://doi.org/10.1016/j.bpsgos.2022.04.002>

2. SLCMA scripts 
R scripts used to analyze each adversity for each of the main analyses using the Structured Life Course Modeling Approach (SLCMA)
	- 1. LARS-noimpute-function-20190516 - SLCMA function used for these analyses
	- 2. runSLCMA.adversity.noever.20200513.R - primary script, used for abuse, mompsych, oneadult, and nbhqual
	- 3. 2021-04-24_r_faminstability_slcma_noever.R - same as primary, but using r_faminstability 
	- 4. 2021-06-01_fscore_slcma_noever.R - same as primary, but using revised Fscore variable
	- 5. 2022-07-13_parcruelty_revised_SLCMA_noever_SI.R - same as primary, but using revised parcruelty variable

3. Analysis scripts
Scripts to analyze the results from the SLCMA and generate the figures
	- 1. Results_structure_2022-08-18.Rmd - concatenate results from the SLCMA scripts 
	- 2. Adversity_summary_2022-08-18.Rmd - summary of childhood adversity measures (Fig 1)
	- 3. Top_hits_biology_2022-08-18.Rmd - analysis of main age 15 findings 
	- 4. Top_hits_stability_2022-12-05.Rmd - stability of DNAm changes from age 7 and 15
	- 5. Top_hits_bootstrap_2022-08-24.Rmd - internal validation using non-parametric bootstrap
	- 6. Top_hits_confounders_2022-12-03.Rmd - sensitivity analyses of confounders and mediators
	- 7. Top_hits_MAD_models_2022-08-26.Rmd - sensitivity analyses adjusting for exposures to other adversities
	- 8. Top_hits_trajectories_2022-09-02.Rmd - analysis of the types of DNAm trajectories
	- 9. Additional_analyses_2022-12-09.Rmd - from revisions (QQ plots; stability of DNAm without adversity; EWAS catalog lookup; missMethyl gene ontology)



Please reach out to alussier[at]mgh.harvard.edu or edunn2[at]mgh.harvard.edu for summary statistics or further details of these analyses. 
PS - these documents can also be found here github.com/alussier17/alussier_scripts/adversity-DNAm/
