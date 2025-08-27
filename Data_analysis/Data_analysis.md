his folder contains the HTML file “Reuter_et_al_2025_analysis_manuscript.html” of the R code used for data analysis and visualisation for the manuscript Reuter et al. “Single capsid mutations modulating phage adsorption, persistence, and plaque morphology shape evolutionary trajectories in PhiX174”.

Packages required to perform the analysis are found in an R chunk at the beginning of the script.

Data from seven different experiments or assays were acquired and analysed for this manuscript. The experiments as well as the data sets used for the analysis are specified in the code as well as in this file below. 
All data files are available on Zenodo: (REF)

In the following, the individual experiments and the corresponding raw data sets are specified.

1) Experiment 1: 3-hour transfer regime
Phage plaque and bacteria colony counting data for five independent selection lines of phage PhiX174 evolved during 3-hour transfers. Counting data was recorded for each transfer until the first heritable small plaques appeared when plated. Data was used to calculate the estimated frequency of small plaque mutants in the phage population at the first transfer; heritable small plaques were isolated. Figure Fig. S2A was generated based on this data. 
Two data files were used for the analysis of this experiment
- data set 1: “experiment_1_3_h_transfers.csv”; containing raw phage and bacteria counts for each transfer. Data was used to monitor the MOI input for each transfer. 
- data set 2: “experiment_1_3_h_transfers_small_plaque_frequency.csv”; containing frequencies of heritable small plaques isolated compared to the total number of plaques at the last transfer. This data was used to calculate the estimated frequency of small plaque mutants in the phage population at the first transfer, at which small plaques were observed. The data set was used to generate Figure Fig. S2A


2) Experiment 2: 30-min transfer regime
Phage plaque and bacteria colony counting data for five independent selection lines of phage PhiX174 evolved during 30-min transfers. Each passage day consisted of 5 consecutive transfers with four 30-min transfers followed by one 3-h transfer. For transfers 1-4 of each passage day, phages were quantified by spot tests; phages after each 5th transfer were quantified by plaque assay after chloroform extraction. Data was used to monitor the MOI input for each transfer. We did not make a figure in the manuscript for this data set. Data was taken for 10 passage days, corresponding to 40 30-min transfers and 10 3-h transfers.
- data set: “experiment_2_30_min_transfer.csv”, containing raw phage and bacteria counts for each transfer. 

3) Experiment 3: Adsorption constant analysis
Phage plaque and bacteria counting data during a 6-minute adsorption period to determine adsorption dynamics and calculate adsorption constants for the ancestral PhiX174 and mutants F(T101A), F(G322D), and F(S427*). We performed statistical analyses on the calculated adsorption constants to test for statistically significant differences between genotypes. Data was used to generate Figure Fig. 3A.
- data set: “experiment_3_adsorption_assay.csv”, containing raw phage and bacteria counts of multiple experimental replicates for each phage genotype

4) Experiment 4: Image analysis of plaque assays
Data for plaque area and turbidity extracted from the image analysis found in the folder “Plaque_analysis” to visualise morphology differences between phage genotypes. Four independent plaque assay plates were imaged for each phage genotype ancestor, F(T101A), F(G322D), and F(S427*) for image analysis. The data set was used to generate figures Fig. 3B and Fig. S4
-data set: “experiment_4_plaque_image_analysis.csv”; containing area and turbidity values extracted from the image analysis described in folder “Plaque_analysis”

5) Experiment 5: Free phages transferred during 30-minute transfer experiment 
Phage and bacterial counting data for multiple 30-min transfers for phage genotypes ancestor and F(T101A). Free phages were quantified at the end of each 30-min transfer to test the simulation prediction of slow-adsorbers (here F(T101A)) being present at higher rates as free phages at the moment of transfer. Data was used to generate figure Fig. 4D.
- data set: “experiment_5_30_min_free_phages.csv”; containing raw free phage and bacteria counts

6) Experiment 6: Decay analysis
Phage and bacteria counts during a long, host-limiting infection period (< 3 hours) to investigate phage population dynamics of phage genotypes ancestor, F(T101A), F(G322D), and F(S427*) with special emphasis on the decay period. Data was used to generate figures Fig. 5B and Fig. S6.
-data set: “experiment_6_decay_analysis.csv”; containing raw phage and bacteria counts

7) Experiment 7: One-step growth curve of ancestral phage
Phage and bacteria counts during the eclipse and rise period of one infection cycle of ancestral phage PhiX174. This data was used to estimate the burst size and lysis time of the ancestral phage. Data was used to generate Figure Fig. S3.
- data set: “experiment_7_OSGC_WT.csv”; containing raw phage and bacteria counts
