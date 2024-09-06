These R source files and raw data files contain all code to reproduce the figures 1-7. and S1-S3 for the paper 'Group integrative dynamic factor models with application to multiple subject brain connectivity' by Y. Kim, Z. F. Fisher, and V. Pipiras. The version that contains processed files including simulation and data application results from these R codes can be found at https://zenodo.org/records/13711280.

# ---------------------------------------------------------------------------------- #
The version of R used in the simulation and data application is 4.0.5 (2021-03-31). All necessary packages are listed below. The versions of the package may not be the latest, however, the reported versions are based on the latest execution by the authors. The information provided is from the result of sessionInfo(). Before execution, PLEASE MAKE SURE that all required packages are installed.

R version 4.0.5 (2021-03-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 9.0 (Blue Onyx)

Matrix products: default
BLAS:   /programs/R-4.0.5-r9/lib/libRblas.so
LAPACK: /programs/R-4.0.5-r9/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] gridExtra_2.3      reshape2_1.4.4     RColorBrewer_1.1-3 ggpubr_0.5.0
 [5] latex2exp_0.9.6    ggplot2_3.4.0      readxl_1.4.1       ica_1.0-3
 [9] multiway_1.0-6     CMLS_1.0-1         quadprog_1.5-8     combinat_0.0-8
[13] Matrix_1.5-1       mvtnorm_1.2-5

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.12      plyr_1.8.7       cellranger_1.1.0 pillar_1.8.1
 [5] compiler_4.0.5   tools_4.0.5      lifecycle_1.0.3  tibble_3.1.8
 [9] gtable_0.3.1     lattice_0.20-45  pkgconfig_2.0.3  rlang_1.0.6
[13] cli_3.4.1        DBI_1.1.3        withr_2.5.0      dplyr_1.0.10
[17] stringr_1.4.1    generics_0.1.3   vctrs_0.5.1      grid_4.0.5
[21] tidyselect_1.2.0 glue_1.6.2       R6_2.5.1         rstatix_0.7.1
[25] fansi_1.0.3      carData_3.0-5    car_3.1-1        tidyr_1.2.1
[29] purrr_0.3.5      magrittr_2.0.3   backports_1.4.1  scales_1.2.1
[33] abind_1.4-5      assertthat_0.2.1 colorspace_2.0-3 ggsignif_0.6.4
[37] utf8_1.2.2       stringi_1.7.8    munsell_0.5.0    broom_1.0.1

# ---------------------------------------------------------------------------------- #
This directory structure outlines the stored R source files and raw data for the data application provided by the authors. Please visit https://zenodo.org/records/13711280 to get the processed files from the simulation and data application codes. And make sure that the downloaded files are placed in the right subdirectories. The parent folder contains R source codes and two subfolders, and their subdirectories also contain the subfolders and different R source codes.

README.txt: this file.
AJIVE_retrieve.R : Since many years has passed from the last update from the maintainers, 
                   some command lines in ajive are not working as desired. To align with the current R version, 
                   the main function of ajive package is extracted. See Iain Carmichael (2020) for the ajive package.
		   [citation: r jive: First Github release (2020), GitHub Repository, https://github.com/idc9/r_jive]
gica.R : this function is used to implement group ICA as a benchmark, 
         which is extracted from online supplementary file of Helwig, N.E. and Snodgress, M.A. (2019, NeuroImage)
	 [citation: Helwig, N.E. and Snodgress, M.A. (2019). "Exploring individual and group differences in latent brain networks 
                    using cross-validated simultaneous component analysis", NeuroImage, 201:116019.]
Rotational_bootstrap.R : this function is created to implement principal angle-based estimation methods for 
                         initial rank of multi-view data by using bootstrap by Prothero et al. (2024, TEST).
                         [citation: Prothero, J., Jiang, M., Hannig, J., Tran-Dinh, Q., Ackerman, A., 
                                    and Marron, J.S. (2024). Data integration via analysis of subspaces (DIVAS). TEST, 14, 1-42.]
library_simulation.R : This source code contains packages required, sim_model (data generating function), 
                       factor_regression (estimating factor series), and YW_compute (computing VAR transition matrices 
                       and covariance matrices via Yule-Walker equations) used for the Illustrative_examples section.
library_application.R : This source code contains packages required, factor_regression (estimating factor series), 
                        and YW_compute (computing VAR transition matrices and covariance matrices via Yule-Walker equations) 
                        used for the Illustrative_examples section.

./Illustrative_examples/:
A package of R source codes that reproduce the result in the Illustrative examples section.
 
   ./figures/ : a subfolder where Figures 1 - 5 are stored.
   ./result/ : a subfolder where the simulation result RData files SHOULD BE stored. 

   Simulation_first.R : code for producing results for the first simulation setting in the Illustrative example section.
   Simulation_second.R : code for producing results for the second simulation setting in the Illustrative example section.
   Simulation_third.R : code for producing results for the third simulation setting in the Illustrative example section.
   Figure1.R : code for producing results for plotting Figure 1 in the Illustrative example section
   Figure2.R : code for producing results for plotting Figure 2 in the Illustrative example section
   Figure3.R : code for producing results for plotting Figure 3 in the Illustrative example section
   Figure4.R : code for producing results for plotting Figure 4 in the Illustrative example section
   Figure5.R : code for producing results for plotting Figure 5 in the Illustrative example section

./Data_application/:
A package of raw data and R source codes that reproduce the result in the Data application section.

   ./data/ : under this folder, there are three types of Excel files that the authors collected and stored.
    - cvs files : the fMRI data from each subject collected at each site. 
                  The data can be accessible through the ABIDE Preprocessed website 
                  http://preprocessed-connectomes-project.org/abide/download.html and 
                  how to download is exactly explained on the page.
		  [citation: Dosenbach, Nico UF, Binyam Nardos, Alexander L. Cohen, Damien A. Fair, Jonathan D. Power, 
                             Jessica A. Church, Steven M. Nelson et al. (2010) "Prediction of individual brain maturity using fMRI.", 
                             Science, 329(5997), 1358-1361.]

    - dos160_labels.xlsx : information of ROI numbers, ROI labels, names of networks that each ROI belongs to, 
                           abbreviations of the networks, and pairs of network abbreviations and ROI numbers.
                           The atlas used in this study follows Dosenbach et al. (2010, Science).
    - file_name_list.xlsx : information of the collecting sites, subject (file) number, and group indicators; 
                            1 for the subject with autism spectrum disorder and 2 for the normal subjects under control
   ./figures/ : a subfolder where Figures 6,7, and S1 - S3 are stored.
   ./result/ : a subfolder where the data application result RData files SHOULD BE stored. 

   Application.R : code for producing results in the Data application section. 
   Figure6.R : code for producing results for plotting Figure 6 in the Data application section.
   Figure7.R : code for producing results for plotting Figure 7 in the Data application section.
   FigureS1.R : code for producing results for plotting Figure S1 in Supplemental material
   FigureS2.R : code for producing results for plotting Figure S2 in Supplemental material
   FigureS3.R : code for producing results for plotting Figure S3 in Supplemental material


# ---------------------------------------------------------------------------------- #
The required directories and data files are included in the folder "GRIDY". Below are a couple of ways to reproduce all figures in the manuscript. On a computer with a single core, the whole computation takes about "3 DAYS". Therefore, we recommend choosing Way 1 and Way 2 below, depending on whether the computation time is limited, for example, 10 hours. Way 3 is the case where the ready-simulated results are used for merely producing the figures.

############################################
# Way 1: On a terminal, use the following command lines:
# Note: Please MAKE SURE that you can spend enough computation time for Simulation_first.R and Simulation_second.R.

cd ~/GRIDY/Illustrative_example
Rscript Simulation_first.R
Rscript Simulation_second.R
Rscript Simulation_third.R

Rscript Figure1.R
Rscript Figure2.R
Rscript Figure3.R
Rscript Figure4.R
Rscript Figure5.R

cd ~/GRIDY/Data_application
Rscript Application.R

Rscript FigureS1.R
Rscript FigureS2.R
Rscript FigureS3.R
Rscript Figure6.R
Rscript Figure7.R

############################################
# Way 2: If using RStudio on Windows or Mac:

## For the figures in Illustrative example section: 
1. Open "Simulation_first.R" and make sure the directory is "./GRIDY/Illustrative_example".
2. Execute all lines. The simulated files named "sim1_d-T-K-_type-.RData" are stored under "./GRIDY/Illustrative_example/result".
3. Open "Simulation_second.R". 
4. Execute all lines. The simulated files named "sim2_snr-_type-.RData" are stored under "./GRIDY/Illustrative_example/result".
5. Open "Simulation_third.R". 
6. Execute all lines. The simulated files named "sim3_set-_d-T-K-.RData" are stored under "./GRIDY/Illustrative_example/result".
7. Open "Figure1.R" and execute all lines. "Figure1.pdf" is stored under "./GRIDY/Illustrative_example/figures". 
This R file will load "sim1_d-T-K-_type1.RData".
8. Open "Figure2.R" and execute all lines. "Figure2.pdf" is stored under "./GRIDY/Illustrative_example/figures". 
This R file will load "sim1_d-T-K-_type2.RData".
9. Open "Figure3.R" and execute all lines. "Figure3.pdf" is stored under "./GRIDY/Illustrative_example/figures". 
This R file will load "sim2_snr-_type1.RData".
10. Open "Figure4.R" and execute all lines. "Figure4.pdf" is stored under "./GRIDY/Illustrative_example/figures". 
This R file will load "sim2_snr-_type2.RData".
11. Open "Figure5.R" and execute all lines. "Figure5.pdf" is stored under "./GRIDY/Illustrative_example/figures". 
This R file will load "sim3_set-_d-T-K-.RData".

## For the figures in Data application section: 
1. Open "Application.R" and make sure the directory is "./GRIDY/Data_application".
2. execute all lines. Separate files named "result_intermediate.RData", "result_final_X.RData",  "result_final_GRIDY.RData", "result_final_SCA_P.RData", and "result_final_GICA.RData" are stored under "./GRIDY/Data_application/result"
3. Open "FigureS1.R" and execute all lines. "FigureS1.pdf" is stored under "./GRIDY/Data_application/figures". 
This R file will load "result_intermediate.RData".
4. Open "FigureS2.R" and execute all lines. "FigureS2.pdf" is stored under "./GRIDY/Data_application/figures". 
This R file will load "result_final_X.RData","result_final_GRIDY.RData","result_final_SCA_P.RData", and "result_final_GICA.RData".
5. Open "FigureS3.R" and execute all lines. "FigureS3.pdf" is stored under "./GRIDY/Data_application/figures". 
This R file will load "result_final_X.RData","result_final_GRIDY.RData","result_final_SCA_P.RData", and "result_final_GICA.RData".
6. Open "Figure6.R" and execute all lines. "Figure6.pdf" is stored under "./GRIDY/Data_application/figures". 
This R file will load "result_final_X.RData","result_final_GRIDY.RData","result_final_SCA_P.RData", and "result_final_GICA.RData".
7. Open "Figure7.R" and execute all lines. "Figure7.pdf" is stored under "./GRIDY/Data_application/figures". 
This R file will load "result_final_X.RData","result_final_GRIDY.RData","result_final_SCA_P.RData", and "result_final_GICA.RData".


############################################
# Way 3: If the ready-simulated results are used to produce the figures:

## 1) On the command:
### For the figures in Illustrative example section: 
cd ~/GRIDY/Illustrative_example
Rscript Figure1.R
Rscript Figure2.R
Rscript Figure3.R
Rscript Figure4.R
Rscript Figure5.R

### For the figures in Data application section: 
cd ~/GRIDY/Data_application
Rscript FigureS1.R
Rscript FigureS2.R
Rscript FigureS3.R
Rscript Figure6.R
Rscript Figure7.R



## 2) On RStudio:
### For the figures in Illustrative example section: 
1. Open "Figure1.R" and execute all lines. 
2. Open "Figure2.R" and execute all lines. 
3. Open "Figure3.R" and execute all lines. 
4. Open "Figure4.R" and execute all lines. 
5. Open "Figure5.R" and execute all lines. 

### For the figures in Data application section: 
1. Open "FigureS1.R" and execute all lines. 
2. Open "FigureS2.R" and execute all lines. 
3. Open "FigureS3.R" and execute all lines. 
4. Open "Figure6.R" and execute all lines. 
5. Open "Figure7.R" and execute all lines. 
