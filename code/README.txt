These R source files and shell-script contain all code to reproduce the figures 1-7. and S1-S3 for the paper 'Group integrative dynamic factor models with application to multiple subject brain connectivity' by Y. Kim, Z. F. Fisher, and V. Pipiras.

The version of R used in the simulation and data application is 4.0.5 (2021-03-31). All necessary packages are listed below. The versions of the package may not be the latest, however, the reported versions are based on the latest execution by the authors. Before execution, PLEASE MAKE SURE that all required packages are installed.

mvtnorm 1.2.5
Matrix	1.5.1
combinat 0.0.8
multiway 1.0.6
ica 1.0.3
ggplot2 3.4.0
latex2exp 0.9.6
ggpubr 0.5.0
RColorBrewer 1.1.3
reshape2 1.4.4
gridExtra 2.3


The required directories and data files are included in the folder "GRIDY". Below are a couple of ways to reproduce all figures in the manuscript. On a computer with 80 cores, the whole computation takes about "3 DAYS" (Two days for the Illustrative example section and one day for the Data application section), in particular, the simulation for the first 4 figures (by Simulation_first.R and Simulation_second.R) takes most of the time. Therefore, we recommend choosing Way 1 and Way 2 below, depending on whether the computation time is limited, for example, 10 hours. Way 3 is the case where the ready-simulated results are used for merely producing the figures.

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
2. Execute all lines. The simulated files named "sim1_d...RData" are stored under "./GRIDY/Illustrative_example/result".
3. Open "Simulation_second.R". 
4. Execute all lines. The simulated files named "sim2_snr...RData" are stored under "./GRIDY/Illustrative_example/result".
5. Open "Simulation_third.R". 
6. Execute all lines. The simulated files named "sim3_set...RData" are stored under "./GRIDY/Illustrative_example/result".
7. Open "Figure1.R" and execute all lines. This file will load and use "sim1_d...RData".
8. Open "Figure2.R" and execute all lines. This file will load and use "sim1_d...RData".
9. Open "Figure3.R" and execute all lines. This file will load and use "sim2_snr...RData".
10. Open "Figure4.R" and execute all lines. This file will load and use "sim2_snr...RData".
11. Open "Figure5.R" and execute all lines. This file will load and use "sim3_set...RData".

## For the figures in Data application section: 
1. Open "Application.R" and make sure the directory is "./GRIDY/Data_application".
2. execute all lines. Two separate files named "result_intermediate.RData" and "result_final.RData" are stored under "./GRIDY/Data_application/result"
3. Open "FigureS1.R" and execute all lines. This file will load and use "result_intermediate.RData".
4. Open "FigureS2.R" and execute all lines. This file will load and use "result_final.RData".
5. Open "FigureS3.R" and execute all lines. This file will load and use "result_final.RData".
6. Open "Figure6.R" and execute all lines. This file will load and use "result_final.RData".
7. Open "Figure7.R" and execute all lines. This file will load and use "result_final.RData".

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
