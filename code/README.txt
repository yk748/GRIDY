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


The required directories and data files are included in the folder "GRIDY". Below are three ways to reproduce all figures in the manuscript. On a computer with 80 cores, the whole computation takes about "3 DAYS" (Two days for Illustrative example section and one day for Data application section). Therefore, we recommend following Way2 or Way3, if the time used for computation is limited, for example, 10 hours.


############################################
# Way1: When using shell script file "Computation_ALL.sh" on a terminal 
# Note: Please follow this only when there is no time limit on your computation environment:

cd ~/GRIDY/
chmod +x /path/to/Computation_ALL.sh
./Computation_ALL.sh

############################################
# Way2: On a terminal, use the following command lines:

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
# Way3: If using RStudio on Windows or Mac:

## For the figures in Illustrative example section: 
1. Open "Simulation_first.R" and make sure the directory is "./GRIDY/Illustrative_example".
2. execute all lines.
3. Open "Simulation_second.R" and execute all lines.
4. Open "Simulation_third.R" and execute all lines.
5. Open "Figure1.R" and execute all lines.
6. Open "Figure2.R" and execute all lines.
6. Open "Figure3.R" and execute all lines.
7. Open "Figure4.R" and execute all lines.
8. Open "Figure5.R" and execute all lines.

## For the figures in Data application section: 
1. Open "Application.R" and make sure the directory is "./GRIDY/Data_application".
2. execute all lines.
3. Open "FigureS1.R" and execute all lines.
4. Open "FigureS2.R" and execute all lines.
5. Open "FigureS3.R" and execute all lines.
6. Open "Figure6.R" and execute all lines.
7. Open "Figure7.R" and execute all lines.

############################################


