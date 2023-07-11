#!/bin/bash
# Brenda Pardo: 2023-06-08
# run analysis pipeline
# pipeline requires
# gofigure: https://gitlab.com/evogenlab/GO-Figure
# go figure has several python library requirements, we recommend installing this
# in a conda environment and activating it below.

# org R library: org.Pcanaliculata.eg.db
# install.packages("org.Pcanaliculata.eg.db", repos=NULL)

# the following public R libraries
#clusterProfiler
#corrplot
#cowplot
#dplyr
#edgeR
#enrichplot
#ggplot2
#ggpubr
#ggrepel
#GOSemSim
#here
#matrixStats
#openxlsx
#purrr
#RColorBrewer
#readr
#stringr
#tibble
#tidyr
#tidyverse
#xlsx


###############################################
# Set environment for R and gofigure
###############################################

# if you need to activate a conda environment for gofigure to run do it here.
source /home/ejr/anaconda3/etc/profile.d/conda.sh
conda activate snail-pipeline
GOFIGURE="/n/projects/ejr/aai/brenda_pipeline/GO-Figure/gofigure.py"
###############################################

###############################################
# Expression Analysis using 1dpa as reference
###############################################
if [ ! -e 1dp1-reference.finished ]
then
DIR=1dpa-reference

#Run Differential Expression Analysis
Rscript ${DIR}/02_dea/01_dea-1dpa-ref.R
Rscript ${DIR}/02_dea/02_create-logFC-table.R
Rscript ${DIR}/02_dea/03_create-logFC-logCPM-pval-table.R
Rscript ${DIR}/02_dea/04_dea-plots.R

#Run GO enrichment Analysis
Rscript ${DIR}/03_go/01_go-input.R
Rscript ${DIR}/03_go/02_clusterProfiler.R
Rscript ${DIR}/03_go/03_go-degenes-search.R
Rscript ${DIR}/03_go/04_go-plot.R

#GO-figure
Rscript ${DIR}/04_go-figure/01_go-figure-input.R
Rscript ${DIR}/04_go-figure/02_run-go-figure.R ${GOFIGURE}
touch 1dp1-reference.finished
fi

###############################################
# Expression Analysis using intact eyestalk  as reference
###############################################
if [ ! -e intact-reference.finished ]
then
DIR=intact-reference

#Run Differential Expression Analysis
Rscript ${DIR}/02_dea/01_dea-intact-ref.R
Rscript ${DIR}/02_dea/02_create-logFC-table.R
Rscript ${DIR}/02_dea/03_create-logFC-logCPM-pval-table.R
Rscript ${DIR}/02_dea/04_dea-plots.R

#Run GO enrichment Analysis
Rscript ${DIR}/03_go/01_go-input.R
Rscript ${DIR}/03_go/02_clusterProfiler.R
Rscript ${DIR}/03_go/03_go-degenes-search.R
Rscript ${DIR}/03_go/04_go-plot.R

#GO-figure
Rscript ${DIR}/04_go-figure/01_go-figure-input.R
Rscript ${DIR}/04_go-figure/02_run-go-figure.R ${GOFIGURE}
touch intact-reference.finished
fi
###############################################
# Expression Analysis using previous timepoints (Sequentially) as reference
###############################################
if [ ! -e sequential-reference.finished ]
then
DIR=sequential-reference

#Run Differential Expression Analysis
Rscript ${DIR}/02_dea/01_dea-seq-ref.R
Rscript ${DIR}/02_dea/02_create-logFC-table.R
Rscript ${DIR}/02_dea/03_create-logFC-logCPM-pval-table.R
Rscript ${DIR}/02_dea/04_dea-plots.R

#Run GO enrichment Analysis
Rscript ${DIR}/03_go/01_go-input.R
Rscript ${DIR}/03_go/02_clusterProfiler.R
Rscript ${DIR}/03_go/03_go-degenes-search.R
Rscript ${DIR}/03_go/04_go-plot.R

#GO-figure
Rscript ${DIR}/04_go-figure/01_go-figure-input.R
Rscript ${DIR}/04_go-figure/02_run-go-figure.R ${GOFIGURE}
touch sequential-reference.finished
fi
###############################################
# analysis combining results from more than one previous analysis
###############################################
if [ ! -e mixed.finished ]
then
DIR=mixed-analysis
#Make go plots
Rscript ${DIR}/02_go/01_go-plot-mixed-analysis.R
Rscript ${DIR}/02_go/02_go-genes-plot-mixed-analysis.R
touch mixed.finished
fi
