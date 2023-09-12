#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 2 -n 1 --mem 16gb --out logs/plot_figures.log

module load R

#Rscript -e 'library(rmarkdown); rmarkdown::render("scripts/plot_ampliconSummary.Rmd", "html_document")'
Rscript scripts/plot_amplicon_simple.R
