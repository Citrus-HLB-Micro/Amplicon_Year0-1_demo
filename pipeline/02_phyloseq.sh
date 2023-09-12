#!/usr/bin/bash -l
#SBATCH -p short -C xeon -N 1 -n 2 --mem 2gb

module load R

Rscript -e 'library(rmarkdown); rmarkdown::render("scripts/plot_ampliconSummary.Rmd", "html_document")'
