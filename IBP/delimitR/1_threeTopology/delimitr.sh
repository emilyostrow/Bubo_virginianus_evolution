#!/bin/bash
#SBATCH --job-name=delimitr             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=12               # CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH -o
#SBATCH --mail-user=emily.ostrow@ku.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --chdir=/panfs/pfs.local/home/e001o910/../../work/bi/e001o910/delimitrBubo/threeModel       # Set working d$
#SBATCH --mem-per-cpu=10gb            # memory requested
#SBATCH --time=1440:00:00
#SBATCH --output=delimitr_%j.log


module load R
R -e "Sys.setenv(RSTUDIO_PANDOC='/panfs/pfs.local/work/bi/bin/pandoc/bin');  rmarkdown::render('clusterDelimitr.Rmd',output_file='clusterDelimitr.html')"