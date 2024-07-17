#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="dada2"
#SBATCH --output=dada2_workflow.dada2_core.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate r-dada2_new

cd /mnt/home/garnas/ewj4/EDRR_patho

#Run first R scripts
#mkdir /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_dada2_core

#Rscript ~/repo/BIL_pathogens/sequence_processing_and_dada2/3_dada2_core.r /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_out /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_dada2_core

#mkdir /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core

Rscript ~/repo/BIL_pathogens/sequence_processing_and_dada2/3_dada2_core.r /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_out /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core


