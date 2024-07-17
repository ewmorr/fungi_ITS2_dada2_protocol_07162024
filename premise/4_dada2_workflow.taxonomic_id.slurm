#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="dada2"
#SBATCH --output=dada2_workflow.taxonomic_id.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate r-dada2_new

cd /mnt/home/garnas/ewj4/EDRR_patho


Rscript ~/repo/BIL_pathogens/sequence_processing_and_dada2/4_UNITE_taxonomic_classification.r /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core/dada2_seq_table_no_chim.rds /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core /mnt/home/garnas/ewj4/blast_dbs/unite_07252023/sh_general_allEuk_dynamic_singletons_25072023.ITS2.fasta


