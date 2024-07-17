#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="r-dada2"
#SBATCH --output=dada2.out
#SBATCH --partition=shared
##SBATCH --cpus-per-task=8
##SBATCH --mem=128000

module purge
module load anaconda/colsa

###
#This script runs the primary dada2 algorithm, i.e., assuming that quality filtering and ITS extraction has already been performed

conda activate dada2-check

# first run dada2 for the concatenated files
cd ~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_files_cat
mkdir dada2_processing_tables_figs
mkdir intermediate_RDS

srun Rscript ~/repo/neonectria_barcoding_012220/dada2-slurm.r
srun Rscript ~/repo/neonectria_barcoding_012220/UNITE_taxonomic_classification-slurm.r

# next running the non-concatenated files (this is pooled across all run1/run2 samples; the pooling within run only is already done)

cd ~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_files_sep
mkdir dada2_processing_tables_figs
mkdir intermediate_RDS

srun Rscript ~/repo/neonectria_barcoding_012220/dada2-slurm_files_sep.r
srun Rscript ~/repo/neonectria_barcoding_012220/UNITE_taxonomic_classification-slurm.r

#Then running run1 and run2 separately

cd ~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_files_sep
mkdir dada2_processing_tables_figs
mkdir intermediate_RDS

srun Rscript ~/repo/neonectria_barcoding_012220/dada2-slurm_files_sep.r
srun Rscript ~/repo/neonectria_barcoding_012220/UNITE_taxonomic_classification-slurm.r

cd ~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run2_files_sep
mkdir dada2_processing_tables_figs
mkdir intermediate_RDS

srun Rscript ~/repo/neonectria_barcoding_012220/dada2-slurm_files_sep.r
srun Rscript ~/repo/neonectria_barcoding_012220/UNITE_taxonomic_classification-slurm.r

conda deactivate

