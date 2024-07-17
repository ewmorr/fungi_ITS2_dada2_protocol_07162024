#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="dada2"
#SBATCH --output=dada2_workflow.dada2_core.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate r-dada2_new

workDir=$HOME/$1
cd $workDir

Rscript ~/repo/fungi_ITS2_dada2_protocol_07162024/sequence_processing_and_dada2/3_dada2_core.r $workDir


