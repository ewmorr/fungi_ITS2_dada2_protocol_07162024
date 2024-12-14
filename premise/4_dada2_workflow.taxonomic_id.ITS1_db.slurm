#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="dada2"
#SBATCH --output=dada2_workflow.taxonomic_id.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa
conda activate r-dada2_new

workDir=$HOME/$1
cd $workDir

indir=$workdir/itsxpress_out
outdir=$workdir/dada2_core


Rscript ~/repo/fungi_ITS2_dada2_protocol_07162024/sequence_processing_and_dada2/4_UNITE_taxonomic_classification.r $workDir/dada2_core/dada2_seq_table_no_chim.rds $workDir/dada2_core /mnt/home/garnas/ewj4/blast_dbs/unite_07252023/sh_general_allEuk_dynamic_singletons_25072023.ITS1.fasta


