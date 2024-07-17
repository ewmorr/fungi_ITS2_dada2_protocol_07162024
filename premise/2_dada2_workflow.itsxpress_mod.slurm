#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="itsxpress"
#SBATCH --output=dada2_workflow.itsxpress_mod.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
conda activate trim_3p

cd /mnt/home/garnas/ewj4/EDRR_patho

mkdir itsxpress_mod_out
outdir=itsxpress_mod_out

for i in cutadapt/*R1*
do(

    dir=${i%/*}
    r1File=${i##*/}
    pre=${r1File%R1*}
    post=${r1File##*R1}
    r2File=${pre}R2${post}

    #if [ ! -f itsxpress_mod_out/$r2File ]; then
        ~/.local/bin/itsxpress --fastq $dir/$r1File --fastq2 $dir/$r2File \
            --reversed_primers --threads 24 --cluster_id 1 --taxa Fungi --region ITS2 \
            --outfile $outdir/$r1File --outfile2 $outdir/$r2File
    #fi
)
done

