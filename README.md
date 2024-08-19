# Workflow for processing ITS2 reads to ASVs on the permise server using dada2 and ITSxpress

## software setup
If conda environments are not already installed on the user premise account install them as follows
```
module load anaconda/colsa

#itsxpress env
mamba create -n itsxpressenv -c bioconda -c conda-forge itsxpress
mamba activate itsxpressenv
itsxpress -h
conda deactivate

#dada2 env
conda create --name r-dada2_new --clone template
conda activate r-dada2_new
conda install bioconda::bioconductor-dada2
conda install r-dplyr
```

Then open R and check that dada2 is installed and loads. Open an R interpreter by typing a capital R on the command line
```
R
```
Once R is open you can run R commands as normal
```
library(dada2) #load dada2 
q(save="no") #exit R without saving
```
Once you have the environment set up you can run dada2 via slurm by setting up a slurm script to first activate the dada2 conda environment and then calling an R script using `Rscript` command. 

### Finally...
Also clone this git repo to your local directory
```
mkdir repo
cd repo
git clone https://github.com/ewmorr/fungi_ITS2_dada2_protocol_07162024
```

## Now we are ready to process sequences. 
You should set up a parent directory with some meaningful name to hold all of the work that will be done. Within this directory put the sequences in a subdirectory called 'reads' (you could call the reads dir whatever you want, but this part is hardcoded in the script so you would have to change that manually :) . For example: `test_run/reads`
```
mkdir test_run
cd test_run
mkdir reads
cp path_to/seqs/*fastq.gz ./reads
```
You can check the number of reads per sample.
```
for i in reads/*R1*fastq*
do(
    r1File=${i##*/}
    sample=${r1File%_S*} #you can change this to pick out a different part of the file
    numReads=$(echo $(zcat $i | wc -l) / 4 | bc)
    echo -e "$sample\t$numReads" >> reads_per_sample.txt
)
done
```

### Step 1. check primer/adapter orientation, remove adapters, and check quality
Run the first script. You pass the script to sbatch, and pass the name of your parent directory or "working directory" to the script as an argument
```
sbatch ~/repo/fungi_ITS2_dada2_protocol_07162024/premise/1_dada2_workflow.primer_and_qual_checks.slurm test_run
```
The script has written some initial quality checks to the `dada2_processing_tables_figs` directory. Download these and take a look

### Step 2. Run itsxpress
Make sure to pass your working directory as the last argument
```
sbatch ~/repo/fungi_ITS2_dada2_protocol_07162024/premise/2_dada2_workflow.itsxpress.slurm test_run
```
Check the output
```
less dada2_workflow.itsxpress.out
```

### Step 3. Run the core dada2 algorithm

```
sbatch ~/repo/fungi_ITS2_dada2_protocol_07162024/premise/3_dada2_workflow.dada2.slurm test_run
```
Check `dada2_processing_tables_figs/read_processing_tracking.csv` for a breakdown of the number of reads remaining after each prcessing step.

### Step 4. Taxonopmic classification
```
sbatch ~/repo/fungi_ITS2_dada2_protocol_07162024/premise/4_dada2_workflow.taxonomic_id.slurm test_run
```
