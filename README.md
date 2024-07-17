# Workflow for processing ITS2 reads through dada2 on the permise server using dada2 and ITSxpress

If conda environments are not already installed on the user premise account install them as follows


Also clone this git repo to your local directory
```
mkdir repo
cd repo
clone ...
```

## Now we are ready to process sequences. 
You should set up a parent directory with some meaningful name to hold all of the work that will be done. Within this directory put the sequences in a subdirectory called 'reads' (you could call it whatever you want, but this part is hardcoded in the script :) . For example: `meaningful_name/reads`

Run the first script. You pass the script to sbatch, and pass the name of your parent directory or "working directory" to the script as an argument
```
sbatch ~/repo/fungi_ITS2_dada2_protocol_07162024/premise/1_dada2_workflow.primer_and_qual_checks.slurm meaningful_name
```

