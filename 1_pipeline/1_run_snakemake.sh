#snakemake -np --cluster "sbatch -t 00:02:00 --cores 1" --jobs 1 --printshellcmds --unlock
#snakemake -p --cluster "sbatch -t 24:00:00 --cores 4" --jobs 6 --printshellcmds
snakemake -p --cluster "sbatch --cores 4" --jobs 6 --printshellcmds

