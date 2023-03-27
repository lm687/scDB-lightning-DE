snakemake -p --cluster "sbatch --cores 4 -t 00:20:00 --error=logs/%j.err --output=logs/%j.out" --jobs 20 --printshellcmds



