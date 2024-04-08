snakemake -p --cluster "sbatch -t 10:30:00 --cores 1 --ntasks 7 " --jobs 20 --printshellcmds

#snakemake -p --cluster "sbatch --job-name=snakeDB -t 10:30:00 --cores 1 --ntasks 7 --partition=short --out=slurm/%j.out --err=slurm/%j.err" --jobs 20 --printshellcmds 
#snakemake -p --cluster "sbatch -t 20:00:00 --job-name=snakeDB --cores 1 --ntasks 7 --partition=long --out=slurm/%j.out --err=slurm/%j.err" --jobs 30 --printshellcmds 

## tiny jobs
#snakemake -p --cluster "sbatch --job-name=snakeDB -t 00:30:00 --cores 1 --ntasks 7 --partition=short --out=slurm/%j.out --err=slurm/%j.err" --jobs 20 --printshellcmds 


#snakemake -p --cluster "sbatch --job-name=snakeDB -t 3-23:30:00 --cores 12 --mem=450G --ntasks 7 --partition=long --out=slurm/%j.out --err=slurm/%j.err" --jobs 20 --printshellcmds 
#snakemake -p --cluster "sbatch --job-name=snakeDB -t 23:30:00 --cores 12 --mem=250G --ntasks 7 --partition=short --out=slurm/%j.out --err=slurm/%j.err" --jobs 20 --printshellcmds


