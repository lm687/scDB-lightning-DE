#!/bin/sh
# properties = {"type": "single", "rule": "kallisto_alignment", "local": false, "input": ["../data/lightning/about_files/about-condition-G4_rep1.txt"], "output": ["/project/sims-lab/lmorrill/data/lightning/alignment/KallistoOUT/dmel649Chrimson/G4_rep1"], "wildcards": {"sample": "G4_rep1"}, "params": {"cells": "/project/cncb/shared/proj002/backup/CELLRANGER_OUT/636200_17 /project/cncb/shared/proj002/backup/CELLRANGER_OUT/654654_17 /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660727_17 /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660728_17 /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660729_17 /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660730_17 /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660731_17 /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660732_17"}, "log": ["logs/kallisto_count_G4_rep1"], "threads": 1, "resources": {"mem_mb": 1000, "mem_mib": 954, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>"}, "jobid": 7, "cluster": {}}
cd '/Filers/home/l/lmorrill/projects/lighting/1_pipeline' && /home/l/lmorrill/t1data/conda/rufus_env/envs/py310/bin/python3.10 -m snakemake --snakefile '/Filers/home/l/lmorrill/projects/lighting/1_pipeline/Snakefile' --target-jobs 'kallisto_alignment:sample=G4_rep1' --allowed-rules 'kallisto_alignment' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=1000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files '/Filers/home/l/lmorrill/projects/lighting/1_pipeline/.snakemake/tmp.fhn7faz3' '../data/lightning/about_files/about-condition-G4_rep1.txt' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'mtime' 'params' 'code' 'input' 'software-env' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 5 --scheduler 'ilp' --scheduler-solver-path '/home/l/lmorrill/t1data/conda/rufus_env/envs/py310/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/Filers/home/l/lmorrill/projects/lighting/1_pipeline/.snakemake/tmp.fhn7faz3/7.jobfinished' || (touch '/Filers/home/l/lmorrill/projects/lighting/1_pipeline/.snakemake/tmp.fhn7faz3/7.jobfailed'; exit 1)

