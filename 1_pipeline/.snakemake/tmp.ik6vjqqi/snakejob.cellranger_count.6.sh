#!/bin/sh
# properties = {"type": "single", "rule": "cellranger_count", "local": false, "input": ["../data/lightning/about_files/about-condition-G3_rep2.txt"], "output": ["/home/l/lmorrill/t1data/data/lightning/alignment/cellrangerOUT/G3_rep2/filtered_feature_bc_matrix/features.tsv.gz"], "wildcards": {"sample": "G3_rep2"}, "params": {"cells": "636200_20,654654_20,660727_20,660728_20,660729_20,660730_20,660731_20,660732_20"}, "log": ["logs/cellranger_count_G3_rep2"], "threads": 1, "resources": {"mem_mb": 1000, "mem_mib": 954, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>"}, "jobid": 6, "cluster": {}}
cd '/Filers/home/l/lmorrill/projects/lighting/1_pipeline' && /home/l/lmorrill/t1data/conda/rufus_env/envs/snakemake_scRNA/bin/python3.11 -m snakemake --snakefile '/Filers/home/l/lmorrill/projects/lighting/1_pipeline/Snakefile' --target-jobs 'cellranger_count:sample=G3_rep2' --allowed-rules 'cellranger_count' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=1000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files '/Filers/home/l/lmorrill/projects/lighting/1_pipeline/.snakemake/tmp.ik6vjqqi' '../data/lightning/about_files/about-condition-G3_rep2.txt' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'software-env' 'mtime' 'code' 'input' 'params' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 5 --scheduler 'ilp' --scheduler-solver-path '/home/l/lmorrill/t1data/conda/rufus_env/envs/snakemake_scRNA/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/Filers/home/l/lmorrill/projects/lighting/1_pipeline/.snakemake/tmp.ik6vjqqi/6.jobfinished' || (touch '/Filers/home/l/lmorrill/projects/lighting/1_pipeline/.snakemake/tmp.ik6vjqqi/6.jobfailed'; exit 1)

