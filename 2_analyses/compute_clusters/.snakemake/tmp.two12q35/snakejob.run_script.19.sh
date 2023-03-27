#!/bin/sh
# properties = {"type": "single", "rule": "run_script", "local": false, "input": ["ABOUT"], "output": ["/home/l/lmorrill/projects/lighting/data/robjects/KB-dmel649ChrimsonV2/clustering_nPCs/28.RDS"], "wildcards": {"nPC": "28"}, "params": {}, "log": ["logs/log_compute_clusters_PCs_resolution_nPC28"], "threads": 1, "resources": {"mem_mb": 1000, "mem_mib": 954, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>"}, "jobid": 19, "cluster": {}}
cd '/Filers/home/l/lmorrill/projects/lighting/2_analyses/compute_clusters' && /home/l/lmorrill/t1data/conda/rufus_env/envs/snakemake_scRNA/bin/python3.11 -m snakemake --snakefile '/Filers/home/l/lmorrill/projects/lighting/2_analyses/compute_clusters/Snakefile' --target-jobs 'run_script:nPC=28' --allowed-rules 'run_script' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=1000' 'mem_mib=954' 'disk_mb=1000' 'disk_mib=954' --wait-for-files '/Filers/home/l/lmorrill/projects/lighting/2_analyses/compute_clusters/.snakemake/tmp.two12q35' 'ABOUT' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'code' 'mtime' 'params' 'input' 'software-env' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 5 --scheduler 'ilp' --scheduler-solver-path '/home/l/lmorrill/t1data/conda/rufus_env/envs/snakemake_scRNA/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/Filers/home/l/lmorrill/projects/lighting/2_analyses/compute_clusters/.snakemake/tmp.two12q35/19.jobfinished' || (touch '/Filers/home/l/lmorrill/projects/lighting/2_analyses/compute_clusters/.snakemake/tmp.two12q35/19.jobfailed'; exit 1)

