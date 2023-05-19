#module add fastq_screen; module add bowtie2; module add bwa; fastq_screen --aligner 'bwa' --conf fastq_screen.conf /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660730_18/660730_18_S5_L006_R2_001.fastq.gz --outdir /project/sims-lab/lmorrill/data/lightning/fastQCScreen/

source activate bioinformatics_software; module add bwa; /home/l/lmorrill/t1data/conda/rufus_env/envs/bioinformatics_software/bin/fastq_screen  --force --aligner 'bwa' --conf fastq_screen.conf /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660730_18/660730_18_S5_L006_R2_001.fastq.gz --outdir /project/sims-lab/lmorrill/data/lightning/fastQCScreen/



