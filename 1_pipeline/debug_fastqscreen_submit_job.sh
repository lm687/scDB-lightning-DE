#!/bin/bash -l
#SBATCH --time=0-00:25:00
#SBATCH --cores=4
#SBATCH --mem=50GB

#module add bwa; /home/l/lmorrill/t1data/conda/rufus_env/envs/bioinformatics_software/bin/fastq_screen  --force --aligner 'bwa' --conf fastq_screen.conf /project/cncb/shared/proj002/backup/CELLRANGER_OUT/660730_18/660730_18_S5_L006_R2_001.fastq.gz --outdir /project/sims-lab/lmorrill/data/lightning/fastQCScreen/

#conda activate snakemake_scRNA
module add perl


#module add bowtie2; /home/l/lmorrill/t1data/conda/rufus_env/envs/bioinformatics_software/bin/fastq_screen  --force --aligner 'bowtie2' --conf fastq_screen.conf /project/cncb/imorgan/3_analysis/sanity_check/links_fastq_water_mem/690901_GX35_S3_L001_R2_001.fastq.gz --outdir /project/sims-lab/lmorrill/data/lightning/fastQCScreen/

module add bowtie2; /home/l/lmorrill/t1data/conda/rufus_env/envs/bioinformatics_software/bin/fastq_screen  --force --aligner 'bowtie2' --conf fastq_screen.conf /project/cncb/imorgan/1_fastq/FASTQ/973017_TY33_S7_L002_R2_001.fastq.gz --outdir /project/sims-lab/lmorrill/data/lightning/fastQCScreen/ --filter 333----------- --tag 



#source deactivate


