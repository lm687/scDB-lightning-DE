#!/bin/bash
#SBATCH --time=2-20:00:00
#SBATCH --partition=long
cd ../data/1_alignment/cellrangerOUT//dmel649ChrimsonV2/G1_rep2;

/project/sims-lab/lmorrill/software/cellranger/cellranger-7.1.0/bin/cellranger count --id=G1_rep2 --fastqs=/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_18,654654_18,660727_18,660728_18,660729_18,660730_18,660731_18,660732_18 --transcriptome=/home/l/lmorrill/projects/general/genomes/dmel649ChrimsonV2 --chemistry=SC3Pv3 --force-cells=10000 --nosecondary

