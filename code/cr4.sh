#!/bin/bash
#SBATCH --time=2-20:00:00
#SBATCH --partition=long

cd ../data/1_alignment/cellrangerOUT//dmel649ChrimsonV2/G1_rep1;

/project/sims-lab/lmorrill/software/cellranger/cellranger-7.1.0/bin/cellranger count --id=G1_rep1 --fastqs=/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_14,654654_14,660727_14,660728_14,660729_14,660730_14,660731_14,660732_14 --transcriptome=/home/l/lmorrill/projects/general/genomes/dmel649ChrimsonV2 --chemistry=SC3Pv3 --force-cells=10000 --nosecondary


