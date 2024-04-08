#!/bin/bash
#SBATCH --time=2-20:00:00
#SBATCH --partition=long

cd ../data/1_alignment/cellrangerOUT//dmel649ChrimsonV2/G2_rep2;

/project/sims-lab/lmorrill/software/cellranger/cellranger-7.1.0/bin/cellranger count --id=G2_rep2 --fastqs=/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_19,654654_19,660727_19,660728_19,660729_19,660730_19,660731_19,660732_19 --transcriptome=/home/l/lmorrill/projects/general/genomes/dmel649ChrimsonV2 --chemistry=SC3Pv3 --force-cells=10000 --nosecondary


