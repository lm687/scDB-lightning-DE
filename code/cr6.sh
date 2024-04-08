#!/bin/bash
#SBATCH --time=2-20:00:00
#SBATCH --partition=long

cd ../data/1_alignment/cellrangerOUT//dmel649ChrimsonV2/G3_rep2;

/project/sims-lab/lmorrill/software/cellranger/cellranger-7.1.0/bin/cellranger count --id=G3_rep2 --fastqs=/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_20,654654_20,660727_20,660728_20,660729_20,660730_20,660731_20,660732_20 --transcriptome=/home/l/lmorrill/projects/general/genomes/dmel649ChrimsonV2 --chemistry=SC3Pv3 --force-cells=10000 --nosecondary


