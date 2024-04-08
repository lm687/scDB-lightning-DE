#!/bin/bash
#SBATCH --time=2-20:00:00
#SBATCH --partition=long

cd ../data/1_alignment/cellrangerOUT//dmel649ChrimsonV2/G3_rep1;

/project/sims-lab/lmorrill/software/cellranger/cellranger-7.1.0/bin/cellranger count --id=G3_rep1 --fastqs=/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_16,654654_16,660727_16,660728_16,660729_16,660730_16,660731_16,660732_16 --transcriptome=/home/l/lmorrill/projects/general/genomes/dmel649ChrimsonV2 --chemistry=SC3Pv3 --force-cells=10000 --nosecondary


