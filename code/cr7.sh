#!/bin/bash
#SBATCH --time=2-20:00:00
#SBATCH --partition=long

cd ../data/1_alignment/cellrangerOUT//dmel649ChrimsonV2/G4_rep1;

/project/sims-lab/lmorrill/software/cellranger/cellranger-7.1.0/bin/cellranger count --id=G4_rep1 --fastqs=/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_17,654654_17,660727_17,660728_17,660729_17,660730_17,660731_17,660732_17 --transcriptome=/home/l/lmorrill/projects/general/genomes/dmel649ChrimsonV2 --chemistry=SC3Pv3 --force-cells=10000 --nosecondary


