#!/bin/bash
#SBATCH --time=2-20:00:00
#SBATCH --partition=long

cd ../data/1_alignment/cellrangerOUT//dmel649ChrimsonV2/G4_rep2;

/project/sims-lab/lmorrill/software/cellranger/cellranger-7.1.0/bin/cellranger count --id=G4_rep2 --fastqs=/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_21,654654_21,660727_21,660728_21,660729_21,660730_21,660731_21,660732_21 --transcriptome=/home/l/lmorrill/projects/general/genomes/dmel649ChrimsonV2 --chemistry=SC3Pv3 --force-cells=10000 --nosecondary


