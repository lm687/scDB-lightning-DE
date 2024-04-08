#!/bin/bash
#SBATCH --time=2-20:00:00
#SBATCH --partition=long
cd ../data/1_alignment/cellrangerOUT//dmel649ChrimsonV2/G2_rep1;
/project/sims-lab/lmorrill/software/cellranger/cellranger-7.1.0/bin/cellranger count --id=G2_rep1 --fastqs=/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_15,654654_15,660727_15,660728_15,660729_15,660730_15,660731_15,660732_15  --transcriptome=/home/l/lmorrill/projects/general/genomes/dmel649ChrimsonV2 --chemistry=SC3Pv3 --force-cells=10000 --nosecondary


