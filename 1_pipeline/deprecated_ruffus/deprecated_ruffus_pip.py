from ruffus import *
import pandas as pd
import os
import sys


## this the folder where the results from cellranger are going to go
folder_out_cellranger="/home/l/lmorrill/t1data/data/lightning/alignment/cellrangerOUT/"


SAMPLES = pd.read_csv("samples.csv", sep='\t')
print(SAMPLES)
## this fle looks like this:
#     sample                                              cells
# 0  G1_rep1  636200_14,654654_14,660727_14,660728_14,660729...
# 1  G1_rep2  636200_18,654654_18,660727_18,660728_18,660729...
# 2  G2_rep1  636200_15,654654_15,660727_15,660728_15,660729...


DATETIME = "date +'%Y-%m-%dT%H:%M:%S'"

## example of one of the scripts that Charly ran
## change path to the output path for cellranger (output will be created only in the current folder)
#/home/l/lmorrill/t1data/cellranger/cellranger-7.1.0/bin/cellranger count --id=G1_rep1 --fastqs=/t1-data/project/cncb/shared/proj002/backup/CELLRANGER_OUT/ --sample=636200_14,654654_14,660727_14,660728_14,660729_14,660730_14,660731_14,660732_14 --transcriptome=/home/l/lmorrill/t1data/data/genomes/dmel649CORRECTED --chemistry=SC3Pv3 --localcores=$LC --localmem=$LM --force-cells=10000 --nosecondary

## fixed parameters
transcriptome="/home/l/lmorrill/t1data/data/genomes/dmel649CORRECTED"
chemistry="SC3Pv3"
job_threads=2
local_memory="1GB"

input_files = os.listdir("../data/lightning/about_files/")
print(input_files)

@follows("cd %{folder_out_cellranger}s")
@follows(mkdir("count"))
@transform("../data/lightning/about_files/about-condition-*.txt",
           regex(r"../data/lightning/about_files/about-condition-([A-Za-z0-9_]*).txt"),
           r"/home/l/lmorrill/t1data/data/lightning/alignment/cellrangerOUT/\1/filtered_feature_bc_matrix/features.tsv.gz")
def cellranger_count(infile, outfile):
    '''Docstring'''

    sample = re.search('data/([A-Za-z0-9_]*)/.sample', infile).group(1)

    #fastqs = SAMPLES['fastqs'][sample]
    #cells = SAMPLES['cells'][sample]
    cells = SAMPLES[SAMPLES['sample'] == wildcards.sample]['cells'].iloc[0]
    #chemistry = SAMPLES['chemistry'][sample]

    #transcriptome = PARAMS["transcriptome"]

    datetime = DATETIME

    #job_threads = PARAMS["cellranger"]["count"]["threads"]
    #job_memory = PARAMS["cellranger"]["count"]["memory"]

    #local_memory = int(job_memory.replace("G", "")) * job_threads

    statement = """
    %(datetime)s > count/%(sample)s.time &&
    cellranger count
        --id %(sample)s
	--sample %{cells}s
        --transcriptome %(transcriptome)s
        --fastqs %(fastqs)s
        --expect-cells %(cells)s
        --chemistry %(chemistry)s
        --localcores %(job_threads)s
        --localmem %(local_memory)s &&
        mv %(sample)s count/ &&
        touch %(outfile)s &&
    %(datetime)s >> count/%(sample)s.time
    """

    P.run(statement)




def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


