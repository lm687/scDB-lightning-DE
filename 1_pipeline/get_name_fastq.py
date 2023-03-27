import pandas as pd
import os
import re


flder='/project/cncb/shared/proj002/backup/CELLRANGER_OUT/'

SAMPLES = pd.read_csv("samples.csv", sep='\t')

print(SAMPLES)



def give_fastq_for_sample(sample_it):
	append_all = []
	subfolder=([flder + x for x in SAMPLES[SAMPLES['sample'] == sample_it]['cells'].iloc[0].split(',')])
	subfolder = [x.strip()+'/' for x in subfolder]
	#print(subfolder)
	#print('/////////////')
	for subfolder_it in subfolder:
		#print(subfolder_it)
		subfiles=([x for x in os.listdir(subfolder_it)])
		#print(subfiles)
		#print(re.search('_R2_', subfiles))
		append_all.append([subfolder_it+x for x in subfiles if '_R1_001.fastq.gz' in x][0])		
		append_all.append([subfolder_it+x for x in subfiles if '_R2_001.fastq.gz' in x][0])
	return(' '.join(append_all))


for sample_it in SAMPLES['sample'] :
	print(give_fastq_for_sample(sample_it))


#for subfolder_it in subfolders:	
#	print(subfolder_it)
#	subfiles=([ os.listdir(x.strip()+'/') for x in subfolder_it])
#	print(subfiles)
#	print('/////////////')
#	#print(re.match('_R2_', subfiles))


