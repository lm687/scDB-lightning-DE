# lightning-DE

- Conda environment
```
module load R-cbrg
conda activate snakemake_scRNA
```

## Computational requirements

Step | Time (Day-Hours:Minutes:Seconds) | CPU | Mmeory (GB) | 
--- | --- | --- | --- |
Cellranger alignment | 1-00:40:00 | 1 | 15 | 

## Alignment
### Kallisto

```
kb ref ../genomes/dmel-all-chromosome-r6.49_Chrimson.fasta ../genomes/dmel-all-r6.49_CORRECTED_Chrimson.gtf -i dmel-all-r6.49_CORRECTED_Chrimson_KallistoINDEX -g dmel-all-r6.49_CORRECTED_Chrimson_KallistoT2G -f1 dmel-all-r6.49_CORRECTED_ChrimsonFASTA
```
I assume that, if the chemistry if SC3Pv3, then the technology is 10XV3.


