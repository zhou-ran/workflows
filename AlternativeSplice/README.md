
# Alternative splice analysis using Exonic Part


`Note: Bedtools2.23`



```shell
snakemake -s Snakefile.smk
```

## RNA SEQ (main)

### QC

Include:
 - FastQC
 - picard
 - multiQC
 - GeneBodyCoverge

### Trim (optional)

Inlucde:

 - fastp
 - Trimmomatic

### STAR

All config information was in `config` file, yaml format.

### Quantificaiton

#### RSEM

All config information was in `config` file, yaml format.

#### Htseq

R script

### Alternative splice

[ExonicPart](https://github.com/LuChenLab/MeDAS)

### MultiQC

generate the html report


