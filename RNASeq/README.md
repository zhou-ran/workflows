
# NGS pipeline

All analysis pipeline was based the snakemake

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

TODO

### MultiQC

generate the html report


### TODO

[] Alternative splice, rMATS_4.1
