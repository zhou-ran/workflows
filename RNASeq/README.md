
# NGS pipeline

All analysis pipeline was based the snakemake

## RNA SEQ (main)

### QC

Include:
 - FastQC
 - picard
 - multiQC

### Trim (optional)

Inlucde:
 - fastp (TODO)
 - Trimmomatic (TODO)

### STAR

All config information was in `config` file, yaml format.

### Quantificaiton


#### RSEM

All config information was in `config` file, yaml format.

#### Htseq

R script

### MultiQC

generate the html report


### Alternative splice

Include:
 - DSU