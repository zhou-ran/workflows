
# DropSeq pipeline

All analysis pipeline was based the snakemake.

How to run:
```shell

snakemake -s Snakefile.smk

```

## Prepare

### 1. sample information

| names          | type | cb                              |
| -------------- | ---- | ------------------------------- |
| Adult-Adipose1 | PE   | example_data/Adult-Adipose1.tsv |
| Adult-Liver1   | PE   | example_data/Adult-Liver1.tsv   |


### 2. config information

All config information was stored in `config.yaml`.

## Run

### 1. QC

Include:
 - FastQC

### 2. Trim (optional)

Inlucde:

 - fastp
 - Trimmomatic

### 3. Fetch barcode and umi information by DropTools
All commands were stored in `rule/pre_process/process_fq.rule`

Include:
 - FastqToSam (rule.fq2bam)
 - TagBamWithReadSequenceExtended (rule.addBarcode)
 - AddPolyTInfo (rule.addToanother)
 - TagBamWithReadSequenceExtended (rule.addUMI)

### 4. Alignment

All commands were stored in `rule/alignment/alignment.rule`


### 5. Run SCAPE

All commands were stored in `rule/apa/scape.rule`
