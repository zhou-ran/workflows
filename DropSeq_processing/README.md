
# DropSeq pipeline

All analysis pipeline was based the snakemake.


How to run:
```shell

snakemake -s Snakefile.smk

```

If you are not familiar with snakemake, I have listed all the commands, you can run it directly.

## Prepare

### 1. sample information

| names          | type | cb                              |
| -------------- | ---- | ------------------------------- |
| Adult-Adipose1 | PE   | example_data/Adult-Adipose1.tsv |
| Adult-Liver1   | PE   | example_data/Adult-Liver1.tsv   |

The fastq files were download from [CNP0000325](https://db.cngb.org/search/project/CNP0000325/).

### 2. config information

All config information was stored in `config.yaml`.
including:
 - The format and location of fastq file format.
 - The location and version of soft and script.
 - Genome, genome index and annotation files.
 - Output dir and temporary dir.


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

```shell

# Convert fastq file into unmapped bam file.

java -Xmx20G -jar picard.jar\
FastqToSam \
F1=sample_R1.fastq.gz \
F2=sample_R2.fastq.gz \
SM=DS \
O=sample_unmapped.bam

```

 - TagBamWithReadSequenceExtended (rule.addBarcode)

```shell

# Add barcode sequence into tags of unmapped bam.
# In Microwell-seq, the barcode sequence located in 1-6:22-27:43-48 of R1 reads.
# In 10X genomics, the barcode sequence located in 1:16 of R1 reads.

droptools/TagBamWithReadSequenceExtended \
SUMMARY=addBarcode.log \
BASE_RANGE=1-6:22-27:43-48 \
DISCARD_READ=false \
BARCODED_READ=1 \
TAG_NAME=CB \
NUM_BASES_BELOW_QUALITY=10 \
INPUT=sample_unmapped.bam \
OUTPUT=sample_BC.bam

```

 - AddPolyTInfo (rule.addToanother)

```shell

# Add polyT length information into bam tags. 
# The polyT length was used to infer the accurate polyadenylation sites.

python script/manibam.py \
-i sample_BC.bam \
-o sample_BC_dT.bam \
-b 54

```

 - TagBamWithReadSequenceExtended (rule.addUMI)

```shell

# Add UMI sequence into tags of unmapped bam.
# In Microwell-seq, the UMI sequence located in 49-54 of R1 reads.
# In 10X genomics, the UMI sequence located in 17:26 of R1 reads.

droptools/TagBamWithReadSequenceExtended \
SUMMARY=addBarcode.log \
BASE_RANGE=49-54 \
DISCARD_READ=true \
BARCODED_READ=1 \
TAG_NAME=CB \
NUM_BASES_BELOW_QUALITY=10 \
INPUT=sample_BC_dT.bam \
OUTPUT=sample_BC_dT_UMI.bam

```

### 4. Alignment

All commands were stored in `rule/alignment/alignment.rule`

```shell

# Convert the unmmaped bam into fq and align to reference genome

java -Djava.io.tmpdir=tmp \
-Xmx20G -jar picard.jar SamToFastq \
INPUT=sample_BC_dT_UMI.bam FASTQ=/dev/stdout | STAR \
--genomeDir star_index \
--alignMatesGapMax 5000 --runThreadN 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes All \
--readFilesIn /dev/stdin \
--sjdbGTFfile gtf \
--outFileNamePrefix alignment/sample_ \
--outFilterScoreMinOverLread 0.4 \
--outFilterMatchNminOverLread 0.4 \
--limitBAMsortRAM 70000000000

# Merge unmapped and mapped file to add the barcode, umi and dT tag into mapped file.

java -Djava.io.tmpdir=tmp \
-Xmx20G -jar picard.jar \
MergeBamAlignment \
REFERENCE_SEQUENCE=genome.fa \
UNMAPPED_BAM=sample_BC_dT_UMI.bam \
ALIGNED_BAM=alignment/sample_Aligned.sortedByCoord.out.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false OUTPUT=sample.bam

# index the bam file

samtools index sample.bam

```

### 5. Run SCAPE

All commands were stored in `rule/apa/scape.rule`

```shell

python scape/main.py apamix \
--bed utr.bed \
--bam sample.bam \
--out sample/ \
--cores 12 \
--cb sample.tsv

```