

## WES

All pipelien followed the GATK best practice

All sample must be paired, (normal vs ca).

### QC

Include:
 - FastQC
 - picard
 - multiQC

### Trim (optional)

Inlucde:
 - fastp (NOT used)
 - Trimmomatic (NOT used)

### BWA

All config information was in `config` file, yaml format.

### GATK

#### Germline mutation

- BQSR

- Filter mutaiton with hard criterion

- Mutation annotation (VEP)

- PCA

#### Somatic mutation

- Mutect2

- Mutation annotaiton (VEP)

- vcf2maf and Maftools analysis

### CNV analysis

- CNVkitv2

- gistic2

- maftools analysis

### MultiQC

generate the html report

