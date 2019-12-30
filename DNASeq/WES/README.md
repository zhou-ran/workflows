

## WES (not complete)

All pipelien followed the GATK best practice

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
- BQSR
- Mutect2
- Annotation

### MultiQC

generate the html report
