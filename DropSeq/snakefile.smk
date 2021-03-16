import os
import glob
from itertools import islice


samples = [os.path.basename(i).replace('.bam', '') for i in glob.glob('data/bam/*.bam')]


BARCODE = '1-12'
UMI = '13-20'

GTF = '/mnt/raid61/Personal_data/zhouran/reference/STAR/Homo_sapiens/Homo_sapiens.GRCh38.93.sorted.gtf'
GENOME = '/mnt/raid61/Personal_data/zhouran/reference/STAR/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
INDEX = '/mnt/raid61/Personal_data/zhouran/reference/STAR/Homo_sapiens/release93/STAR_index_2.7.1a'



picard = '/home/zhouran/data/soft/DropSeqTools/Drop-seq_tools-1.13/3rdParty/picard/picard.jar'
droptools= '/home/zhouran/data/soft/DropSeqTools/Drop-seq_tools-1.13'
Tmp =  './tmp'

rule all:
	input:
		expand('data/{sample}/{sample}.merge.bam.bai', sample = samples)


rule addBarcode:
	input:'data/bam/{sample}.bam'
	output:
		temp('data/{sample}/{sample}.barcode.bam')
	log:'data/{sample}/{sample}.barcode.log'
	shell:
		"{droptools}/TagBamWithReadSequenceExtended SUMMARY={log} \
		BASE_RANGE={BARCODE} DISCARD_READ=false BARCODED_READ=1 TAG_NAME=CB NUM_BASES_BELOW_QUALITY=10 INPUT={input} \
		OUTPUT={output}"

rule addUMI:
	input:'data/{sample}/{sample}.barcode.bam'
	output:
		temp('data/{sample}/{sample}.barcode.UMI.bam')
	log:
		'data/{sample}/{sample}.UMI.log'
	shell:
		"{droptools}/TagBamWithReadSequenceExtended SUMMARY={log} \
		BASE_RANGE={UMI} DISCARD_READ=true BARCODED_READ=1 TAG_NAME=UB NUM_BASES_BELOW_QUALITY=10 INPUT={input} \
		OUTPUT={output}"


rule SamToFastq_align:
	input:
		'data/{sample}/{sample}.barcode.UMI.bam'
	output:
		'data/{sample}/STAR/{sample}.Aligned.sortedByCoord.out.bam'
	params:
		'data/{sample}/STAR/{sample}.'
	log:
		'data/{sample}/{sample}.SamToFastq.log'
	shell:
		"java -Djava.io.tmpdir={Tmp} -Xmx20G -jar {picard} SamToFastq \
		INPUT={input} FASTQ=/dev/stdout | STAR --genomeDir {INDEX} \
            --alignMatesGapMax 5000 --runThreadN 4 --outSAMtype BAM SortedByCoordinate --outSAMattributes All \
            --readFilesIn /dev/stdin\
            --sjdbGTFfile {GTF}\
            --outFileNamePrefix {params} --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 \
            --limitBAMsortRAM 70000000000"



rule MergeBam:
	input:
		'data/{sample}/{sample}.barcode.UMI.bam', 'data/{sample}/STAR/{sample}.Aligned.sortedByCoord.out.bam'
	output:
		bam = 'data/{sample}/{sample}.merge.bam'
	log:'data/{sample}/{sample}.MergeBam.log'
	shell:
		"java -Djava.io.tmpdir={Tmp} -Xmx20G -jar {picard} \
		MergeBamAlignment REFERENCE_SEQUENCE={GENOME} UNMAPPED_BAM={input[0]} \
		ALIGNED_BAM={input[1]} INCLUDE_SECONDARY_ALIGNMENTS=false OUTPUT={output.bam} 2>{log}"

rule indexbam:
	input:
		rules.MergeBam.output.bam
	output:
		rules.MergeBam.output.bam + '.bai'
	shell:
		"samtools index {input}"

