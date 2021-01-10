import os
import sys
import re
import pandas as pd

## Configuration file
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Read metadata
if not os.path.isfile(config["metatxt"]):
	sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

samples = pd.read_csv(config["metatxt"], sep='\t')

try:
	_ = samples.version
except AttributeError:
	print("There was no strand parameter, please add it")
	sys.exit(0)

## Sanitize provided input and output directories
def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

output_dir = getpath(config["output"])
FASTQdir = getpath(config["FASTQ"])


rule STARsolo:
	"""
	STARsolo to align single cell or single nuceli data
	specify --soloFeatures Gene for single-cell  data
	specify --soloFeatures GeneFull for single-nuclei data
	specify --soloFeatures Gene GeneFull for getting both counts in exons level and exon + intron level (velocity)
	"""
	input: 
		mapindex = config["genome"]["mapindex"],
		gtf = config["genome"]["gtf"],
		whitelist = config["barcode"]["whitelist"]
	output: 
		bam = output_dir + "STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		bai = output_dir + "STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
		rawmtx = output_dir + "STAR/{sample}/{sample}Solo.out/Gene/raw/matrix.mtx",
		feature = output_dir + "STAR/{sample}/{sample}Solo.out/Gene/raw/features.tsv",
		barcode = output_dir + "STAR/{sample}/{sample}Solo.out/Gene/raw/barcodes.tsv"
	params:
		star_custom = config.get("STARsolo_custom", ""),
		outprefix = output_dir +  "STAR/{sample}/{sample}",
		transcript = lambda wildcards: ','.join(FILES[wildcards.sample]["R2"]),
		barcode = lambda wildcards: ','.join(FILES[wildcards.sample]["R1"]),
		barcodestart = config["barcode"]["barcodestart"],
		barcodelength = config["barcode"]["barcodelength"],
		umistart = config["barcode"]["umistart"],
		umilength = config["barcode"]["umilength"]
	version: STAR_VERSION        
	log:
		"Result/Log/{sample}_STAR.log" 
	benchmark:
		"Result/Benchmark/{sample}_STAR.benchmark"
	threads:
		config.get("STARsolo_threads", "")

	shell:
		"""
		STAR \
			--runMode alignReads \
			--genomeDir {input.mapindex} \
			--sjdbGTFfile {input.gtf} \
			--runThreadN {threads} \
			--outFileNamePrefix {params.outprefix} \
			--outSAMtype BAM SortedByCoordinate \
			--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
			--soloType CB_UMI_Simple \
			{params.star_custom} \
			--soloCBwhitelist {input.whitelist} \
			--soloCBstart {params.barcodestart} \
			--soloCBlen {params.barcodelength} \
			--soloUMIstart {params.umistart} \
			--soloUMIlen {params.umilength} \
			--soloCBmatchWLtype 1MM_multi_pseudocounts \
			--soloUMIfiltering MultiGeneUMI \
			--readFilesIn {params.transcript} {params.barcode} --readFilesCommand zcat \
			> {log} 2>&1
			
		samtools index -b -@ {threads} {output.bam} >> {log} 2>&1
				
		"""
