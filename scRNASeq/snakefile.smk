import os
import sys
import re
import pandas as pd
import numpy as np

## Configuration file
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")

## Read metadata
if not os.path.isfile(config["metatxt"]):
	sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

samples = pd.read_csv(config["metatxt"], sep="\t")

try:
	_ = samples.version
except AttributeError:
	print("There was no library version, please add it")
	sys.exit(0)


version_check = samples.groupby("label")["version"].apply(set).to_dict()

# check version value whether unique 
if np.any(np.array(list(map(len, version_check.values()))) != 1):
	sys.exit("10X version value not unique, please chect it")

# check version value whether correct.
pool_version = np.array(['v2','v3'])
if not np.all(np.isin(np.array([i for set_info in version_check.values() for i in set_info]), pool_version)):
	sys.exit("10X version value was not v2 or v3, please chect it")

version_check = pd.Series(samples.version.values.tolist(), index = samples.label).to_dict()

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


rule all:
	input:
		expand(output_dir + "STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam", sample = set(samples.label.values.tolist()))


rule STARsolo:
	"""
	STARsolo to align single cell or single nuceli data
	specify --soloFeatures Gene for single-cell  data
	specify --soloFeatures GeneFull for single-nuclei data
	specify --soloFeatures Gene GeneFull for getting both counts in exons level and exon + intron level (velocity)
	"""
	output: 
		bam = output_dir + "STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
		bai = output_dir + "STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
		rawmtx = output_dir + "STAR/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx",
		feature = output_dir + "STAR/{sample}/{sample}.Solo.out/Gene/raw/features.tsv",
		barcode = output_dir + "STAR/{sample}/{sample}.Solo.out/Gene/raw/barcodes.tsv"

	params:
		r2 = lambda wildcards: ','.join(FASTQdir + f"{sra}_" + str(config["fqext2"]) + "." + str(config["fqsuffix"]) + ".gz" for sra in samples[samples.label == wildcards.sample].sra.values.tolist()),
		r1 = lambda wildcards: ','.join(FASTQdir + f"{sra}_" + str(config["fqext1"]) + "." + str(config["fqsuffix"]) + ".gz" for sra in samples[samples.label == wildcards.sample].sra.values.tolist()),
		# STAR solo
		index = config["STAR"]["index"],
		gtf = config["STAR"]["gtf"],
		STAR = config["soft"]["STAR"],
		threads = config["STAR"]["cpus"],
		star_custom = config["soloFeatures"],
		outprefix = output_dir +  "STAR/{sample}/{sample}.",
		whitelist = lambda wildcards: config["barcode_v2"]["whitelist"] if version_check[wildcards.sample]=='v2' else config["barcode_v3"]["whitelist"],
		barcodestart = lambda wildcards: config["barcode_v2"]["barcodestart"] if version_check[wildcards.sample]=='v2' else config["barcode_v3"]["barcodestart"],
		barcodelength = lambda wildcards: config["barcode_v2"]["barcodelength"] if version_check[wildcards.sample]=='v2' else config["barcode_v3"]["barcodelength"],
		umistart = lambda wildcards: config["barcode_v2"]["umistart"] if version_check[wildcards.sample]=='v2' else config["barcode_v3"]["umistart"],
		umilength = lambda wildcards: config["barcode_v2"]["umilength"] if version_check[wildcards.sample]=='v2' else config["barcode_v3"]["umilength"],
		samtools = config["soft"]["samtools"],
		soloBarcodeReadLength = config.get('soloBarcodeReadLength', '')

	log:
		output_dir +  "Log/{sample}_STAR.log" 
	benchmark:
		output_dir +  "Benchmark/{sample}_STAR.benchmark"

	shell:
		"""
		{params.STAR} \
			--runMode alignReads \
			--genomeDir {params.index} \
			--sjdbGTFfile {params.gtf} \
			--runThreadN {params.threads} \
			--outFileNamePrefix {params.outprefix} \
			--outSAMtype BAM SortedByCoordinate \
			--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
			--soloType CB_UMI_Simple \
			{params.star_custom} \
			{params.soloBarcodeReadLength} \
			--soloCBwhitelist {params.whitelist} \
			--soloCBstart {params.barcodestart} \
			--soloCBlen {params.barcodelength} \
			--soloUMIstart {params.umistart} \
			--soloUMIlen {params.umilength} \
			--soloCBmatchWLtype 1MM_multi_pseudocounts \
			--soloUMIfiltering MultiGeneUMI \
			--readFilesIn {params.r2} {params.r1} --readFilesCommand zcat \
			> {log} 2>&1

		{params.samtools} index -b -@ {params.threads} {output.bam} >> {log} 2>&1
				
		"""
