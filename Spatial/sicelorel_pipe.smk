import os
import sys
import re
import pandas as pd

# sicelorel pipeline

## Configuration file
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("Make sure there is a config.yaml file in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")


if not os.path.isfile(config["metatxt"]):
	sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

samples = pd.read_csv(config["metatxt"], sep='\t').label.values.tolist()


def getpath(str):
	if str in ['', '.', './']:
		return ''
	if str.startswith('./'):
		regex = re.compile('^\./?')
		str = regex.sub('', str)
	if not str.endswith('/'):
		str += '/'
	return str

output_dir = getpath(config['output'])
FASTQdir = getpath(config['FASTQ'])



rule all:
	input:
		expand(output_dir + 'sicelorel/{sample}/{sample}_cellmetrics.txt', sample = samples)


rule barcodeProcess:
	input:
		output_dir + 'ranger/{sample}/barcodes.tsv.gz'
	output:
		output_dir + 'ranger/{sample}/barcodes.tsv'
	shell:
		"""

		zcat {input} | cut -d "-" -f 1 > {output}

		"""

rule IlluminaParser:
	input:
		bam = output_dir + 'ranger/{sample}/possorted_genome_bam.bam',
		bc = output_dir + 'ranger/{sample}/barcodes.tsv'
	output:
		output_dir + 'sicelorel/{sample}/{sample}.illumina.bam.obj'
	params:
		java = config['java']['java'],
		parser_jar = config['soft']['parser'],
		bc = config['bc'],
		gn = config['gn'],
		umi = config['umi'],
		xmx = config['java']['xmx'],
		xms = config['java']['xms']

	shell:
		"""

		{params.java} -Xmx{params.xmx} -Xms{params.xms} \
		-jar {params.parser_jar} \
		-i {input.bam} \
		-o {output} \
		-t {input.bc} -b {params.bc} \
		-g {params.gn} -u {params.umi}

		"""

rule NanoScan:
	input:
		FASTQdir + '{sample}' + '.' + str(config['fqsuffix'])
	output:
		output_dir + 'sicelorel/{sample}/{sample}FWD.fastq'
	params:
		outdir = output_dir + "sicelorel/{sample}/",
		java = config['java']['java'],
		nanoscan = config['soft']['nanoscan']
	shell:
		"""

		{params.java}  -jar {params.nanoscan} \
		-i {input} \
		-o {params.outdir}

		"""

rule minimap:
	input:
		rules.NanoScan.output
	output:
		output_dir + 'sicelorel/{sample}/{sample}.bam'
	params:
		minimap2 = config['soft']['minimap2'],
		samtools = config['soft']['samtools'],
		junc_bed = config['genome']['juncbed'],
		fa = config['genome']['fa']

	shell:
		"""
		
		{params.minimap2} -ax splice -uf --MD --sam-hit-only \
		-t 4 --junc-bed {params.junc_bed} {params.fa} \
		{input} | {params.samtools} sort - -b -o {output}

		"""

rule minimap_index:
	input:
		rules.minimap.output
	output:
		output_dir + 'sicelorel/{sample}/{sample}.bam.bai'
	params:
		samtools = config['soft']['samtools']

	shell:
		"""

		{params.samtools} index {input}

		"""

rule addGeneTag:
	input:
		bai = rules.minimap_index.output,
		bam = rules.minimap.output
	output:
		output_dir + 'sicelore/{sample}/GE.bam'
	params:
		sicelore = config['soft']['sicelore'],
		java = config['java']['java'],
		xmx = config['java']['mostxms'],
		refflat = config['genome']['refFlat'],
		samtools = config['soft']['samtools']

	shell:
		"""

		{params.java} -jar -Xmx{params.xmx} {params.sicelore} \
		AddGeneNameTag \
		I={input.bam} \
		O={output} \
		REFFLAT={params.refflat} \
		GENETAG=GE ALLOW_MULTI_GENE_READS=true \
		USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT

		"""

rule addGeneTag_index:
	input:
		rules.addGeneTag.output
	output:
		output_dir + 'sicelore/{sample}/GE.bam.bai'
	params:
		samtools = config['soft']['samtools']
	shell:
		"""

		{params.samtools} index {input}

		"""

rule AddBamReadSequenceTag:
	input:
		bam = rules.addGeneTag.output,
		bai = rules.addGeneTag_index.output,
		fq = rules.NanoScan.output
	output:
		output_dir + "sicelore/{sample}/GEUS.bam"
	params:
		sicelore = config['soft']['sicelore'],
		java = config['java']['java'],
		xmx = config['java']['mostxms'],
		refflat = config['genome']['refFlat'],
		samtools = config['soft']['samtools']

	shell:
		"""

		{params.java} -jar -Xmx{params.xmx} {params.sicelore} \
		AddBamReadSequenceTag I={input.bam} \
		O={output} FASTQ={input.fq} \
		VALIDATION_STRINGENCY=SILENT

		"""

rule AddBamReadSequenceTag_index:
	input:
		rules.AddBamReadSequenceTag.output
	output:
		output_dir + "sicelore/{sample}/GEUS.bam.bai"
	params:
		samtools = config['soft']['samtools']
	shell:
		"""

		{params.samtools} index {input}

		"""


rule NanoUmi:
	input:
		bai = rules.AddBamReadSequenceTag_index.output,
		bam = rules.AddBamReadSequenceTag.output,
		obj = rules.IlluminaParser.output

	output:
		bam_1 = output_dir + "sicelore/{sample}/GEUS10xAttributes.bam",
		bam_2 = output_dir + "sicelore/{sample}/GEUS10xAttributes_umifound_.bam"
	params:
		log = output_dir + 'sicelore/{sample}/out.log',
		bu_finder = config['soft']['bu_finder'],
		cores = config['cores'],
		java = config['java']['java']

	shell:
		"""
		{params.java} -jar -Xmx100g -Xms30g \
		{params.bu_finder} \
		-i {input.bam} \
		-o {output.bam_1} \
		-k {input.obj} \
		--ncpu {params.cores} --maxUMIfalseMatchPercent 1 \
		--maxBCfalseMatchPercent 5 --logFile {params.log}

		"""

rule NanoUmi_index:
	input:
		bam_1 = rules.NanoUmi.output.bam_1,
		bam_2 = rules.NanoUmi.output.bam_2
	output:
		bam_1 = output_dir + "sicelore/{sample}/GEUS10xAttributes.bam.bai",
		bam_2 = output_dir + "sicelore/{sample}/GEUS10xAttributes_umifound_.bam.bai"
	params:
		samtools = config['soft']['samtools']
	shell:
		"""

		{params.samtools} index {input.bam_1} | {params.samtools} index {input.bam_2}

		"""

rule computeConsensus:
	input:
		bai = rules.NanoUmi_index.output.bam_2,
		bam = rules.NanoUmi.output.bam_2,
	output:
		output_dir + "sicelore/{sample}/consensus.fq"
	params:
		java = config['java']['java'],
		xmx =config['java']['mostxms'],
		sicelore = config['soft']['sicelore'],
		cores = config['cores'],
		tmp = config['tmpdir']

	shell:
		"""
		{params.java} -jar -Xmx{params.xmx} \
		{params.sicelore} ComputeConsensus \
		T={params.cores} I={input} \
		O={output} TMPDIR={params.tmp}

		"""

rule moleculesMap:
	input:
		rules.computeConsensus.output
	output:
		output_dir + 'sicelorel/{sample}/molecule.bam'
	params:
		minimap2 = config['soft']['minimap2'],
		samtools = config['soft']['samtools'],
		junc_bed = config['genome']['juncbed'],
		fa = config['genome']['fa'],
		cores = config['cores']

	shell:
		"""
		{params.minimap2} -ax splice \
		-uf --MD --sam-hit-only -t {params.cores} \
		--junc-bed {params.junc_bed} \
		{params.fa} \
		{input} | {params.samtools} sort - -b -o {output}

		"""

rule moleculesMap_index:
	input:
		rules.moleculesMap.output
	output:
		output_dir + 'sicelorel/{sample}/molecule.bam.bai'
	params:
		samtools = config['soft']['samtools']
	shell:
		"""

		{params.samtools} index {input}

		"""

rule addBC_UMI:
	input:
		bai = rules.moleculesMap_index.output,
		bam = rules.moleculesMap.output

	output:
		output_dir + 'sicelorel/{sample}/molecule.tags.bam'
	params:
		java = config['java']['java'],
		xmx =config['java']['mostxms'],
		sicelore = config['soft']['sicelore']

	shell:
		"""

		{params.java} -jar -Xmx{params.xmx} {params.sicelore} \
		AddBamMoleculeTags I={input.bam} \
		O={output}

		"""

rule addBC_UMI_index:
	input:
		rules.addBC_UMI.output
	output:
		output_dir + 'sicelorel/{sample}/molecule.tags.bam.bai'
	params:
		samtools = config['soft']['samtools']
	shell:
		"""

		{params.samtools} index {input}

		"""

rule addGeNameTag:
	input:
		bai = rules.addBC_UMI_index.output,
		bam = rules.addBC_UMI.output
	output:
		output_dir + 'sicelorel/{sample}/molecule.tags.GE.bam'
	params:
		java = config['java']['java'],
		xmx =config['java']['mostxms'],
		sicelore = config['soft']['sicelore'],
		refFlat = config['genome']['refFlat'],
		gn = config['gn']

	shell:
		"""

		{params.java} -jar -Xmx{params.xmx} {params.sicelore} \
		AddGeneNameTag I={input.bam} \
		O={output} \
		REFFLAT={params.refFlat} \
		GENETAG={params.gn} ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true \
		VALIDATION_STRINGENCY=SILENT

		"""

rule addGeNameTag_index:
	input:
		rules.addGeNameTag.output
	output:
		output_dir + 'sicelorel/{sample}/molecule.tags.GE.bam.bai'
	params:
		samtools = config['soft']['samtools']
	shell:
		"""

		{params.samtools} index {input}

		"""

rule final:
	input:
		bai = rules.addGeNameTag_index.output,
		bam = rules.addGeNameTag.output,
		barcode = output_dir + 'ranger/{sample}/barcodes.tsv'
	output:
		output_dir + 'sicelorel/{sample}/{sample}_cellmetrics.txt'
	params:
		java = config['java']['java'],
		xmx =config['java']['mostxms'],
		sicelore = config['soft']['sicelore'],
		gn = config['gn'],
		refFlat = config['genome']['refFlat'],
		outdir = output_dir + 'sicelorel/{sample}/',
		prefix = '{sample}'

	shell:
		"""
		{params.java} -jar -Xmx{params.xmx} \
		{params.sicelore} IsoformMatrix \
		DELTA=2 METHOD=STRICT ISOBAM=true GENETAG={params.gn} \
		I={input.bam} \
		REFFLAT={params.refFlat} \
		CSV={input.barcode} OUTDIR={params.outdir} \
		PREFIX={params.prefix} VALIDATION_STRINGENCY=SILENT
		"""