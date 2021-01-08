import os
import sys
import re
import pandas as pd

## Configuration file
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("Make sure there is a config.yaml file in " + \
			os.getcwd() + \
			" or specify one with the --configfile commandline parameter.")

## Read metadata
if not os.path.isfile(config["metatxt"]):
	sys.exit("Metadata file " + config["metatxt"] + " does not exist.")

samples = pd.read_csv(config["metatxt"], sep="\t")


## Sanitize provided input and output directories
def getpath(str):
	if str in ["", ".", "./"]:
		return ""
	if str.startswith("./"):
		regex = re.compile("^\./?")
		str = regex.sub("", str)
	if not str.endswith("/"):
		str += "/"
	return str

output_dir = getpath(config["output"])

rule all:
	input:
		expand(output_dir + "{sample}_1.fastq.gz", sample = samples.Accession.values.tolist())


rule download:
	output:
		temp(output_dir + "{sample}.sra")
	params:
		prefix = output_dir,
		prefetch = config['soft']['prefetch'],
		ascp = config['soft']['ascp']
	priority: 50
	shell:
		"""
		{params.prefetch} -t ascp -a "{params.ascp}" {wildcards.sample} --max-size 200G -O {params.prefix}
		"""
rule unzip:
	input:
		rules.download.output
	output:
		output_dir + "{sample}_1.fastq.gz"
	params:
		fq_dump = config['soft']['fq_dump'],
		prefix = output_dir
	shell:
		"""

		{params.fq_dump} --split-files --origfmt --gzip {input} -O {params.prefix}

		"""



