import os
import sys
import re
import pandas as pd

## Configuration file
if len(config) == 0:
	if os.path.isfile("./config.yaml"):
		configfile: "./config.yaml"
	else:
		sys.exit("\033[91m Make sure there is a config.yaml file in " + \
		os.getcwd() + \
		" or specify one with the --configfile commandline parameter. \033[0m")

## Read metadata
if not os.path.isfile(config["metatxt"]):
	sys.exit("\033[91m Metadata file " + config["metatxt"] + " does not exist. \033[0m")

samples = pd.read_csv(config["metatxt"], sep='\t')

try:
	_ = samples.Strand
	_ = samples.Length
except AttributeError:
	sys.exit("\033[91mThere was no Strand or Length column in " + 
	"{}, please add it. \033[0m".format(config["metatxt"]))

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
		output_dir + "MultiQC/multiqc_report.html"

include: "rule/QC/FastQC.rule"
include: "rule/alignment/STAR_PE.rule"
include: "rule/QC/MultiQC.rule"
include: "rule/quant/RSEM.rule"
include: "rule/quant/GeneCounts.rule"
include: "rule/quant/ExonPart_cal.rule"
include: "rule/quant/rMATs.rule"
include: "rule/quant/qapa.rule"
include: "rule/quant/DEXseq_count.rule"
include: "rule/QC/GeneBodyQC.rule"
