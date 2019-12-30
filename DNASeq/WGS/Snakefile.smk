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

samples = pd.read_csv(config["metatxt"])


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
include: "rule/QC/MultiQC.rule"
include: "rule/alignment/BWA.rule"
include: "rule/BamProc/MergeAndMark.rule"
include: "rule/BamProc/BQSR.rule"
include: "rule/BamProc/VQSR.rule"
