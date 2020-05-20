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


pair_com_dic = {}
normal_use = []
for _, df in samples.groupby(['Patient']):
    for i in df.Label[df.Type!='N']:
        normal_label = df.Label[df.Type=='N']
        if normal_label.empty:
            pair_com_dic[i] = i
        else:
            pair_com_dic[i] = normal_label.values.tolist()[0]
            normal_use.append(i)


output_dir = getpath(config["output"])
FASTQdir = getpath(config["FASTQ"])

rule all:
	input:
		output_dir + "MultiQC/multiqc_report.html",
		output_dir + "mutect2/PoN/PoN.vcf.gz"


include: "rule/QC/FastQC.rule"
include: "rule/QC/MultiQC.rule"
include: "rule/Alignment/BWA.rule"
include: "rule/BamProc/MergeAndMark.rule"
include: "rule/BamProc/BQSR.rule"
include: "rule/BamProc/VQSR.rule"
include: "rule/Mutect2/mutect.rule"