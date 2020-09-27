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

meta = pd.read_csv(config["metatxt"], sep='\t')

samples = meta.SampleLabel.values.tolist()

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

pairs = ['R1', 'R2']


rule all:
    input:
        expand(output_dir + "Metrics/{sample}.RNA_metrics.txt", sample = samples),
        expand(output_dir + "rseqc/{sample}.txt", sample = samples),
        expand(output_dir + 'genebody/{sample}.geneBodyCoverage.txt', sample = samples),
        expand(output_dir + "fastqc_filter/{sample}_{pair}.txt", sample=samples, pair=pairs),
        expand(output_dir + "fastqc_raw/{sample}_{pair}.txt", sample=samples, pair=pairs),
        expand(output_dir + 'q30/{sample}_{pair}.txt', sample=samples, pair=pairs),
        expand(output_dir + 'TIN/{sample}.Aligned.sortedByCoord.out.summary.txt', sample=samples),
        output_dir + "fastqc_filter/fastqc.Rds",
        output_dir + "fastqc_raw/fastqc.Rds"

rule fastqc_raw:
    input:
        fastq = "fq_raw/" + "{sample}_" + "{pair}." + str(config["fqsuffix"]) + ".gz"
    output:
        output_dir + "fastqc_raw/{sample}_{pair}_fastqc.zip"
    params:
        outdir = output_dir + "fastqc_raw",
        FastQC = config["soft"]["FastQC"]
    log:
        output_dir + "logs/fastqc_raw_{sample}_{pair}.log"
    shell:
        """
        {params.FastQC} -o {params.outdir} {input.fastq} > {log}
        """

rule fetch_fastqc_raw:
    input:
        output_dir + "fastqc_raw/{sample}_{pair}_fastqc.zip"
    output:
        output_dir + "fastqc_raw/{sample}_{pair}.txt"
    params:
        output_dir + "fastqc_raw/{sample}_{pair}_fastqc",
        output_dir + "fastqc_raw/"
    shell:
        """

        unzip -o -q {input} -d {params[1]} && sed -n "13,51p" {params[0]}/fastqc_data.txt > {output}

        """

rule fastqc_raw_collpase:
    input:
        expand(output_dir + "fastqc_raw/{sample}_{pair}.txt", sample=samples, pair=pairs)
    output:
        output_dir + "fastqc_raw/fastqc.Rds"
    params:
        expand(output_dir + "fastqc_raw/{sample}_{pair}", sample=samples, pair=pairs)
    shell:
        """
        Rscript bin/collapse.R {input} {output}
        """

rule fastqc:
    input:
        fastq = "fq/" + "{sample}_" + "{pair}." + str(config["fqsuffix"]) + ".gz"
    output:
        output_dir + "fastqc_filter/{sample}_{pair}_fastqc.zip"
    params:
        outdir = output_dir + "fastqc_filter",
        FastQC = config["soft"]["FastQC"]
    log:
        output_dir + "logs/fastqc_filter_{sample}_{pair}.log"
    shell:
        """

        {params.FastQC} -o {params.outdir} {input.fastq} > {log}

        """

rule fetch_fastqc:
    input:
        output_dir + "fastqc_filter/{sample}_{pair}_fastqc.zip"
    output:
        output_dir + "fastqc_filter/{sample}_{pair}.txt"
    params:
        output_dir + "fastqc_raw/{sample}_{pair}_fastqc",
        output_dir + "fastqc_raw/"
    shell:
        """


        unzip -o -q {input} -d {params[1]} && sed -n "13,51p" {params[0]}/fastqc_data.txt > {output}

        """

rule fastqc_filter_collpase:
    input:
        expand(output_dir + "fastqc_filter/{sample}_{pair}.txt", sample=samples, pair=pairs)
    output:
        output_dir + "fastqc_filter/fastqc.Rds"

    shell:
        """
        Rscript bin/collapse.R {input} {output}
        """



rule q30:
    input:
        fastq = "fq/" + "{sample}_" + "{pair}." + str(config["fqsuffix"]) + ".gz"
    output:
        output_dir + 'q30/{sample}_{pair}.txt'

    shell:
        """

        python2 bin/q30.py {input} > {output}


        """






rule GeneBodyQC:
    input:
        output_dir + 'STAR/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        output_dir + 'Metrics/{sample}.RNA_metrics.txt'
    params:
        picard = config['soft']['picard'],
        refFlat = config['refFlat'],
        strand = 'SECOND_READ_TRANSCRIPTION_STRAND' if config['strand'] else 'NONE'
    shell:
        """

        java -jar {params.picard} CollectRnaSeqMetrics \
        I={input} \
        O={output} \
        REF_FLAT={params.refFlat} \
        STRAND={params.strand}

        """

rule read_dist:
    input:
        output_dir + 'STAR/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        output_dir + 'rseqc/{sample}.txt'
    params:
        read_dist = config['soft']['read_dist'],
        bed12 = config['bed12']

    shell:
        """

        {params.read_dist} -i {input} -r {params.bed12} > {output}


        """


rule coverage_plot:
    input:
        output_dir + 'STAR/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        output_dir + 'genebody/{sample}.geneBodyCoverage.txt'
    params:
        genebody = config['soft']['genebody'],
        bed12 = config['bed12'],
        prefix = output_dir + 'genebody/{sample}'

    shell:
        """

        {params.genebody} -i {input} -r {params.bed12} -o {params.prefix}


        """
rule tin_calc:
    input:
        output_dir + 'STAR/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        output_dir + 'TIN/{sample}.Aligned.sortedByCoord.out.summary.txt'
    params:
        tin = config['soft']['tin']
        bed12 = config['bed12'],
        prefix = output_dir + 'TIN/'

    shell:
        """
        {params.tin} -i {input} -r {params.bed12} -o {output}
        """
