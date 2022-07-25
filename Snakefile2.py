import os
import pandas as pd
configfile: "./config.yaml"

###########################
# define sample
###########################

df = pd.read_csv(config["sra_list"],sep="\t")
SAMPLES = df["Run"].tolist()


WORKING_DIR = config["result_dir"]
RESULT_DIR = config["result_dir"]
###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function checks if raeds are single or paired end"""
    if os.stat(sample).st_size < 100:
        return True
    else:
        return False



rule all:
    input:
        #fw  = expand(WORKING_DIR + "fastq/{sample}_1.fastq",sample = SAMPLES),
        #qcfile = expand(RESULT_DIR + "fastp/{sample}.html",sample=SAMPLES),
        fq1 = expand(WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",sample = SAMPLES),
        fq2 = expand(WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",sample = SAMPLES),
        bam = expand(WORKING_DIR + "mapped/{sample}.bam", sample = SAMPLES),
    message:
        "Job done! Removing temporary directory"

rule index_hisat:
    conda:
        "mapping.yaml"
    input:
        config["working_dir"] +"MorexV3/"+ config["ref"]
    output:
        [config["working_dir"]+ "MorexV3/" + config["ref"].replace("fasta","") + str(i) + ".ht2" for i in range(1,9)]
    message:
        "indexing genome {output}"
    params:
        index  = config["working_dir"] + "MorexV3/" + config["ref"].replace("fasta",""),
        genome = config["working_dir"] + "MorexV3/" + config["ref"]
    threads: 1
    shell:
        "hisat2-build -p {threads} {params.genome} {params.index} --quiet"


# rule hisat2_align:
#     input:
#       reads=[WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz", WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"]
#     output:
#       WORKING_DIR + "mapped/{sample}.bam"
#     log:
#         "logs/hisat2_align_{sample}.log"
#     params:
#       extra="",
#       idx=config["working_dir"]+ "MorexV3/",
#     threads: 2
#     wrapper:
#       "0.77.0/bio/hisat2/align"



rule hisat_mapping:
    input:
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        indexFiles = [config["working_dir"]+ "MorexV3/" + config["ref"].replace("fasta","") + str(i) + ".ht2" for i in range(1,9)]
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        met   = RESULT_DIR + "logs/{sample}_met.txt"
    params:
        indexName = config["working_dir"]+ "MorexV3/" + config["ref"].replace(".fasta",""),
        sampleName = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    conda:
        "mapping.yaml"
    message:
        "mapping reads to genome to bam files {params.sampleName}."
    threads: 2
    run:
        if sample_is_single_end(params.sampleName):
            shell("hisat2 -p {threads} --met-file {output.met} -x {params.indexName} \
            -U {input.fq1} > {output.sams}") #| samtools view -Sb -F 4 -o {output.bams}
        else:
            shell("hisat2 -p {threads} --met-file {output.met} -x {params.indexName} \
            -1 {input.fq1} -2 {input.fq2} | samtools view -Sb -F 4 -o {output.bams}")
    # shell:
    #     "hisat2 -p {threads} --met-file {output.met} -x {params.indexName} \
    #     -1 {input.fq1} -2 {input.fq2} | samtools view -Sb -F 4 -o {output.bams}"
