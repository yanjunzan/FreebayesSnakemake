configfile: "./Tobacco-freebayes.yaml"

###########################
# define sample snakemake -npr --use-conda  snakemake -s Snakefile --jobs 2
###########################

#df = pd.read_csv(config["sra_list"],sep="\t")
#df["Run"].tolist()


WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]
REF_dir = config["refdir"]
genome = config["ref"]


SAMPLES = [f for f in os.listdir(WORKING_DIR+"fastq") if re.search("FDSW",f)] 

###########################
# Input functions for rules
###########################

rule all:
    input:
        vcf = WORKING_DIR + "vcf/all.vcf",
        reffai = REF_dir + genome + ".fasta.fai",
        bamlist= WORKING_DIR + "mapped/BAMList.txt",
        outbam = expand(WORKING_DIR + "mapped/{sample}_sort_index_rg_nodup.bam",sample=SAMPLES),
        samstats = expand(WORKING_DIR + "bamstats/{sample}.bamstats.stat",sample=SAMPLES),
        #sort_index_rgbam = expand(WORKING_DIR + "mapped/{sample}_sort_index_rg.bam",sample=SAMPLES),
        genomeF = REF_dir + genome + ".bwt",
        qcfile = expand(RESULT_DIR + "fastp/{sample}.html",sample=SAMPLES),
        fq1 = expand(WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",sample = SAMPLES),
        fq2 = expand(WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",sample = SAMPLES),
        #bam = expand(WORKING_DIR + "mapped/{sample}_sort_index.bam", sample = SAMPLES)
    message:
        "Job done! Removing temporary directory"

rule bwa_index:
    #conda:
    #    "./Tobsoft.yaml"
    input:
        REF_dir + genome + ".fasta"
    output:
        REF_dir + genome + ".amb",
        REF_dir + genome + ".ann",
        REF_dir + genome + ".bwt",
        REF_dir + genome + ".pac",
        REF_dir + genome + ".sa"
    log:
        "./logs/bwa_index/" + genome + ".log"
    params:
        prefix = REF_dir + genome,
        algorithm="bwtsw"
    shell:
        '''
        bwa index -p {params.prefix} -a {params.algorithm} {input} {log}
        '''
rule fastp:
    input:
        fw = WORKING_DIR + "fastq/{sample}/{sample}_DS_1.clean.fq.gz",
        rev= WORKING_DIR + "fastq/{sample}/{sample}_DS_1.clean.fq.gz"
    output:
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads to {output.fq1}"
    threads: 1
    log:
        WORKING_DIR + "logs/fastp/{sample}.log.txt"
    params:
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    shell:
        "fastp --thread {threads}  --html {output.html} \
        --qualified_quality_phred {params.qualified_quality_phred} \
        --in1 {input.fw} --in2 {input.rev} --out1 {output.fq1} --out2 {output.fq2}; \
        2> {log}"
        
rule bwamem_map:
    input:
        #ref = REF_dir + genome + ".fasta",
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    params:
        threads = 2,
        ref = REF_dir + genome
    output:
        #sort_index_bam = WORKING_DIR + "mapped/{sample}_sort_index.bam",
        sort_bam = temp(WORKING_DIR + "mapped/{sample}_sort_index.bam")
    shell:
        "bwa mem  -t {params.threads} {params.ref} {input.fq1} {input.fq2} | samtools sort -O BAM -@ {params.threads}  -o \
        {output.sort_bam} "
        " && samtools index -@ {params.threads} {output.sort_bam}"

rule addReadGroup:
    input:
        picard = config["picard"],
        sort_index_bam = WORKING_DIR + "mapped/{sample}_sort_index.bam"
    output:
        sort_index_rgbam = temp(WORKING_DIR + "mapped/{sample}_sort_index_rg.bam")
    params:
        sample = "{sample}"
    shell:
        "java -jar {input.picard} AddOrReplaceReadGroups \
        I={input.sort_index_bam}  O={output.sort_index_rgbam} RGID={params.sample} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM={params.sample}"
        
rule samtools_stat:
    input:
        bam = WORKING_DIR + "mapped/{sample}_sort_index_rg.bam"
    output:
        samstats = WORKING_DIR + "bamstats/{sample}.bamstats.stat"
    params:
        threads=2
    shell:
        "samtools stats -@ {params.threads} {input.bam} > {output.samstats}"
        
rule picard_remove_duplicates:
    input:
        picard = config["picard"],
        bam = WORKING_DIR + "mapped/{sample}_sort_index_rg.bam"
    output:
        outbam = protected(WORKING_DIR + "mapped/{sample}_sort_index_rg_nodup.bam"),
        metrics = WORKING_DIR + "dupstats/{sample}.picard.marked_dup_metrics.txt"
    shell:
        "java -jar {input.picard}  MarkDuplicates -I {input.bam} -O {output.outbam} -M {output.metrics} --REMOVE_DUPLICATES true"
        

rule generateBAMList:
    params:
        bam = WORKING_DIR + "mapped/",
        suff= "_sort_index_rg_nodup.bam"
    output:
        bamlist= WORKING_DIR + "mapped/BAMList.txt"
    #shell:
        #"ls {params.bam} | grep {params.suff} > {output.bamlist}"
        #" && for i in `cat {output.bamlist}`; do echo {params.bam}"$i"; done"
    run:
        import os
        import re
        F = [f for f in os.listdir(params.bam) if re.search(params.suff,f)]
        out = open(output[0],"w")
        for f in F:
            print(params.bam+f,file=out)
        out.close()
        
rule samtools_faidx:
    output:
        reffai = REF_dir + genome + ".fasta.fai"
    input:
        ref = REF_dir + genome + ".fasta"
    shell:
        "samtools faidx {input.ref}"
#        
#rule freeb:
#    input:
#        bamlist= WORKING_DIR + "mapped/BAMList.txt",
#        ref = REF_dir + genome + ".fasta"
#    output:
#        vcf = WORKING_DIR + "vcf/all.vcf"
#    shell:
#        "freebayes -f {input.ref} -L {input.bamlist} > {output.vcf}"
#        
rule freeb_pa:
    input:
        bamlist= WORKING_DIR + "mapped/BAMList.txt",
        ref = REF_dir + genome + ".fasta",
        fai = REF_dir + genome + ".fasta.fai"
    output:
        vcf = WORKING_DIR + "vcf/all.vcf"
    params:
        threads=4
    shell:
        "freebayes-parallel <(fasta_generate_regions.py {input.fai} 100000) {params.threads} -f {input.ref}\
        -L {input.bamlist} > {output.vcf}"


