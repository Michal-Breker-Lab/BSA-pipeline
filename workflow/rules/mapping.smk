rule bwa_mem:
    input:
        reads=rules.fastp_pe.output.trimmed,
        idx=rules.bwa_index_ref.output.idx,
    output:
        bam="results/mapping/{sample}.bam",
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",
        sort_order="coordinate", 
        sort_extra="",
    threads: 8
    wrapper:
        "v7.2.0/bio/bwa/mem"


rule removeduplicates_bam:
    input:
        bams=rules.bwa_mem.output.bam
    output:
        bam="results/mapping/{sample}.dedup.bam",
        metrics="results/qc/picard_mark_dup/{sample}.metrics.txt"
    params:
        extra="--REMOVE_DUPLICATES true"
    log:
        "logs/removeduplicates_bam/{sample}.log"
    wrapper:
        "v7.2.0/bio/picard/markduplicates"

rule samtools_index:
    input:
        rules.bwa_mem.output.bam
    output:
        "results/mapping/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",
    threads: 4
    wrapper:
        "v7.2.0/bio/samtools/index"

rule samtools_index_markdup:
    input:
        rules.removeduplicates_bam.output.bam
    output:
        "results/mapping/{sample}.dedup.bam.bai"
    log:
        "logs/samtools_index_markdup/{sample}.log",
    params:
        extra="",
    threads: 4
    wrapper:
        "v7.2.0/bio/samtools/index"
