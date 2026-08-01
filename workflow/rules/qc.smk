rule fastqc_raw:
    input:
        get_fastqc_raw_reads,
    output:
        html="results/qc/fastqc_raw/{sample}_{fq}.html",
        zip="results/qc/fastqc_raw/{sample}_{fq}_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_raw/{sample}_{fq}.log",
    resources:
        mem_mb = 1024,
    wrapper:
        "v9.15.0/bio/fastqc"


rule fastqc_trimmed:
    input:
        "results/trimming/{sample}/{sample}_{fq}.fastq.gz",
    output:
        html="results/qc/fastqc_trimmed/{sample}_{fq}.html",
        zip="results/qc/fastqc_trimmed/{sample}_{fq}_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_trimmed/{sample}_{fq}.log",
    resources:
        mem_mb = 1024,
    wrapper:
        "v9.15.0/bio/fastqc"


rule samtools_stats:
    input:
        bam=get_bam,
    output:
        "results/qc/samtools_stats/{sample}.txt",
    params:
        extra="",
    log:
        "logs/samtools_stats/{sample}.log",
    wrapper:
        "v9.15.0/bio/samtools/stats"


rule samtools_idxstats:
    input:
        bam=get_bam,
        idx=get_bai,
    output:
        "results/qc/samtools_idxstats/{sample}.bam.idxstats",
    log:
        "logs/samtools_idxstats/{sample}.log",
    params:
        extra="",
    wrapper:
        "v9.15.0/bio/samtools/idxstats"


rule samtools_coverage:
    input:
        bam=get_bam,
    output:
        "results/qc/samtools_coverage/{sample}.txt",
    log:
        "logs/samtools_coverage/{sample}.log",
    params:
        extra="",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools coverage {params.extra} -o {output} {input.bam}"


rule insert_size:
    input:
        get_bam,
    output:
        txt="results/qc/picard_insert_size/{sample}.txt",
        pdf=temp("results/qc/picard_insert_size/{sample}.pdf"),
    log:
        "logs/picard_insert_size/{sample}.log",
    params:
        extra="--VALIDATION_STRINGENCY LENIENT --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE",
    wrapper:
        "v9.15.0/bio/picard/collectinsertsizemetrics"


rule multiqc:
    input:
        multiqc_input,
        config=config["multiqc_config"],
    output:
        "results/qc/multiqc_{experiment}.html",
    params:
        extra="--verbose",
        use_input_files_only=True,
    log:
        "logs/multiqc/{experiment}.log",
    wrapper:
        "v9.15.0/bio/multiqc"
