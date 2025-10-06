rule fastp_pe:
    input:
        sample=get_raw_reads,
    output:
        trimmed=temp(
            [
                "results/trimming/{sample}/{sample}_R1.fastq.gz",
                "results/trimming/{sample}/{sample}_R2.fastq.gz",
            ]
        ),
        unpaired=temp("results/trimming/{sample}/{sample}.singletons.fastq.gz"),
        html="results/qc/trimming/{sample}.html",
        json="results/qc/trimming/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    params:
        extra=config["params"]["fastp"],
    threads: 2
    wrapper:
        "v7.2.0/bio/fastp"
