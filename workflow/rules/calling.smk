rule freebayes:
    input:
        alns=get_bams_experiment,
        idsx=get_bai_experiment,
        ref=get_ref_fasta,
        fai=rules.samtools_index_ref.output,
    output:
        vcf="results/calling/{experiment}.vcf.gz",
    log:
        "logs/freebayes/{experiment}.log",
    params:
        chunksize=config["params"]["freebayes"]["chunksize"],
        extra=f"--pooled-continuous {config['params']['freebayes']['extra']}",
    threads: 30
    wrapper:
        "v7.2.0/bio/freebayes"


rule norm_vcf:
    input:
        rules.freebayes.output.vcf,
        ref=get_ref_fasta,
    output:
        temp("results/calling/{experiment}.norm.vcf.gz"),
    log:
        "logs/norm_vcf/{experiment}.log",
    wrapper:
        "v7.2.0/bio/bcftools/norm"

rule mask_reference:
    input:
        get_ref_fasta,
    output:
        "resources/ref/masked_regions_reference.bed",
    log:
        "logs/mask_reference/mask.log",
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/mask_reference.py"

# TODO more fileformats
rule masked_gff2bed:
    input:
        config["ref"].get("masked_regions", []),
    output:
        "resources/masked_any2bed/masked_gff2bed.bed",
    log:
        "logs/masked_any2bed/masked_gff2bed.log",
    conda:
        "../envs/bedops.yaml"
    shell:
        "gff2bed < {input} > {output} 2> {log}"


# TODO add log to the script
rule mask_vcf:
    input:
        vcf=rules.norm_vcf.output,
        bed=get_masked_regions_bed,
    output:
        vcf=temp("results/calling/{experiment}.masked.vcf.gz"),
        tbi=temp("results/calling/{experiment}.masked.vcf.gz.tbi"),
    log:
        "logs/mask_vcf/{experiment}.log",
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/mask_vcf.py"
