rule freebayes:
    input:
        alns=get_bams_experiment,
        idsx=get_bai_experiment,
        ref=config["ref_fasta"],
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
        ref=config["ref_fasta"],
    output:
        temp("results/calling/{experiment}.norm.vcf.gz"),
    log:
        "logs/norm_vcf/{experiment}.log",
    wrapper:
        "v7.2.0/bio/bcftools/norm"


# TODO more fileformats
rule regions2bed:
    input:
        get_regions_to_mask,
    output:
        temp("resources/regions_bed/{regions_file}.bed"),
    log:
        "logs/regions2bed/{regions_file}.log",
    conda:
        "../envs/bedops.yaml"
    shell:
        "gff2bed < {input} > {output} 2> {log}"


# TODO add log to the script
rule mask_vcf:
    input:
        vcf=rules.norm_vcf.output,
        bed=get_bed,
    output:
        vcf=temp("results/calling/{experiment}.overlaps.vcf.gz"),
    log:
        "logs/mask_vcf/{experiment}.log",
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/add_overlaps_vcf.py"
