rule mock_config:
    output:
        temp("resources/snpeff/{experiment}.mock.config"),
    log:
        "logs/snpeff/{experiment}.mock.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/snpEff_mock_conf.py"


rule snpeff:
    input:
        calls=get_snpeff_input,
        db=get_snpeff_db,
        conf=rules.mock_config.output,
    output:
        calls="results/calling/{experiment}.annot.vcf.gz",
    log:
        "logs/snpeff/{experiment}.log",
    params:
        extra=lambda wc, input: f"-c {input.conf} {config['snpEff']['extra']}",
    resources:
        mem_mb=4096,
    wrapper:
        "v7.2.0/bio/snpeff/annotate"
