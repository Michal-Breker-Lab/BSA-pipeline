rule bwa_index_ref:
    input:
        config["ref_fasta"],
    output:
        idx=multiext(
            config["ref_fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    log:
        "logs/bwa_index_ref/reference.log",
    threads: 6
    wrapper:
        "v7.2.0/bio/bwa/index"


rule samtools_index_ref:
    input:
        config["ref_fasta"],
    output:
        f"{config["ref_fasta"]}.fai",
    log:
        "logs/samtools_index_ref/reference.log",
    threads: 6
    wrapper:
        "v7.2.0/bio/samtools/faidx"
