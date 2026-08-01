rule prepare_gff:
    input:
        gff=get_annot_gff
    output:
        gff="resources/vep/annotations.gff.gz",
        tbi="resources/vep/annotations.gff.gz.tbi"
    log:
        "logs/annotation/prepare_gff.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        # I use sed '/^#/d' and not just grep -v "^#" 
        # Because if the input file is empty grep fails
        """
        sed '/^#/d' {input.gff} \
        | sort -t$'\t' -k1,1 -k4,4n -k5,5n \
        | bgzip -c > {output.gff} 2> {log};

        tabix -p gff {output.gff} 2>> {log}
        """

rule plugins_folder:
    output:
        directory("resources/vep/plugins")
    log:
        "logs/annotation/dummy_plugins.log"
    shell:
        "mkdir -p {output}"

rule annotate_variants:
    input:
        calls=get_annotation_input,
        plugins=rules.plugins_folder.output,
        fasta=get_ref_fasta,
        fai=rules.samtools_index_ref.output, # fasta index
        gff=rules.prepare_gff.output.gff, # note: GFF files must be sorted, gzipped and indexed
        csi=rules.prepare_gff.output.tbi, # tabix index
    output:
        calls="results/calling/{experiment}.annot.vcf.gz",  # .vcf, .vcf.gz or .bcf
        stats="results/qc/vep/{experiment}.html",
    params:
        plugins="",
        extra=lambda wildcards: (
            "--allele_number "
            + "--warning_file STDERR "
            + ("--shift_3prime 1 " if config["vep"]["shift_3prime"] else "")
            + ("--hgvs " if config["vep"]["hgvs"] else "")
            + config["vep"]["extra"]
        ).strip(),
    log:
        "logs/annotation/{experiment}_vep_annotate.log",
    threads: 8
    wrapper:
        "v8.0.0/bio/vep/annotate"