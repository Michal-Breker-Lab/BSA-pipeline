rule download_ncbi:
    output:
        temp(
            expand(
                "resources/ref/ncbi/{acc}.{ext}",
                acc=ncbi_assemblies + ncbi_sequences, 
                ext=["fa", "gff"]
            )
        )
    log:
        "logs/ref/download_ncbi.log"
    params:
        ncbi_assemblies=ncbi_assemblies,
        ncbi_sequences=ncbi_sequences,
        email=config["ref"].get("ncbi", {}).get("email", ""),
        out_dir=lambda w, output: str(Path(output[0]).parent)
    conda:
        "../envs/ncbi_download.yaml"
    script:
        "../scripts/download_ncbi.py"

if not config["ref"].get("fasta", None):
    rule concat_ncbi_fasta:
        input:
            expand(
                "resources/ref/ncbi/{acc}.fa",
                acc=ncbi_assemblies + ncbi_sequences, 
            )
        output:
            "resources/ref/ncbi_ref.fa"
        log:
            "logs/ref/concat_ncbi_fasta.log"
        shell:
            "cat {input} > {output} 2> {log} "

    rule concat_ncbi_gff:
        input:
            expand(
                "resources/ref/ncbi/{acc}.gff",
                acc=ncbi_assemblies + ncbi_sequences, 
            )
        output:
            gff="resources/ref/ncbi_annot.gff"
        log:
            "logs/ref/concat_ncbi_gff.log"
        shell:
            "cat {input} > {output.gff} 2> {log}" 

rule bwa_index_ref:
    input:
        get_ref_fasta,
    output:
        idx=multiext(
            get_ref_fasta(""), ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    log:
        "logs/bwa_index_ref/reference.log",
    threads: 6
    wrapper:
        "v7.2.0/bio/bwa/index"


rule samtools_index_ref:
    input:
        get_ref_fasta,
    output:
        f"{get_ref_fasta("")}.fai",
    log:
        "logs/samtools_index_ref/reference.log",
    threads: 6
    wrapper:
        "v7.2.0/bio/samtools/faidx"
