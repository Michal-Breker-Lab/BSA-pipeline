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

# The reference is staged into resources/ref/ so that every index below is
# written beside the staged copy rather than beside the user's own FASTA,
# which may sit on a read-only or shared path.
rule stage_reference:
    input:
        get_ref_source(),
    output:
        REF_FASTA,
    log:
        "logs/ref/stage_reference.log",
    shell:
        "ln -sfr {input} {output} 2> {log}"


rule bwa_index_ref:
    input:
        REF_FASTA,
    output:
        idx=multiext(REF_FASTA, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index_ref/reference.log",
    threads: 6
    wrapper:
        "v7.2.0/bio/bwa/index"


rule samtools_index_ref:
    input:
        REF_FASTA,
    output:
        f"{REF_FASTA}.fai",
    log:
        "logs/samtools_index_ref/reference.log",
    threads: 6
    wrapper:
        "v7.2.0/bio/samtools/faidx"
