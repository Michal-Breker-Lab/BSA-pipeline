import pandas as pd
from pathlib import Path


samples = pd.read_csv(config["samples"], dtype=str, comment="#").set_index(
    "sample_name", drop=False
)
experiments = samples["Experiment"].dropna().unique()
ncbi_assemblies = config["ref"].get("ncbi", {}).get("assemblies", [])
ncbi_sequences = config["ref"].get("ncbi", {}).get("sequences", [])


wildcard_constraints:
    sample="|".join(samples["sample_name"].tolist()),
    fq="R1|R2|single",
    experiment="|".join(experiments),


def get_ref_fasta(wildcards):
    cfg_ref = config["ref"]
    cfg_fasta = cfg_ref.get("fasta", None)

    if cfg_fasta:
        if not Path(cfg_fasta).exists():
            raise ValueError(f"The specified reference FASTA file does not exist: {cfg_fasta}")
        return cfg_fasta
    elif ncbi_assemblies or ncbi_sequences:
        return "resources/ref/ncbi_ref.fa"
    else:
        raise ValueError("No reference FASTA file or NCBI accession was provided in the config.")


def get_annot_gff(wildcards):
    fasta_provided = config["ref"].get("fasta", None)
    # Use gff param only if fasta file provided
    if fasta_provided and config["ref"].get("gff", None):
        gff = config["ref"]["gff"]
        if not Path(gff).exists():
            raise ValueError(f"The specified GFF file does not exist: {gff}")
        return gff
    # Use NCBI GFF file only if there is no "fasta" file
    elif not fasta_provided and (ncbi_assemblies or ncbi_sequences):
        return "resources/ref/ncbi_annot.gff"


def get_masked_regions_bed(wildcards):
    if not config["ref"].get("fasta", None):
        return "resources/ref/masked_regions_reference.bed"

    mask_file = Path(config["ref"].get("masked_regions", "resources/ref/masked_regions_reference.bed"))
    if mask_file.suffix == ".bed":
        return mask_file
    elif mask_file.suffix in [".gff3", ".gff"]:
        return f"resources/masked_any2bed/masked_gff2bed.bed"
    else:
        raise ValueError(
            f"Invalid file extension '{mask_file.suffix}'. "
            "Only .bed, .gff and .gff3 file extensions are supported."
        )


def is_single(sample_name):
    sample = samples.loc[sample_name]

    if pd.isna(sample["R2"]):
        return True

    return False


def get_raw_reads(wildcards):
    sample = samples.loc[wildcards.sample]
    # TEMP solution?
    if isinstance(sample, pd.DataFrame):
        dup = sample.duplicated(subset=["sample_name", "R1", "R2"], keep=False)
        if dup.all():
            sample = sample.iloc[0]
        else:
            raise ValueError(f"{wildcards.sample} has more than one set of reads")
    if pd.isna(sample["R2"]):
        return [sample["R1"]]

    return sample["R1"], sample["R2"]


def get_fastqc_raw_reads(wildcards):
    if wildcards.fq == "R2":
        return get_raw_reads(wildcards)[1]
    if wildcards.fq in ["single", "R1"]:
        return get_raw_reads(wildcards)[0]
    else:
        raise ValueError(f"'fq' wildcard {wildcards.fq} is invalid")


def get_bam(wildcards):
    if config.get("remove_duplicates"):
        return f"results/mapping/{wildcards.sample}.dedup.bam"
    else:
        return f"results/mapping/{wildcards.sample}.bam"


def get_bai(wildcards):
    if config.get("remove_duplicates"):
        return f"results/mapping/{wildcards.sample}.dedup.bam.bai"
    else:
        return f"results/mapping/{wildcards.sample}.bam.bai"


def get_bams_experiment(wildcards):
    experiment_df = samples.loc[samples["Experiment"] == wildcards.experiment]
    dedup = ".dedup" if config.get("remove_duplicates") else ""
    return expand(
        f"results/mapping/{{sample}}{dedup}.bam", sample=experiment_df.index
    )


def get_bai_experiment(wildcards):
    experiment_df = samples.loc[samples["Experiment"] == wildcards.experiment]
    dedup = ".dedup" if config.get("remove_duplicates") else ""
    return expand(
        f"results/mapping/{{sample}}{dedup}.bam.bai", sample=experiment_df.index
    )


def get_annotation_input(wildcards):
    if config["mask_regions"]:
        return "results/calling/{experiment}.masked.vcf.gz"
    return "results/calling/{experiment}.norm.vcf.gz"


def get_candidate_mutations_input(wildcards):
    if config["annotate"] and get_annot_gff(wildcards):
        return "results/calling/{experiment}.annot.vcf.gz"
    return get_annotation_input(wildcards)


def get_wt_samples(wildcards):
    experiment_df = samples.loc[samples["Experiment"] == wildcards.experiment]
    wt_samples = experiment_df.index[experiment_df["Condition"] == "WT"].tolist()
    return wt_samples if wt_samples else None


def multiqc_input(wildcards):
    output = []
    experiment_df = samples.loc[samples["Experiment"] == wildcards.experiment]
    pe_samples = experiment_df[experiment_df["R1"].notna() & experiment_df["R2"].notna()].index
    se_samples = experiment_df[experiment_df["R1"].notna() & experiment_df["R2"].isna()].index

    output.extend(
        expand(
            "results/qc/fastqc_raw/{sample}_single_fastqc.zip",
            sample=se_samples,
        )
    )
    output.extend(
        expand(
            "results/qc/fastqc_raw/{sample}_{fq}_fastqc.zip",
            sample=pe_samples,
            fq=["R1", "R2"],
        )
    )
    output.extend(
        expand(
            "results/qc/trimming/{sample}.json",
            sample=experiment_df.index,
        )
    )
    output.extend(
        expand(
            "results/qc/fastqc_trimmed/{sample}_single_fastqc.zip",
            sample=se_samples,
        )
    )
    output.extend(
        expand(
            "results/qc/fastqc_trimmed/{sample}_{fq}_fastqc.zip",
            sample=pe_samples,
            fq=["R1", "R2"],
        )
    )
    output.extend(
        expand(
            "results/qc/fastqc_trimmed/{sample}_single_fastqc.zip",
            sample=se_samples,
        )
    )
    output.extend(
        expand(
            "results/qc/fastqc_trimmed/{sample}_{fq}_fastqc.zip",
            sample=pe_samples,
            fq=["R1", "R2"],
        )
    )
    output.extend(
        expand(
            "results/qc/samtools_stats/{sample}.txt", sample=experiment_df.index
        )
    )
    output.extend(
        expand(
            "results/qc/samtools_idxstats/{sample}.bam.idxstats",
            sample=experiment_df.index,
        )
    )
    output.extend(
        expand(
            "results/qc/samtools_coverage/{sample}.txt",
            sample=experiment_df.index,
        )
    )
    output.extend(
        expand(
            "results/qc/picard_insert_size/{sample}.txt",
            sample=experiment_df.index,
        )
    )
    if config["remove_duplicates"]:
        output.extend(
            expand(
                "results/qc/picard_mark_dup/{sample}.metrics.txt",
                sample=experiment_df.index,
            )
        )
    
    return output


def get_final_output(wildcards):
    output = []
    
    output.extend(
        expand(
            "results/qc/multiqc_{experiment}.html",
            experiment=experiments
        )
    )

    output.extend(
        expand(
            "results/candidate_mutations/{experiment}.tsv",
            experiment=experiments,
        )
    )

    return output
