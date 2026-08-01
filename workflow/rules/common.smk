import re

import pandas as pd
from pathlib import Path
from snakemake.exceptions import WorkflowError


samples = pd.read_csv(config["samples"], dtype=str, comment="#").set_index(
    "sample_name", drop=False
)
experiments = samples["Experiment"].dropna().unique()
ncbi_assemblies = config["ref"].get("ncbi", {}).get("assemblies", [])
ncbi_sequences = config["ref"].get("ncbi", {}).get("sequences", [])


def validate_samples(samples):
    """Reject sample sheets the workflow cannot actually process.

    Both conditions below used to fail deep inside a wrapper (or, worse,
    silently pass the same BAM to the caller twice), so they are caught here
    while the error can still name the offending row.
    """
    required = ["sample_name", "R1", "R2", "Experiment", "Condition"]
    missing_columns = [c for c in required if c not in samples.columns]
    if missing_columns:
        raise WorkflowError(
            f"Missing column(s) in {config['samples']}: "
            f"{', '.join(missing_columns)}. Required: {', '.join(required)}."
        )

    duplicated = samples.index[samples.index.duplicated()].unique().tolist()
    if duplicated:
        raise WorkflowError(
            f"Duplicate sample_name(s) in {config['samples']}: "
            f"{', '.join(duplicated)}. "
            "Each sample must appear exactly once; multiple sequencing runs "
            "per sample are not supported."
        )

    missing_r2 = samples.index[samples["R2"].isna()].tolist()
    if missing_r2:
        raise WorkflowError(
            f"Sample(s) with no R2 in {config['samples']}: "
            f"{', '.join(missing_r2)}. "
            "This workflow supports paired-end reads only."
        )


validate_samples(samples)


def name_alternatives(names):
    """Regex alternation over literal names, longest first.

    re.escape matters because a sample name containing '.', '+' or '(' would
    otherwise be interpolated into the constraint as a metacharacter.
    """
    return "|".join(re.escape(n) for n in sorted(set(names), key=len, reverse=True))


wildcard_constraints:
    sample=name_alternatives(samples["sample_name"].tolist()),
    fq="R1|R2",
    experiment=name_alternatives(experiments),


# The reference is staged into resources/ref/ and indexed there, so that bwa
# and samtools never write their index files next to the user's own FASTA
# (which may live on a read-only or shared path).
REF_FASTA = "resources/ref/reference.fa"


def get_ref_source():
    """Path the staged reference is built from. Evaluated at parse time."""
    cfg_fasta = config["ref"].get("fasta", None)

    if cfg_fasta:
        if not Path(cfg_fasta).exists():
            raise WorkflowError(
                f"The specified reference FASTA file does not exist: {cfg_fasta}"
            )
        if Path(cfg_fasta) == Path(REF_FASTA):
            raise WorkflowError(
                f"ref: fasta must not point at {REF_FASTA}; that path is "
                "reserved for the staged copy built by the workflow."
            )
        return cfg_fasta
    elif ncbi_assemblies or ncbi_sequences:
        return "resources/ref/ncbi_ref.fa"
    else:
        raise WorkflowError(
            "No reference FASTA file or NCBI accession was provided in the config."
        )


def get_ref_fasta(wildcards):
    return REF_FASTA


def get_annot_gff(wildcards):
    fasta_provided = config["ref"].get("fasta", None)
    # Use gff param only if fasta file provided
    if fasta_provided and config["ref"].get("gff", None):
        gff = config["ref"]["gff"]
        if not Path(gff).exists():
            raise WorkflowError(f"The specified GFF file does not exist: {gff}")
        return gff
    # Use NCBI GFF file only if there is no "fasta" file
    elif not fasta_provided and (ncbi_assemblies or ncbi_sequences):
        return "resources/ref/ncbi_annot.gff"
    elif config["annotate"]:
        # Falling through to None here used to silently drop VEP annotation
        # from the candidate table with no explanation.
        raise WorkflowError(
            "annotate is true but no annotation is available: set ref: gff "
            "(alongside ref: fasta), or use the ref: ncbi options, or set "
            "annotate: false."
        )


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


def get_raw_reads(wildcards):
    # validate_samples guarantees one row per sample, both mates present.
    sample = samples.loc[wildcards.sample]
    return sample["R1"], sample["R2"]


def get_fastqc_raw_reads(wildcards):
    if wildcards.fq == "R1":
        return get_raw_reads(wildcards)[0]
    elif wildcards.fq == "R2":
        return get_raw_reads(wildcards)[1]
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

    output.extend(
        expand(
            "results/qc/fastqc_raw/{sample}_{fq}_fastqc.zip",
            sample=experiment_df.index,
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
            "results/qc/fastqc_trimmed/{sample}_{fq}_fastqc.zip",
            sample=experiment_df.index,
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
