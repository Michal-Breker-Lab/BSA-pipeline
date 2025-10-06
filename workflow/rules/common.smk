import pandas as pd
from pathlib import Path

samples = pd.read_csv(config["samples"], dtype=str, comment="#").set_index(
    "sample_name", drop=False
)
experiments = samples["Experiment"].dropna().unique()

wildcard_constraints:
    sample="|".join(samples["sample_name"].tolist()),
    fq="R1|R2|single",
    experiment="|".join(experiments),


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
    if config["remove_duplicates"]:
        return f"results/mapping/{wildcards.sample}.dedup.bam"
    else:
        return f"results/mapping/{wildcards.sample}.bam"

def get_bai(wildcards):
    if config["remove_duplicates"]:
        return f"results/mapping/{wildcards.sample}.dedup.bam.bai"
    else:
        return f"results/mapping/{wildcards.sample}.bam.bai"


def get_bams_experiment(wildcards):
    experiment_df = samples.loc[samples["Experiment"] == wildcards.experiment]
    dedup = ".dedup" if config["remove_duplicates"] else ""
    return expand(
        f"results/mapping/{{sample}}{dedup}.bam", sample=experiment_df.index
    )


def get_bai_experiment(wildcards):
    experiment_df = samples.loc[samples["Experiment"] == wildcards.experiment]
    dedup = ".dedup" if config["remove_duplicates"] else ""
    return expand(
        f"results/mapping/{{sample}}{dedup}.bam.bai", sample=experiment_df.index
    )

def get_snpeff_input(wildcards):
    if config.get("regions_to_mask"):
        return "results/calling/{experiment}.overlaps.vcf.gz" #rules.mask_vcf.output.vcf

    return "results/calling/{experiment}.norm.vcf.gz" #rules.norm_vcf.output


def get_candidate_mutations_input(wildcards):
    if config.get("snpEff", {}).get("custom_db"):
        return "results/calling/{experiment}.annot.vcf.gz" #rules.snpeff.output.calls

    return get_snpeff_input(wildcards)

def get_regions_to_mask(wildcards):
    if config.get("regions_to_mask"):
        return get_bed(wildcards)

    return None

def get_snpeff_db(wildcards):
    if config.get("snpEff", {}).get("custom_db"):
        return config.get("snpEff", {}).get("custom_db")

    return None

def get_bed(wildcards):
    masking_file = Path(config["regions_to_mask"])
    if masking_file.suffix == ".bed":
        return config["regions_to_mask"]
    elif masking_file.suffix in [".gff3", ".gff"]:
        return f"resources/regions_bed/{masking_file.stem}.bed"
    else:
        raise ValueError(
            f"Invalid file extension '{masking_file.suffix}'. "
            "Only .bed and .gff files are supported."
        )


def get_wt_samples(wildcards):
    experiment_df = samples.loc[samples["Experiment"] == wildcards.experiment]
    wt_samples = experiment_df.index[experiment_df["Condition"] == "WT"].tolist()
    return wt_samples if wt_samples else None


def get_vcf_to_analyze(wildcard):
    if config["snpEff"]["custom_db"]:
        return "results/calling/{experiment}.annot.vcf.gz"

    return "results/calling/{experiment}.norm.vcf.gz"


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
