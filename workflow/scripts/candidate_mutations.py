from snakemake.script import snakemake
import sys

sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")

import re
from cyvcf2 import VCF
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def ann2df(ann, alt):
    vep_severity_order = [
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "initiator_codon_variant",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "intergenic_variant"
    ]
    severity_rank = {term: i for i, term in enumerate(vep_severity_order)}
    data = [record.split('|') for record in ann.split(',')]
    df = pd.DataFrame(data, columns=annot_fields)
    df = df.loc[df["ALLELE_NUM"] == str(alt)]

    def get_severity_rank(consequence_str):
        """
        Handles cases like: 'missense_variant&splice_region_variant'
        Most severe term = smallest rank.
        """
        terms = consequence_str.split('&')
        ranks = [severity_rank.get(t, 9999) for t in terms]  # unknown = very low severity
        return min(ranks)

    df["severity_rank"] = df["Consequence"].apply(get_severity_rank)
    df = df.sort_values("severity_rank").drop(columns="severity_rank")
    
    return df


def delta_no_wt(prop, threshold):
        snp_data = []
        rows, cols = np.where(prop[:,1:] >= threshold)
        for r,c in zip(rows, cols + 1):
            snp_index = prop[r][c]
            delta = snp_index - prop[:,c]
            delta[r] = np.nan # Do not count itself
            delta[delta < 0] = np.nan
            if np.any(delta >= 0) and np.nanmin(delta) >= threshold:
                wt_index = np.nanargmax(delta)
                snp_data.append([r,c,np.nanmin(delta), wt_index])
        return snp_data


def delta_wt(prop, wt_indexes, threshold):
        snp_data = []
        prop_alt = prop[:,1:]
        mut_indexes = [i for i in range(prop.shape[0]) if i not in wt_indexes]
        prop_wt = prop_alt[wt_indexes]
        n_alts = prop_alt.shape[1]

        for mut_index in mut_indexes:
            deltas = prop_alt[mut_index] - prop_wt
            if not (np.abs(deltas) >= threshold).any():
                 continue
            for alt_idx in range(n_alts):
                alt_delta = deltas[:,alt_idx]
                wt_index = np.nanargmin(np.abs(alt_delta))
                snp_delta = alt_delta[wt_index]
                if np.abs(snp_delta) < threshold:
                    continue
                snp_data.append([mut_index,alt_idx + 1,snp_delta, wt_indexes[wt_index]])

        return snp_data


wt_samples = snakemake.params['wt_samples']
report_delta = snakemake.config['candidate_mutations']['report_delta']
skip_masked = snakemake.config['candidate_mutations']['skip_masked']
keep_all_annotations = snakemake.config['candidate_mutations']['keep_all_annotations']
min_dp = snakemake.config['candidate_mutations']['min_dp']
plot_delta = snakemake.config['candidate_mutations']['plot_delta']
delta_seq_name = snakemake.config['candidate_mutations'].get("delta_seq_name", {})

vcf = VCF(snakemake.input.vcf)
sample_names = vcf.samples
# Extract VEP column names from VCF INFO field
annot_description = re.search(r'##INFO=<ID=CSQ[^>]*Format: ([^"]+)"', vcf.raw_header)
annot_fields = annot_description.group(1).split('|') if annot_description else []

records = []
plot_data = []

if plot_delta > report_delta:
    raise ValueError(f"'plot_delta' cannot be larger than 'report_delta'")

for record in vcf:
    if record.FILTER and skip_masked:
        continue
    ad_array = record.format("AD")
    
    passed_genotypes = [i for i,g in enumerate(record.genotypes) if g[0] >= 0]
    passed_ad_dp = np.where(ad_array.sum(axis=1) >= min_dp)[0]
    passed_samples_idx = list(set(passed_genotypes) & set(passed_ad_dp))

    ad_array = ad_array[passed_samples_idx]
    if ad_array.shape[0] < 2:
        continue

    record_plot_delta = delta_seq_name.get(record.CHROM, plot_delta)
    record_report_delta = delta_seq_name.get(record.CHROM, report_delta)

    passed_sample_names = [sample_names[i] for i in passed_samples_idx]
    ad_prop = ad_array / ad_array.sum(axis=1, keepdims=True)

    if wt_samples:
        wt_indexes = [i for i,s in enumerate(passed_sample_names) if s in wt_samples]
        if not wt_indexes:
            continue
        res = delta_wt(ad_prop, wt_indexes, record_plot_delta)
    else: 
        res = delta_no_wt(ad_prop, record_plot_delta)

    for mut in res:
        sample_idx, alt_idx, delta, wt_idx = mut
        sample = passed_sample_names[sample_idx]

        color = "red" if delta > 0 else "blue"
        plot_data.append({
            "CHROM": record.CHROM,
            "POS": record.POS,
            "Sample": sample,
            "delta": abs(delta),
            "color": color,
            "facecolor": "none" if record.FILTER else color
        })

        if abs(delta) < record_report_delta:
                continue

        allele_record = {
            "REC_ID": len(records) + 1,
            "CHROM": record.CHROM,
            "POS": record.POS,
            "REF": record.REF,
            "ALT": "|".join(record.ALT),  # ALT is 0-based
            "FILTER": record.FILTER,
            "WT": passed_sample_names[wt_idx],
            "Sample": sample,
            "AD_wt": "|".join(map(str,ad_array[wt_idx])),
            "AD_sample": "|".join(map(str,ad_array[sample_idx])),
            "ALT_index": alt_idx,
            "delta": round(delta, 3)
        }
        
        allele_df = pd.DataFrame([allele_record])
        ann_record = record.INFO.get('CSQ', '')

        # If CSQ field is missing append without annotation
        if not ann_record:
            records.append(allele_df)
            continue

        annot_df = ann2df(ann_record, alt_idx)
        # Keep only most severe if not specified
        if not keep_all_annotations:
            annot_df = annot_df.iloc[[0]]
        
        # Cartesian join
        allele_df["key"] = 1
        annot_df["key"]   = 1

        merged_df = (
            allele_df
            .merge(annot_df, on="key")
            .drop(columns="key")
        )
        
        records.append(merged_df)

if records:
    df = pd.concat(records, ignore_index=True).dropna(axis=1, how="all")
else:
    df = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "FILTER", "Sample"])

df.to_csv(snakemake.output.tsv, sep="\t", index=False)


### Plotting
if plot_data:
    plot_df = pd.DataFrame(plot_data)
else:
    plot_df = pd.DataFrame(columns=["CHROM", "Sample"])

# Keep only chromosome from the samples that have candidate mutations
mask = plot_df.set_index(['CHROM', 'Sample']).index.isin(
    df.set_index(['CHROM', 'Sample']).index
)
plot_df = plot_df[mask]

info_file = Path(snakemake.output.info)
plot_dir = info_file.parent

for chr, sample in plot_df[['CHROM', 'Sample']].drop_duplicates().values:
    subset = plot_df[(plot_df['CHROM'] == chr) & (plot_df['Sample'] == sample)]
    subset = subset.sort_values(by='POS')

    plt.figure(figsize=(10, 6))
    plt.scatter(
        subset['POS'],
        subset['delta'], 
        edgecolors=subset["color"], 
        facecolors=subset["facecolor"], 
        linewidths=1, 
        label=None
    )
    plt.title(f"{sample}:{chr}")
    plt.axhline(y=delta_seq_name.get(chr, report_delta), color="red", linestyle="--", linewidth=1)
    plt.ylim(delta_seq_name.get(chr, plot_delta) - 0.02, 1.02)
    plt.ylabel("Absolute delta SNP index")
    plt.xlabel("Genomic position")

    legend_elements = [
        Line2D(
            [0], [0], marker='o', color='red', markerfacecolor='red', 
            markersize=10, label='Delta > 0', linestyle='None'
        )
    ]

    if subset["facecolor"].eq("blue").any():
        Line2D(
            [0], [0], marker='o', color='blue', markerfacecolor='blue', 
            markersize=10, label='Delta < 0', linestyle='None'
        )

    if subset['facecolor'].eq('none').any():
        legend_elements.append(
            Line2D(
                [0], [0], marker='o', color='black', markerfacecolor='none', 
                markersize=10, label='Masked', linestyle='None'
                )
        )

    plt.legend(
        handles=legend_elements,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.1), 
        ncol=3
    )
    plt.tight_layout()
    plt.savefig(plot_dir / f"{sample}:{chr}.png", dpi=300)
    plt.savefig(plot_dir / f"{sample}:{chr}.svg", dpi=300)
    plt.close()

with open(info_file, "w") as f:
    f.write(f"Total candidate mutations: {len(df)}\n")
    f.write(f"Total plotted variants: {len(plot_df)}\n")
    f.write("Plots are saved in the same directory as this info file.\n")
