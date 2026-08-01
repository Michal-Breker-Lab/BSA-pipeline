"""
Candidate Mutations Analysis Script

This script processes VCF files to identify candidate mutations in BSA (Bulk Segregant Analysis)
experiments. It compares mutant samples to wild-type samples, filters variants based on
allele frequency differences (delta SNP index), and generates annotated reports with plots.

Key features:
- Filters variants by minimum depth and delta SNP index thresholds
- Annotates variants using VEP (Variant Effect Predictor) annotations
- Generates scatter plots showing mutation deltas across chromosomes
- Supports both wild-type specified and wild-type inferred modes
"""

from snakemake.script import snakemake
import sys

# Redirect stdout to log file for Snakemake logging
sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")

import re
from cyvcf2 import VCF
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless: no display needed, silences Qt plugin warnings

# Keep vector output editable in Adobe Illustrator:
# fonttype 42 embeds TrueType so glyphs stay real, selectable text rather than
# being converted to outlines (the default Type 3 does not open cleanly in
# Illustrator). svg.fonttype "none" leaves text as <text> elements.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

# Arial, with metric-compatible fallbacks. Arial itself is proprietary and
# comes from the 'mscorefonts' package in this rule's conda environment;
# Liberation Sans and Arimo are metrically identical to it, so a fallback
# changes the glyph shapes slightly but not the layout, and text can be
# reset to Arial in Illustrator without reflowing.
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [
    "Arial",
    "Liberation Sans",
    "Arimo",
    "Helvetica",
    "Nimbus Sans",
    "DejaVu Sans",
]
# Use ASCII hyphen rather than U+2212 MINUS: the Unicode minus is a common
# cause of missing glyphs and copy-paste surprises in Adobe applications.
matplotlib.rcParams["axes.unicode_minus"] = False

# Font fallback is silent by default, which is how a whole batch of figures
# ends up in the wrong typeface unnoticed. Resolve it once and say so.
from matplotlib import font_manager

# Pass the family list, not the string "sans-serif": a bare string containing
# a hyphen is parsed as a fontconfig pattern and raises. The list is also what
# we actually want resolved -- findfont returns the first name available.
_font_path = font_manager.findfont(
    font_manager.FontProperties(family=matplotlib.rcParams["font.sans-serif"])
)
_font_name = font_manager.FontProperties(fname=_font_path).get_name()
print(f"Plot font resolved to: {_font_name} ({_font_path})")
if _font_name != "Arial":
    print(
        f"WARNING: Arial was not found in this environment; figures use "
        f"{_font_name} instead. Install the 'mscorefonts' conda package to "
        f"get genuine Arial."
    )

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def ann2df(ann, alt):
    """
    Convert VEP annotation string to DataFrame and filter by allele number.
    
    This function parses the VEP CSQ (Consequence) annotation field, which contains
    pipe-delimited annotation records. It filters annotations for a specific alternate
    allele and sorts them by severity (most severe first).
    
    Args:
        ann (str): VEP annotation string from CSQ INFO field (comma-separated records)
        alt (int): Alternate allele index (1-based) to filter annotations for
    
    Returns:
        pd.DataFrame: DataFrame with VEP annotations, sorted by severity (most severe first)
    """
    # VEP consequence terms ordered by severity (most severe to least severe)
    # This ordering is used to prioritize the most impactful variant effects
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
    
    # Parse annotation string: split by comma (multiple transcripts) then by pipe (fields)
    data = [record.split('|') for record in ann.split(',')]
    df = pd.DataFrame(data, columns=annot_fields)
    
    # Filter annotations for the specific alternate allele
    df = df.loc[df["ALLELE_NUM"] == str(alt)]

    def get_severity_rank(consequence_str):
        """
        Calculate severity rank for a consequence string.
        
        Handles cases where multiple consequences are combined with '&' (e.g., 
        'missense_variant&splice_region_variant'). Returns the rank of the most
        severe consequence (smallest rank number).
        
        Args:
            consequence_str (str): Consequence string, possibly with multiple terms joined by '&'
        
        Returns:
            int: Severity rank (lower = more severe)
        """
        terms = consequence_str.split('&')
        ranks = [severity_rank.get(t, 9999) for t in terms]  # unknown = very low severity
        return min(ranks)

    # Sort by severity (most severe first) to prioritize important annotations
    df["severity_rank"] = df["Consequence"].apply(get_severity_rank)
    df = df.sort_values("severity_rank").drop(columns="severity_rank")
    
    return df


def delta_no_wt(prop, threshold):
    """
    Calculate delta SNP index when wild-type samples are not explicitly specified.
    
    This function identifies candidate mutations by comparing all samples pairwise.
    For each sample with a high alternate allele frequency (>= threshold), it finds
    the sample with the lowest frequency for the same allele. The difference (delta)
    must be at least the threshold to be considered a candidate mutation.
    
    Args:
        prop (np.ndarray): 2D array of allele proportions [samples x alleles]
                          First column (index 0) is REF, subsequent columns are ALT alleles
        threshold (float): Minimum delta SNP index to consider a candidate mutation
    
    Returns:
        list: List of [mutant_sample_idx, alt_allele_idx, delta, wt_sample_idx] tuples
    """
    snp_data = []
    # Find all positions where alternate allele frequency >= threshold
    rows, cols = np.where(prop[:,1:] >= threshold)
    for r,c in zip(rows, cols + 1):
        snp_index = prop[r][c]
        # Calculate delta: difference between this sample and all other samples
        delta = snp_index - prop[:,c]
        delta[r] = np.nan # Do not count itself
        delta[delta < 0] = np.nan  # Only consider cases where mutant > other samples
        # Check if there's a sample with significantly lower frequency (potential WT)
        if np.any(delta >= 0) and np.nanmin(delta) >= threshold:
            # Report the sample that produced the *smallest* delta -- the binding
            # constraint, and the same delta returned below. Using nanargmax here
            # named a different sample than the delta beside it in the output.
            wt_index = np.nanargmin(delta)
            snp_data.append([r, c, delta[wt_index], wt_index])
    return snp_data


def delta_wt(prop, wt_indexes, threshold):
    """
    Calculate delta SNP index when wild-type samples are explicitly specified.
    
    This function compares mutant samples to wild-type samples. For each mutant sample
    and alternate allele, it finds the WT sample with the most similar allele frequency
    and calculates the delta. Only deltas >= threshold are reported.
    
    Args:
        prop (np.ndarray): 2D array of allele proportions [samples x alleles]
                          First column (index 0) is REF, subsequent columns are ALT alleles
        wt_indexes (list): List of sample indices that are wild-type
        threshold (float): Minimum absolute delta SNP index to consider a candidate mutation
    
    Returns:
        list: List of [mutant_sample_idx, alt_allele_idx, delta, wt_sample_idx] tuples
              Delta can be positive (mutant has higher ALT freq) or negative (mutant has lower ALT freq)
    """
    snp_data = []
    prop_alt = prop[:,1:]  # Extract alternate alleles only (exclude REF)
    mut_indexes = [i for i in range(prop.shape[0]) if i not in wt_indexes]
    prop_wt = prop_alt[wt_indexes]  # Allele proportions for WT samples only
    n_alts = prop_alt.shape[1]

    for mut_index in mut_indexes:
        # Calculate deltas between mutant and all WT samples for all alternate alleles
        deltas = prop_alt[mut_index] - prop_wt
        # Skip if no delta meets threshold
        if not (np.abs(deltas) >= threshold).any():
             continue
        # For each alternate allele, find the WT sample with closest frequency
        for alt_idx in range(n_alts):
            alt_delta = deltas[:,alt_idx]
            wt_index = np.nanargmin(np.abs(alt_delta))  # WT sample with closest frequency
            snp_delta = alt_delta[wt_index]
            # Only report if delta meets threshold
            if np.abs(snp_delta) < threshold:
                continue
            # alt_idx + 1 because alt_idx is 0-based in prop_alt, but 1-based in original prop
            snp_data.append([mut_index,alt_idx + 1,snp_delta, wt_indexes[wt_index]])

    return snp_data


# ============================================================================
# Configuration: Load parameters from Snakemake config
# ============================================================================
wt_samples = snakemake.params['wt_samples']  # List of wild-type sample names (None if not specified)
report_delta = snakemake.config['candidate_mutations']['report_delta']  # Minimum delta to report in TSV
skip_masked = snakemake.config['candidate_mutations']['skip_masked']  # Skip variants with FILTER flag
keep_all_annotations = snakemake.config['candidate_mutations']['keep_all_annotations']  # Keep all VEP annotations or just most severe
min_dp = snakemake.config['candidate_mutations']['min_dp']  # Minimum read depth to consider variant
plot_delta = snakemake.config['candidate_mutations']['plot_delta']  # Minimum delta to include in plots
delta_seq_name = snakemake.config['candidate_mutations'].get("delta_seq_name", {})  # Sequence-specific delta thresholds

# ============================================================================
# VCF Processing: Open VCF and extract metadata
# ============================================================================
vcf = VCF(snakemake.input.vcf)
sample_names = vcf.samples

# Extract VEP column names from VCF INFO field header
# VEP annotations are stored in CSQ field with format: Field1|Field2|Field3|...
annot_description = re.search(r'##INFO=<ID=CSQ[^>]*Format: ([^"]+)"', vcf.raw_header)
annot_fields = annot_description.group(1).split('|') if annot_description else []

# Initialize data structures for storing results
records = []  # List of DataFrames, one per candidate mutation
plot_data = []  # List of dicts for plotting (includes all variants above plot_delta)

# Validate configuration
if plot_delta > report_delta:
    raise ValueError(f"'plot_delta' cannot be larger than 'report_delta'")

# ============================================================================
# Main Processing Loop: Iterate through VCF records
# ============================================================================
for record in vcf:
    # Skip masked variants if configured to do so
    if record.FILTER and skip_masked:
        continue
    
    # Extract allele depth (AD) format field: [REF_depth, ALT1_depth, ALT2_depth, ...]
    ad_array = record.format("AD")
    
    # Filter samples: keep only those with valid genotypes and sufficient depth
    passed_genotypes = [i for i,g in enumerate(record.genotypes) if g[0] >= 0]  # Valid genotype (not missing)
    passed_ad_dp = np.where(ad_array.sum(axis=1) >= min_dp)[0]  # Total depth >= min_dp
    # sorted(), not list(): these indices must stay row-aligned with both
    # ad_array and sample_names, and set iteration order is not a contract.
    passed_samples_idx = sorted(set(passed_genotypes) & set(passed_ad_dp))

    # Subset AD array to only passed samples
    ad_array = ad_array[passed_samples_idx]
    if ad_array.shape[0] < 2:
        continue  # Need at least 2 samples to compare

    # Get sequence-specific delta thresholds (if configured) or use defaults
    record_plot_delta = delta_seq_name.get(record.CHROM, plot_delta)
    record_report_delta = delta_seq_name.get(record.CHROM, report_delta)

    passed_sample_names = [sample_names[i] for i in passed_samples_idx]
    # Calculate allele proportions (frequencies) from allele depths
    ad_prop = ad_array / ad_array.sum(axis=1, keepdims=True)

    # Calculate delta SNP index: use different function depending on whether WT is specified
    if wt_samples:
        # Mode 1: WT samples explicitly specified
        wt_indexes = [i for i,s in enumerate(passed_sample_names) if s in wt_samples]
        if not wt_indexes:
            continue  # Skip if no WT samples in passed samples
        res = delta_wt(ad_prop, wt_indexes, record_plot_delta)
    else: 
        # Mode 2: Infer WT by comparing all samples pairwise
        res = delta_no_wt(ad_prop, record_plot_delta)

    # Process each candidate mutation found for this variant
    for mut in res:
        sample_idx, alt_idx, delta, wt_idx = mut
        sample = passed_sample_names[sample_idx]

        # Prepare plot data: color indicates direction of delta
        # Red = positive delta (mutant has higher ALT frequency), Blue = negative delta
        color = "red" if delta > 0 else "blue"
        plot_data.append({
            "CHROM": record.CHROM,
            "POS": record.POS,
            "Sample": sample,
            "delta": abs(delta),
            "color": color,
            "facecolor": "none" if record.FILTER else color  # Masked variants have no fill
        })

        # Only report variants that meet the report_delta threshold
        if abs(delta) < record_report_delta:
                continue

        # Create base record with variant information
        allele_record = {
            "REC_ID": len(records) + 1,
            "CHROM": record.CHROM,
            "POS": record.POS,
            "REF": record.REF,
            "ALT": "|".join(record.ALT),  # All alternate alleles (pipe-separated)
            "FILTER": record.FILTER,
            "WT": passed_sample_names[wt_idx],
            "Sample": sample,
            "AD_wt": "|".join(map(str,ad_array[wt_idx])),  # Allele depths for WT
            "AD_sample": "|".join(map(str,ad_array[sample_idx])),  # Allele depths for mutant
            "ALT_index": alt_idx,  # Which alternate allele this mutation refers to
            "delta": round(delta, 3)
        }
        
        allele_df = pd.DataFrame([allele_record])
        ann_record = record.INFO.get('CSQ', '')  # VEP annotation field

        # If CSQ field is missing, append record without annotation
        if not ann_record:
            records.append(allele_df)
            continue

        # Parse VEP annotations and filter for this alternate allele
        annot_df = ann2df(ann_record, alt_idx)
        
        # Keep only most severe annotation if configured (otherwise keep all)
        if not keep_all_annotations:
            annot_df = annot_df.iloc[[0]]  # Most severe is first after sorting
        
        # Cartesian join: combine variant info with all annotations
        # This creates one row per annotation (e.g., one per transcript)
        allele_df["key"] = 1
        annot_df["key"]   = 1

        merged_df = (
            allele_df
            .merge(annot_df, on="key")
            .drop(columns="key")
        )
        
        records.append(merged_df)

# ============================================================================
# Output: Combine all records and write TSV file
# ============================================================================
if records:
    # Concatenate all DataFrames and remove columns that are all NaN
    df = pd.concat(records, ignore_index=True).dropna(axis=1, how="all")
else:
    # Create empty DataFrame with required columns if no mutations found
    df = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "FILTER", "Sample"])

# Write candidate mutations to TSV file
df.to_csv(snakemake.output.tsv, sep="\t", index=False)


# ============================================================================
# Plotting: Generate scatter plots for each chromosome-sample combination
# ============================================================================
if plot_data:
    plot_df = pd.DataFrame(plot_data)
else:
    plot_df = pd.DataFrame(columns=["CHROM", "Sample"])

# Filter plot data: only keep chromosomes/samples that have candidate mutations
# This avoids creating plots with no significant mutations
mask = plot_df.set_index(['CHROM', 'Sample']).index.isin(
    df.set_index(['CHROM', 'Sample']).index
)
plot_df = plot_df[mask]

plot_dir = Path(snakemake.output.plots)
plot_dir.mkdir(parents=True, exist_ok=True)
info_file = plot_dir / "info.txt"


def safe_filename(*parts):
    """Join name parts into a filesystem-safe filename stem.

    Sequence IDs routinely contain '.' and '|', and the previous
    '{sample}:{chrom}' pattern is not a legal filename on Windows or SMB.
    """
    return "__".join(re.sub(r"[^A-Za-z0-9._-]", "_", str(p)) for p in parts)

# Generate one plot per chromosome-sample combination
for chr, sample in plot_df[['CHROM', 'Sample']].drop_duplicates().values:
    # Extract data for this chromosome-sample combination
    subset = plot_df[(plot_df['CHROM'] == chr) & (plot_df['Sample'] == sample)]
    subset = subset.sort_values(by='POS')  # Sort by position for plotting

    # Create scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(
        subset['POS'],
        subset['delta'], 
        edgecolors=subset["color"],  # Edge color indicates delta direction
        facecolors=subset["facecolor"],  # Face color (none for masked variants)
        linewidths=1, 
        label=None
    )
    plt.title(f"{sample}:{chr}")
    
    # Add horizontal line at report_delta threshold
    plt.axhline(y=delta_seq_name.get(chr, report_delta), color="red", linestyle="--", linewidth=1)
    
    # Set y-axis limits: from plot_delta to 1.02 (slightly above max possible delta of 1.0)
    plt.ylim(delta_seq_name.get(chr, plot_delta) - 0.02, 1.02)
    plt.ylabel("Absolute delta SNP index")
    plt.xlabel("Genomic position")

    # Build legend: show different marker types based on what's in the plot
    legend_elements = [
        Line2D(
            [0], [0], marker='o', color='red', markerfacecolor='red', 
            markersize=10, label='Delta > 0', linestyle='None'
        )
    ]

    # Add blue marker legend if negative deltas are present
    if subset["facecolor"].eq("blue").any():
        legend_elements.append(
            Line2D(
                [0], [0], marker='o', color='blue', markerfacecolor='blue', 
                markersize=10, label='Delta < 0', linestyle='None'
            )
        )

    # Add masked variant legend if any masked variants are present
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
    
    # Save as PNG (raster preview) plus PDF and SVG (fully vector, text kept
    # as text -- see the rcParams above). No artist is rasterized, so the
    # vector outputs contain no embedded bitmaps.
    stem = safe_filename(sample, chr)
    plt.savefig(plot_dir / f"{stem}.png", dpi=300)
    plt.savefig(plot_dir / f"{stem}.pdf")
    plt.savefig(plot_dir / f"{stem}.svg")
    plt.close()

# Write summary info file
with open(info_file, "w") as f:
    f.write(f"Total candidate mutations: {len(df)}\n")
    f.write(f"Total plotted variants: {len(plot_df)}\n")
    f.write("Plots are saved in the same directory as this info file.\n")
