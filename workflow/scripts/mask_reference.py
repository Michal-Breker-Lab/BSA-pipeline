from snakemake.script import snakemake
import sys

sys.stdout = sys.stderr
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

from Bio import SeqIO

def find_repeats(seq):
    """Finds lowercase (repeat) regions in a sequence."""
    regions = []
    in_repeat = False
    start = 0

    for i, base in enumerate(seq):
        if base not in "ATGC":
            if not in_repeat:
                start = i
                in_repeat = True
        else:
            if in_repeat:
                regions.append((start, i))
                in_repeat = False
    if in_repeat:
        regions.append((start, len(seq)))
    return regions

def fasta_to_bed(input_fasta, output_bed):
    with open(output_bed, "w") as bed:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            repeats = find_repeats(seq)
            for start, end in repeats:
                bed.write(f"{record.id}\t{start}\t{end}\n")

fasta_to_bed(snakemake.input[0], snakemake.output[0])