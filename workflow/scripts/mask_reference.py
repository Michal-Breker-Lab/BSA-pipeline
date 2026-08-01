from typing import TYPE_CHECKING
import sys

if TYPE_CHECKING:
    from snakemake.script import Snakemake

# Snakemake binds this object in the module globals through the preamble it
# prepends to script: rules. The bare annotation declares it for type checkers
# and linters without importing or assigning anything at runtime -- importing
# the name works only because that preamble mutates snakemake.script, which
# breaks editors, linting and any standalone run of this file.
snakemake: "Snakemake"

sys.stdout = sys.stderr
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

import re

from Bio import SeqIO

# Anything that is not an uppercase A/T/G/C: soft-masked (lowercase) repeats
# plus ambiguity codes and N runs.
NON_ACGT = re.compile(r"[^ATGC]+")

def find_repeats(seq):
    """Finds masked (lowercase or ambiguous) regions in a sequence.

    Returns 0-based half-open (start, end) intervals -- BED convention.
    Uses one regex scan rather than a per-base Python loop, which took
    minutes to hours on a real assembly.
    """
    return [m.span() for m in NON_ACGT.finditer(seq)]

def fasta_to_bed(input_fasta, output_bed):
    with open(output_bed, "w") as bed:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            repeats = find_repeats(seq)
            for start, end in repeats:
                bed.write(f"{record.id}\t{start}\t{end}\n")

fasta_to_bed(snakemake.input[0], snakemake.output[0])