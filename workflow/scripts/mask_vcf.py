from snakemake.script import snakemake
from snakemake.shell import shell
import tempfile
import pandas as pd

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "(awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, \"MASKED\"}}' {snakemake.input.bed}"
        " | bgzip -c > {tmpdir}/exclude.filt.bed.gz &&"
        " tabix -p bed {tmpdir}/exclude.filt.bed.gz &&"
        " echo '##FILTER=<ID=MASKED,Description=\"Variant positions masked by BED \">' > {tmpdir}/add_overlap.hdr &&"
        " bcftools annotate"
        " -a {tmpdir}/exclude.filt.bed.gz"
        " -h {tmpdir}/add_overlap.hdr"
        " -c CHROM,FROM,TO,FILTER"
        # -O z, not -O v: the output path ends in .vcf.gz and downstream
        # consumers (VEP, cyvcf2) require it to actually be BGZF-compressed
        " -O z -o {snakemake.output.vcf}"
        " {snakemake.input.vcf} &&"
        " bcftools index -t {snakemake.output.vcf}"
        ") {log}"
    )
