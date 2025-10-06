rule find_candidate_mutations:
    input:
        vcf=get_candidate_mutations_input,
    output:
        tsv="results/candidate_mutations/{experiment}.tsv",
        info="results/candidate_mutations/{experiment}_plots/info.txt",
    params:
        wt_samples=get_wt_samples,
    log:
        "logs/candidate_mutations/{experiment}.log",
    conda:
        "../envs/candidate_mutations.yaml"
    script:
        "../scripts/candidate_mutations.py"
