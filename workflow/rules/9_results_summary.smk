# -----------------------------------------------------------------------------
# One-line-per-genome summary: gene counts, BGC regions, InterProScan6
# coverage, OMArk completeness/consistency, locus_tag status, ...
# -----------------------------------------------------------------------------
rule results_summary:
    input:
        gff3 = expand("results/{genome}/06_annotate/{genome}.gff3", genome=GENOMES),
        omark = expand("results/{genome}/07_omark/{genome}.sum", genome=GENOMES),
    output:
        tsv = "results/final_summary.tsv",
    params:
        genomes = GENOMES,
        genome_config = config["genomes"],
        busco_lineages = {g: busco_lineage(g) for g in GENOMES},
    log:
        "logs/results_summary.log",
    resources:
        mem_mb = resources["results_summary"]["mem_mb"],
        runtime = resources["results_summary"]["runtime"],
    script:
        "../scripts/summarize_results.py"
