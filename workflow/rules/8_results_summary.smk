# -----------------------------------------------------------------------------
# One-line-per-genome summary (gene count, BGC count, locus_tag status, ...)
# -----------------------------------------------------------------------------
rule results_summary:
    input:
        expand("results/{genome}/06_annotate/{genome}.annotate.done", genome=GENOMES),
    output:
        tsv = "results/final_summary.tsv",
    params:
        genomes = GENOMES,
        genome_config = config["genomes"],
    log:
        "logs/results_summary.log",
    resources:
        mem_mb = resources["results_summary"]["mem_mb"],
        runtime = resources["results_summary"]["runtime"],
    script:
        "../scripts/summarize_results.py"
