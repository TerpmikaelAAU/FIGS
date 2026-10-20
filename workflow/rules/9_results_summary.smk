# -----------------------------------------------------------------------------
# Per genome, one row per protein with its functional annotation
# (06_annotate/<genome>.annotations.tsv), read from funannotate2's .tbl.
# -----------------------------------------------------------------------------
rule annotation_table:
    input:
        tbl = rules.funannotate2_annotate.output.tbl,
    output:
        tsv = "results/{genome}/06_annotate/{genome}.annotations.tsv",
    log:
        "logs/{genome}/annotation_table.log",
    resources:
        mem_mb = resources["annotation_table"]["mem_mb"],
        runtime = resources["annotation_table"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        python workflow/scripts/annotation_table.py -i {input.tbl} -o {output.tsv} > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# One-line-per-genome summary: gene counts, BGC regions, InterProScan6
# coverage, OMArk completeness/consistency, locus_tag status, ...
# Written for all genomes (results/final_summary.tsv) and, as a single row,
# into each genome's folder (results/<genome>/<genome>.final_summary.tsv).
# -----------------------------------------------------------------------------
rule results_summary:
    input:
        gff3 = expand("results/{genome}/06_annotate/{genome}.gff3", genome=GENOMES),
        omark = expand("results/{genome}/07_omark/{genome}.sum", genome=GENOMES),
    output:
        tsv = "results/final_summary.tsv",
        per_genome = expand("results/{genome}/{genome}.final_summary.tsv", genome=GENOMES),
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
