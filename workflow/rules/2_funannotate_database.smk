# -----------------------------------------------------------------------------
# One-off, shared funannotate database setup (PFAM, InterPro, UniProtKB,
# MEROPS, CAZyme, GO, BUSCO lineages, ...). Every funannotate subcommand
# elsewhere in the workflow points FUNANNOTATE_DB at this same directory
# (set once via shell.prefix() in the Snakefile).
# -----------------------------------------------------------------------------
rule funannotate_database:
    output:
        marker = touch("data/databases/funannotatedb"),
    params:
        db_dir = FUNANNOTATE_DB,
        busco_db = config["funannotate"]["busco_db"],
    log:
        "logs/funannotate_database.log",
    resources:
        mem_mb = resources["funannotate_database"]["mem_mb"],
        runtime = resources["funannotate_database"]["runtime"],
    container:
        FUNANNOTATE_CONTAINER
    conda:
        "../envs/funannotate.yaml"
    shell:
        r"""
        mkdir -p {params.db_dir}
        funannotate setup -d {params.db_dir} --busco_db {params.busco_db} > {log} 2>&1
        """
