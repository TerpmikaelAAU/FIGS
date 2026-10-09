# -----------------------------------------------------------------------------
# One-off, shared antiSMASH database download. All per-genome antismash rules
# depend on the marker file this rule produces. PYTHONUNBUFFERED so the
# log shows how far it got even if the job is killed.
# -----------------------------------------------------------------------------
rule antismash_database:
    output:
        marker = touch("data/databases/antismashdatabase"),
    params:
        db_dir = config["antismash"]["database_dir"],
    log:
        "logs/antismash_database.log",
    resources:
        mem_mb = resources["antismash_database"]["mem_mb"],
        runtime = resources["antismash_database"]["runtime"],
    container:
        ANTISMASH_CONTAINER
    conda:
        "../envs/antismash.yaml"
    shell:
        r"""
        mkdir -p {params.db_dir}
        PYTHONUNBUFFERED=1 download-antismash-databases --database-dir {params.db_dir} > {log} 2>&1
        """
