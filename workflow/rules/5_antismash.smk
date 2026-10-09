# -----------------------------------------------------------------------------
# antiSMASH, seeded with the gene models via --genefinding-gff3 so it doesn't
# re-call genes itself (--genefinding-tool none). Produces a GenBank file
# with BGC annotations, which rule external_annotations turns into
# funannotate2 annotation tables.
# -----------------------------------------------------------------------------
rule antismash:
    input:
        genome = rules.softmask.output,
        gff3 = rules.gene_models.output.gff3,
        db_marker = "data/databases/antismashdatabase",
    output:
        out_dir = directory("results/{genome}/03_antismash"),
        gbk = "results/{genome}/03_antismash/{genome}.gbk",
    params:
        db_dir = config["antismash"]["database_dir"],
        extra = config["antismash"].get("extra", ""),
    threads: 8
    log:
        "logs/{genome}/antismash.log",
    resources:
        mem_mb = resources["antismash"]["mem_mb"],
        runtime = resources["antismash"]["runtime"],
    container:
        ANTISMASH_CONTAINER
    conda:
        "../envs/antismash.yaml"
    shell:
        r"""
        antismash {input.genome} \
            --taxon fungi \
            --cpus {threads} \
            --databases {params.db_dir} \
            --genefinding-tool none \
            --genefinding-gff3 {input.gff3} \
            --output-dir {output.out_dir} \
            --output-basename {wildcards.genome} \
            {params.extra} \
            > {log} 2>&1
        """
