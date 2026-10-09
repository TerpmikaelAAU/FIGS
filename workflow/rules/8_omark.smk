# -----------------------------------------------------------------------------
# OMArk (https://github.com/DessimozLab/OMArk): proteome quality check
# against the OMA hierarchical orthologous groups. Beyond BUSCO-style
# completeness (conserved HOGs that are single/duplicated/missing), it reports
# how many of ALL predicted proteins are taxonomically consistent, fragmented,
# of unknown origin, or likely contamination - i.e. it also flags spurious
# gene models, which completeness alone can't.
#
# OMAmer first places every protein in a HOG, then OMArk scores the result.
# It runs on the gene models (the same models funannotate2 annotates), so it
# doesn't wait for annotate. geneML runs with --max-transcripts 1, so there
# is one protein per gene and no isoform file is needed.
#
# omark_database fetches, once, the OMAmer database (LUCA.h5, ~10 GB, pinned
# to a Zenodo release in the config and checked against its md5) and the
# NCBI taxonomy OMArk reads through ete3. Without -e, ete3 would build that
# taxonomy in ~/.etetoolkit on first use, with parallel jobs racing to do so.
# -----------------------------------------------------------------------------
rule omark_database:
    output:
        db = protected(OMARK_DB),
        taxa = protected(os.path.join(config["omark"]["db_dir"], "taxa.sqlite")),
    params:
        url = config["omark"]["db_url"],
        md5 = config["omark"].get("db_md5", ""),
        db_dir = config["omark"]["db_dir"],
        log = os.path.abspath("logs/omark_database.log"),
    log:
        "logs/omark_database.log",
    resources:
        mem_mb = resources["omark_database"]["mem_mb"],
        runtime = resources["omark_database"]["runtime"],
    container:
        OMARK_CONTAINER
    conda:
        "../envs/omark.yaml"
    shell:
        r"""
        mkdir -p {params.db_dir}
        # Resumes after dropped connections (and across re-runs, from
        # LUCA.h5.part); an empty md5 skips the check.
        python workflow/scripts/download.py "{params.url}" {output.db} \
            --md5 "{params.md5}" > {log} 2>&1

        # ete3 downloads taxdump.tar.gz into the working directory.
        cd {params.db_dir}
        python -c "from ete3 import NCBITaxa
NCBITaxa(dbfile='taxa.sqlite')" >> {params.log} 2>&1
        rm -f taxdump.tar.gz
        """


rule omark:
    input:
        proteins = rules.gene_models.output.proteins,
        db = OMARK_DB,
        taxa = os.path.join(config["omark"]["db_dir"], "taxa.sqlite"),
    output:
        omamer = "results/{genome}/07_omark/{genome}.omamer",
        # OMArk names its files after the OMAmer file: {genome}.sum etc.
        summary = "results/{genome}/07_omark/{genome}.sum",
    params:
        out_dir = "results/{genome}/07_omark",
        extra = config["omark"].get("extra", ""),
    threads: 8
    log:
        "logs/{genome}/omark.log",
    resources:
        mem_mb = resources["omark"]["mem_mb"],
        runtime = resources["omark"]["runtime"],
    container:
        OMARK_CONTAINER
    conda:
        "../envs/omark.yaml"
    shell:
        r"""
        omamer search \
            --db {input.db} \
            --query {input.proteins} \
            --out {output.omamer} \
            --nthreads {threads} \
            > {log} 2>&1

        omark \
            -f {output.omamer} \
            -d {input.db} \
            -e {input.taxa} \
            -of {input.proteins} \
            -o {params.out_dir} \
            {params.extra} \
            >> {log} 2>&1
        """
