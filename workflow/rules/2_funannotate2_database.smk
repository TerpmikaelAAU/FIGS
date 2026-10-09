# -----------------------------------------------------------------------------
# One-off, shared funannotate2 database setup (Pfam, dbCAN, MEROPS,
# UniProtKB/Swiss-Prot, GO, MIBiG, InterPro, gene2product). Every
# funannotate2 subcommand elsewhere in the workflow points FUNANNOTATE2_DB at
# this same directory (set once via shell.prefix() in the Snakefile).
#
# `install -d all` would also build the RefSeq mitochondrial database, which
# needs minimap2 (not in the funannotate2 container) and is only used by
# `funannotate2 predict` to set mitochondrial contigs aside. FIGS predicts
# genes with geneML instead, so that database is left out.
#
# funannotate2 annotate downloads its BUSCO lineage into FUNANNOTATE2_DB the
# first time it needs it. Every lineage the config uses (BUSCO_LINEAGES) is
# fetched here instead, so parallel annotate jobs never unpack the same
# lineage into the shared directory at the same time.
# -----------------------------------------------------------------------------
rule funannotate2_database:
    output:
        marker = touch("data/databases/funannotate2db"),
    params:
        db_dir = FUNANNOTATE2_DB,
        lineages = " ".join(BUSCO_LINEAGES),
    log:
        "logs/funannotate2_database.log",
    resources:
        mem_mb = resources["funannotate2_database"]["mem_mb"],
        runtime = resources["funannotate2_database"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        mkdir -p {params.db_dir}
        funannotate2 install -d merops uniprot dbCAN pfam go mibig interpro gene2product > {log} 2>&1

        for lineage in {params.lineages}; do
            python -c "import logging, sys
logging.basicConfig(level=logging.INFO)
from funannotate2.utilities import ensure_busco_lineage
print(ensure_busco_lineage(sys.argv[1], logging.getLogger('busco')))" "$lineage" >> {log} 2>&1
        done
        """
