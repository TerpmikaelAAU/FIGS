# -----------------------------------------------------------------------------
# eggNOG-mapper (https://github.com/eggnogdb/eggnog-mapper): orthology-based
# gene names and product descriptions. Next to UniProtKB/Swiss-Prot (which
# funannotate2 annotate searches itself, with strict cut-offs) it is where
# most gene names come from, as it was in funannotate 1. Its results reach
# funannotate2 through rule external_annotations.
#
# eggnog_database downloads the eggNOG 5.0.2 data once (~12 GB, several
# times that unpacked). eggNOG-mapper's own download_eggnog_data.py is not used:
# up to 2.1.15 it fetches from eggnogdb.embl.de, which no longer exists.
# The files have no published checksums; download.py checks each one is
# complete, and gunzip checks its CRC on unpacking. Files already unpacked
# are skipped on a re-run.
# -----------------------------------------------------------------------------
rule eggnog_database:
    output:
        marker = touch("data/databases/eggnogdb"),
    params:
        db_dir = config["eggnog"]["db_dir"],
        url = config["eggnog"]["db_url"],
    log:
        "logs/eggnog_database.log",
    resources:
        mem_mb = resources["eggnog_database"]["mem_mb"],
        runtime = resources["eggnog_database"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        mkdir -p {params.db_dir}
        : > {log}
        # Unpack to a temporary name, so a killed job never leaves a
        # half-unpacked file that the next run would take as done.
        for f in eggnog.db eggnog_proteins.dmnd; do
            [ -s {params.db_dir}/$f ] && continue
            [ -f {params.db_dir}/$f.gz ] || \
                python workflow/scripts/download.py "{params.url}/$f.gz" {params.db_dir}/$f.gz >> {log} 2>&1
            gunzip -c {params.db_dir}/$f.gz > {params.db_dir}/$f.tmp
            mv {params.db_dir}/$f.tmp {params.db_dir}/$f
            rm {params.db_dir}/$f.gz
        done
        if [ ! -s {params.db_dir}/eggnog.taxa.db ]; then
            python workflow/scripts/download.py "{params.url}/eggnog.taxa.tar.gz" \
                {params.db_dir}/eggnog.taxa.tar.gz >> {log} 2>&1
            tar -xzf {params.db_dir}/eggnog.taxa.tar.gz -C {params.db_dir}
            rm {params.db_dir}/eggnog.taxa.tar.gz
        fi
        ls -l {params.db_dir} >> {log}
        """


rule eggnog_mapper:
    input:
        proteins = rules.gene_models.output.proteins,
        db = rules.eggnog_database.output.marker,
    output:
        annotations = "results/{genome}/04_eggnog/{genome}.emapper.annotations",
    params:
        db_dir = config["eggnog"]["db_dir"],
        out_dir = "results/{genome}/04_eggnog",
        tmpdir = "results/{genome}/04_eggnog/tmp",
        extra = config["eggnog"].get("extra", ""),
    threads: 16
    log:
        "logs/{genome}/eggnog_mapper.log",
    resources:
        mem_mb = resources["eggnog_mapper"]["mem_mb"],
        runtime = resources["eggnog_mapper"]["runtime"],
    container:
        EGGNOG_CONTAINER
    conda:
        "../envs/eggnog_mapper.yaml"
    shell:
        r"""
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}
        emapper.py -i {input.proteins} --itype proteins -m diamond \
            --data_dir {params.db_dir} \
            --output {wildcards.genome} --output_dir {params.out_dir} \
            --temp_dir {params.tmpdir} \
            --cpu {threads} --override \
            {params.extra} \
            > {log} 2>&1
        rm -rf {params.tmpdir}
        """
