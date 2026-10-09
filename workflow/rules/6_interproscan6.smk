# -----------------------------------------------------------------------------
# InterProScan6 (https://github.com/ebi-pf-team/interproscan6), a Nextflow
# pipeline that starts its own container for every analysis step. That is
# why this is the one rule WITHOUT a `container:`: wrapping it in an
# Apptainer container would mean Apptainer-in-Apptainer. It runs on the
# host instead, with Nextflow + Java from the Snakemake env
# (Snakemake_env.yml), and Nextflow pulls InterProScan6's images itself
# with `-profile apptainer` (interproscan6.profile in the config; use
# `docker` on a machine with a Docker daemon).
#
# interproscan6_setup runs once before any genome: it fetches the pipeline,
# its images and its databases, so parallel genome jobs don't all try to
# download them into the same place at the same time. Every genome then
# gets its own Nextflow launch directory, so concurrent jobs
# don't share one .nextflow/ cache and work/ dir. The images are pulled
# once into interproscan6.container_cache and shared by every genome.
#
# NOTE: funannotate's --iprscan flag was written against InterProScan5's XML
# schema. InterProScan6 is a very recent rewrite and its XML output has not
# been exhaustively verified against funannotate's parser here - if
# `funannotate annotate` (rule 7) errors while reading {genome}.xml, compare
# it against an InterProScan5 XML example and adjust/convert as needed.
# -----------------------------------------------------------------------------
rule interproscan6_setup:
    output:
        marker = touch("data/databases/interproscan6_ready"),
    params:
        version = config["interproscan6"]["version"],
        profile = config["interproscan6"].get("profile", "apptainer"),
        datadir = os.path.abspath(config["interproscan6"]["datadir"]),
        cache = os.path.abspath(config["interproscan6"].get(
            "container_cache", "data/databases/interproscan6_containers")),
        launch_dir = os.path.abspath("data/interproscan6_runs/setup"),
        log = os.path.abspath("logs/interproscan6_setup.log"),
    threads: 4
    log:
        "logs/interproscan6_setup.log",
    resources:
        mem_mb = resources["interproscan6_setup"]["mem_mb"],
        runtime = resources["interproscan6_setup"]["runtime"],
    shell:
        r"""
        export PATH={SNAKEMAKE_ENV_BIN}:$PATH
        export NXF_APPTAINER_CACHEDIR={params.cache}
        export NXF_SINGULARITY_CACHEDIR={params.cache}
        mkdir -p {params.datadir} {params.cache}
        rm -rf {params.launch_dir}
        mkdir -p {params.launch_dir}
        cd {params.launch_dir}

        nextflow pull ebi-pf-team/interproscan6 -r {params.version} > {params.log} 2>&1
        # InterProScan6's own quick test: annotates its bundled test proteins,
        # which downloads the databases into datadir and pulls the images.
        nextflow run ebi-pf-team/interproscan6 \
            -r {params.version} \
            -profile {params.profile},test \
            --datadir {params.datadir} \
            --outdir {params.launch_dir}/test_output \
            >> {params.log} 2>&1

        cd - > /dev/null
        rm -rf {params.launch_dir}
        """


rule interproscan6:
    input:
        proteins = rules.geneml.output.proteins,
        ready = "data/databases/interproscan6_ready",
    output:
        out_dir = directory("results/{genome}/04_interproscan"),
        xml = "results/{genome}/04_interproscan/{genome}.xml",
    params:
        version = config["interproscan6"]["version"],
        profile = config["interproscan6"].get("profile", "apptainer"),
        # Absolute: Nextflow runs from its own launch directory (below).
        proteins = lambda wc: os.path.abspath(f"results/{wc.genome}/02_geneml/{wc.genome}.geneml.proteins.faa"),
        out_dir = lambda wc: os.path.abspath(f"results/{wc.genome}/04_interproscan"),
        log = lambda wc: os.path.abspath(f"logs/{wc.genome}/interproscan6.log"),
        launch_dir = lambda wc: os.path.abspath(f"data/interproscan6_runs/{wc.genome}"),
        datadir = os.path.abspath(config["interproscan6"]["datadir"]),
        cache = os.path.abspath(config["interproscan6"].get(
            "container_cache", "data/databases/interproscan6_containers")),
        formats = config["interproscan6"].get("formats", "xml,tsv,json,gff3"),
        extra = config["interproscan6"].get("extra", ""),
    threads: 8
    log:
        "logs/{genome}/interproscan6.log",
    resources:
        mem_mb = resources["interproscan6"]["mem_mb"],
        runtime = resources["interproscan6"]["runtime"],
    shell:
        r"""
        export PATH={SNAKEMAKE_ENV_BIN}:$PATH
        mkdir -p {params.datadir} {params.cache} {params.out_dir}
        rm -rf {params.launch_dir}
        mkdir -p {params.launch_dir}
        cd {params.launch_dir}

        export NXF_APPTAINER_CACHEDIR={params.cache}
        export NXF_SINGULARITY_CACHEDIR={params.cache}

        nextflow run ebi-pf-team/interproscan6 \
            -r {params.version} \
            -profile {params.profile} \
            --datadir {params.datadir} \
            --input {params.proteins} \
            --outdir {params.out_dir} \
            --outprefix {wildcards.genome} \
            --formats {params.formats} \
            {params.extra} \
            > {params.log} 2>&1

        # The work/ dir holds a copy of every intermediate; the results are
        # already in out_dir.
        cd - > /dev/null
        rm -rf {params.launch_dir}
        """
