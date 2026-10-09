# -----------------------------------------------------------------------------
# InterProScan6 (https://github.com/ebi-pf-team/interproscan6), a Nextflow
# pipeline that starts its own container for every analysis step. That is
# why this is the one rule WITHOUT a `container:`: wrapping it in an
# Apptainer container would mean Apptainer-in-Apptainer. It runs on the
# host instead, with the pinned Nextflow + Java that get_nextflow installs
# into resources/, and Nextflow pulls InterProScan6's images itself
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
# Nextflow runs every task with its local executor inside this one SLURM
# job, so --max-workers caps the parallel tasks at the job's threads.
# The InterPro data release is pinned (interproscan6.interpro_version) and is
# part of the setup marker's name, so changing it re-runs the setup.
#
# GO terms and pathways are OFF by default in InterProScan6, so they are
# switched on here (interproscan6.goterms / .pathways in the config).
# funannotate2 reads the TSV (via rule external_annotations), not the XML:
# the XML groups identical protein sequences under one entry, and the
# parser only picks up the first ID of such a group.
# -----------------------------------------------------------------------------
IPS6_READY = (f"data/databases/interproscan6_{config['interproscan6']['version']}"
              f"_interpro{config['interproscan6']['interpro_version']}.ready")

# Put the pinned Java and Nextflow first on PATH in the InterProScan6 rules.
NEXTFLOW_ENV = f"export JAVA_HOME={NEXTFLOW_DIR}/jdk; export PATH={NEXTFLOW_DIR}/jdk/bin:{NEXTFLOW_DIR}:$PATH"


rule get_nextflow:
    output:
        nextflow = protected(f"{NEXTFLOW_DIR}/nextflow"),
        java = protected(f"{NEXTFLOW_DIR}/jdk/bin/java"),
    params:
        dir = NEXTFLOW_DIR,
        url = config["nextflow"]["url"],
        sha256 = config["nextflow"]["sha256"],
        java_url = config["nextflow"]["java_url"],
        java_sha256 = config["nextflow"]["java_sha256"],
    log:
        "logs/get_nextflow.log",
    resources:
        mem_mb = resources["get_nextflow"]["mem_mb"],
        runtime = resources["get_nextflow"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        rm -rf {params.dir}/jdk
        mkdir -p {params.dir}/jdk
        python workflow/scripts/download.py "{params.java_url}" {params.dir}/jdk.tar.gz \
            --sha256 {params.java_sha256} > {log} 2>&1
        tar -xzf {params.dir}/jdk.tar.gz -C {params.dir}/jdk --strip-components 1
        rm {params.dir}/jdk.tar.gz
        python workflow/scripts/download.py "{params.url}" {output.nextflow} \
            --sha256 {params.sha256} >> {log} 2>&1
        chmod +x {output.nextflow}
        {output.java} -version >> {log} 2>&1
        """


rule interproscan6_setup:
    input:
        nextflow = rules.get_nextflow.output.nextflow,
        java = rules.get_nextflow.output.java,
    output:
        marker = touch(IPS6_READY),
    params:
        version = config["interproscan6"]["version"],
        interpro = config["interproscan6"]["interpro_version"],
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
        {NEXTFLOW_ENV}
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
            --interpro {params.interpro} \
            --max-workers {threads} \
            --datadir {params.datadir} \
            --outdir {params.launch_dir}/test_output \
            >> {params.log} 2>&1

        cd - > /dev/null
        rm -rf {params.launch_dir}
        """


rule interproscan6:
    input:
        proteins = rules.gene_models.output.proteins,
        ready = IPS6_READY,
        nextflow = rules.get_nextflow.output.nextflow,
    output:
        out_dir = directory("results/{genome}/04_interproscan"),
        tsv = "results/{genome}/04_interproscan/{genome}.tsv",
    params:
        version = config["interproscan6"]["version"],
        interpro = config["interproscan6"]["interpro_version"],
        profile = config["interproscan6"].get("profile", "apptainer"),
        # Absolute: Nextflow runs from its own launch directory (below).
        proteins = lambda wc: os.path.abspath(f"results/{wc.genome}/02_geneml/{wc.genome}.models.proteins.faa"),
        out_dir = lambda wc: os.path.abspath(f"results/{wc.genome}/04_interproscan"),
        log = lambda wc: os.path.abspath(f"logs/{wc.genome}/interproscan6.log"),
        launch_dir = lambda wc: os.path.abspath(f"data/interproscan6_runs/{wc.genome}"),
        datadir = os.path.abspath(config["interproscan6"]["datadir"]),
        cache = os.path.abspath(config["interproscan6"].get(
            "container_cache", "data/databases/interproscan6_containers")),
        # tsv is what the rest of the workflow reads, so it is always written.
        formats = ",".join(dict.fromkeys(
            ["tsv"] + [f.strip() for f in config["interproscan6"].get("formats", "tsv").split(",") if f.strip()]
        )),
        go_pathways = " ".join(
            flag for flag, key in (("--goterms", "goterms"), ("--pathways", "pathways"))
            if config["interproscan6"].get(key, True)
        ),
        extra = config["interproscan6"].get("extra", ""),
    threads: 8
    log:
        "logs/{genome}/interproscan6.log",
    resources:
        mem_mb = resources["interproscan6"]["mem_mb"],
        runtime = resources["interproscan6"]["runtime"],
    shell:
        r"""
        {NEXTFLOW_ENV}
        mkdir -p {params.datadir} {params.cache} {params.out_dir}
        rm -rf {params.launch_dir}
        mkdir -p {params.launch_dir}
        cd {params.launch_dir}

        export NXF_APPTAINER_CACHEDIR={params.cache}
        export NXF_SINGULARITY_CACHEDIR={params.cache}

        nextflow run ebi-pf-team/interproscan6 \
            -r {params.version} \
            -profile {params.profile} \
            --interpro {params.interpro} \
            --max-workers {threads} \
            --datadir {params.datadir} \
            --input {params.proteins} \
            --outdir {params.out_dir} \
            --outprefix {wildcards.genome} \
            --formats {params.formats} \
            {params.go_pathways} \
            {params.extra} \
            > {params.log} 2>&1

        # The work/ dir holds a copy of every intermediate; the results are
        # already in out_dir.
        cd - > /dev/null
        rm -rf {params.launch_dir}
        """
