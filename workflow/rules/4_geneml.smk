# -----------------------------------------------------------------------------
# geneML: deep-learning fungal gene prediction (https://github.com/hexagonbio/geneML)
# Runs on the soft-masked genome and emits GFF3 + gene/protein FASTA.
# gene_models then turns geneML's GFF3 into the gene set every later step
# shares (antiSMASH, InterProScan6, funannotate2 annotate, OMArk), so their
# results all key on the same transcript IDs.
#
# geneML is only on PyPI, so get_geneml installs the pinned version
# (geneml.version in the config) once into resources/geneml-<version>/ inside
# the official python image, and geneml runs from there in that same image --
# the same approach as the T2T pipeline's BAGS rules.
# -----------------------------------------------------------------------------
rule get_geneml:
    output:
        geneml = protected(f"{GENEML_ENV}/bin/geneml"),
    params:
        env = GENEML_ENV,
        version = GENEML_VERSION,
    threads: 2
    log:
        "logs/get_geneml.log",
    resources:
        mem_mb = resources["get_geneml"]["mem_mb"],
        runtime = resources["get_geneml"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        rm -rf {params.env}
        python -m venv {params.env} > {log} 2>&1
        {params.env}/bin/pip install --no-cache-dir "geneml=={params.version}" >> {log} 2>&1
        {output.geneml} --version >> {log} 2>&1
        """


rule geneml:
    input:
        geneml = f"{GENEML_ENV}/bin/geneml",
        genome = rules.softmask.output,
    output:
        gff3 = "results/{genome}/02_geneml/{genome}.geneml.gff3",
        genes = "results/{genome}/02_geneml/{genome}.geneml.genes.fna",
        proteins = "results/{genome}/02_geneml/{genome}.geneml.proteins.faa",
    params:
        # gene-id-prefix must be alphanumeric-ish; sanitize accessions like
        # GCA_052058355.1 -> GCA_052058355_1
        prefix = lambda wc: "".join(c if c.isalnum() else "_" for c in wc.genome),
        min_gene_score = config["geneml"]["min_gene_score"],
        extra = config["geneml"].get("extra", ""),
    threads: 8
    log:
        "logs/{genome}/geneml.log",
    resources:
        mem_mb = resources["geneml"]["mem_mb"],
        runtime = resources["geneml"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        {input.geneml} {input.genome} \
            -o {output.gff3} \
            -g {output.genes} \
            -p {output.proteins} \
            --gene-id-prefix {params.prefix} \
            --min-gene-score {params.min_gene_score} \
            --max-transcripts 1 \
            --cpu-only \
            -c {threads} \
            {params.extra} \
            > {log} 2>&1
        """


# With a locus_tag in the config, gfftk rename gives genes NCBI-style IDs
# (<locus_tag>_000001, transcripts <locus_tag>_000001-T1), numbered by
# position; this replaces funannotate 1's `annotate --rename`. Without one,
# gfftk sanitize keeps geneML's own IDs. Either way the protein FASTA is
# written from the final GFF3, so its headers are the final transcript IDs.
rule gene_models:
    input:
        genome = rules.softmask.output,
        gff3 = rules.geneml.output.gff3,
    output:
        gff3 = "results/{genome}/02_geneml/{genome}.models.gff3",
        proteins = "results/{genome}/02_geneml/{genome}.models.proteins.faa",
    params:
        locus_tag = lambda wc: (config["genomes"][wc.genome].get("locus_tag") or "").strip(),
    log:
        "logs/{genome}/gene_models.log",
    resources:
        mem_mb = resources["gene_models"]["mem_mb"],
        runtime = resources["gene_models"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        if [ -n "{params.locus_tag}" ]; then
            gfftk rename -f {input.genome} -g {input.gff3} -l {params.locus_tag} -o {output.gff3} > {log} 2>&1
        else
            gfftk sanitize -f {input.genome} -g {input.gff3} -o {output.gff3} > {log} 2>&1
        fi
        gfftk convert -f {input.genome} -i {output.gff3} --output-format proteins -n -o {output.proteins} >> {log} 2>&1
        """
