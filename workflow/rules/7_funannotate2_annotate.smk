# -----------------------------------------------------------------------------
# funannotate2 annotate: Pfam, dbCAN (CAZymes), MEROPS, UniProtKB/Swiss-Prot
# and BUSCO searches on the gene models, merged with the antiSMASH,
# InterProScan6 and eggNOG-mapper results, then written out as GFF3, GenBank
# (table2asn), TBL and protein/transcript FASTA, with the tRNA and rRNA genes.
#
# funannotate2 takes outside results as 3-column tables (transcript ID,
# annotation type, value). funannotate2-addons (f2a) converts the
# InterProScan6 TSV and the antiSMASH GenBank into that format, and
# fix_annotations.py then fixes what f2a 26.3.7 gets wrong for this setup
# (see the script):
#   * antiSMASH results are keyed by gene ID; funannotate2 wants transcript IDs
#   * GO terms are written as `go_term`, which annotate ignores (`go_terms`)
#   * InterProScan6 GO terms keep their (source) suffix
#   * InterPro's Reactome pathways for every species (millions of notes,
#     nearly all from animals) are passed on; they are dropped
# eggNOG-mapper's table is written by emapper_annotations.py instead of f2a,
# whose eggNOG descriptions would never become products (see the script).
# -----------------------------------------------------------------------------
rule get_funannotate2_addons:
    output:
        f2a = protected(f"{F2A_ENV}/bin/f2a"),
    params:
        env = F2A_ENV,
        version = F2A_VERSION,
    threads: 2
    log:
        "logs/get_funannotate2_addons.log",
    resources:
        mem_mb = resources["get_funannotate2_addons"]["mem_mb"],
        runtime = resources["get_funannotate2_addons"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        rm -rf {params.env}
        python -m venv {params.env} > {log} 2>&1
        {params.env}/bin/pip install --no-cache-dir "funannotate2-addons=={params.version}" >> {log} 2>&1
        {output.f2a} --version >> {log} 2>&1
        """


rule external_annotations:
    input:
        f2a = f"{F2A_ENV}/bin/f2a",
        gff3 = rules.gene_models.output.gff3,
        ips_tsv = rules.interproscan6.output.tsv,
        antismash_gbk = rules.antismash.output.gbk,
        emapper = rules.eggnog_mapper.output.annotations,
    output:
        iprscan = "results/{genome}/05_annotations/{genome}.iprscan.annotations.txt",
        antismash = "results/{genome}/05_annotations/{genome}.antismash.annotations.txt",
        clusters = "results/{genome}/05_annotations/{genome}.antismash.clusters.txt",
        emapper = "results/{genome}/05_annotations/{genome}.emapper.annotations.txt",
    params:
        out_dir = "results/{genome}/05_annotations",
        raw_dir = "results/{genome}/05_annotations/f2a",
    log:
        "logs/{genome}/external_annotations.log",
    resources:
        mem_mb = resources["external_annotations"]["mem_mb"],
        runtime = resources["external_annotations"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        rm -rf {params.raw_dir}
        mkdir -p {params.raw_dir}

        # f2a logs parsing errors but still exits 0, so check for its output
        # (the clusters file is legitimately empty if antiSMASH found no BGCs).
        {input.f2a} iprscan --parse {input.ips_tsv} -o {params.raw_dir} > {log} 2>&1
        {input.f2a} antismash --parse {input.antismash_gbk} -o {params.raw_dir} >> {log} 2>&1
        for f in iprscan.annotations.txt antismash.annotations.txt antismash.clusters.txt; do
            if [ ! -f {params.raw_dir}/$f ]; then
                echo "ERROR: f2a did not write {params.raw_dir}/$f" >> {log}
                exit 1
            fi
        done

        python workflow/scripts/fix_annotations.py --gff3 {input.gff3} \
            -i {params.raw_dir}/iprscan.annotations.txt -o {output.iprscan} >> {log} 2>&1
        python workflow/scripts/fix_annotations.py --gff3 {input.gff3} \
            -i {params.raw_dir}/antismash.annotations.txt -o {output.antismash} >> {log} 2>&1
        cp {params.raw_dir}/antismash.clusters.txt {output.clusters}

        python workflow/scripts/emapper_annotations.py -i {input.emapper} \
            -o {params.raw_dir}/emapper.annotations.txt >> {log} 2>&1
        python workflow/scripts/fix_annotations.py --gff3 {input.gff3} \
            -i {params.raw_dir}/emapper.annotations.txt -o {output.emapper} >> {log} 2>&1
        """


def _annotate_slug(genome):
    # funannotate2's naming_slug(): how annotate names its result files.
    meta = config["genomes"][genome]
    slug = meta["species"].replace(" ", "_").capitalize()
    strain = (meta.get("strain") or "").replace(" ", "")
    return f"{slug}_{strain}" if strain else slug


rule funannotate2_annotate:
    input:
        genome = rules.softmask.output,
        # the gene models plus the tRNA/rRNA genes (4b_ncrna.smk)
        gff3 = rules.ncrna_models.output.gff3,
        iprscan = rules.external_annotations.output.iprscan,
        antismash = rules.external_annotations.output.antismash,
        emapper = rules.external_annotations.output.emapper,
        db_marker = "data/databases/funannotate2db",
    output:
        gff3 = "results/{genome}/06_annotate/{genome}.gff3",
        gbk = "results/{genome}/06_annotate/{genome}.gbk",
        tbl = "results/{genome}/06_annotate/{genome}.tbl",
        proteins = "results/{genome}/06_annotate/{genome}.proteins.fa",
        transcripts = "results/{genome}/06_annotate/{genome}.transcripts.fa",
        summary = "results/{genome}/06_annotate/{genome}.summary.json",
        # Gene names/products funannotate2's name cleaner could not make
        # NCBI-safe; worth a look before an NCBI submission.
        curate = "results/{genome}/06_annotate/{genome}.need-curating.txt",
    params:
        out_dir = "results/{genome}/06_annotate",
        tmpdir = "results/{genome}/06_annotate/tmp",
        species = lambda wc: config["genomes"][wc.genome]["species"],
        strain = lambda wc: config["genomes"][wc.genome].get("strain") or "",
        slug = lambda wc: _annotate_slug(wc.genome),
        busco_lineage = lambda wc: busco_lineage(wc.genome),
        extra = config["funannotate2"].get("extra", ""),
    threads: 8
    log:
        "logs/{genome}/funannotate2_annotate.log",
    resources:
        mem_mb = resources["funannotate2_annotate"]["mem_mb"],
        runtime = resources["funannotate2_annotate"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        # annotate reuses any search results already in annotate_misc, which
        # would be stale if the gene models changed since the last run.
        rm -rf {params.out_dir}/annotate_misc {params.out_dir}/annotate_results {params.tmpdir}
        mkdir -p {params.tmpdir}

        strain_flag=()
        [ -n "{params.strain}" ] && strain_flag=(--strain "{params.strain}")

        funannotate2 annotate \
            -f {input.genome} \
            -g {input.gff3} \
            -o {params.out_dir} \
            -a {input.iprscan} {input.antismash} {input.emapper} \
            -s "{params.species}" \
            "${{strain_flag[@]}}" \
            --busco-lineage {params.busco_lineage} \
            --cpus {threads} \
            --tmpdir {params.tmpdir} \
            {params.extra} \
            > {log} 2>&1

        res={params.out_dir}/annotate_results/{params.slug}
        cp "$res.gff3"           {output.gff3}
        cp "$res.gbk"            {output.gbk}
        cp "$res.tbl"            {output.tbl}
        cp "$res.proteins.fa"    {output.proteins}
        cp "$res.transcripts.fa" {output.transcripts}
        cp "$res.summary.json"   {output.summary}
        curate={params.out_dir}/annotate_results/Gene2Products.need-curating.txt
        if [ -f "$curate" ]; then cp "$curate" {output.curate}; else : > {output.curate}; fi
        rm -rf {params.tmpdir}
        """
