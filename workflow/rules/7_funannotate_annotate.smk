# -----------------------------------------------------------------------------
# funannotate annotate: pulls PFAM/InterPro/EggNog/UniProtKB/MEROPS/CAZyme/GO
# annotation onto geneML's gene models, folding in the antiSMASH BGC calls
# and InterProScan6 domain hits. If config["genomes"][genome]["locus_tag"] is
# set, gene models are renamed to NCBI-style locus_tags via --rename; leave
# it blank to keep geneML's native GML###### IDs.
# -----------------------------------------------------------------------------
rule funannotate_annotate:
    input:
        masked = rules.funannotate_mask.output,
        gff3 = rules.geneml.output.gff3,
        antismash_gbk = rules.antismash.output.gbk,
        ips_xml = rules.interproscan6.output.xml,
        db_marker = "data/databases/funannotatedb",
    output:
        out_dir = directory("results/{genome}/06_annotate"),
        done = touch("results/{genome}/06_annotate/{genome}.annotate.done"),
    params:
        species = lambda wc: config["genomes"][wc.genome]["species"],
        strain = lambda wc: config["genomes"][wc.genome].get("strain", ""),
        locus_tag = lambda wc: config["genomes"][wc.genome].get("locus_tag", ""),
        sbt = config.get("sbt_template", ""),
        busco_db = config["funannotate"]["busco_db"],
        extra = config["funannotate"].get("extra", ""),
    threads: 8
    log:
        "logs/{genome}/funannotate_annotate.log",
    resources:
        mem_mb = resources["funannotate_annotate"]["mem_mb"],
        runtime = resources["funannotate_annotate"]["runtime"],
    container:
        FUNANNOTATE_CONTAINER
    conda:
        "../envs/funannotate.yaml"
    shell:
        r"""
        strain_flag=()
        [ -n "{params.strain}" ] && strain_flag=(--strain "{params.strain}")

        rename_flag=()
        [ -n "{params.locus_tag}" ] && rename_flag=(--rename {params.locus_tag})

        sbt_flag=()
        [ -n "{params.sbt}" ] && [ -f "{params.sbt}" ] && sbt_flag=(--sbt {params.sbt})

        funannotate annotate \
            --gff {input.gff3} \
            --fasta {input.masked} \
            -s "{params.species}" \
            "${{strain_flag[@]}}" \
            --antismash {input.antismash_gbk} \
            --iprscan {input.ips_xml} \
            --busco_db {params.busco_db} \
            "${{rename_flag[@]}}" \
            "${{sbt_flag[@]}}" \
            --cpus {threads} \
            -o {output.out_dir} \
            {params.extra} \
            > {log} 2>&1
        """
