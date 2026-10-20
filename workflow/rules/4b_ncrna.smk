# -----------------------------------------------------------------------------
# tRNA and rRNA genes, which geneML (protein-coding only) does not predict.
#
#   trnascan     tRNAscan-SE 2 in its default eukaryotic mode, as funannotate
#                1 and funannotate2 predict run it.
#   rrna_search  Infernal cmsearch with Rfam's eukaryotic rRNA models (5S,
#                5.8S, SSU/18S, LSU/28S) and Rfam's gathering thresholds, the
#                approach of Rfam's own genome annotation and of NCBI's
#                eukaryotic annotation pipeline.
#   ncrna_models adds them to the protein-coding gene models for funannotate2
#                annotate (see merge_ncrna.py).
#
# Both tools are in the funannotate2 container (tRNAscan-SE needs Infernal),
# so there is no extra image. The RNA genes are added after the protein
# searches, numbered on from the last protein-coding gene, so adding them
# does not change any protein's ID or re-run antiSMASH, InterProScan6,
# eggNOG-mapper or OMArk. ncrna_gff3.py explains which predictions are kept.
# -----------------------------------------------------------------------------
RFAM_RRNA_MODELS = "5S_rRNA 5_8S_rRNA SSU_rRNA_eukarya LSU_rRNA_eukarya"


rule rfam_rrna_models:
    output:
        cm = f"{config['ncrna']['rfam_dir']}/Rfam{config['ncrna']['rfam_version']}.rRNA.cm",
    params:
        rfam = lambda wc, output: output.cm.replace(".rRNA.cm", ".cm"),
        url = config["ncrna"]["rfam_url"],
        sha256 = config["ncrna"]["rfam_sha256"],
        models = RFAM_RRNA_MODELS,
    log:
        "logs/rfam_rrna_models.log",
    resources:
        mem_mb = resources["rfam_rrna_models"]["mem_mb"],
        runtime = resources["rfam_rrna_models"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output.cm})
        python workflow/scripts/download.py "{params.url}" {params.rfam}.gz \
            --sha256 {params.sha256} > {log} 2>&1
        gunzip -c {params.rfam}.gz > {params.rfam}
        cmfetch --index {params.rfam} >> {log} 2>&1
        for m in {params.models}; do cmfetch {params.rfam} $m; done > {output.cm}.tmp 2>> {log}
        mv {output.cm}.tmp {output.cm}
        rm -f {params.rfam} {params.rfam}.gz {params.rfam}.ssi
        grep "^NAME" {output.cm} >> {log}
        """


rule trnascan:
    input:
        genome = rules.softmask.output,
    output:
        raw = "results/{genome}/02_ncrna/{genome}.tRNAscan.out",
        gff3 = "results/{genome}/02_ncrna/{genome}.tRNA.gff3",
    params:
        extra = config["ncrna"].get("trnascan_extra", ""),
    threads: 4
    log:
        "logs/{genome}/trnascan.log",
    resources:
        mem_mb = resources["trnascan"]["mem_mb"],
        runtime = resources["trnascan"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        rm -f {output.raw}
        tRNAscan-SE -o {output.raw} --thread {threads} -Q {params.extra} \
            {input.genome} > {log} 2>&1
        python workflow/scripts/ncrna_gff3.py trna -i {output.raw} -o {output.gff3} >> {log} 2>&1
        """


rule rrna_search:
    input:
        genome = rules.softmask.output,
        cm = rules.rfam_rrna_models.output.cm,
    output:
        tblout = "results/{genome}/02_ncrna/{genome}.rRNA.tblout",
        gff3 = "results/{genome}/02_ncrna/{genome}.rRNA.gff3",
    params:
        min_coverage = config["ncrna"]["rrna_min_coverage"],
        extra = config["ncrna"].get("cmsearch_extra", ""),
    threads: 8
    log:
        "logs/{genome}/rrna_search.log",
    resources:
        mem_mb = resources["rrna_search"]["mem_mb"],
        runtime = resources["rrna_search"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        cmsearch --cpu {threads} --cut_ga --rfam --nohmmonly --noali \
            --tblout {output.tblout} -o /dev/null {params.extra} \
            {input.cm} {input.genome} > {log} 2>&1
        python workflow/scripts/ncrna_gff3.py rrna -i {output.tblout} --cm {input.cm} \
            --min-coverage {params.min_coverage} -o {output.gff3} >> {log} 2>&1
        """


rule ncrna_models:
    input:
        genome = rules.softmask.output,
        models = rules.gene_models.output.gff3,
        rrna = rules.rrna_search.output.gff3,
        trna = rules.trnascan.output.gff3,
        seqinfo = rules.sequence_info.output,
    output:
        gff3 = "results/{genome}/02_ncrna/{genome}.models_with_rna.gff3",
    params:
        # the same ID prefix geneml uses
        prefix = lambda wc: "".join(c if c.isalnum() else "_" for c in wc.genome),
        locus_tag = lambda wc: (config["genomes"][wc.genome].get("locus_tag") or "").strip(),
        merged = "results/{genome}/02_ncrna/{genome}.merged.gff3",
    log:
        "logs/{genome}/ncrna_models.log",
    resources:
        mem_mb = resources["ncrna_models"]["mem_mb"],
        runtime = resources["ncrna_models"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        python workflow/scripts/merge_ncrna.py --models {input.models} \
            --rna {input.rrna} {input.trna} --seqinfo {input.seqinfo} \
            --prefix {params.prefix} --locus-tag "{params.locus_tag}" \
            -o {params.merged} > {log} 2>&1
        gfftk sanitize -f {input.genome} -g {params.merged} -o {output.gff3} >> {log} 2>&1
        rm -f {params.merged}
        """
