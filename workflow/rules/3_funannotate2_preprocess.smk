# -----------------------------------------------------------------------------
# Assembly preprocessing ahead of gene prediction: funannotate2 clean drops
# short and duplicated contigs, sorts them largest first and renames them
# <genome>_1, <genome>_2, ... (what funannotate clean + sort did). The
# assembly is then soft-masked with pytantan through funannotate2's own
# softmask_fasta, the same tantan masking `funannotate mask` used by default
# (funannotate2 has no mask command). geneML and antiSMASH both run on the
# soft-masked assembly.
#
# Genomes with keep_contigs (always the case for an NCBI "update", whose
# sequences must match the deposited ones) skip clean: the assembly is used
# as staged, with only the FASTA descriptions dropped.
# -----------------------------------------------------------------------------
rule funannotate2_clean:
    input:
        "data/genomes/{genome}.fna",
    output:
        "results/{genome}/01_preprocess/{genome}.clean.fna",
    params:
        minlen = config["funannotate2"]["clean_minlen"],
        # clean deletes its tmpdir when it finishes, so each genome gets its own.
        tmpdir = "results/{genome}/01_preprocess/clean_tmp",
        keep_contigs = lambda wc: "yes" if keep_contigs(wc.genome) else "",
    threads: 4
    log:
        "logs/{genome}/funannotate2_clean.log",
    resources:
        mem_mb = resources["funannotate2_clean"]["mem_mb"],
        runtime = resources["funannotate2_clean"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        if [ -n "{params.keep_contigs}" ]; then
            awk '/^>/ {{ print $1; next }} {{ print }}' {input} > {output}
            echo "keep_contigs: {input} used as is, funannotate2 clean skipped" > {log}
            exit 0
        fi
        funannotate2 clean \
            -f {input} \
            -o {output} \
            --minlen {params.minlen} \
            --rename {wildcards.genome}_ \
            --cpus {threads} \
            --tmpdir {params.tmpdir} \
            > {log} 2>&1
        """


rule softmask:
    input:
        rules.funannotate2_clean.output,
    output:
        "results/{genome}/01_preprocess/{genome}.masked.fna",
    log:
        "logs/{genome}/softmask.log",
    resources:
        mem_mb = resources["softmask"]["mem_mb"],
        runtime = resources["softmask"]["runtime"],
    container:
        FUNANNOTATE2_CONTAINER
    conda:
        "../envs/funannotate2.yaml"
    shell:
        r"""
        python -c "import sys
from funannotate2.fastx import softmask_fasta
softmask_fasta(sys.argv[1], sys.argv[2])" {input} {output} > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# One row per sequence of the cleaned assembly: its ID, the ID it gets in an
# NCBI submission, its chromosome name and, for organelles, its location.
# gene_models leaves organelle sequences out of the gene models, and
# ncbi_submission uses the rest. For an NCBI "update" this comes from the
# GenBank records (needs network access); otherwise IDs are kept and
# organelles come from the config.
# -----------------------------------------------------------------------------
rule sequence_info:
    input:
        rules.funannotate2_clean.output,
    output:
        "results/{genome}/01_preprocess/{genome}.seqinfo.tsv",
    params:
        ncbi_update = lambda wc: f"--ncbi-update {wc.genome}" if ncbi_mode(wc.genome) == "update" else "",
        organelles = lambda wc: " ".join(
            f"--organelle {seqid}={location}"
            for seqid, location in (config["genomes"][wc.genome].get("organelles") or {}).items()),
    log:
        "logs/{genome}/sequence_info.log",
    resources:
        mem_mb = resources["sequence_info"]["mem_mb"],
        runtime = resources["sequence_info"]["runtime"],
    container:
        PYTHON_CONTAINER
    shell:
        r"""
        python workflow/scripts/sequence_info.py \
            --fasta {input} \
            {params.ncbi_update} \
            {params.organelles} \
            -o {output} \
            > {log} 2>&1
        """
