# -----------------------------------------------------------------------------
# Assembly preprocessing ahead of gene prediction: funannotate2 clean drops
# short and duplicated contigs, sorts them largest first and renames them
# <genome>_1, <genome>_2, ... (what funannotate clean + sort did). The
# assembly is then soft-masked with pytantan through funannotate2's own
# softmask_fasta, the same tantan masking `funannotate mask` used by default
# (funannotate2 has no mask command). geneML and antiSMASH both run on the
# soft-masked assembly.
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
