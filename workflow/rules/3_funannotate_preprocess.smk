# -----------------------------------------------------------------------------
# Standard funannotate preprocessing (clean duplicated/short contigs, sort +
# rename headers, soft-mask repeats) ahead of gene prediction. geneML and
# antiSMASH both then run on the soft-masked assembly.
# -----------------------------------------------------------------------------
rule funannotate_clean:
    input:
        "data/genomes/{genome}.fna",
    output:
        "results/{genome}/01_preprocess/{genome}.clean.fna",
    params:
        minlen = config["funannotate"]["clean_minlen"],
    log:
        "logs/{genome}/funannotate_clean.log",
    resources:
        mem_mb = resources["funannotate_clean"]["mem_mb"],
        runtime = resources["funannotate_clean"]["runtime"],
    container:
        FUNANNOTATE_CONTAINER
    conda:
        "../envs/funannotate.yaml"
    shell:
        "funannotate clean -i {input} -o {output} --minlen {params.minlen} > {log} 2>&1"


rule funannotate_sort:
    input:
        rules.funannotate_clean.output,
    output:
        "results/{genome}/01_preprocess/{genome}.sorted.fna",
    log:
        "logs/{genome}/funannotate_sort.log",
    resources:
        mem_mb = resources["funannotate_sort"]["mem_mb"],
        runtime = resources["funannotate_sort"]["runtime"],
    container:
        FUNANNOTATE_CONTAINER
    conda:
        "../envs/funannotate.yaml"
    shell:
        "funannotate sort -i {input} -o {output} -b {wildcards.genome} > {log} 2>&1"


rule funannotate_mask:
    input:
        rules.funannotate_sort.output,
    output:
        "results/{genome}/01_preprocess/{genome}.masked.fna",
    threads: 8
    log:
        "logs/{genome}/funannotate_mask.log",
    resources:
        mem_mb = resources["funannotate_mask"]["mem_mb"],
        runtime = resources["funannotate_mask"]["runtime"],
    container:
        FUNANNOTATE_CONTAINER
    conda:
        "../envs/funannotate.yaml"
    shell:
        "funannotate mask -i {input} -o {output} --cpus {threads} > {log} 2>&1"
