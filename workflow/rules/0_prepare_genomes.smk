# -----------------------------------------------------------------------------
# Stage one genome into data/genomes/{genome}.fna: a local assembly is copied
# in, an NCBI accession is downloaded with `datasets download genome
# accession`. The list of genomes is worked out when the workflow is parsed
# (GENOMES in the Snakefile), so every genome is its own job and a dry run
# shows the whole DAG.
# -----------------------------------------------------------------------------
rule prepare_genome:
    input:
        # The local assembly, so a changed file re-stages it; nothing for an
        # NCBI accession.
        lambda wc: LOCAL_GENOMES.get(wc.genome, []),
    output:
        "data/genomes/{genome}.fna",
    params:
        tmp = "data/downloads/ncbi/{genome}",
    log:
        "logs/{genome}/prepare_genome.log",
    resources:
        mem_mb = resources["prepare_genome"]["mem_mb"],
        runtime = resources["prepare_genome"]["runtime"],
    container:
        NCBI_DATASETS_CONTAINER
    conda:
        "../envs/ncbi_datasets.yaml"
    shell:
        r"""
        if [ -n "{input}" ]; then
            cp {input} {output}
            echo "[local] staged {input} -> {output}" > {log}
        else
            rm -rf {params.tmp}
            mkdir -p {params.tmp}
            datasets download genome accession {wildcards.genome} \
                --include genome \
                --filename {params.tmp}/{wildcards.genome}.zip > {log} 2>&1
            unzip -q -o {params.tmp}/{wildcards.genome}.zip -d {params.tmp}
            fna=$(find {params.tmp}/ncbi_dataset/data/{wildcards.genome} -name '*.fna' | sort | head -n 1)
            if [ -z "$fna" ]; then
                echo "[NCBI] no FASTA found for {wildcards.genome}" >> {log}
                exit 1
            fi
            mv "$fna" {output}
            rm -rf {params.tmp}
            echo "[NCBI] staged {wildcards.genome} -> {output}" >> {log}
        fi
        """
