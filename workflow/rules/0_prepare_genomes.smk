# -----------------------------------------------------------------------------
# Stage one genome into data/genomes/{genome}.fna: stage_genome copies a
# local assembly (a genome with `fasta`) in, download_genome downloads an
# NCBI accession with `datasets download genome accession`. The list of
# genomes is worked out when the workflow is parsed (GENOMES in the
# Snakefile), so every genome is its own job and a dry run shows the whole
# DAG.
#
# stage_genome runs without a container: the assembly can be anywhere, and
# a container only sees the working directory unless other paths are bound
# into it. Every later rule reads the staged copy.
# -----------------------------------------------------------------------------
NCBI_GENOME_NAMES = [_g for _g in GENOMES if _g not in LOCAL_GENOMES]


def _one_of(names):
    """A wildcard constraint matching exactly these names (or nothing)."""
    return "|".join(re.escape(n) for n in names) or "(?!)"


rule stage_genome:
    input:
        lambda wc: LOCAL_GENOMES[wc.genome],
    output:
        "data/genomes/{genome}.fna",
    wildcard_constraints:
        genome = _one_of(LOCAL_GENOMES),
    log:
        "logs/{genome}/stage_genome.log",
    resources:
        mem_mb = resources["stage_genome"]["mem_mb"],
        runtime = resources["stage_genome"]["runtime"],
    shell:
        r"""
        cp {input} {output} 2> {log}
        echo "staged {input} -> {output}" >> {log}
        """


rule download_genome:
    output:
        "data/genomes/{genome}.fna",
    wildcard_constraints:
        genome = _one_of(NCBI_GENOME_NAMES),
    params:
        tmp = "data/downloads/ncbi/{genome}",
    log:
        "logs/{genome}/download_genome.log",
    resources:
        mem_mb = resources["download_genome"]["mem_mb"],
        runtime = resources["download_genome"]["runtime"],
    container:
        NCBI_DATASETS_CONTAINER
    conda:
        "../envs/ncbi_datasets.yaml"
    shell:
        r"""
        rm -rf {params.tmp}
        mkdir -p {params.tmp}
        datasets download genome accession {wildcards.genome} \
            --include genome \
            --filename {params.tmp}/{wildcards.genome}.zip > {log} 2>&1
        unzip -q -o {params.tmp}/{wildcards.genome}.zip -d {params.tmp}
        fna=$(find {params.tmp}/ncbi_dataset/data/{wildcards.genome} -name '*.fna' | sort | head -n 1)
        if [ -z "$fna" ]; then
            echo "no FASTA found for {wildcards.genome}" >> {log}
            exit 1
        fi
        mv "$fna" {output}
        rm -rf {params.tmp}
        echo "downloaded {wildcards.genome} -> {output}" >> {log}
        """
