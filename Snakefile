import glob
import os
import sys
from snakemake.utils import min_version

min_version("9.0")

configfile: "config/config.yaml"

# ============================================================================
#  FIGS -- Funannotate + InterProScan6 + GeneML + antiSMASH (Secondary
#  metabolites). Fungal genome annotation, the counterpart to BAGS.
# ============================================================================

# -----------------------------------------------------------------------------
# PER-RULE RESOURCES
# runtime is an INTEGER NUMBER OF MINUTES (Snakemake >=8, SLURM executor
# plugin). Partitions are not set: BioCloud assigns them from the
# mem-per-CPU ratio.
# -----------------------------------------------------------------------------
GB = 1024

resources = {
    "prepare_genome":       {"mem_mb": 4  * GB, "runtime": 120},
    "antismash_database":   {"mem_mb": 4  * GB, "runtime": 240},
    "funannotate_database": {"mem_mb": 8  * GB, "runtime": 360},
    "funannotate_clean":    {"mem_mb": 8  * GB, "runtime": 60},
    "funannotate_sort":     {"mem_mb": 4  * GB, "runtime": 30},
    "funannotate_mask":     {"mem_mb": 16 * GB, "runtime": 240},
    "get_geneml":           {"mem_mb": 8  * GB, "runtime": 120},
    "geneml":               {"mem_mb": 32 * GB, "runtime": 180},
    "antismash":            {"mem_mb": 32 * GB, "runtime": 360},
    "interproscan6_setup":  {"mem_mb": 8  * GB, "runtime": 360},
    "interproscan6":        {"mem_mb": 16 * GB, "runtime": 360},
    "funannotate_annotate": {"mem_mb": 16 * GB, "runtime": 240},
    "results_summary":      {"mem_mb": 2  * GB, "runtime": 15},
}

# -----------------------------------------------------------------------------
# CONTAINERS
# Every rule runs in a pinned image (profile: software-deployment-method:
# apptainer). The conda envs under workflow/envs/ list the same tool+version,
# for reference or for running with --software-deployment-method conda.
# InterProScan6 is the exception: it is a Nextflow pipeline that starts its
# own containers, so it runs on the host with Nextflow from the Snakemake
# env (see 6_interproscan6.smk).
# -----------------------------------------------------------------------------
NCBI_DATASETS_CONTAINER = "docker://staphb/ncbi-datasets:18.37.0"
FUNANNOTATE_CONTAINER   = "docker://quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5"
ANTISMASH_CONTAINER     = "docker://quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_1"
PYTHON_CONTAINER        = "docker://python:3.12-slim"

# The Snakemake env's bin/ (Nextflow + Java live there), put on PATH by the
# InterProScan6 rules so they don't rely on SLURM passing the login shell's
# PATH on to the job.
SNAKEMAKE_ENV_BIN = os.path.dirname(sys.executable)

# geneML is only on PyPI (no container), so rule get_geneml installs it once
# into resources/ inside the python image, and geneml runs from there.
GENEML_VERSION = config["geneml"].get("version", "1.2.0")
GENEML_ENV = f"resources/geneml-{GENEML_VERSION}"

# All funannotate subcommands (clean/sort/mask/annotate) need FUNANNOTATE_DB
# set. Rather than repeating an `export` in every rule, set it once for every
# shell command Snakemake runs. Absolute, so it means the same thing inside
# the container and from any directory a tool changes into.
FUNANNOTATE_DB = os.path.abspath(config["funannotate"]["db_dir"])
shell.prefix(f"export FUNANNOTATE_DB={FUNANNOTATE_DB}; ")

# -----------------------------------------------------------------------------
# GENOMES
# Known when the workflow is parsed: every FASTA in local_genome_dir (named
# after the file, extension stripped) plus every NCBI accession in the config.
# Each is staged as data/genomes/{genome}.fna by rule prepare_genome.
# -----------------------------------------------------------------------------
LOCAL_GENOMES = {}
_local_dir = config.get("local_genome_dir", "")
if _local_dir and os.path.isdir(_local_dir):
    for _ext in ("fna", "fa", "fasta"):
        for _f in sorted(glob.glob(os.path.join(_local_dir, f"*.{_ext}"))):
            _name = os.path.splitext(os.path.basename(_f))[0]
            if _name in LOCAL_GENOMES:
                raise WorkflowError(
                    f"Two local assemblies are both named '{_name}': "
                    f"{LOCAL_GENOMES[_name]} and {_f}. Rename one of them.")
            LOCAL_GENOMES[_name] = _f

NCBI_ACCESSIONS = list(config.get("ncbi_accessions") or [])
for _acc in NCBI_ACCESSIONS:
    if _acc in LOCAL_GENOMES:
        raise WorkflowError(
            f"'{_acc}' is both a local assembly ({LOCAL_GENOMES[_acc]}) and an "
            f"NCBI accession in ncbi_accessions. Keep only one.")

GENOMES = list(LOCAL_GENOMES) + NCBI_ACCESSIONS
if not GENOMES:
    raise WorkflowError(
        "No genomes to annotate: put assemblies in local_genome_dir and/or "
        "list accessions under ncbi_accessions in config/config.yaml.")

for _g in GENOMES:
    if _g not in config.get("genomes", {}):
        raise WorkflowError(
            f"Genome '{_g}' has no entry under 'genomes:' in config/config.yaml "
            f"(funannotate annotate needs at least its species).")

# Keep {genome} from swallowing path separators, e.g. matching
# "GCA_1/01_preprocess/GCA_1" in results/{genome}/... patterns.
wildcard_constraints:
    genome = r"[^/]+",

# -----------------------------------------------------------------------------
# TARGET RULE
# -----------------------------------------------------------------------------
rule all:
    input:
        expand("results/{genome}/06_annotate/{genome}.annotate.done", genome=GENOMES),
        "results/final_summary.tsv",

# -----------------------------------------------------------------------------
# INCLUDES
# -----------------------------------------------------------------------------
include: "workflow/rules/0_prepare_genomes.smk"
include: "workflow/rules/1_antismash_database.smk"
include: "workflow/rules/2_funannotate_database.smk"
include: "workflow/rules/3_funannotate_preprocess.smk"
include: "workflow/rules/4_geneml.smk"
include: "workflow/rules/5_antismash.smk"
include: "workflow/rules/6_interproscan6.smk"
include: "workflow/rules/7_funannotate_annotate.smk"
include: "workflow/rules/8_results_summary.smk"
