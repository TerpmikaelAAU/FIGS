import os
import re
from snakemake.utils import min_version

min_version("9.0")

# Tool settings have defaults in workflow/defaults.yaml; config/config.yaml
# (the genomes, plus any setting to change) is merged over them.
configfile: "workflow/defaults.yaml"
configfile: "config/config.yaml"

# ============================================================================
#  FIGS -- Funannotate2 + InterProScan6 + GeneML + antiSMASH (Secondary
#  metabolites). Fungal genome annotation, the counterpart to BAGS, with
#  OMArk as a quality check of the final gene set.
# ============================================================================

# -----------------------------------------------------------------------------
# PER-RULE RESOURCES
# runtime is an INTEGER NUMBER OF MINUTES (Snakemake >=8, SLURM executor
# plugin). Partitions are not set: BioCloud assigns them from the
# mem-per-CPU ratio.
# -----------------------------------------------------------------------------
GB = 1024

resources = {
    "stage_genome":            {"mem_mb": 1  * GB, "runtime": 30},
    "download_genome":         {"mem_mb": 4  * GB, "runtime": 120},
    "sequence_info":           {"mem_mb": 2  * GB, "runtime": 30},
    "antismash_database":      {"mem_mb": 16 * GB, "runtime": 240},
    "funannotate2_database":   {"mem_mb": 8  * GB, "runtime": 360},
    "funannotate2_clean":      {"mem_mb": 8  * GB, "runtime": 60},
    "softmask":                {"mem_mb": 8  * GB, "runtime": 120},
    "get_geneml":              {"mem_mb": 8  * GB, "runtime": 120},
    "geneml":                  {"mem_mb": 32 * GB, "runtime": 180},
    "gene_models":             {"mem_mb": 4  * GB, "runtime": 30},
    "antismash":               {"mem_mb": 32 * GB, "runtime": 360},
    "get_nextflow":            {"mem_mb": 2  * GB, "runtime": 30},
    "interproscan6_setup":     {"mem_mb": 8  * GB, "runtime": 360},
    "interproscan6":           {"mem_mb": 32 * GB, "runtime": 360},
    "get_funannotate2_addons": {"mem_mb": 4  * GB, "runtime": 60},
    "external_annotations":    {"mem_mb": 4  * GB, "runtime": 30},
    "funannotate2_annotate":   {"mem_mb": 16 * GB, "runtime": 240},
    "omark_database":          {"mem_mb": 4  * GB, "runtime": 1440},
    "omark":                   {"mem_mb": 32 * GB, "runtime": 120},
    "get_table2asn":           {"mem_mb": 2  * GB, "runtime": 30},
    "ncbi_submission":         {"mem_mb": 16 * GB, "runtime": 120},
    "results_summary":         {"mem_mb": 2  * GB, "runtime": 15},
}

# -----------------------------------------------------------------------------
# CONTAINERS
# Every rule runs in a pinned image (profile: software-deployment-method:
# apptainer). The conda envs under workflow/envs/ list the same tool+version,
# for reference or for running with --software-deployment-method conda.
# InterProScan6 is the exception: it is a Nextflow pipeline that starts its
# own containers, so it runs on the host with a pinned Nextflow + Java
# (NEXTFLOW_DIR below; see 6_interproscan6.smk).
# -----------------------------------------------------------------------------
NCBI_DATASETS_CONTAINER = "docker://staphb/ncbi-datasets:18.37.0"
FUNANNOTATE2_CONTAINER  = "docker://quay.io/biocontainers/funannotate2:26.6.21--pyhdfd78af_0"
ANTISMASH_CONTAINER     = "docker://quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_1"
OMARK_CONTAINER         = "docker://quay.io/biocontainers/omark:0.5.0--pyhdfd78af_0"
PYTHON_CONTAINER        = "docker://python:3.12-slim"

# Nextflow + Java (Eclipse Temurin) for InterProScan6, pinned releases
# installed once into resources/ by rule get_nextflow, so InterProScan6
# doesn't depend on the env Snakemake is started from.
NEXTFLOW_DIR = os.path.abspath(f"resources/nextflow-{config['nextflow']['version']}")

# geneML is only on PyPI (no container), so rule get_geneml installs it once
# into resources/ inside the python image, and geneml runs from there.
GENEML_VERSION = config["geneml"].get("version", "1.2.0")
GENEML_ENV = f"resources/geneml-{GENEML_VERSION}"

# funannotate2-addons (f2a) parses the InterProScan6 and antiSMASH results
# into funannotate2's annotation format. Also PyPI-only, so it gets the same
# treatment as geneML: rule get_funannotate2_addons installs it into resources/.
F2A_VERSION = config["funannotate2"].get("addons_version", "26.3.7")
F2A_ENV = f"resources/funannotate2-addons-{F2A_VERSION}"

# All funannotate2 subcommands need FUNANNOTATE2_DB set. Rather than repeating
# an `export` in every rule, set it once for every shell command Snakemake
# runs. Absolute, so it means the same thing inside the container and from
# any directory a tool changes into. shell.prefix() replaces Snakemake's
# default "set -euo pipefail; " prefix, so that is repeated here.
FUNANNOTATE2_DB = os.path.abspath(config["funannotate2"]["db_dir"])
shell.prefix(f"set -euo pipefail; export FUNANNOTATE2_DB={FUNANNOTATE2_DB}; ")

# BUSCO lineage per genome: genomes.<genome>.busco_lineage if set, else
# funannotate2.busco_lineage. Rule funannotate2_database downloads each one
# up front, so parallel annotate jobs don't race to unpack the same lineage
# into the shared database directory.
def busco_lineage(genome):
    return (config["genomes"][genome].get("busco_lineage")
            or config["funannotate2"].get("busco_lineage", "fungi"))

# OMArk's OMAmer database (LUCA.h5, ~10 GB), shared by every genome.
OMARK_DB = os.path.join(config["omark"]["db_dir"], "LUCA.h5")

# table2asn for the NCBI submission: NCBI's own pinned release, installed
# once into resources/ by rule get_table2asn (bioconda lags a release behind).
TABLE2ASN_VERSION = config["ncbi"]["table2asn_version"]
TABLE2ASN = f"resources/table2asn-{TABLE2ASN_VERSION}/table2asn"

# NCBI submission mode per genome: "" (none), "update" or "new" (see config).
def ncbi_mode(genome):
    return (config["genomes"][genome].get("ncbi_submission") or "").strip()

# Genomes whose contigs are used exactly as staged (funannotate2 clean is
# skipped): asked for with keep_contigs, and always for an NCBI update,
# whose sequences must match the deposited ones.
def keep_contigs(genome):
    return ncbi_mode(genome) == "update" or bool(config["genomes"][genome].get("keep_contigs"))

# -----------------------------------------------------------------------------
# GENOMES
# One per entry under genomes: in the config. An entry with `fasta` is a
# local assembly; one without is downloaded from NCBI, with its key as the
# assembly accession. Each is staged as data/genomes/{genome}.fna by rule
# stage_genome (local) or download_genome (NCBI).
# -----------------------------------------------------------------------------
GENOMES = list(config.get("genomes") or {})
if not GENOMES:
    raise WorkflowError(
        "No genomes to annotate: add an entry per genome under 'genomes:' "
        "in config/config.yaml.")

LOCAL_GENOMES = {}
for _g in GENOMES:
    _meta = config["genomes"][_g] or {}
    if not (_meta.get("species") or "").strip():
        raise WorkflowError(f"Genome '{_g}' has no species in config/config.yaml.")
    _fasta = (_meta.get("fasta") or "").strip()
    if _fasta:
        if not os.path.isfile(_fasta):
            raise WorkflowError(f"Genome '{_g}': fasta {_fasta} does not exist.")
        LOCAL_GENOMES[_g] = _fasta
    elif not re.fullmatch(r"GC[AF]_\d{9}\.\d+", _g):
        raise WorkflowError(
            f"Genome '{_g}' has no fasta, and '{_g}' is not an NCBI assembly "
            f"accession (GCA_/GCF_...) to download.")

BUSCO_LINEAGES = sorted({busco_lineage(_g) for _g in GENOMES})

# NCBI submission settings are checked here, so mistakes stop the run at
# start-up rather than after hours of annotation.
for _g in GENOMES:
    _meta = config["genomes"][_g]
    _mode = ncbi_mode(_g)
    if _mode not in ("", "update", "new"):
        raise WorkflowError(
            f"Genome '{_g}': ncbi_submission must be \"update\", \"new\" or "
            f"\"\", not \"{_mode}\".")
    if _mode == "update" and _g in LOCAL_GENOMES:
        raise WorkflowError(
            f"Genome '{_g}': ncbi_submission \"update\" adds annotation to an "
            f"assembly already in GenBank, so it must be downloaded from NCBI: "
            f"remove fasta and use the assembly accession as the genome's name.")
    _lt = (_meta.get("locus_tag") or "").strip()
    if _mode and not _lt:
        raise WorkflowError(
            f"Genome '{_g}' is set up for NCBI submission but has no locus_tag; "
            f"use the prefix registered to its BioProject/BioSample.")
    for _key in ("bioproject", "biosample"):
        if _mode and not (_meta.get(_key) or "").strip():
            raise WorkflowError(
                f"Genome '{_g}' is set up for NCBI submission but has no {_key}.")
    if _lt and not re.fullmatch(r"[A-Za-z][A-Za-z0-9]{2,11}", _lt):
        raise WorkflowError(
            f"Genome '{_g}': locus_tag '{_lt}' is not a valid NCBI prefix "
            f"(3-12 letters/digits, starting with a letter).")
    if _meta.get("organelles") and not keep_contigs(_g):
        raise WorkflowError(
            f"Genome '{_g}': organelles needs keep_contigs: true, since "
            f"funannotate2 clean renames the contigs.")

NCBI_GENOMES = [_g for _g in GENOMES if ncbi_mode(_g)]
# A warning rather than an error, so dry runs (and CI) work without one;
# only from the main process, not again in every SLURM job's log.
if (workflow.is_main_process and NCBI_GENOMES
        and not os.path.isfile(config["ncbi"]["sbt_template"])):
    logger.warning(
        f"No NCBI submission template at {config['ncbi']['sbt_template']}; "
        f"ncbi_submission will fail for {', '.join(NCBI_GENOMES)} until you make "
        f"one at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/")

# Keep {genome} from swallowing path separators, e.g. matching
# "GCA_1/01_preprocess/GCA_1" in results/{genome}/... patterns.
wildcard_constraints:
    genome = r"[^/]+",

# -----------------------------------------------------------------------------
# TARGET RULE
# -----------------------------------------------------------------------------
rule all:
    input:
        expand("results/{genome}/06_annotate/{genome}.gff3", genome=GENOMES),
        expand("results/{genome}/07_omark/{genome}.sum", genome=GENOMES),
        "results/final_summary.tsv",
        expand("results/{genome}/08_ncbi/{genome}.sqn", genome=NCBI_GENOMES),

# -----------------------------------------------------------------------------
# INCLUDES
# -----------------------------------------------------------------------------
include: "workflow/rules/0_prepare_genomes.smk"
include: "workflow/rules/1_antismash_database.smk"
include: "workflow/rules/2_funannotate2_database.smk"
include: "workflow/rules/3_funannotate2_preprocess.smk"
include: "workflow/rules/4_geneml.smk"
include: "workflow/rules/5_antismash.smk"
include: "workflow/rules/6_interproscan6.smk"
include: "workflow/rules/7_funannotate2_annotate.smk"
include: "workflow/rules/8_omark.smk"
include: "workflow/rules/9_results_summary.smk"
include: "workflow/rules/10_ncbi_submission.smk"
