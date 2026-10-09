# FIGS: fungal genome annotation

**F**unannotate + **I**nterProScan6 + **G**eneML + anti**S**MASH, as a
Snakemake 9 pipeline. The annotation counterpart to
[BAGS](https://github.com/TerpmikaelAAU/BAGS), set up the same way as the
[T2T pipeline](https://github.com/TerpmikaelAAU/T2T_fungal_genome_pipeline)
(Snakemake 9, SLURM executor plugin, one Apptainer container per rule).

It takes your own assembled fungal genomes and/or NCBI accessions and runs:

```
prepare genome (local copy or NCBI datasets download)
        │
funannotate clean → sort → mask   (repeat-soft-masked assembly)
        │
     geneML                       (deep-learning gene prediction → GFF3 + proteins)
        │            │
 antiSMASH   InterProScan6 (Nextflow)
        │            │
        └── funannotate annotate ──┘   (functional annotation, optional locus_tag rename)
        │
  final_summary.tsv
```

## Directory layout

```
Snakefile
Snakemake_env.yml      # Snakemake 9.26.1 + SLURM plugin + Nextflow
slurm_submit.sbatch    # runs the whole workflow as one small SLURM job
config/
  config.yaml          # edit this
  template.sbt         # optional, for NCBI --rename/--sbt submission prep
profile/
  config.yaml          # BioCloud SLURM executor profile
workflow/
  rules/               # 0_ ... 8_, one stage per file
  envs/                # tool+version reference (see "Containers")
  scripts/             # summarize_results.py
data/
  local_assemblies/    # drop your own *.fna/*.fa/*.fasta here (optional)
  genomes/             # staged, unified input (created by the pipeline)
  databases/           # antiSMASH / funannotate / InterProScan6 DBs
resources/             # geneML install (created by the pipeline)
results/
  <genome>/
    01_preprocess/
    02_geneml/
    03_antismash/
    04_interproscan/
    06_annotate/
  final_summary.tsv
```

## Installation

**1. Create the Snakemake conda environment.**
```bash
conda env create -f Snakemake_env.yml
conda activate snakemake_figs
```
This pins Snakemake 9.26.1 (the same as the T2T pipeline), the native SLURM
executor plugin, and Nextflow + Java for InterProScan6. The pipeline needs
Snakemake >= 9; it will not run under Snakemake 7.

**2. Apptainer** must be available on the compute nodes (it is on BioCloud).
Every rule's container is pulled automatically the first time it's needed.

## Configure

Edit `config/config.yaml`:

- `local_genome_dir` / `ncbi_accessions` - your genome inputs. Both can be
  used together. The genome list is read when Snakemake starts, so a dry run
  (`-n`) shows every genome's jobs.
- `genomes:` - per-genome `species`, `strain`, and **optional** `locus_tag`.
  Every genome needs an entry here (at least `species`); Snakemake stops at
  start-up and names the genome if one is missing. Set `locus_tag` to get
  `funannotate annotate --rename LOCUSTAG` (NCBI-style gene renaming); leave
  it blank to keep geneML's native IDs.
- `sbt_template` - only needed alongside `locus_tag` if you intend to
  submit to NCBI (`--sbt`); get one from NCBI's
  [template generator](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/).
- `geneml.version` - geneML release installed from PyPI (default `1.2.0`).
- `antismash.database_dir` / `interproscan6.datadir` / `funannotate.db_dir` -
  where the (large, shared) reference databases live. Downloaded once by
  rules `antismash_database` / `funannotate_database`; InterProScan6 fetches
  its own data into `datadir` on first run.
- `interproscan6.profile` - the Nextflow profile InterProScan6 runs its own
  containers with: `apptainer` (default, for the cluster), `singularity`, or
  `docker` on a machine with a Docker daemon.
- `interproscan6.container_cache` - where those images are pulled to, once,
  shared by every genome.

## Containers

Every rule runs in one pinned image, set under `container:` in its rule file
(image names at the top of the `Snakefile`):

| Rules | Image |
|---|---|
| `prepare_genome` | `staphb/ncbi-datasets:18.37.0` (NCBI `datasets` is not on bioconda) |
| `funannotate_*` | `quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5` |
| `antismash*` | `quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_1` |
| `get_geneml`, `geneml` | `python:3.12-slim`; geneML is only on PyPI, so `get_geneml` installs it once into `resources/geneml-<version>/` |
| `interproscan6` | none, see below |
| `results_summary` | none (plain Python, runs in the Snakemake env) |

**InterProScan6** is a Nextflow pipeline that starts its own container for
every analysis step, so it is the one rule without a `container:` (wrapping
it would mean Apptainer-in-Apptainer). It runs on the compute node with the
Nextflow from `Snakemake_env.yml`, which the SLURM jobs inherit from the
environment Snakemake was started in. Each genome gets its own Nextflow
launch directory under `data/interproscan6_runs/`, removed after a
successful run. Before the first genome, rule `interproscan6_setup` runs
InterProScan6's own test once, which fetches the pipeline, pulls its images
into `container_cache` and downloads its databases into `datadir`, so
parallel genome jobs don't all download them at the same time.

`workflow/envs/*.yaml` lists the same tool versions as conda environments,
as a reference only: the profile uses Apptainer, and geneML (PyPI) and
InterProScan6 (Nextflow) have no conda env.

Databases and outputs are inside the working directory, which Snakemake binds
into every container. Anything you point at outside it -- `local_genome_dir`,
`funannotate.db_dir`, `antismash.database_dir` -- must be bound too, e.g.
`--apptainer-args "--bind /shared/databases,/path/to/assemblies"` (or
`apptainer-args:` in `profile/config.yaml`).

## Run

**Dry run first**, to see what will run:
```bash
snakemake --profile profile -n
```

**Run the workflow.** `slurm_submit.sbatch` wraps it in a single lightweight
SLURM job (1 CPU, 4G); each *rule* still gets its own appropriately sized job
via the executor plugin:
```bash
sbatch slurm_submit.sbatch
squeue --me
```
Or run it on the login node inside `tmux`:
```bash
tmux new -s figs
conda activate snakemake_figs
snakemake --profile profile
```

**If a run dies and leaves the working directory locked:**
```bash
snakemake --profile profile --unlock
```

**Resources** are set per rule in the `resources` dict at the top of the
`Snakefile`. `runtime` is an integer number of minutes (Snakemake >= 8).

## Notes / caveats

- **InterProScan6 + funannotate compatibility**: `funannotate annotate
  --iprscan` was written against InterProScan5's XML output. InterProScan6
  is a recent rewrite; its XML schema is expected to be compatible but this
  hasn't been exhaustively verified here. If `funannotate annotate` (rule
  `funannotate_annotate`) chokes on `{genome}.xml`, diff it against an
  InterProScan5 example and add a small conversion step in
  `workflow/rules/6_interproscan6.smk` if needed.
- **geneML gene-ID prefixes** are derived automatically from the `{genome}`
  wildcard (non-alphanumeric characters replaced with `_`), so NCBI
  accessions like `GCA_052058355.1` become a valid `--gene-id-prefix`.
- `funannotate mask` defaults to `tantan` for soft-masking. If you want
  RepeatMasker/RepeatModeler instead, add the relevant flags via
  `funannotate.extra` in the config, or edit
  `workflow/rules/3_funannotate_preprocess.smk` directly.
- The antiSMASH and funannotate database rules download large, shared
  reference data (several GB+) - run those once and reuse `database_dir`
  across projects rather than re-downloading per genome.

## Automatic checks

`.github/workflows/dry-run.yml` byte-compiles `workflow/scripts/*.py` and
runs `snakemake -n` with Snakemake 9.26.1 on every push and pull request.
