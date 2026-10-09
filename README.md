# FIGS: fungal genome annotation

**F**unannotate2 + **I**nterProScan6 + **G**eneML + anti**S**MASH, with
[OMArk](https://github.com/DessimozLab/OMArk) quality control, as a
Snakemake 9 pipeline. The annotation counterpart to
[BAGS](https://github.com/TerpmikaelAAU/BAGS), set up the same way as the
[T2T pipeline](https://github.com/TerpmikaelAAU/T2T_fungal_genome_pipeline)
(Snakemake 9, SLURM executor plugin, one Apptainer container per rule).

It takes your own assembled fungal genomes and/or NCBI accessions and runs:

```
prepare genome (local copy or NCBI datasets download)
        │
funannotate2 clean → soft-mask (tantan)
        │
     geneML               (deep-learning gene prediction)
        │
  gene models             (gfftk: optional locus_tag rename, protein FASTA)
        │
   ┌────┼───────────────┬─────────────┐
antiSMASH  InterProScan6 (Nextflow)   OMArk (completeness + consistency)
   │            │                     │
   └─ f2a + fix_annotations ─┘        │
        │                             │
funannotate2 annotate                 │
(Pfam, dbCAN, MEROPS, Swiss-Prot, BUSCO, GO + the above)
        │                             │
        │                             │
        ├──────── final_summary.tsv ──┘
        │
NCBI submission (table2asn, validated)   only for genomes with ncbi_submission set
```

## Directory layout

```
Snakefile
Snakemake_env.yml      # Snakemake 9.26.1 + SLURM plugin + Nextflow
slurm_submit.sbatch    # runs the whole workflow as one small SLURM job
config/
  config.yaml          # edit this
profile/
  config.yaml          # BioCloud SLURM executor profile
workflow/
  rules/               # 0_ ... 10_, one stage per file
  envs/                # tool+version reference (see "Containers")
  scripts/             # fix_annotations.py, summarize_results.py, NCBI scripts
data/
  local_assemblies/    # drop your own *.fna/*.fa/*.fasta here (optional)
  genomes/             # staged, unified input (created by the pipeline)
  databases/           # antiSMASH / funannotate2 / InterProScan6 / OMArk DBs
resources/             # geneML, funannotate2-addons, table2asn (created by the pipeline)
results/
  <genome>/
    01_preprocess/     # cleaned + soft-masked assembly, <genome>.seqinfo.tsv
    02_geneml/         # geneML output + the final gene models (*.models.*)
    03_antismash/
    04_interproscan/
    05_annotations/    # antiSMASH + InterProScan6 as funannotate2 tables
    06_annotate/       # <genome>.gff3/.gbk/.tbl/.proteins.fa/... (final)
    07_omark/          # <genome>.sum + OMArk plots/lists
    08_ncbi/           # <genome>.sqn to upload, .gbf, .validation.txt (NCBI genomes only)
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
- `genomes:` - per-genome `species`, `strain`, and **optional** `locus_tag`
  and `busco_lineage`. Every genome needs an entry here (at least
  `species`); Snakemake stops at start-up and names the genome if one is
  missing. Set `locus_tag` to have `gfftk rename` give the gene models
  NCBI-style IDs (`LOCUSTAG_000001`, transcripts `LOCUSTAG_000001-T1`);
  leave it blank to keep geneML's native IDs. `busco_lineage` overrides
  `funannotate2.busco_lineage` (default `fungi`) for that genome, e.g.
  `sordariomycetes`, `eurotiomycetes`, `basidiomycota`.
- `genomes:` NCBI keys - `ncbi_submission`, `bioproject`, `biosample`,
  `keep_contigs`, `organelles`, `mito_gcode`; see
  [NCBI submission](#ncbi-submission).
- `ncbi:` - the submission template and table2asn settings, see
  [NCBI submission](#ncbi-submission).
- `geneml.version` - geneML release installed from PyPI (default `1.2.0`).
- `funannotate2.addons_version` - funannotate2-addons release installed from
  PyPI (default `26.3.7`).
- `antismash.database_dir` / `interproscan6.datadir` / `funannotate2.db_dir`
  / `omark.db_dir` - where the (large, shared) reference databases live.
  Downloaded once by rules `antismash_database` / `funannotate2_database` /
  `omark_database`; InterProScan6 fetches its own data into `datadir` on
  first run.
- `omark.db_url` / `omark.db_md5` - the OMAmer database (`LUCA.h5`, ~10 GB),
  pinned to a Zenodo release (OMA May 2026) and checked against its md5.
- `interproscan6.interpro_version` - the InterPro data release (default
  `110.0`, the newest compatible with InterProScan 6.0.x). It is pinned so
  results are reproducible, and so a new InterPro release doesn't make every
  genome job download it at once. Changing it re-runs
  `interproscan6_setup`, which downloads the new release once.
- `interproscan6.goterms` / `.pathways` - GO terms and pathways, which
  InterProScan6 leaves off by default; on here, since funannotate2 takes its
  GO annotation from them.
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
| `funannotate2_*`, `softmask`, `gene_models` | `quay.io/biocontainers/funannotate2:26.6.21--pyhdfd78af_0` (includes gfftk) |
| `antismash*` | `quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_1` |
| `omark*` | `quay.io/biocontainers/omark:0.5.0--pyhdfd78af_0` (OMAmer 2.1.2) |
| `get_geneml`, `geneml` | `python:3.12-slim`; geneML is only on PyPI, so `get_geneml` installs it once into `resources/geneml-<version>/` |
| `get_funannotate2_addons`, `external_annotations` | `python:3.12-slim`, same approach for funannotate2-addons (`resources/funannotate2-addons-<version>/`) |
| `sequence_info`, `get_table2asn`, `ncbi_submission` | `python:3.12-slim`; `get_table2asn` downloads NCBI's table2asn binary once into `resources/table2asn-<version>/` and checks its sha256 |
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
as a reference only: the profile uses Apptainer, and geneML and
funannotate2-addons (PyPI) and InterProScan6 (Nextflow) have no conda env.

The funannotate2 biocontainer is used rather than the `nextgenusfs/funannotate2`
Docker image, which bakes a read-only `FUNANNOTATE2_DB` into the image;
funannotate2 needs to write BUSCO lineages into that directory.

Databases and outputs are inside the working directory, which Snakemake binds
into every container. Anything you point at outside it -- `local_genome_dir`,
`funannotate2.db_dir`, `antismash.database_dir`, `omark.db_dir` -- must be bound too, e.g.
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

## NCBI submission

For genomes with `ncbi_submission` set, rule `ncbi_submission` builds
`results/<genome>/08_ncbi/<genome>.sqn`, the file you upload to GenBank. It
takes funannotate2 annotate's `.tbl` and the unmasked assembly, and runs
NCBI's [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/) in
genome mode (`-M n`). This gives the validator messages and the
discrepancy report. The rule **fails** if either has a problem that NCBI
requires you to fix: an ERROR or REJECT validator message, or a FATAL
discrepancy category. `08_ncbi/<genome>.validation.txt` summarises them,
and table2asn's own `.val`/`.dr`/`.stats` files stay in
`08_ncbi/table2asn/out/`. Set `ncbi.allow_errors: true` to get the `.sqn`
anyway, for example to look at it before fixing things.

**One-time setup.**

1. Make a submission template (`.sbt`) at
   <https://submit.ncbi.nlm.nih.gov/genbank/template/submission/> with
   your contact details and author list. Save it as `config/template.sbt`
   (`ncbi.sbt_template`). `config/*.sbt` is git-ignored because it holds
   personal details.
2. Every genome needs a BioProject, a BioSample and a registered
   `locus_tag` prefix. NCBI assigns the prefix when you register the
   BioSample under the BioProject.

**Per genome** (under `genomes:`):

- `ncbi_submission`:
  - `"update"` adds this annotation to an existing, unannotated GenBank WGS
    assembly. The genome must come from `ncbi_accessions`. The sequences
    are named as NCBI's
    [WGS update instructions](https://www.ncbi.nlm.nih.gov/genbank/wgs_update/)
    ask: `gnl|WGS:<prefix>|<SeqID>|gb|<contig accession>`. Protein and
    transcript IDs become `gnl|WGS:<prefix>|...`. Rule `sequence_info`
    looks up these names from NCBI (Datasets API + E-utilities), so it
    needs internet access. The assembly is not reordered or renamed, so
    the coordinates match the GenBank records.
  - `"new"` is for a genome that is not in GenBank yet. The `.sqn` is part
    of a new WGS submission.
- `bioproject`, `biosample`, `locus_tag` are required for NCBI genomes.
  Snakemake checks this at start-up.
- `keep_contigs: true` keeps the assembly's own sequence order and names.
  This is implied by `"update"`. Without it, funannotate2 clean sorts the
  contigs, drops duplicates and filters short ones.
- `organelles: {contig_id: mitochondrion}` lists organelle sequences for
  `"new"` genomes; `"update"` takes them from NCBI. Genes on organelle
  sequences are dropped, because geneML predicts with the nuclear genetic
  code. Those sequences are submitted unannotated, with
  `[location=mitochondrion]` and `mgcode` = `mito_gcode` (default 4, the
  mould/protozoan mitochondrial code; use 3 for yeasts).

**Under `ncbi:`:**

- `structured_comment` is an optional table2asn `-w` file, for example a
  Genome-Assembly-Data comment for a `"new"` genome.
- `gaps_min` / `linkage_evidence`: runs of at least this many Ns become
  `assembly_gap` features with this evidence.
- `extra` passes any other table2asn flags.

table2asn runs with `-T`, so it looks the organism up in NCBI Taxonomy.
This also needs internet access. Without it, table2asn applies bacterial
checks.

**Uploading.** Read the `.validation.txt` and the warnings in the `.val`
file. Then also check
`06_annotate/<genome>.need-curating.txt`, the gene names/products that
funannotate2 flagged. Then submit the `.sqn`:

- For `"update"`: the update route is through the
  [Genome Submission Portal](https://submit.ncbi.nlm.nih.gov/subs/genome/)
  ("update an existing submission"), or by emailing the `.sqn` to
  genomes@ncbi.nlm.nih.gov with the assembly accession. NCBI's update page
  covers changing the annotation of genomes that are already annotated.
  Adding annotation for the first time is not described there, so ask
  genomes@ncbi.nlm.nih.gov which route they want before uploading.
- For `"new"`: create a new genome submission in the portal and give the
  `.sqn` as the annotated sequence file.

NCBI also screens new genomes for contamination with
[FCS-GX](https://github.com/ncbi/fcs). It is not part of the pipeline
because it needs about 500 GB of RAM. For a `"new"` genome, consider
running it first.

## Notes / caveats

- **Network access from compute nodes.** Besides the database rules,
  `funannotate2 annotate` looks the species up in an online taxonomy
  service on every run (failures are logged and ignored), and InterProScan6
  can look up precalculated matches online. The database rules
  (`antismash_database`, `funannotate2_database`, `omark_database`,
  `interproscan6_setup`, `get_table2asn`) need internet access, and so do
  `sequence_info` for `"update"` genomes and `ncbi_submission`
  (taxonomy lookup).
- **funannotate2-addons fixes.** `workflow/scripts/fix_annotations.py`
  post-processes the funannotate2-addons (f2a 26.3.7) tables before
  `funannotate2 annotate` reads them: antiSMASH results are keyed by gene ID
  and are mapped to transcript IDs; f2a's `go_term` key is renamed to the
  `go_terms` key annotate actually reads (otherwise every InterProScan6 GO
  term is silently dropped); InterProScan6's `GO:...(InterPro)` source
  suffixes are stripped. The InterProScan6 TSV is parsed rather than the
  XML, because the XML groups identical proteins and f2a only reads the
  first ID of each group.
- **06_annotate/ vs 08_ncbi/.** funannotate2 annotate runs table2asn
  itself, with a placeholder submitter and without genome mode, so the
  `.gbk` in `06_annotate/` is for reading, not for submission. Submit the
  `.sqn` from `08_ncbi/`.
- **Soft-masking** is tantan (funannotate2's own `softmask_fasta`, the same
  default as `funannotate mask`), which only masks low-complexity/tandem
  repeats. For repeat-rich genomes, consider masking with
  RepeatModeler + RepeatMasker before the pipeline.
- **geneML gene-ID prefixes** are derived automatically from the `{genome}`
  wildcard (non-alphanumeric characters replaced with `_`), so NCBI
  accessions like `GCA_052058355.1` become a valid `--gene-id-prefix`.
- **OMArk** runs on the gene models (one protein per gene, since geneML
  runs with `--max-transcripts 1`). The `S:` line in `07_omark/<genome>.sum`
  is completeness against conserved HOGs; the `A:` line shows how much of
  the whole proteome is taxonomically consistent vs. inconsistent,
  contaminant or unknown -- a high unknown/fragmented share flags spurious
  gene models.
- The antiSMASH, funannotate2 and OMArk database rules download large,
  shared reference data (several GB+ each) - run those once and reuse the
  database directories across projects rather than re-downloading.

## Automatic checks

`.github/workflows/dry-run.yml` byte-compiles `workflow/scripts/*.py` and
runs `snakemake -n` with Snakemake 9.26.1 on every push and pull request.
