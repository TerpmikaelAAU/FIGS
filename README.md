# FIGS: fungal genome annotation

**F**unannotate2 + **I**nterProScan6 + **G**eneML + anti**S**MASH, with
[OMArk](https://github.com/DessimozLab/OMArk) quality control, as a
Snakemake 9 pipeline (SLURM executor plugin, one Apptainer container per
rule). The annotation counterpart to
[BAGS](https://github.com/TerpmikaelAAU/BAGS).

```
assembly (local FASTA or NCBI accession)
  → funannotate2 clean → soft-mask (tantan)
  → geneML gene models (optional locus_tag rename)
  → antiSMASH, InterProScan6, eggNOG-mapper, OMArk
  → + tRNA/rRNA genes (tRNAscan-SE, Rfam)
  → funannotate2 annotate (Pfam, dbCAN, MEROPS, Swiss-Prot, BUSCO + the above)
  → annotations.tsv, final_summary.tsv
  → NCBI .sqn (table2asn, validated)   only for genomes with ncbi_submission
```

## Install

```bash
conda env create -f Snakemake_env.yml     # Snakemake 9.26.1 + SLURM plugin
conda activate snakemake_figs
```

Apptainer must be on the compute nodes. Every tool comes from a container
or is installed by the pipeline; the first run downloads the databases.

## Configure

`config/config.yaml`, one entry per genome:

```yaml
genomes:
  my_strain:                          # results/my_strain/
    fasta: "/path/to/my_strain.fna"
    species: "Aspergillus niger"
    strain: "ABC123"

  GCA_052058355.1:                    # no fasta: downloaded from NCBI
    species: "Trichoderma asperellum"
    strain: "TA1"
```

| Key | |
|---|---|
| `fasta` | The assembly. Leave it out to download it from NCBI; the name must then be the accession (`GCA_`/`GCF_`). |
| `species` | Required. Species and strain are not fetched from NCBI. |
| `strain` | Optional. |
| `locus_tag` | Optional; required for NCBI. Renames genes to `LOCUSTAG_000001`. Must be the prefix registered to the BioSample. |
| `busco_lineage` | Optional, default `fungi`, e.g. `eurotiomycetes`, `sordariomycetes`. |
| `ncbi_submission`, `bioproject`, `biosample`, `keep_contigs`, `organelles`, `mito_gcode` | See [NCBI submission](#ncbi-submission). |

Tool settings are in `workflow/defaults.yaml`; repeat one under its
section in `config/config.yaml` to change it.

## Run

```bash
snakemake --profile profile -n     # dry run
sbatch slurm_submit.sbatch         # or: snakemake --profile profile, in tmux
snakemake --profile profile --unlock   # after a crash left the directory locked
```

Per-rule resources are in the `resources` dict at the top of the
`Snakefile` (`runtime` in minutes).

Results go to `results/<genome>/`. The final annotation is in
`06_annotate/` (`.gff3`, `.gbk`, `.tbl`, proteins), with
`<genome>.annotations.tsv` (one row per protein, tRNA and rRNA, a column per
annotation source). OMArk is in `07_omark/` and the NCBI files in
`08_ncbi/`. `<genome>.final_summary.tsv` has the genome's counts, and
`results/final_summary.tsv` has one row per genome.

## NCBI submission

For genomes with `ncbi_submission`, rule `ncbi_submission` runs table2asn in
genome mode on the annotation and writes `08_ncbi/<genome>.sqn`. It fails
on validator ERROR/REJECT messages or FATAL discrepancies, summarised in
`08_ncbi/<genome>.validation.txt`. `ncbi.allow_errors: true` writes the
`.sqn` anyway.

Setup:
- A submission template from
  <https://submit.ncbi.nlm.nih.gov/genbank/template/submission/>, saved as
  `config/template.sbt`. It is git-ignored, since it holds personal details.
- Per genome: `bioproject`, `biosample` and a registered `locus_tag`.
  These are checked at start-up.
- `ncbi_submission: "update"` adds annotation to an existing GenBank WGS
  assembly. The genome must be downloaded by accession. Contigs keep their
  order and are named the way NCBI's WGS update asks
  (`gnl|WGS:<prefix>|<SeqID>|gb|<accession>`), looked up online.
- `ncbi_submission: "new"` is for a genome not yet in GenBank.
  `keep_contigs: true` keeps your contig order and names. `organelles:
  {contig: mitochondrion}` marks organelle contigs (needs `keep_contigs`);
  `mito_gcode` defaults to 4 (moulds); use 3 for yeasts.

Before uploading:
- Read the `.val` warnings and `06_annotate/<genome>.need-curating.txt`.
- Submit the `.sqn` from `08_ncbi/`, not the `.gbk` in `06_annotate/`.
  funannotate2 builds that `.gbk` with a placeholder submitter and outside
  genome mode.
- For `"update"`: NCBI documents changing existing annotation, not adding
  it for the first time. Ask genomes@ncbi.nlm.nih.gov which route they want.
- NCBI screens new genomes with [FCS-GX](https://github.com/ncbi/fcs).
  FIGS doesn't run it, because it needs ~500 GB RAM.

## Caveats and lessons learned

- **Network access** is needed by the database rules, `download_genome`,
  `sequence_info` (for `"update"`), funannotate2 annotate (taxonomy lookup,
  failures ignored) and `ncbi_submission`. table2asn runs with `-T`; without
  internet it applies bacterial checks.
- **Databases** are large (eggNOG ~12 GB download, OMArk `LUCA.h5` ~10 GB,
  slow from Zenodo). Downloads resume after a dropped connection
  (`<file>.part`) and are checksum-verified. A database directory outside
  the working directory must be bound into the containers, e.g.
  `--apptainer-args "--bind /shared/databases"`.
- **Containers.**
  - antiSMASH uses the official `standalone-lite` image. The biocontainer
    lacks the pre-built data and can't write it into its read-only install.
  - funannotate2 uses the biocontainer, not the `nextgenusfs` Docker image,
    which has a read-only `FUNANNOTATE2_DB` that BUSCO lineages can't be
    written into.
  - InterProScan6 has no container: it is Nextflow starting its own
    containers, and wrapping it would mean Apptainer-in-Apptainer.
    `interproscan6_setup` runs once first, so parallel genome jobs don't all
    download its data.
- **funannotate2-addons (f2a 26.3.7) needs fixing** (`fix_annotations.py`):
  - It writes GO terms as `go_term`, which annotate silently ignores.
  - antiSMASH results are keyed by gene, not transcript.
  - It only reads the first protein of each identical-protein group in the
    InterProScan XML, so the TSV is used instead.
  - Reactome pathways are dropped. InterPro links each entry to every
    species' pathways, which for one genome came to 5.1 million notes,
    nearly all from animals.
  - InterPro entries also become `db_xref`s, because gfftk's
    `summary.json` only counts InterPro from `db_xref`s.
- **eggNOG** gives most gene names and products. Its table is converted by
  `emapper_annotations.py`, because via f2a its descriptions never become
  products. EC numbers end up as `EC:` notes (funannotate2 26.6.21 has no
  EC qualifier). eggNOG GO/KEGG terms are not added; they are in
  `04_eggnog/`.
- **tRNA/rRNA genes** (geneML predicts protein-coding genes only):
  - tRNAscan-SE pseudogenes and `Undet`/`Sup` tRNAs are dropped, and
    products are written as `tRNA-Met`, not `iMet`. NCBI rejects those.
  - rRNAs must cover 80% of the Rfam model and not be cut by a contig end,
    since gfftk can't mark them partial.
  - An RNA gene overlapping a protein-coding CDS is dropped, not the
    protein; on a test genome those proteins were real. The dropped genes
    are listed in `logs/<genome>/ncrna_models.log`.
  - RNA genes are numbered after the last protein-coding gene, so adding
    them doesn't change protein IDs or rerun the protein searches.
- **Soft-masking** is tantan only (low-complexity/tandem repeats). For
  repeat-rich genomes, mask with RepeatModeler + RepeatMasker first.
- **Organelle contigs** are left unannotated, because geneML predicts with
  the nuclear genetic code.

## License

[MIT](LICENSE). The tools and databases FIGS downloads keep their own
licenses.
