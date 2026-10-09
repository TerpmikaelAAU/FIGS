"""
Aggregate a one-line-per-genome summary: gene counts (geneML, final
funannotate2 annotation), antiSMASH region count, proteins with an
InterProScan6 hit, OMArk completeness and consistency, the BUSCO lineage
funannotate2 used, and whether locus_tags were applied.
"""
import csv
import os

genomes = snakemake.params.genomes  # list of genome names
config_genomes = snakemake.params.genome_config  # config["genomes"] dict
busco_lineages = snakemake.params.busco_lineages  # genome -> lineage


def count_genes(gff3):
    n = 0
    if os.path.exists(gff3):
        with open(gff3) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3 and fields[2] == "gene":
                    n += 1
    return n


def omark_scores(sum_file):
    """Return the completeness (S/D/M) and consistency (A/I/C/U) lines with
    percentages from an OMArk .sum file, e.g. 'S:95.10%,D:1.20%[...],M:3.70%'."""
    completeness = consistency = ""
    clade = ""
    if os.path.exists(sum_file):
        with open(sum_file) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#The selected clade was"):
                    clade = line.split("was", 1)[1].strip()
                elif line.startswith("S:") and "%" in line:
                    completeness = line
                elif line.startswith("A:") and "%" in line:
                    consistency = line
    return clade, completeness, consistency


rows = []
for genome in genomes:
    antismash_gbk = f"results/{genome}/03_antismash/{genome}.gbk"
    ips_tsv = f"results/{genome}/04_interproscan/{genome}.tsv"

    n_regions = 0
    if os.path.exists(antismash_gbk):
        with open(antismash_gbk) as fh:
            for line in fh:
                if "/region_number=" in line:
                    n_regions += 1

    ips_proteins = set()
    if os.path.exists(ips_tsv):
        with open(ips_tsv) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    ips_proteins.add(line.split("\t", 1)[0])

    clade, completeness, consistency = omark_scores(
        f"results/{genome}/07_omark/{genome}.sum")

    meta = config_genomes.get(genome, {})
    rows.append({
        "genome": genome,
        "species": meta.get("species", ""),
        "strain": meta.get("strain", "") or "",
        "genes_geneml": count_genes(f"results/{genome}/02_geneml/{genome}.geneml.gff3"),
        "genes_final": count_genes(f"results/{genome}/06_annotate/{genome}.gff3"),
        "antismash_regions": n_regions,
        "proteins_with_interpro_hit": len(ips_proteins),
        "omark_clade": clade,
        "omark_completeness": completeness,
        "omark_consistency": consistency,
        "busco_lineage": busco_lineages.get(genome, ""),
        "locus_tag": (meta.get("locus_tag") or "").strip(),
        "annotation": f"results/{genome}/06_annotate/{genome}.gbk",
    })

fieldnames = ["genome", "species", "strain", "genes_geneml", "genes_final",
              "antismash_regions", "proteins_with_interpro_hit", "omark_clade",
              "omark_completeness", "omark_consistency", "busco_lineage",
              "locus_tag", "annotation"]
with open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
