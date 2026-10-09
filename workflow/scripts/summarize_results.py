"""
Aggregate a one-line-per-genome summary: predicted gene count (geneML),
antiSMASH region count, whether InterProScan6 results are present, whether
a locus_tag rename was applied, and the funannotate output folder.
"""
import csv
import os

genomes = snakemake.params.genomes  # list of genome names
config_genomes = snakemake.params.genome_config  # config["genomes"] dict

rows = []
for genome in genomes:
    geneml_gff3 = f"results/{genome}/02_geneml/{genome}.geneml.gff3"
    antismash_gbk = f"results/{genome}/03_antismash/{genome}.gbk"
    ips_xml = f"results/{genome}/04_interproscan/{genome}.xml"
    annotate_dir = f"results/{genome}/06_annotate"

    n_genes = 0
    if os.path.exists(geneml_gff3):
        with open(geneml_gff3) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3 and fields[2] == "gene":
                    n_genes += 1

    n_regions = 0
    if os.path.exists(antismash_gbk):
        with open(antismash_gbk) as fh:
            for line in fh:
                if line.startswith("FEATURES"):
                    continue
                if "/region_number=" in line:
                    n_regions += 1

    meta = config_genomes.get(genome, {})
    rows.append({
        "genome": genome,
        "species": meta.get("species", ""),
        "strain": meta.get("strain", ""),
        "predicted_genes_geneml": n_genes,
        "antismash_regions": n_regions,
        "interproscan6_run": os.path.exists(ips_xml),
        "locus_tag_applied": bool(meta.get("locus_tag", "")),
        "annotate_out_dir": annotate_dir,
    })

with open(snakemake.output.tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=list(rows[0].keys()) if rows else
                             ["genome", "species", "strain",
                              "predicted_genes_geneml", "antismash_regions",
                              "interproscan6_run", "locus_tag_applied",
                              "annotate_out_dir"],
                             delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
