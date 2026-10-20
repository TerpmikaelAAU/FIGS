#!/usr/bin/env python3
"""
One row per protein with its functional annotation, from funannotate2
annotate's .tbl: gene and transcript ID, location, gene name, product, and
one column per source (InterPro, Pfam, GO, eggNOG, CAZymes, MEROPS,
antiSMASH, ...), followed by one row per tRNA and rRNA gene (type column).
The counterpart of funannotate 1's annotations.txt.

The .tbl is read rather than the GFF3, because funannotate2 joins a gene's
notes with unescaped commas in the GFF3, which splits descriptions that
contain a comma; the .tbl has one qualifier per line. Several values in a
column are separated by "; ".

Usage: annotation_table.py -i genome.tbl -o genome.annotations.tsv
"""
import argparse
import csv
import re
import sys

COLUMNS = ["gene_id", "transcript_id", "type", "contig", "start", "end", "strand",
           "name", "product", "ec_number", "interpro", "pfam", "go_terms",
           "eggnog", "cog", "cazyme", "protease", "antismash", "busco",
           "pathways", "swissprot", "other_notes"]

# note prefix -> column (the prefix itself is left out of the column)
NOTE_PREFIXES = [
    ("InterPro:", "interpro"),
    ("EggNog:", "eggnog"),
    ("COG:", "cog"),
    ("EC:", "ec_number"),
    ("CAZy:", "cazyme"),
    ("MEROPS:", "protease"),
    ("antiSMASH:", "antismash"),
    ("BUSCO:", "busco"),
    ("Pathway:", "pathways"),
]
SMCOG_RE = re.compile(r"^SMCOG\d+")
FEATURE_TYPES = ("CDS", "tRNA", "rRNA")
GO_ASPECTS = {"go_function", "go_process", "go_component"}


def coord(text):
    return int(re.sub(r"\D", "", text))  # drop partial markers (< >)


def parse_tbl(path):
    """Yield one dict per CDS, tRNA and rRNA: type, contig, intervals, gene
    qualifiers, and its own qualifiers as a list of (key, value)."""
    contig = None
    gene = {}
    feature = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Feature"):
                if feature and feature["type"] in FEATURE_TYPES:
                    yield feature
                contig = line.split(None, 1)[1].strip()
                feature = None
                continue
            cols = line.split("\t")
            if line.startswith("\t\t\t"):
                if feature is None or len(cols) < 4:
                    continue
                key, value = cols[3], "\t".join(cols[4:]).strip()
                if feature["type"] == "gene":
                    gene.setdefault(key, []).append(value)
                else:
                    feature["quals"].append((key, value))
            elif len(cols) >= 3 and cols[2]:
                if feature and feature["type"] in FEATURE_TYPES:
                    yield feature
                if cols[2] == "gene":
                    gene = {}
                feature = {"type": cols[2], "contig": contig, "gene": gene,
                           "intervals": [(coord(cols[0]), coord(cols[1]))],
                           "quals": []}
            elif len(cols) == 2 and feature is not None:
                feature["intervals"].append((coord(cols[0]), coord(cols[1])))
    if feature and feature["type"] in FEATURE_TYPES:
        yield feature


def row_for(cds):
    cells = {c: [] for c in COLUMNS}
    for key, value in cds["quals"]:
        if key == "product":
            cells["product"].append(value)
        elif key == "protein_id":
            cells["transcript_id"].append(value.rsplit("|", 1)[-1])
        elif key == "EC_number":
            cells["ec_number"].append(value)
        elif key in GO_ASPECTS:
            desc, _, rest = value.partition("|")
            go_id = rest.split("|", 1)[0]
            cells["go_terms"].append(f"GO:{go_id} {desc}" if go_id else desc)
        elif key == "db_xref":
            if value.startswith("PFAM:"):
                cells["pfam"].append(value[len("PFAM:"):])
            elif value.startswith("UniProtKB/Swiss-Prot:"):
                cells["swissprot"].append(value[len("UniProtKB/Swiss-Prot:"):])
            elif value.startswith("InterPro:"):
                pass  # the same entries as the InterPro notes, which have descriptions
            else:
                cells["other_notes"].append(value)
        elif key == "note":
            for prefix, column in NOTE_PREFIXES:
                if value.startswith(prefix):
                    cells[column].append(value[len(prefix):].strip())
                    break
            else:
                if SMCOG_RE.match(value):
                    cells["antismash"].append(value)
                elif "% identical to " in value:
                    cells["swissprot"].append(value)
                else:
                    cells["other_notes"].append(value)

    starts = [min(a, b) for a, b in cds["intervals"]]
    ends = [max(a, b) for a, b in cds["intervals"]]
    first = cds["intervals"][0]
    gene = cds["gene"]
    cells["gene_id"] = gene.get("locus_tag", [])[:1]
    cells["name"] = gene.get("gene", [])[:1]
    cells["type"] = [cds["type"]]
    cells["contig"] = [cds["contig"]]
    cells["start"] = [str(min(starts))]
    cells["end"] = [str(max(ends))]
    cells["strand"] = ["-" if first[0] > first[1] else "+"]
    # De-duplicate, keeping the order they came in.
    return {c: "; ".join(dict.fromkeys(v)) for c, v in cells.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-i", "--tbl", required=True, help="funannotate2 annotate .tbl")
    ap.add_argument("-o", "--out", required=True, help="output TSV")
    args = ap.parse_args()

    with open(args.out, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        # proteins first, then the RNA genes
        rows = sorted((row_for(f) for f in parse_tbl(args.tbl)),
                      key=lambda r: FEATURE_TYPES.index(r["type"]))
        writer.writerows(rows)
    n = {t: sum(r["type"] == t for r in rows) for t in FEATURE_TYPES}
    print(f"{args.tbl}: {n['CDS']} proteins, {n['tRNA']} tRNAs and {n['rRNA']} "
          f"rRNAs written to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
