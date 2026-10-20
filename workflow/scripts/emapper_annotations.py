#!/usr/bin/env python3
"""
Turn eggNOG-mapper's .emapper.annotations into a funannotate2 annotation
table (3 columns: transcript ID, type, value) for `funannotate2 annotate -a`.

funannotate2-addons' `f2a emapper` parser is not used, because funannotate2
26.6.21 would ignore most of what it writes: it puts the eggNOG description
in a `note`, so it never becomes a product, and its `EC_number` key is not
read by annotate. Here, per protein:

  * name + product: eggNOG's Preferred_name and Description (first
    sentence), as funannotate 1 did, only when the name looks like a gene
    name. funannotate2's name cleaner then swaps in its curated product for
    known names and fixes or replaces descriptions NCBI would reject.
    "Belongs to the ... family" descriptions are not used as products.
  * note `EggNog:<OG> <description>`: the orthologous group at eggNOG's
    annotation level, as InterPro notes are written.
  * note `COG:<categories>` and one note `EC:<number>` per EC number
    (funannotate2 26.6.21 cannot carry EC numbers as their own qualifier).

Usage: emapper_annotations.py -i genome.emapper.annotations -o out.txt
"""
import argparse
import re
import sys

EC_RE = re.compile(r"^\d+\.(\d+|-)\.(\d+|-)\.(\d+|-)$")


def gene_name_ok(name):
    # The same test funannotate2-addons applies to eggNOG names.
    return (len(name) > 2 and not name.startswith("_") and "_" not in name
            and "." not in name and re.search("[A-Za-z]", name) is not None)


def annotation_og(ogs, level):
    """The OG at eggNOG's annotation level (max_annot_lvl), e.g.
    3NU6E@4751|Fungi -> ENOG503NU6E; COG/KOG IDs are kept as they are."""
    for og in ogs.split(","):
        acc, _, og_level = og.partition("@")
        if og_level == level and acc:
            return acc if acc.startswith(("COG", "KOG")) else f"ENOG50{acc}"
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-i", "--input", required=True, help=".emapper.annotations")
    ap.add_argument("-o", "--out", required=True, help="annotation table")
    args = ap.parse_args()

    header = None
    rows = []
    n_hits = n_named = 0
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("##") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                header = [c.lstrip("#") for c in cols]
                continue
            if header is None:
                sys.exit(f"ERROR: no #query header line in {args.input}")
            rec = {k: ("" if v == "-" else v.strip()) for k, v in zip(header, cols)}
            query = rec.get("query", "")
            if not query:
                continue
            n_hits += 1

            name = rec.get("Preferred_name", "")
            desc = rec.get("Description", "").split(". ")[0].rstrip(".").strip()
            if gene_name_ok(name):
                rows.append((query, "name", name))
                if desc and not desc.startswith("Belongs to"):
                    rows.append((query, "product", desc))
                n_named += 1

            og = annotation_og(rec.get("eggNOG_OGs", ""), rec.get("max_annot_lvl", ""))
            if og or desc:
                rows.append((query, "note", f"EggNog:{og} {desc}".strip()))
            cog = rec.get("COG_category", "").replace(" ", "")
            if cog:
                rows.append((query, "note", f"COG:{','.join(cog)}"))
            for ec in rec.get("EC", "").split(","):
                if EC_RE.match(ec.strip()):
                    rows.append((query, "note", f"EC:{ec.strip()}"))

    with open(args.out, "w") as out:
        out.write("#gene_id\tannotation_type\tannotation_value\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    print(f"{args.input}: {n_hits} proteins with an eggNOG hit, {n_named} with a "
          f"gene name, {len(rows)} lines written", file=sys.stderr)


if __name__ == "__main__":
    main()
