#!/usr/bin/env python3
"""
Make a funannotate2-addons annotation table (3 columns: id, type, value)
usable by `funannotate2 annotate -a`:

  * IDs: funannotate2 matches annotations to transcript (mRNA) IDs. antiSMASH,
    run with --genefinding-gff3, names its CDS features after the gene
    instead, so gene IDs (and gene/transcript Aliases) are mapped to the
    gene's transcript IDs using the gene models GFF3. Lines whose ID is not
    in the GFF3 are dropped and counted.
  * GO terms: funannotate2-addons 26.3.7 writes `go_term`, but annotate only
    picks up `go_terms`, so the key is renamed. InterProScan6 also tags each
    GO term with its source, e.g. GO:0005515(InterPro); only the GO ID is kept.
  * Placeholder values ("-", or a trailing " -" where InterProScan6 had no
    description) are removed, and duplicate lines are dropped.

Usage: fix_annotations.py --gff3 models.gff3 -i in.annotations.txt -o out.txt
"""
import argparse
import re
import sys

TRANSCRIPT_TYPES = {"mRNA", "tRNA", "rRNA", "ncRNA", "transcript"}
GO_RE = re.compile(r"^(GO:\d{7})")


def parse_attributes(col):
    attrs = {}
    for part in col.strip().strip(";").split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key.strip()] = value.strip().split(",")
    return attrs


def id_map(gff3):
    """Return {any known ID or Alias: [transcript IDs]}."""
    gene_aliases = {}   # gene ID -> its aliases
    transcripts = {}    # gene ID -> [transcript IDs]
    lookup = {}
    with open(gff3) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            attrs = parse_attributes(cols[8])
            fid = attrs.get("ID", [None])[0]
            if fid is None:
                continue
            if cols[2] == "gene":
                gene_aliases[fid] = attrs.get("Alias", [])
                transcripts.setdefault(fid, [])
            elif cols[2] in TRANSCRIPT_TYPES:
                lookup[fid] = [fid]
                for alias in attrs.get("Alias", []):
                    lookup.setdefault(alias, []).append(fid)
                for parent in attrs.get("Parent", []):
                    transcripts.setdefault(parent, []).append(fid)
    for gene, tx in transcripts.items():
        if not tx:
            continue
        for key in [gene] + gene_aliases.get(gene, []):
            lookup.setdefault(key, [])
            for t in tx:
                if t not in lookup[key]:
                    lookup[key].append(t)
    return lookup


def resolve(fid, lookup):
    if fid in lookup:
        return lookup[fid]
    # antiSMASH names a gene's CDSs <gene>_0, <gene>_1, ... if it has several.
    trimmed = re.sub(r"_\d+$", "", fid)
    return lookup.get(trimmed, [])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--gff3", required=True, help="gene models GFF3")
    ap.add_argument("-i", "--input", required=True, help="f2a annotations table")
    ap.add_argument("-o", "--out", required=True, help="fixed annotations table")
    args = ap.parse_args()

    lookup = id_map(args.gff3)
    seen = set()
    rows = []
    n_in = n_unmapped = n_malformed = n_dropped = 0
    unmapped_examples = []
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            n_in += 1
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 3:
                n_malformed += 1
                continue
            fid, key, value = (c.strip() for c in cols)
            if key == "go_term":
                key = "go_terms"
            if key == "go_terms":
                m = GO_RE.match(value)
                value = m.group(1) if m else ""
            value = re.sub(r"\s+-$", "", value).strip()
            if not key or value in ("", "-"):
                n_dropped += 1
                continue
            targets = resolve(fid, lookup)
            if not targets:
                n_unmapped += 1
                if len(unmapped_examples) < 5:
                    unmapped_examples.append(fid)
                continue
            for t in targets:
                row = (t, key, value)
                if row not in seen:
                    seen.add(row)
                    rows.append(row)

    with open(args.out, "w") as out:
        out.write("#gene_id\tannotation_type\tannotation_value\n")
        for row in rows:
            out.write("\t".join(row) + "\n")

    print(f"{args.input}: {n_in} lines in, {len(rows)} written, "
          f"{n_dropped} empty values dropped, {n_malformed} malformed, "
          f"{n_unmapped} with an ID not in {args.gff3}", file=sys.stderr)
    if unmapped_examples:
        print(f"  unmapped IDs, e.g.: {', '.join(unmapped_examples)}", file=sys.stderr)
    if n_in and not rows:
        sys.exit(f"ERROR: none of the annotations in {args.input} matched a "
                 f"gene or transcript in {args.gff3}")


if __name__ == "__main__":
    main()
