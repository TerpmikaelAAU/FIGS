#!/usr/bin/env python3
"""
Add the tRNA and rRNA genes (ncrna_gff3.py) to the protein-coding gene
models, for funannotate2 annotate:

  * RNA genes on organelle sequences (a location in sequence_info) are
    dropped, as gene_models drops geneML's genes there.
  * An RNA gene overlapping one already taken is dropped (files are taken
    in the order given, so list rRNAs first).
  * An RNA gene overlapping a protein-coding CDS is dropped, not the
    protein: NCBI's validator flags CDS/RNA overlaps, and in a test genome
    the proteins involved were real (a glutamyl-tRNA synthetase, a METTL16
    methyltransferase), each overlapping a single 5S rRNA or tRNA hit.
  * RNA genes get final IDs. With a locus_tag, they are numbered on from the
    highest protein-coding gene (<locus_tag>_000001 style, as gfftk rename
    numbers), so no protein-coding gene changes ID and nothing searched on
    the proteins has to be redone. Without one: <prefix>_tRNA00001, ...

The output is not sorted; gfftk sanitize sorts it.

Usage: merge_ncrna.py --models models.gff3 --rna rRNA.gff3 tRNA.gff3
           --seqinfo seqinfo.tsv --prefix P [--locus-tag T] -o out.gff3
"""
import argparse
import re
import sys
from collections import Counter


def parse_attrs(col):
    return dict(p.split("=", 1) for p in col.strip().strip(";").split(";") if "=" in p)


def read_genes(path):
    """Return [(gene ID, [GFF3 lines])] in file order, children with their gene."""
    rows = []
    parent_of = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            attrs = parse_attrs(cols[8])
            fid = attrs.get("ID")
            parent = attrs.get("Parent", "").split(",")[0] or None
            if fid:
                parent_of[fid] = parent
            rows.append((fid, parent, cols))

    def root(fid):
        seen = set()
        while parent_of.get(fid) and fid not in seen:
            seen.add(fid)
            fid = parent_of[fid]
        return fid

    genes = {}
    for fid, parent, cols in rows:
        gid = root(parent) if parent else fid
        genes.setdefault(gid, []).append(cols)
    return list(genes.items())


def span(lines, ftypes):
    return [(int(c[3]), int(c[4])) for c in lines if c[2] in ftypes]


def overlaps(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--models", required=True, help="protein-coding gene models GFF3")
    ap.add_argument("--rna", nargs="+", required=True, help="RNA gene GFF3s, rRNA first")
    ap.add_argument("--seqinfo", required=True, help="<genome>.seqinfo.tsv")
    ap.add_argument("--prefix", required=True, help="ID prefix without a locus_tag")
    ap.add_argument("--locus-tag", default="", help="the genome's locus_tag, if any")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    order = {}
    organelle = set()
    with open(args.seqinfo) as fh:
        next(fh)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            order[cols[0]] = len(order)
            if len(cols) > 3 and cols[3]:
                organelle.add(cols[0])

    # RNA genes: (contig, (start, end), type, provisional ID, lines)
    rnas = []
    dropped = Counter()
    for path in args.rna:
        for gid, lines in read_genes(path):
            contig = lines[0][0]
            rtype = next((c[2] for c in lines if c[2] in ("tRNA", "rRNA")), "RNA")
            loc = (int(lines[0][3]), int(lines[0][4]))
            if contig in organelle:
                dropped[f"{rtype} on an organelle sequence"] += 1
            elif any(r[0] == contig and overlaps(r[1], loc) for r in rnas):
                dropped[f"{rtype} overlapping another RNA gene"] += 1
            else:
                rnas.append((contig, loc, rtype, gid, lines))

    models = read_genes(args.models)
    cds = {}
    for gid, lines in models:
        for c in span(lines, {"CDS"}):
            cds.setdefault(lines[0][0], []).append((c, gid))
    kept = []
    removed = []
    for r in rnas:
        hit = next((gid for c, gid in cds.get(r[0], []) if overlaps(r[1], c)), None)
        if hit:
            removed.append(f"{r[2]} at {r[0]}:{r[1][0]}-{r[1][1]} (CDS of {hit})")
        else:
            kept.append(r)
    rnas = kept

    # Final RNA IDs, in genome order.
    rnas.sort(key=lambda r: (order.get(r[0], len(order)), r[1][0]))
    tag = args.locus_tag.rstrip("_")
    if tag:
        pat = re.compile(rf"^{re.escape(tag)}_(\d+)$")
        nums = [int(m.group(1)) for gid, _ in models if (m := pat.match(gid))]
        next_num = max(nums, default=0) + 1
    per_type = Counter()
    out_lines = []
    for contig, _, rtype, old, lines in rnas:
        per_type[rtype] += 1
        if tag:
            new = f"{tag}_{next_num:06d}"
            next_num += 1
        else:
            new = f"{args.prefix}_{rtype}{per_type[rtype]:05d}"
        for cols in lines:
            attrs = parse_attrs(cols[8])
            for key in ("ID", "Parent"):
                if key in attrs and attrs[key].startswith(old):
                    attrs[key] = new + attrs[key][len(old):]
            out_lines.append(cols[:8] + [";".join(f"{k}={v}" for k, v in attrs.items())])

    with open(args.out, "w") as out:
        out.write("##gff-version 3\n")
        for _, lines in models:
            for cols in lines:
                out.write("\t".join(cols) + "\n")
        for cols in out_lines:
            out.write("\t".join(cols) + "\n")

    print(f"{len(models)} protein-coding genes, {dict(per_type)} RNA genes "
          f"written to {args.out}", file=sys.stderr)
    for reason, n in dropped.items():
        print(f"  dropped {n} {reason}", file=sys.stderr)
    print(f"  dropped {len(removed)} RNA genes overlapping a protein-coding CDS",
          file=sys.stderr)
    for r in removed:
        print(f"    {r}", file=sys.stderr)


if __name__ == "__main__":
    main()
