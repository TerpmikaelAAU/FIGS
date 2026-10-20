#!/usr/bin/env python3
"""
Turn tRNAscan-SE and Infernal cmsearch (Rfam rRNA models) results into
GFF3 gene models (gene -> tRNA/rRNA -> exon) that gfftk and funannotate2
read. IDs are provisional; merge_ncrna.py gives the final ones.

trna: tRNAscan-SE's tabular output (-o).
  * Products are written the way NCBI expects (tRNA-Met for iMet/fMet,
    tRNA-Sec for SeC). funannotate2 predict's own conversion writes
    tRNA-iMet, tRNA-SeC and tRNA-Undet, which NCBI rejects.
  * Pseudogenes, undetermined (Undet), suppressor (Sup) and unknown
    isotypes are dropped, and so are tRNAs outside 50-150 nt (as
    funannotate2 predict does, for NCBI's length rules).
  * Introns become a two-exon tRNA.

rrna: cmsearch --tblout with --cut_ga (Rfam's gathering thresholds).
  * Only hits covering at least --min-coverage of the model are kept;
    shorter fragments (degenerate copies) are dropped and counted. So are
    rRNAs cut off by a contig end (cmsearch's trunc column): gfftk cannot
    mark an rRNA as partial, and NCBI would take it as complete.
  * Where hits of different models overlap, the higher-scoring one wins.

Usage: ncrna_gff3.py trna -i tRNAscan.out -o tRNA.gff3
       ncrna_gff3.py rrna -i rRNA.tblout --cm rRNA.cm -o rRNA.gff3
"""
import argparse
import sys
from collections import Counter

AMINO_ACIDS = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His",
               "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
               "Tyr", "Val", "Pyl"}
ISOTYPE_ALIASES = {"iMet": "Met", "fMet": "Met", "SeC": "Sec", "Ile2": "Ile"}

# Rfam model name -> product, as NCBI names fungal rRNAs.
RRNA_PRODUCTS = {
    "5S_rRNA": "5S ribosomal RNA",
    "5_8S_rRNA": "5.8S ribosomal RNA",
    "SSU_rRNA_eukarya": "small subunit ribosomal RNA",
    "LSU_rRNA_eukarya": "large subunit ribosomal RNA",
}


def gff_lines(contig, source, rtype, gid, exons, strand, score, attrs):
    """gene -> RNA -> exon lines; exons as (start, end), start <= end."""
    start = min(s for s, _ in exons)
    end = max(e for _, e in exons)
    tid = f"{gid}-T1"
    attr = ";".join(f"{k}={v}" for k, v in attrs.items())
    out = [f"{contig}\t{source}\tgene\t{start}\t{end}\t{score}\t{strand}\t.\tID={gid}",
           f"{contig}\t{source}\t{rtype}\t{start}\t{end}\t{score}\t{strand}\t.\t"
           f"ID={tid};Parent={gid};{attr}"]
    ordered = sorted(exons, reverse=(strand == "-"))
    for n, (s, e) in enumerate(ordered, 1):
        out.append(f"{contig}\t{source}\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                   f"ID={tid}.exon{n};Parent={tid}")
    return out


def trna(args):
    lines = []
    dropped = Counter()
    n = 0
    with open(args.input) as fh:
        for line in fh:
            cols = [c.strip() for c in line.rstrip("\n").split("\t")]
            if len(cols) < 9 or not cols[1].isdigit():
                continue  # the three header lines
            contig, _, begin, end, isotype, anticodon, ib, ie, score = cols[:9]
            note = " ".join(cols[9:])
            begin, end, ib, ie = int(begin), int(end), int(ib), int(ie)
            aa = ISOTYPE_ALIASES.get(isotype, isotype)
            if "pseudo" in note:
                dropped["pseudogene"] += 1
                continue
            if aa not in AMINO_ACIDS:
                dropped[f"isotype {isotype}"] += 1
                continue
            length = abs(end - begin) + 1
            if length < 50 or length > 150:
                dropped["length outside 50-150 nt"] += 1
                continue
            strand = "+" if begin < end else "-"
            lo, hi = min(begin, end), max(begin, end)
            if ib and ie:
                ilo, ihi = min(ib, ie), max(ib, ie)
                exons = [(lo, ilo - 1), (ihi + 1, hi)]
            else:
                exons = [(lo, hi)]
            n += 1
            lines += gff_lines(contig, "tRNAscan-SE", "tRNA", f"tRNA_{n}", exons,
                               strand, score,
                               {"product": f"tRNA-{aa}",
                                "note": f"anticodon {anticodon}"})
    write(args.out, lines)
    print(f"{args.input}: {n} tRNAs written; dropped: "
          f"{dict(dropped) or 'none'}", file=sys.stderr)


def model_lengths(cm_file):
    lengths = {}
    name = None
    with open(cm_file) as fh:
        for line in fh:
            if line.startswith("NAME"):
                name = line.split()[1]
            elif line.startswith("CLEN") and name:
                lengths.setdefault(name, int(line.split()[1]))
    return lengths


def rrna(args):
    clen = model_lengths(args.cm)
    hits = []
    short = Counter()
    truncated = Counter()
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            contig, model, acc = f[0], f[2], f[3]
            mfrom, mto, sfrom, sto = int(f[5]), int(f[6]), int(f[7]), int(f[8])
            strand, trunc, score, inc = f[9], f[10], float(f[14]), f[16]
            if inc != "!" or model not in RRNA_PRODUCTS:
                continue
            if trunc != "no":
                truncated[model] += 1
                continue
            cover = (mto - mfrom + 1) / clen[model]
            if cover < args.min_coverage:
                short[model] += 1
                continue
            hits.append((score, contig, min(sfrom, sto), max(sfrom, sto),
                         strand, model, acc, cover))

    # Best-scoring first; drop hits overlapping one already kept.
    kept = []
    n_overlap = 0
    for hit in sorted(hits, reverse=True):
        _, contig, start, end = hit[:4]
        if any(k[1] == contig and k[2] <= end and start <= k[3] for k in kept):
            n_overlap += 1
            continue
        kept.append(hit)

    lines = []
    per_model = Counter()
    for n, (score, contig, start, end, strand, model, acc, _) in enumerate(
            sorted(kept, key=lambda h: (h[1], h[2])), 1):
        per_model[model] += 1
        lines += gff_lines(contig, "Infernal", "rRNA", f"rRNA_{n}", [(start, end)],
                           strand, f"{score:g}",
                           {"product": RRNA_PRODUCTS[model],
                            "note": f"Rfam:{acc} {model}"})
    write(args.out, lines)
    print(f"{args.input}: {len(kept)} rRNAs written {dict(per_model)}; dropped "
          f"{sum(short.values())} covering < {args.min_coverage:.0%} of the model "
          f"{dict(short)}, {sum(truncated.values())} cut off by a contig end "
          f"{dict(truncated)}, {n_overlap} overlapping a better hit", file=sys.stderr)


def write(path, lines):
    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for line in lines:
            out.write(line + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("trna", help="tRNAscan-SE output to GFF3")
    t.add_argument("-i", "--input", required=True, help="tRNAscan-SE -o output")
    t.add_argument("-o", "--out", required=True, help="GFF3")
    r = sub.add_parser("rrna", help="cmsearch --tblout to GFF3")
    r.add_argument("-i", "--input", required=True, help="cmsearch --tblout")
    r.add_argument("--cm", required=True, help="the CM file searched (for model lengths)")
    r.add_argument("--min-coverage", type=float, default=0.8,
                   help="minimum fraction of the model a hit must cover (default 0.8)")
    r.add_argument("-o", "--out", required=True, help="GFF3")
    args = ap.parse_args()
    trna(args) if args.cmd == "trna" else rrna(args)


if __name__ == "__main__":
    main()
