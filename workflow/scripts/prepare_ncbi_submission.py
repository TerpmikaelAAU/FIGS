#!/usr/bin/env python3
"""
Turn funannotate2's feature table and the assembly into table2asn input:
<name>.fsa, <name>.tbl, plus a source table of per-sequence qualifiers.

  * Sequences are renamed to their submission IDs from sequence_info (for an
    NCBI update, gnl|WGS:<prefix>|<SeqID>|gb|<contig accession>), and the
    sequence is written in upper case.
  * The .tbl's per-sequence REFERENCE block (gfftk writes a placeholder
    "CFMR 12345" reference, which table2asn rejects as an unknown
    qualifier) is removed.
  * For an update, protein_id/transcript_id move from gnl|ncbi|... to the
    project's gnl|WGS:<prefix>|..., as NCBI's WGS update instructions ask.
  * Features on organelle sequences are dropped (gene_models should already
    have left them out).
  * The source table gives chromosome names (from NCBI, for an update) and
    marks organelle sequences with their location and genetic code.

Usage: prepare_ncbi_submission.py --fasta genome.fna --tbl annotation.tbl
           --seqinfo seqinfo.tsv --indir DIR --name NAME --src source.tsv
           [--mito-gcode 4]
"""
import argparse
import csv
import os
import re
import sys


def read_seqinfo(path):
    with open(path) as fh:
        return {r["seqid"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fasta", required=True, help="assembly FASTA (unmasked)")
    ap.add_argument("--tbl", required=True, help="funannotate2 feature table")
    ap.add_argument("--seqinfo", required=True, help="sequence_info TSV")
    ap.add_argument("--indir", required=True, help="table2asn input directory")
    ap.add_argument("--name", required=True, help="basename of the .fsa/.tbl")
    ap.add_argument("--src", required=True, help="source qualifier table to write")
    ap.add_argument("--mito-gcode", default="4", help="mitochondrial genetic code")
    args = ap.parse_args()

    info = read_seqinfo(args.seqinfo)
    os.makedirs(args.indir, exist_ok=True)

    # The WGS dbname for protein/transcript IDs, when this is an update.
    wgs_db = ""
    for r in info.values():
        m = re.match(r"gnl\|(WGS:[A-Z0-9]+)\|", r["submission_id"])
        if m:
            wgs_db = m.group(1)
            break

    seen = set()
    with open(args.fasta) as fin, open(os.path.join(args.indir, f"{args.name}.fsa"), "w") as fout:
        for line in fin:
            if line.startswith(">"):
                seqid = line[1:].split()[0]
                if seqid not in info:
                    sys.exit(f"ERROR: {seqid} from {args.fasta} is not in {args.seqinfo}")
                seen.add(seqid)
                fout.write(f">{info[seqid]['submission_id']}\n")
            else:
                fout.write(line.upper())
    if seen != set(info):
        sys.exit(f"ERROR: {args.fasta} and {args.seqinfo} list different sequences")

    n_features = n_dropped = 0
    with open(args.tbl) as fin, open(os.path.join(args.indir, f"{args.name}.tbl"), "w") as fout:
        skip_seq = in_reference = False
        for line in fin:
            if line.startswith(">Feature"):
                seqid = line.split()[1]
                if seqid not in info:
                    sys.exit(f"ERROR: {args.tbl} has features on {seqid}, "
                             f"which is not in {args.fasta}")
                skip_seq = bool(info[seqid]["location"])
                if skip_seq:
                    print(f"WARNING: dropping features on organelle sequence {seqid}",
                          file=sys.stderr)
                    continue
                fout.write(f">Feature {info[seqid]['submission_id']}\n")
                continue
            if skip_seq:
                cols = line.split("\t")
                if line[0] != "\t" and len(cols) >= 3 and cols[2].strip() not in ("", "REFERENCE"):
                    n_dropped += 1
                continue
            cols = line.rstrip("\n").split("\t")
            # A feature line has a location in columns 1-2 and a key in 3;
            # qualifier lines start with three empty columns.
            if line[0] != "\t":
                in_reference = len(cols) >= 3 and cols[2] == "REFERENCE"
                if in_reference:
                    continue
                if len(cols) >= 3 and cols[2]:
                    n_features += 1
            elif in_reference:
                continue
            if wgs_db:
                line = line.replace("\tgnl|ncbi|", f"\tgnl|{wgs_db}|")
            fout.write(line)

    rows = []
    for seqid, r in info.items():
        if r["location"]:
            rows.append((r["submission_id"], "", r["location"],
                         args.mito_gcode if r["location"] == "mitochondrion" else ""))
        elif r["chromosome"]:
            rows.append((r["submission_id"], r["chromosome"], "", ""))
    with open(args.src, "w") as out:
        if rows:
            out.write("Sequence_ID\tchromosome\tlocation\tmgcode\n")
            for row in rows:
                out.write("\t".join(row) + "\n")

    print(f"{len(info)} sequences, {n_features} features written"
          + (f", {n_dropped} features on organelles dropped" if n_dropped else "")
          + (f", protein/transcript IDs in gnl|{wgs_db}|" if wgs_db else ""),
          file=sys.stderr)


if __name__ == "__main__":
    main()
