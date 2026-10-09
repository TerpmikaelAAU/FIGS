#!/usr/bin/env python3
"""
Write one row per sequence of an assembly: its ID, the ID it gets in an NCBI
submission, its chromosome name and, for organelle sequences, its location
(e.g. mitochondrion). gene_models leaves organelle sequences out of the gene
models and ncbi_submission labels them.

Without --ncbi-update, IDs are kept and organelles are taken from
--organelle ID=location.

With --ncbi-update <assembly accession>, the assembly is an unannotated
GenBank WGS genome that is getting annotation as an update. NCBI wants each
sequence named gnl|WGS:<prefix>|<SeqID>|gb|<contig accession>, where SeqID
is the name it had in the original submission. Those come from NCBI:
  * the assembly's WGS prefix (e.g. JBQFLM01 -> JBQFLM) and sequence report
    (GenBank accession, submitted name, chromosome, location) from the
    Datasets API
  * the contig accession behind each chromosome record (CM...), from the
    CONTIG line of that record via E-utilities.
The FASTA staged from NCBI is named by GenBank accession, so every FASTA ID
must be in the sequence report.

Usage: sequence_info.py --fasta genome.fna [--ncbi-update GCA_...]
                        [--organelle ID=mitochondrion ...] -o seqinfo.tsv
"""
import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/{}"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def fetch(url, params=None, tries=4):
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception as e:  # network hiccups, NCBI rate limits
            if attempt == tries - 1:
                sys.exit(f"ERROR: could not fetch {url}: {e}")
            time.sleep(5 * (attempt + 1))


def sequence_report(accession):
    reports, token = [], None
    while True:
        params = {"page_size": 1000}
        if token:
            params["page_token"] = token
        data = json.loads(fetch(DATASETS.format(accession, "sequence_reports"), params))
        reports += data.get("reports", [])
        token = data.get("next_page_token")
        if not token:
            return reports


def contig_accessions(accessions):
    """{record accession.version: [WGS contig accession.version, ...]} from
    the CONTIG lines of CON records (e.g. chromosomes)."""
    result = {}
    for i in range(0, len(accessions), 50):
        batch = accessions[i:i + 50]
        text = fetch(EFETCH, {"db": "nuccore", "id": ",".join(batch),
                              "rettype": "gb", "retmode": "text"})
        for record in text.split("\n//"):
            version = re.search(r"^VERSION\s+(\S+)", record, re.M)
            contig = re.search(r"^CONTIG\s+(.*?)(?=^\S)", record + "\nEND", re.M | re.S)
            if version and contig:
                result[version.group(1)] = re.findall(
                    r"([A-Z]{4,6}\d{8,}\.\d+):", contig.group(1))
        time.sleep(0.4)  # stay under E-utilities' 3 requests/second
    return result


def ncbi_update_rows(accession, ids):
    report = json.loads(fetch(DATASETS.format(accession, "dataset_report")))["reports"][0]
    wgs = (report.get("wgs_info") or {}).get("wgs_project_accession", "")
    if not wgs:
        sys.exit(f"ERROR: {accession} is not a WGS assembly; the WGS update "
                 f"route (ncbi_submission: update) does not apply to it.")
    prefix = re.sub(r"\d+$", "", wgs)  # JBQFLM01 -> JBQFLM
    print(f"{accession}: WGS project {wgs}, prefix {prefix}", file=sys.stderr)

    by_acc = {r["genbank_accession"]: r for r in sequence_report(accession)
              if r.get("genbank_accession")}
    missing = [i for i in ids if i not in by_acc]
    if missing:
        sys.exit(f"ERROR: {len(missing)} FASTA IDs are not GenBank accessions "
                 f"of {accession}, e.g. {', '.join(missing[:3])}")

    # Records that are not WGS contigs themselves (chromosomes, CM...) point
    # to their contig(s) through their CONTIG line.
    con = [i for i in ids if not i.startswith(prefix)]
    contigs = contig_accessions(con) if con else {}

    rows = []
    for seqid in ids:
        r = by_acc[seqid]
        if seqid.startswith(prefix):
            contig = seqid
        else:
            parts = contigs.get(seqid, [])
            if len(parts) != 1:
                sys.exit(f"ERROR: {seqid} is built from {len(parts)} WGS contigs "
                         f"({', '.join(parts) or 'none found'}); an update must "
                         f"name each contig, which this pipeline does not do.")
            contig = parts[0]
        name = r.get("sequence_name") or seqid
        submission_id = f"gnl|WGS:{prefix}|{name}|gb|{contig.split('.')[0]}"
        location_type = r.get("assigned_molecule_location_type", "")
        location = "" if location_type in ("", "Chromosome") else location_type.lower()
        chromosome = (r.get("chr_name", "") if location_type == "Chromosome"
                      and r.get("role") == "assembled-molecule" else "")
        rows.append((seqid, submission_id, chromosome, location))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fasta", required=True, help="assembly FASTA")
    ap.add_argument("--ncbi-update", metavar="ACCESSION",
                    help="GenBank assembly the annotation is an update of")
    ap.add_argument("--organelle", action="append", default=[], metavar="ID=LOCATION",
                    help="organelle sequence, e.g. contig_9=mitochondrion")
    ap.add_argument("-o", "--out", required=True, help="output TSV")
    args = ap.parse_args()

    ids = fasta_ids(args.fasta)
    if args.ncbi_update:
        rows = ncbi_update_rows(args.ncbi_update, ids)
    else:
        rows = [(i, i, "", "") for i in ids]

    organelles = dict(o.split("=", 1) for o in args.organelle)
    unknown = set(organelles) - set(ids)
    if unknown:
        sys.exit(f"ERROR: organelles not in {args.fasta}: {', '.join(sorted(unknown))}")
    rows = [(s, sub, "" if s in organelles else chrom, organelles.get(s, loc))
            for s, sub, chrom, loc in rows]

    with open(args.out, "w") as out:
        out.write("seqid\tsubmission_id\tchromosome\tlocation\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    n_org = sum(1 for r in rows if r[3])
    print(f"{len(rows)} sequences, {n_org} organelle", file=sys.stderr)
    for r in rows:
        print("\t".join(r), file=sys.stderr)


if __name__ == "__main__":
    main()
