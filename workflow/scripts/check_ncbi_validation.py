#!/usr/bin/env python3
"""
Summarise table2asn's checks of an NCBI submission and decide whether it is
ready to upload:

  * validator messages (*.val): ERROR and REJECT level must be fixed
    ("All Errors and Rejects need to be fixed", NCBI genome submission guide)
  * discrepancy report (*.dr): categories marked FATAL "are nearly always
    unacceptable and must be fixed". BACTERIAL_* categories are left out:
    NCBI says those do not apply to eukaryotes (they only show up when
    table2asn could not look up the taxonomy).

Writes a plain-text report and exits 1 if anything blocking is found,
unless --allow-errors is given.

Usage: check_ncbi_validation.py --outdir TABLE2ASN_OUT --report report.txt
           [--curate need-curating.txt] [--allow-errors]
"""
import argparse
import collections
import glob
import os
import re
import sys

BLOCKING = ("ERROR", "REJECT", "FATAL")
VAL_RE = re.compile(r"^(Info|Warning|Error|Reject|Fatal): valid \[([^\]]+)\] (.*)")


def read_val(outdir):
    counts = collections.Counter()      # (level, code) -> n
    examples = collections.defaultdict(list)
    for path in sorted(glob.glob(os.path.join(outdir, "*.val"))):
        with open(path) as fh:
            for line in fh:
                m = VAL_RE.match(line.strip())
                if not m:
                    continue
                key = (m.group(1).upper(), m.group(2))
                counts[key] += 1
                if len(examples[key]) < 3:
                    examples[key].append(m.group(3)[:300])
    return counts, examples


def read_dr(outdir):
    """FATAL lines from the summary part of the discrepancy report(s)."""
    fatal, ignored = [], []
    for path in sorted(glob.glob(os.path.join(outdir, "*.dr"))):
        with open(path) as fh:
            for line in fh:
                if line.startswith("Detailed Report"):
                    break
                line = line.strip()
                if line.startswith("FATAL:"):
                    target = ignored if line.split()[1].startswith("BACTERIAL_") else fatal
                    if line not in target:
                        target.append(line)
    return fatal, ignored


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--outdir", required=True, help="table2asn output directory")
    ap.add_argument("--report", required=True, help="report to write")
    ap.add_argument("--curate", help="funannotate2 need-curating list")
    ap.add_argument("--allow-errors", action="store_true",
                    help="exit 0 even if blocking problems are found")
    args = ap.parse_args()

    if not glob.glob(os.path.join(args.outdir, "*.sqn")):
        sys.exit(f"ERROR: table2asn wrote no .sqn in {args.outdir}; see the log")

    counts, examples = read_val(args.outdir)
    fatal, ignored = read_dr(args.outdir)
    n_blocking = sum(n for (level, _), n in counts.items() if level in BLOCKING) + len(fatal)

    lines = [f"NCBI submission check: {args.outdir}", ""]
    by_level = collections.Counter()
    for (level, _), n in counts.items():
        by_level[level] += n
    lines.append("Validator messages: " + (", ".join(
        f"{by_level[lv]} {lv}" for lv in ("REJECT", "FATAL", "ERROR", "WARNING", "INFO")
        if by_level[lv]) or "none"))
    for (level, code), n in sorted(counts.items(),
                                   key=lambda kv: (kv[0][0] not in BLOCKING, -kv[1])):
        lines.append(f"  {level:7} {code}: {n}")
        for ex in examples[(level, code)]:
            lines.append(f"            e.g. {ex}")
    lines.append("")
    lines.append(f"Discrepancy report FATAL categories: {len(fatal) or 'none'}")
    lines += [f"  {f}" for f in fatal]
    if ignored:
        lines.append("  (ignored, bacteria only: " + "; ".join(ignored) + ")")
    if args.curate and os.path.exists(args.curate):
        with open(args.curate) as fh:
            n_curate = sum(1 for l in fh if l.strip() and not l.startswith("#"))
        if n_curate:
            lines += ["", f"{n_curate} gene names/products flagged for manual curation "
                          f"in {args.curate}"]
    lines.append("")
    if n_blocking:
        lines.append(f"NOT READY: {n_blocking} problems that NCBI requires to be fixed. "
                     f"Details in the .val and .dr files in {args.outdir}.")
    else:
        lines.append("READY: no validator errors and no FATAL discrepancies. "
                     "Review the warnings and the .dr file before uploading the .sqn.")

    text = "\n".join(lines) + "\n"
    with open(args.report, "w") as out:
        out.write(text)
    print(text)
    if n_blocking and not args.allow_errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
