#!/usr/bin/env python3
"""
Download a URL to a file, resuming where it left off when the connection
drops, and check its checksum. Used for the one-off downloads (OMArk's
~10 GB LUCA.h5, table2asn, Nextflow, Java), where urllib's urlretrieve gives
up at the first interruption.

The partial download is kept as OUT.part, so a re-run of the rule resumes
it too. A checksum mismatch deletes it.

Usage: download.py URL OUT [--sha256 HEX | --md5 HEX] [--tries 20]
"""
import argparse
import hashlib
import http.client
import os
import sys
import time
import urllib.error
import urllib.request

CHUNK = 1 << 20


def fetch(url, part, tries):
    for attempt in range(1, tries + 1):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers),
                                        timeout=120) as r:
                if have and r.status != 206:  # range ignored: start again
                    have = 0
                length = r.headers.get("Content-Length")
                total = have + int(length) if length is not None else None
                with open(part, "ab" if have else "wb") as out:
                    while chunk := r.read(CHUNK):
                        out.write(chunk)
            size = os.path.getsize(part)
            if total is None or size == total:
                return
            print(f"attempt {attempt}/{tries}: connection closed at {size} of "
                  f"{total} bytes", file=sys.stderr)
        except urllib.error.HTTPError as e:
            if e.code == 416 and have:  # nothing left to fetch
                return
            if e.code < 500 and e.code != 429:
                sys.exit(f"ERROR: {url}: {e}")
            print(f"attempt {attempt}/{tries}: {e}", file=sys.stderr)
        except (urllib.error.URLError, http.client.HTTPException, OSError) as e:
            print(f"attempt {attempt}/{tries}: {e}", file=sys.stderr)
        time.sleep(min(60, 5 * attempt))
    sys.exit(f"ERROR: could not download {url} in {tries} attempts")


def file_hash(path, algorithm):
    h = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("url")
    ap.add_argument("out")
    check = ap.add_mutually_exclusive_group()
    check.add_argument("--sha256", default="")
    check.add_argument("--md5", default="")
    ap.add_argument("--tries", type=int, default=20)
    args = ap.parse_args()

    part = args.out + ".part"
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fetch(args.url, part, args.tries)

    algorithm, expected = ("sha256", args.sha256) if args.sha256 else ("md5", args.md5)
    if expected:
        got = file_hash(part, algorithm)
        if got != expected.lower():
            os.remove(part)
            sys.exit(f"ERROR: {args.url}: {algorithm} {got}, expected {expected}")
        print(f"{algorithm} OK: {got}", file=sys.stderr)
    os.replace(part, args.out)
    print(f"{args.url} -> {args.out} ({os.path.getsize(args.out)} bytes)", file=sys.stderr)


if __name__ == "__main__":
    main()
