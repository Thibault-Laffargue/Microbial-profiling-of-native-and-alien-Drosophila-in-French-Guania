#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path
from collections import defaultdict

# 1) Regex robuste : cluster_id + (nom avec espaces) + identity + type
# Ex: "12   my seq name with spaces   0.97   C"
LINE_RE = re.compile(r"^\s*(\S+)\s+(.*?)\s+(\S+)\s+(\S+)\s*$")

TYPE_LETTER_RE = re.compile(r"[A-Za-z]")

def normalize_type(raw: str) -> str:
    raw = (raw or "").strip()
    if not raw:
        return ""
    m = TYPE_LETTER_RE.search(raw)
    return m.group(0).upper() if m else ""

def iter_id_candidates(s: str) -> list[str]:
    """
    Given a raw sequence name (from clstr or fasta header), return a list of
    candidate IDs that might match between files.

    Keeps order from most-specific to most-likely.
    """
    s = (s or "").strip()
    if not s:
        return []

    # remove leading '>' if present
    if s.startswith(">"):
        s = s[1:].strip()

    cands = []

    def add(x: str):
        x = (x or "").strip()
        if not x:
            return
        if x not in cands:
            cands.append(x)

    # full string as-is (can include spaces)
    add(s)

    # first token
    tok = s.split()[0]
    add(tok)

    # common trims on token
    for sep in ["|", ";", ",", ":", "/"]:
        if sep in tok:
            add(tok.split(sep)[0])

    # strip common trailing punctuation
    add(tok.rstrip(";,.:"))

    # sometimes IDs are like "id/1" or "id/2"
    add(tok.split("/")[0])

    # sometimes IDs are like "id;size=123"
    add(tok.split(";")[0])

    # also try removing quotes
    add(tok.strip('"').strip("'"))

    return cands

def read_fasta_multiindex(fasta_path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """
    Read FASTA and build:
      - seq_by_key: maps MANY possible keys -> sequence
      - primary_by_key: maps key -> primary_id (for traceability)

    Primary id = first token of header.
    """
    seq_by_key: dict[str, str] = {}
    primary_by_key: dict[str, str] = {}

    header_raw = None
    primary = None
    chunks: list[str] = []

    def flush():
        nonlocal header_raw, primary, chunks
        if primary is None:
            return
        seq = "".join(chunks)

        # index by multiple candidates derived from:
        # - full header (without leading >)
        # - primary token
        # - variants of primary token
        for key in iter_id_candidates(header_raw):
            seq_by_key.setdefault(key, seq)
            primary_by_key.setdefault(key, primary)

        header_raw = None
        primary = None
        chunks = []

    with fasta_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                flush()
                header_raw = line[1:].strip()
                primary = header_raw.split()[0] if header_raw else None
            else:
                chunks.append(line.strip().replace(" ", "").replace("\t", ""))
        flush()

    return seq_by_key, primary_by_key

def parse_clstr_centroids(clstr_path: Path):
    """
    Returns:
      counts[cluster_id] = number of entries (lines) in cluster
      centroid_names[cluster_id] = list of raw centroid seq_name strings (may include spaces)
    """
    counts = defaultdict(int)
    centroid_names = defaultdict(list)

    with clstr_path.open("r", encoding="utf-8", errors="replace") as f:
        for lineno, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line.strip():
                continue

            m = LINE_RE.match(line)
            if not m:
                # fallback: very lenient split, but keep moving
                parts = line.split()
                if len(parts) < 4:
                    continue
                cluster_id = parts[0]
                raw_type = parts[-1]
                seq_name = " ".join(parts[1:-2])
            else:
                cluster_id = m.group(1)
                seq_name = m.group(2)
                raw_type = m.group(4)

            counts[cluster_id] += 1
            if normalize_type(raw_type) == "C":
                centroid_names[cluster_id].append(seq_name)

    return counts, centroid_names

def write_fasta_record(out_handle, header: str, seq: str, width: int = 80):
    out_handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        out_handle.write(seq[i:i + width] + "\n")

def main():
    ap = argparse.ArgumentParser(
        description="Extract centroid sequences for clusters >10 from MeshClust output and write to a new FASTA."
    )
    ap.add_argument("--clstr", default="output.clstr", help="MeshClust output file (default: output.clstr)")
    ap.add_argument("--fasta", default="all.sample.min1000.fas", help="Input FASTA used for clustering")
    ap.add_argument("--min-size", type=int, default=11, help="Minimum cluster size to keep (default: 11 => >10)")
    ap.add_argument("--out", default="centroids_gt10.fasta", help="Output FASTA (default: centroids_gt10.fasta)")
    ap.add_argument("--debug-missing", default="missing_centroids.tsv",
                    help="Write missing centroid lookup info here (default: missing_centroids.tsv)")
    args = ap.parse_args()

    clstr_path = Path(args.clstr)
    fasta_path = Path(args.fasta)
    out_path = Path(args.out)
    debug_path = Path(args.debug_missing)

    if not clstr_path.exists():
        raise SystemExit(f"Missing clstr file: {clstr_path}")
    if not fasta_path.exists():
        raise SystemExit(f"Missing fasta file: {fasta_path}")

    counts, centroid_names = parse_clstr_centroids(clstr_path)
    seq_by_key, primary_by_key = read_fasta_multiindex(fasta_path)

    kept_clusters = [cid for cid, n in counts.items() if n >= args.min_size]
    kept_clusters.sort(key=lambda x: int(x) if x.isdigit() else x)

    written = 0
    skipped_no_centroid = 0
    skipped_missing_seq = 0
    warnings_multi_centroid = 0

    missing_rows = []
    # header: cluster_id, raw_centroid_name, tried_keys (pipe-separated)
    # plus an indicator if any key existed (should be "NA" if missing)
    with out_path.open("w", encoding="utf-8") as out:
        for cid in kept_clusters:
            c_list = centroid_names.get(cid, [])
            if not c_list:
                skipped_no_centroid += 1
                continue

            if len(c_list) > 1:
                warnings_multi_centroid += 1

            raw_centroid = c_list[0]
            candidates = iter_id_candidates(raw_centroid)

            seq = None
            matched_key = None
            for k in candidates:
                if k in seq_by_key:
                    seq = seq_by_key[k]
                    matched_key = k
                    break

            if seq is None:
                skipped_missing_seq += 1
                missing_rows.append((cid, raw_centroid, "|".join(candidates)))
                continue

            # Rename header: "C" + cluster_id (as requested)
            write_fasta_record(out, f"C{cid}", seq)
            written += 1

    # Write debug file if missing
    if skipped_missing_seq:
        with debug_path.open("w", encoding="utf-8") as dbg:
            dbg.write("cluster_id\traw_centroid_name\tcandidate_keys_tried\n")
            for row in missing_rows[:20000]:  # safety cap
                dbg.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")

    print(f"[DONE] clusters total: {len(counts)}")
    print(f"[DONE] clusters kept (size >= {args.min_size}): {len(kept_clusters)}")
    print(f"[DONE] centroid sequences written: {written} -> {out_path.resolve()}")
    if skipped_no_centroid:
        print(f"[WARN] clusters kept but no centroid row detected (type=C): {skipped_no_centroid}")
    if skipped_missing_seq:
        print(f"[WARN] centroid key not found in FASTA headers: {skipped_missing_seq}")
        print(f"[WARN] wrote debug file: {debug_path.resolve()}")
    if warnings_multi_centroid:
        print(f"[WARN] clusters with multiple centroid rows (used first): {warnings_multi_centroid}")

if __name__ == "__main__":
    main()
