#!/usr/bin/env python3
from __future__ import annotations
import sys
from collections import defaultdict

SPECIAL_SAMPLE = "DNA_RNA_Shield"

def parse_sample(header_with_gt: str) -> str:
    h = header_with_gt.strip()
    if h.startswith(">"):
        h = h[1:]
    if h.startswith(SPECIAL_SAMPLE + "_") or h == SPECIAL_SAMPLE:
        return SPECIAL_SAMPLE
    if "_" in h:
        return h.split("_", 1)[0]
    return h

def main() -> int:
    if len(sys.argv) != 3:
        print(f"Usage: python3 {sys.argv[0]} output.clstr 16s_OTU_table.tsv", file=sys.stderr)
        return 1

    clstr_path, out_tsv = sys.argv[1], sys.argv[2]

    counts = defaultdict(lambda: defaultdict(int))
    samples = set()
    clusters = set()

    nonempty = 0
    used = 0
    malformed = 0

    with open(clstr_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            nonempty += 1

            # split on ANY whitespace (tabs or spaces)
            parts = line.strip().split()
            if len(parts) < 2:
                malformed += 1
                continue

            try:
                cluster = int(parts[0])
            except ValueError:
                malformed += 1
                continue

            sample = parse_sample(parts[1])
            counts[cluster][sample] += 1
            samples.add(sample)
            clusters.add(cluster)
            used += 1

    samples = sorted(samples)
    clusters = sorted(clusters)

    with open(out_tsv, "w", encoding="utf-8") as out:
        out.write("OTU\t" + "\t".join(samples) + "\n")
        for c in clusters:
            row = [f"C{c}"]
            for s in samples:
                row.append(str(counts[c].get(s, 0)))
            out.write("\t".join(row) + "\n")

    total_table = sum(sum(d.values()) for d in counts.values())
    print(f"[OK] wrote {out_tsv}", file=sys.stderr)
    print(f"[INFO] non-empty lines in output.clstr: {nonempty}", file=sys.stderr)
    print(f"[INFO] lines used as reads: {used}", file=sys.stderr)
    print(f"[INFO] malformed non-empty lines skipped: {malformed}", file=sys.stderr)
    print(f"[INFO] total reads counted into table: {total_table}", file=sys.stderr)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
