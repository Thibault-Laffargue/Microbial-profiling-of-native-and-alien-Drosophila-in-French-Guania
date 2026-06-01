#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

def load_primer_map(assoc_path: Path) -> dict[str, str]:
    """
    Read TSV: col1=Primer association like 'F1-R1' ; col2=sample name.
    Returns mapping with normalized key 'F1_R1' -> sample
    """
    mapping: dict[str, str] = {}
    with assoc_path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 2:
                continue
            primer_raw = row[0].strip()
            sample = row[1].strip()
            if not primer_raw:
                continue
            key = primer_raw.replace("-", "_")  # F1-R2 -> F1_R2
            mapping[key] = sample
    return mapping

def extract_primer_from_filename(name: str) -> str | None:
    """
    From 'F1_R2_cleaned.fastq' -> 'F1_R2'
    Works even if extra stuff exists, as long as it contains F<d+> and R<d+>.
    """
    mF = re.search(r"(F\d+)", name)
    mR = re.search(r"(R\d+)", name)
    if not (mF and mR):
        return None
    return f"{mF.group(1)}_{mR.group(1)}"

def safe_rename(src: Path, dst: Path) -> tuple[bool, str]:
    """
    Rename src->dst if possible.
    If dst exists, do NOT overwrite; return False with reason.
    """
    if dst.exists():
        return False, f"target_exists:{dst.name}"
    src.rename(dst)
    return True, "renamed"

def main():
    ap = argparse.ArgumentParser(
        description="Rename *_cleaned.fastq using primer->sample map and add sample_name column to quality_filter_stats.tsv"
    )
    ap.add_argument("--assoc", default="16s_PCR_plate_plan_Primer_Association.fixed",
                    help="TSV mapping primer association to sample name")
    ap.add_argument("--stats", default="quality_filter_stats.tsv",
                    help="TSV stats file to annotate")
    ap.add_argument("--dry-run", action="store_true",
                    help="Show what would happen without renaming/modifying files")
    args = ap.parse_args()

    assoc_path = Path(args.assoc)
    stats_path = Path(args.stats)

    if not assoc_path.exists():
        raise SystemExit(f"Mapping file not found: {assoc_path}")
    if not stats_path.exists():
        raise SystemExit(f"Stats file not found: {stats_path}")

    primer_map = load_primer_map(assoc_path)

    # --- Step 1: rename cleaned fastq files (only when a mapping exists)
    cleaned_files = sorted(Path(".").glob("*_cleaned.fastq"))
    rename_log = []

    for f in cleaned_files:
        primer = extract_primer_from_filename(f.name)
        if primer is None or primer not in primer_map:
            rename_log.append((f.name, "no_match", "NA"))
            continue

        sample = primer_map[primer]
        new_name = f"{sample}_cleaned.fastq"
        dst = f.with_name(new_name)

        if args.dry_run:
            status = "would_rename" if not dst.exists() else "would_skip_target_exists"
            rename_log.append((f.name, status, sample))
        else:
            ok, reason = safe_rename(f, dst)
            rename_log.append((f.name, reason if ok else "skip", sample if ok else sample))

    # --- Step 2: add sample_name column to stats TSV
    tmp_out = stats_path.with_suffix(stats_path.suffix + ".tmp")
    # We will match using output_file column (2nd column per your description),
    # but we also read by header name "output_file" if present.
    with stats_path.open("r", encoding="utf-8", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("Stats file seems empty or has no header.")
        fieldnames = list(reader.fieldnames)
        if "sample_name" not in fieldnames:
            fieldnames.append("sample_name")

        rows = []
        for row in reader:
            out_file = row.get("output_file", "") or ""
            primer = extract_primer_from_filename(out_file)
            if primer is not None and primer in primer_map:
                row["sample_name"] = primer_map[primer]
            else:
                row["sample_name"] = "NA"
            rows.append(row)

    if args.dry_run:
        print("[DRY RUN] Would update stats file by adding column 'sample_name'.")
    else:
        with tmp_out.open("w", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        tmp_out.replace(stats_path)

    # --- Print summary
    renamed = sum(1 for _, st, _ in rename_log if st in ("renamed", "would_rename"))
    skipped = len(rename_log) - renamed
    print(f"[DONE] processed {len(cleaned_files)} *_cleaned.fastq files | renamed={renamed} skipped={skipped}")
    print(f"[DONE] stats annotated: {stats_path}")

    # Optional: write a log file for traceability
    log_path = Path("rename_log.tsv")
    header = ["original_file", "status", "sample_name_or_NA"]
    if args.dry_run:
        print("[DRY RUN] Would write rename_log.tsv")
    else:
        with log_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(header)
            w.writerows(rename_log)
        print(f"[DONE] rename log: {log_path}")

if __name__ == "__main__":
    main()
