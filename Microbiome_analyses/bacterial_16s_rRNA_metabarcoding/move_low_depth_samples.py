#!/usr/bin/env python3
import csv
import shutil
from pathlib import Path

STATS_FILE = Path("quality_filter_stats.tsv")
OUT_DIR = Path("filtered_out")
CLEAN_SUFFIX = "_cleaned.fastq"

EXCLUDE_NAMES = {"NA", "Negative", "negative", "DNA_RNA_Shield"}

THRESHOLD = 1000  # cleaned_reads < 1000

def main():
    if not STATS_FILE.exists():
        raise SystemExit(f"Stats file not found: {STATS_FILE}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) Move all F*_cleaned.fastq to filtered_out
    moved_F = 0
    for f in sorted(Path(".").glob(f"F*{CLEAN_SUFFIX}")):
        if f.is_file():
            shutil.move(str(f), str(OUT_DIR / f.name))
            moved_F += 1
    print(f"[INFO] Moved {moved_F} files matching F*{CLEAN_SUFFIX} -> {OUT_DIR}/")

    # 2) Parse stats and collect sample_name where cleaned_reads < 1000
    low_samples = set()

    with STATS_FILE.open("r", encoding="utf-8", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit("quality_filter_stats.tsv has no header or is empty.")

        # Required columns
        required = {"sample_name", "cleaned_reads"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise SystemExit(f"Missing columns in stats file: {', '.join(sorted(missing))}")

        for row in reader:
            sample = (row.get("sample_name") or "").strip()
            if not sample or sample in EXCLUDE_NAMES:
                continue

            cr = (row.get("cleaned_reads") or "").strip()
            try:
                cleaned_reads = int(cr)
            except ValueError:
                continue

            if cleaned_reads < THRESHOLD:
                low_samples.add(sample)

    print(f"[INFO] Found {len(low_samples)} sample(s) with cleaned_reads < {THRESHOLD} (excluding NA/Negative/DNA_RNA_Shield).")

    # 3) Move SAMPLENAME_cleaned.fastq if SAMPLENAME in low_samples
    moved_samples = 0
    skipped_missing = 0

    for sample in sorted(low_samples):
        f = Path(f"{sample}{CLEAN_SUFFIX}")
        if f.exists() and f.is_file():
            shutil.move(str(f), str(OUT_DIR / f.name))
            moved_samples += 1
        else:
            skipped_missing += 1

    print(f"[DONE] Moved {moved_samples} sample file(s) -> {OUT_DIR}/")
    if skipped_missing:
        print(f"[WARN] {skipped_missing} sample(s) listed in stats had no corresponding file '{sample}{CLEAN_SUFFIX}' in this folder.")

if __name__ == "__main__":
    main()
