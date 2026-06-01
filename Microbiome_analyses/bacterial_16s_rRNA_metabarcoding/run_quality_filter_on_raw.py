#!/usr/bin/env python3
import argparse
import subprocess
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

BARCODE_RE = re.compile(r"(?:^|_)F(\d+).*?R(\d+)(?:_|$)")

def count_reads_fastq(path: Path) -> int:
    """Count FASTQ reads (4 lines per read). Uses binary mode for speed/robustness."""
    n_lines = 0
    with path.open("rb") as f:
        for _ in f:
            n_lines += 1
    return n_lines // 4

def run_quality_filter(in_path: Path, out_path: Path, q: int, p: int, Q: int):
    cmd = [
        "fastq_quality_filter",
        "-q", str(q),
        "-p", str(p),
        "-Q", str(Q),
        "-i", str(in_path),
        "-o", str(out_path),
    ]
    subprocess.run(cmd, check=True)

def extract_barcode(name: str):
    """
    Extract (F#, R#) from filenames like:
    Mbiome_16s__FILT_iBact_F1iBact_R8_SP.fastq
    """
    m = BARCODE_RE.search(name)
    if not m:
        return None
    fnum, rnum = m.group(1), m.group(2)
    return f"F{fnum}", f"R{rnum}"

def process_one_file(fpath_str: str, q: int, p: int, Q: int, overwrite: bool, out_dir_str: str):
    in_path = Path(fpath_str)
    out_dir = Path(out_dir_str)

    # Safety: ignore empty files (size 0) immediately
    try:
        if in_path.stat().st_size == 0:
            return {"status": "skip", "input": in_path.name, "reason": "empty_file_size0"}
    except FileNotFoundError:
        return {"status": "skip", "input": in_path.name, "reason": "missing_file"}

    bc = extract_barcode(in_path.name)
    if not bc:
        return {"status": "skip", "input": in_path.name, "reason": "barcode_not_found"}

    F, R = bc
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{F}_{R}_cleaned.fastq"

    # Count reads; if 0 reads, skip (handles "empty but not size0" edge cases)
    in_reads = count_reads_fastq(in_path)
    if in_reads == 0:
        return {"status": "skip", "input": in_path.name, "reason": "empty_file_0reads"}

    # If output exists and not overwriting, compute stats and skip running tool
    if out_path.exists() and (not overwrite):
        out_reads = count_reads_fastq(out_path)
        removed = in_reads - out_reads
        prop_removed = (removed / in_reads) if in_reads else 0.0
        return {
            "status": "ok",
            "input": in_path.name,
            "output": str(out_path),
            "raw_reads": in_reads,
            "cleaned_reads": out_reads,
            "removed_reads": removed,
            "prop_removed": prop_removed,
            "note": "output_exists_skipped_run",
        }

    # Run
    run_quality_filter(in_path, out_path, q, p, Q)
    out_reads = count_reads_fastq(out_path)

    removed = in_reads - out_reads
    prop_removed = (removed / in_reads) if in_reads else 0.0

    return {
        "status": "ok",
        "input": in_path.name,
        "output": str(out_path),
        "raw_reads": in_reads,
        "cleaned_reads": out_reads,
        "removed_reads": removed,
        "prop_removed": prop_removed,
    }

def main():
    ap = argparse.ArgumentParser(
        description="Batch fastq_quality_filter on Mbiome*SP.fastq -> cleaned_seq2/Fx_Ry_cleaned.fastq + removal stats"
    )
    ap.add_argument("--q", type=int, default=10, help="Quality threshold (-q)")
    ap.add_argument("--p", type=int, default=80, help="Percent of bases that must have Q>=q (-p)")
    ap.add_argument("--Q", type=int, default=33, help="Phred offset (-Q)")
    ap.add_argument("--jobs", type=int, default=4, help="Parallel workers")
    ap.add_argument("--glob", default="Mbiome*SP.fastq", help="Input pattern (default: Mbiome*SP.fastq)")
    ap.add_argument("--outdir", default="cleaned_seq2", help="Output directory (default: cleaned_seq2)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    args = ap.parse_args()

    files = sorted([p for p in Path(".").glob(args.glob) if p.is_file()])
    if not files:
        raise SystemExit(f"Aucun fichier trouvé avec le pattern '{args.glob}' dans le dossier courant.")

    results = []
    skipped = 0

    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = [
            ex.submit(
                process_one_file,
                str(p),
                args.q,
                args.p,
                args.Q,
                args.overwrite,
                args.outdir,
            )
            for p in files
        ]
        for fut in as_completed(futs):
            r = fut.result()
            if r["status"] != "ok":
                skipped += 1
                print(f"[WARN] skip {r['input']} ({r.get('reason','')})")
                continue

            note = r.get("note", "")
            extra = f" [{note}]" if note else ""
            print(
                f"[INFO] {r['input']} -> {r['output']} | removed {r['removed_reads']}/{r['raw_reads']} "
                f"({r['prop_removed']:.4f}){extra}"
            )
            results.append(r)

    results.sort(key=lambda x: x["input"])

    out_tsv = Path(args.outdir) / "quality_filter_stats.tsv"
    with out_tsv.open("w", encoding="utf-8") as f:
        f.write("input_file\toutput_file\traw_reads\tcleaned_reads\tremoved_reads\tprop_removed\n")
        for r in results:
            f.write(
                f"{r['input']}\t{r['output']}\t{r['raw_reads']}\t{r['cleaned_reads']}\t"
                f"{r['removed_reads']}\t{r['prop_removed']:.6f}\n"
            )

    print(f"[DONE] Stats: {out_tsv.resolve()}")
    if skipped:
        print(f"[DONE] Fichiers ignorés: {skipped}")

if __name__ == "__main__":
    main()
