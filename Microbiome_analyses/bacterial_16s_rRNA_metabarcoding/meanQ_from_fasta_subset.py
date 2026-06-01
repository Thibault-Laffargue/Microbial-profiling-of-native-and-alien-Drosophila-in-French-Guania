#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse
import math
import statistics

def iter_fasta_headers(fasta_path: Path):
    with fasta_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                yield line[1:].rstrip("\n")

def sample_from_fasta_header(header: str) -> str:
    sample = header.split("_", 1)[0]
    # Special case: sample "DNA_RNA_Shield" got truncated to "DNA" in FASTA
    if sample == "DNA":
        return "DNA_RNA_Shield"
    return sample

def key_from_fasta_header(header: str) -> str:
    """
    FASTA: D6310_<FASTQ_FIRST_TOKEN> rest...
    We use <FASTQ_FIRST_TOKEN> (before first space) as matching key.
    """
    # Special case: FASTA may contain sample prefix "DNA_RNA_Shield_"
    if header.startswith("DNA_RNA_Shield_"):
        after = header[len("DNA_RNA_Shield_"):]
    else:
        after = header.split("_", 1)[1] if "_" in header else header
    return after.split(" ", 1)[0]

def key_from_fastq_header_line(line: str) -> str:
    """
    FASTQ header line starts with '@'.
    Key = first token after '@' (before first space).
    """
    s = line[1:].rstrip("\n")
    return s.split(" ", 1)[0]

def mean_phred_from_qual(qual: str, phred_offset: int) -> float:
    if not qual:
        return float("nan")
    total = 0
    n = len(qual.rstrip("\n"))
    for c in qual.rstrip("\n"):
        total += (ord(c) - phred_offset)
    return total / n if n else float("nan")

def expected_errors_from_qual(qual: str, phred_offset: int) -> float:
    """
    Expected errors (EE) = sum_i 10^(-Qi/10)
    where Qi are per-base Phred scores.
    """
    if not qual:
        return float("nan")
    ee = 0.0
    qstr = qual.rstrip("\n")
    for c in qstr:
        q = ord(c) - phred_offset
        ee += 10 ** (-q / 10.0)
    return ee

def quantile_type7(sorted_vals, q: float) -> float:
    n = len(sorted_vals)
    if n == 0:
        return float("nan")
    if n == 1:
        return float(sorted_vals[0])
    h = (n - 1) * q
    lo = int(math.floor(h))
    hi = int(math.ceil(h))
    if lo == hi:
        return float(sorted_vals[lo])
    frac = h - lo
    return float(sorted_vals[lo] + frac * (sorted_vals[hi] - sorted_vals[lo]))

def scan_fastq_for_keys(fastq_path: Path, needed_keys: set[str], phred_offset: int):
    """
    Read FASTQ in strict 4-line records. Safe even if '@' appears in quality lines.
    Returns dict key -> (meanQ, expected_errors, length) for keys found.
    """
    found = {}
    if not fastq_path.exists():
        return found

    with fastq_path.open("r", encoding="utf-8", errors="replace") as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                break

            if not h.startswith("@"):
                # skip malformed record
                continue

            key = key_from_fastq_header_line(h)
            if key in needed_keys and key not in found:
                qstr = qual.rstrip("\n")
                meanq = mean_phred_from_qual(qual, phred_offset)
                ee = expected_errors_from_qual(qual, phred_offset)
                found[key] = (meanq, ee, len(qstr))

    return found

def main():
    ap = argparse.ArgumentParser(description="Compute per-read meanQ and expected errors for reads in FASTA using per-sample FASTQs.")
    ap.add_argument("--fasta", default="all_samples.min1000.fas")
    ap.add_argument("--fastq_dir", default=".")
    ap.add_argument("--suffix", default="_cleaned.fastq")
    ap.add_argument("--phred_offset", type=int, default=33)
    ap.add_argument("--out_reads_tsv", default="read_quality.tsv")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    fastq_dir = Path(args.fastq_dir)

    # Build needed keys per sample, keep mapping for output
    needed_by_sample = defaultdict(set)
    records = []  # (sample, fasta_header, key)

    for hdr in iter_fasta_headers(fasta_path):
        sample = sample_from_fasta_header(hdr)
        key = key_from_fasta_header(hdr)
        needed_by_sample[sample].add(key)
        records.append((sample, hdr, key))

    if not records:
        raise SystemExit("FASTA vide ou sans headers.")

    # Scan each sample FASTQ once
    info_by_key = {}
    missing_fastq = 0
    for sample, keys in needed_by_sample.items():
        fq = fastq_dir / f"{sample}{args.suffix}"
        if not fq.exists():
            missing_fastq += 1
            continue
        info_by_key.update(scan_fastq_for_keys(fq, keys, args.phred_offset))

    # Write per-read TSV + stats values (EE-based)
    out_tsv = Path(args.out_reads_tsv)
    ee_values = []
    ee_rate_values = []
    q_values = []
    missing_reads = 0

    with out_tsv.open("w", encoding="utf-8") as w:
        w.write("sample\tkey\tmeanQ\texpected_errors\terror_rate\tlength\tfasta_header\n")
        for sample, hdr, key in records:
            info = info_by_key.get(key)
            if info is None:
                missing_reads += 1
                continue
            meanq, ee, length = info
            erate = (ee / length) if length else float("nan")
            q_values.append(meanq)
            ee_values.append(ee)
            ee_rate_values.append(erate)
            w.write(f"{sample}\t{key}\t{meanq:.4f}\t{ee:.6f}\t{erate:.8f}\t{length}\t{hdr}\n")

    if not ee_values:
        raise SystemExit(
            "Aucune qualité retrouvée. Vérifie:\n"
            "- que les FASTQ sont bien nommés SAMPLENAME_cleaned.fastq\n"
            "- que le FASTA contient bien SAMPLENAME_<fastq_first_token> ...\n"
            "- et que --fastq_dir pointe au bon dossier"
        )

    # Stats on expected_errors and error_rate
    ee_values.sort()
    ee_rate_values.sort()
    q_values.sort()

    total_reads = len(records)
    found_reads = len(ee_values)

    def summarize(vals, label):
        mean_val = sum(vals) / len(vals)
        median_val = statistics.median(vals)
        q1 = quantile_type7(vals, 0.25)
        q2 = quantile_type7(vals, 0.50)
        q3 = quantile_type7(vals, 0.75)
        print(label)
        print(f"  Mean:   {mean_val:.6f}")
        print(f"  Median: {median_val:.6f}")
        print(f"  Q1:     {q1:.6f}")
        print(f"  Q2:     {q2:.6f}")
        print(f"  Q3:     {q3:.6f}")
        print(f"  Min:    {vals[0]:.6f}")
        print(f"  Max:    {vals[-1]:.6f}\n")
        return mean_val, median_val, q1, q2, q3

    print(f"Reads dans FASTA: {total_reads}")
    print(f"Reads retrouvés dans FASTQ: {found_reads}")
    print(f"Reads manquants (non retrouvés): {missing_reads}")
    print(f"FASTQ manquants (samples sans fichier): {missing_fastq}\n")

    summarize(ee_values, "Expected errors per read (EE)")
    summarize(ee_rate_values, "Error rate per base (EE/length)")

    # ---- Added summary TSV (no other changes) ----
    L = 1350.0

    def p_from_q(q):
        return 10 ** (-q / 10.0)

    def expected_identity_percent(p):
        # I ≈ (1-p)^2 + p^2/3
        I = (1.0 - p) ** 2 + (p ** 2) / 3.0
        return 100.0 * I

    # Compute the five summary points for Quality and Error_rate
    q_mean, q_median, q1_q, q2_q, q3_q = (
        sum(q_values) / len(q_values),
        statistics.median(q_values),
        quantile_type7(q_values, 0.25),
        quantile_type7(q_values, 0.50),
        quantile_type7(q_values, 0.75),
    )

    er_mean, er_median, er1, er2, er3 = (
        sum(ee_rate_values) / len(ee_rate_values),
        statistics.median(ee_rate_values),
        quantile_type7(ee_rate_values, 0.25),
        quantile_type7(ee_rate_values, 0.50),
        quantile_type7(ee_rate_values, 0.75),
    )

    rows = [
        ("Mean",   q_mean,   er_mean),
        ("Median", q_median, er_median),
        ("Q1",     q1_q,     er1),
        ("Q2",     q2_q,     er2),
        ("Q3",     q3_q,     er3),
    ]

    summary_path = Path("summary_stats.tsv")
    with summary_path.open("w", encoding="utf-8") as s:
        s.write("stat\tQuality\tExpected_nt_dif_qual\tExpected_id_qual\tError_rate\tExpected_nt_dif_ER\tExpected_id_ER\n")
        for name, qv, erv in rows:
            p_q = p_from_q(qv)
            exp_nt_q = L * p_q
            exp_id_q = expected_identity_percent(p_q)

            p_er = erv
            exp_nt_er = L * p_er
            exp_id_er = expected_identity_percent(p_er)

            s.write(
                f"{name}\t"
                f"{qv:.6f}\t"
                f"{exp_nt_q:.6f}\t"
                f"{exp_id_q:.6f}\t"
                f"{erv:.6f}\t"
                f"{exp_nt_er:.6f}\t"
                f"{exp_id_er:.6f}\n"
            )
    # ---- End added summary TSV ----

    print(f"[DONE] TSV per-read: {out_tsv.resolve()}")

if __name__ == "__main__":
    main()
