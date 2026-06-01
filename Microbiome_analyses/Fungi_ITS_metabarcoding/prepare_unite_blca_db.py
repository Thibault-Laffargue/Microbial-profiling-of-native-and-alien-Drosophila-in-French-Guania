#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path
from typing import Dict, Optional, Tuple, List

# UNITE/QIIME style levels: k__Fungi; p__Ascomycota; ...
PREFIX2RANK = {
    "d": "kingdom",   # sometimes d__ used for domain; map to kingdom
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

TAX_START_RE = re.compile(r"(?:^|[|])([dkpcofgs]__.*)$", flags=re.IGNORECASE)

def short_id_from_original(original_id: str, prefix: str, n_hex: int) -> str:
    # stable, deterministic short id
    h = hashlib.sha1(original_id.encode("utf-8", errors="ignore")).hexdigest()[:n_hex]
    return f"{prefix}{h}"

def extract_unite_tax_from_header(header: str) -> Optional[str]:
    """
    UNITE headers like:
      Abrothallus_subhalei|MT153946|SH...|refs|k__Fungi;p__Ascomycota;...
    There is no whitespace; taxonomy starts at k__/d__/...
    We extract substring starting at the first occurrence of [dpkcofgs]__.
    """
    m = TAX_START_RE.search(header)
    if not m:
        return None
    tax = m.group(1).strip()
    return tax if tax else None

def normalize_token(token: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Convert a taxonomy token to (rank, taxon).

    Accepts:
      - k__Fungi / p__Ascomycota / ...
      - kingdom:Fungi
      - bare taxon (rank None; not expected for UNITE but kept as fallback)
    """
    t = token.strip()
    if not t:
        return None, None

    # already rank:taxon
    if ":" in t and not re.match(r"^[dpkcofgs]__", t, flags=re.IGNORECASE):
        r, tax = t.split(":", 1)
        r = r.strip().lower()
        tax = tax.strip()
        return (r if r else None), (tax if tax else None)

    # prefix style like k__Fungi
    m = re.match(r"^([dpkcofgs])__\s*(.*)$", t, flags=re.IGNORECASE)
    if m:
        pre = m.group(1).lower()
        tax = m.group(2).strip()
        rank = PREFIX2RANK.get(pre)
        return (rank if rank else None), (tax if tax else None)

    # bare
    return None, t

def build_blca_tax_string(raw_tax: str) -> str:
    """
    Convert UNITE raw taxonomy string (k__/p__/...) into BLCA format:
      kingdom:Fungi;phylum:Ascomycota;...
    """
    parts = [p.strip() for p in raw_tax.split(";") if p.strip()]

    explicit: Dict[str, str] = {}
    bare: List[str] = []

    for p in parts:
        r, tax = normalize_token(p)
        if tax is None or tax == "":
            continue
        if r is None:
            bare.append(tax)
        else:
            explicit[r] = tax

    out_items: List[str] = []
    if explicit:
        for r in RANK_ORDER:
            if r in explicit and explicit[r]:
                out_items.append(f"{r}:{explicit[r]}")
    else:
        # fallback: assume bare tokens follow standard order
        for r, tax in zip(RANK_ORDER, bare):
            out_items.append(f"{r}:{tax}")

    return ";".join(out_items)

def main():
    ap = argparse.ArgumentParser(
        description="Prepare UNITE FASTA for BLAST/BLCA: short-id FASTA + map + BLCA taxonomy extracted from UNITE headers."
    )
    ap.add_argument("-i", "--in", dest="inp", required=True, help="Input UNITE FASTA (original headers)")
    ap.add_argument("--out-fasta", default="db/UNITE_1seqBYtax.shortid.fasta", help="Output FASTA with short IDs")
    ap.add_argument("--out-map", default="db/UNITE_1seqBYtax.shortid.map.tsv", help="Output mapping TSV")
    ap.add_argument("--out-tax", default="db/UNITE_1seqBYtax.shortid.BLCA.taxonomy", help="Output BLCA taxonomy TSV")
    ap.add_argument("--prefix", default="U", help="Prefix for short IDs (default: U)")
    ap.add_argument("--hex", type=int, default=16, help="Hex length from SHA1 (default: 16)")
    args = ap.parse_args()

    inp = Path(args.inp)
    out_fa = Path(args.out_fasta)
    out_map = Path(args.out_map)
    out_tax = Path(args.out_tax)

    if not inp.exists():
        raise SystemExit(f"Input FASTA not found: {inp}")

    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_map.parent.mkdir(parents=True, exist_ok=True)
    out_tax.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    n_headers = 0
    n_tax_written = 0
    n_skipped_no_tax = 0
    n_collisions = 0

    with inp.open("r", encoding="utf-8", errors="replace") as fin, \
         out_fa.open("w", encoding="utf-8") as fout, \
         out_map.open("w", encoding="utf-8") as fmap, \
         out_tax.open("w", encoding="utf-8") as ftax:

        fmap.write("short_id\toriginal_id\toriginal_header\n")

        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                n_headers += 1

                # original_id = whole header up to whitespace (in UNITE there is none)
                original_id = header.split()[0]

                sid = short_id_from_original(original_id, args.prefix, args.hex)

                # collision guard
                base = sid
                k = 1
                while sid in seen:
                    k += 1
                    sid = f"{base}_{k}"
                    n_collisions += 1
                seen.add(sid)

                # write FASTA header (short id only)
                fout.write(f">{sid}\n")

                # map
                fmap.write(f"{sid}\t{original_id}\t{header}\n")

                # taxonomy extraction from UNITE header
                raw_tax = extract_unite_tax_from_header(header)
                if raw_tax is None:
                    n_skipped_no_tax += 1
                else:
                    tax_str = build_blca_tax_string(raw_tax)
                    if tax_str:
                        ftax.write(f"{sid}\t{tax_str}\n")
                        n_tax_written += 1
                    else:
                        n_skipped_no_tax += 1
            else:
                fout.write(line)

    print(f"[DONE] headers processed: {n_headers}")
    print(f"[DONE] collisions resolved: {n_collisions}")
    print(f"[DONE] wrote shortid FASTA: {out_fa.resolve()}")
    print(f"[DONE] wrote mapping TSV:  {out_map.resolve()}")
    print(f"[DONE] wrote BLCA taxonomy:{out_tax.resolve()}")
    print(f"[INFO] skipped headers without detectable taxonomy: {n_skipped_no_tax}")

if __name__ == "__main__":
    main()
