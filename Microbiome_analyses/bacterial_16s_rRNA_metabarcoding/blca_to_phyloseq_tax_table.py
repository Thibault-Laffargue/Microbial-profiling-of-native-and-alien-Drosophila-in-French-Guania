#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional, List

# Rangs dans ton output BLCA (ordre du plus haut au plus bas)
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

# Noms de colonnes compatibles phyloseq (tu peux renommer Kingdom -> Superkingdom si tu préfères)
COLNAMES = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

def parse_blca_pairs(tax_part: str) -> Dict[str, Tuple[str, Optional[float]]]:
    """
    Parse format:
      rank:taxon;conf;rank:taxon;conf;...
    Returns dict rank -> (taxon, conf_float)
    """
    items = [x.strip() for x in tax_part.strip().split(";")]
    # enlève les vides (à cause du ; final)
    items = [x for x in items if x != ""]

    out: Dict[str, Tuple[str, Optional[float]]] = {}

    i = 0
    while i < len(items):
        token = items[i]

        # attend "rank:taxon"
        if ":" not in token:
            i += 1
            continue

        rank, tax = token.split(":", 1)
        rank = rank.strip().lower()
        tax = tax.strip()

        conf: Optional[float] = None
        if i + 1 < len(items):
            # l'élément suivant est normalement la confiance
            try:
                conf = float(items[i + 1])
                i += 2
            except ValueError:
                # pas une confiance -> on avance d'un
                i += 1
        else:
            i += 1

        # garde seulement les ranks attendus
        if rank in RANKS:
            out[rank] = (tax, conf)

    return out

def apply_confidence_rule(
    rank2tc: Dict[str, Tuple[str, Optional[float]]],
    threshold_0_1: float = 0.8,
) -> List[str]:
    """
    Apply rule:
      - confidences in file are 0..100 -> threshold = threshold_0_1*100
      - walk ranks from top to bottom
      - if rank missing OR conf missing OR conf < threshold: stop
          fill this rank and all below with "<last_ok_tax> sp."
      - if stop happens at top rank: all Unknown
    """
    thr = threshold_0_1 * 100.0

    values = ["Unknown"] * len(RANKS)
    last_ok_tax: Optional[str] = None

    for idx, rank in enumerate(RANKS):
        tax, conf = rank2tc.get(rank, ("", None))

        # absent taxon -> stop
        if not tax:
            if idx == 0:
                return ["Unknown"] * len(RANKS)
            filler = f"{last_ok_tax} sp." if last_ok_tax else "Unknown"
            for j in range(idx, len(RANKS)):
                values[j] = filler
            return values

        # conf missing or too low -> stop
        if conf is None or conf < thr:
            if idx == 0:
                return ["Unknown"] * len(RANKS)
            filler = f"{last_ok_tax} sp." if last_ok_tax else "Unknown"
            for j in range(idx, len(RANKS)):
                values[j] = filler
            return values

        # OK
        values[idx] = tax
        last_ok_tax = tax

    return values

def main():
    ap = argparse.ArgumentParser(
        description="Convert BLCA output (rank:tax;conf;...) to phyloseq-ready tax_table with confidence cutoff."
    )
    ap.add_argument("--in", dest="inp", default="centroids.SILVA.blca.out", help="BLCA output file")
    ap.add_argument("--out", dest="out", default="16s_tax_table.tsv", help="Output phyloseq tax_table TSV")
    ap.add_argument("--threshold", type=float, default=0.8, help="Confidence threshold in [0,1] (default 0.8)")
    args = ap.parse_args()

    inp = Path(args.inp)
    outp = Path(args.out)
    if not inp.exists():
        raise SystemExit(f"Input not found: {inp}")

    with inp.open("r", encoding="utf-8", errors="replace") as fin, outp.open("w", encoding="utf-8") as fout:
        fout.write("Cluster\t" + "\t".join(COLNAMES) + "\n")

        n = 0
        for line in fin:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue

            # format: ClusterID<TAB>taxonomy_string
            if "\t" not in line:
                continue
            cluster, tax_part = line.split("\t", 1)
            cluster = cluster.strip()
            tax_part = tax_part.strip()

            rank2tc = parse_blca_pairs(tax_part)
            vals = apply_confidence_rule(rank2tc, threshold_0_1=args.threshold)

            fout.write(cluster + "\t" + "\t".join(vals) + "\n")
            n += 1

    print(f"[DONE] wrote {n} rows -> {outp.resolve()}")

if __name__ == "__main__":
    main()
