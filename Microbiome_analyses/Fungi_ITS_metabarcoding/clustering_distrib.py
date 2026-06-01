#!/usr/bin/env python3
from pathlib import Path
from collections import Counter, defaultdict
import argparse

def cluster_size_distribution(path: Path):
    # cluster_id -> nb de séquences
    cluster_counts = defaultdict(int)

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # ton fichier est tab-separated, cluster_id = premier champ
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue

            cid = parts[0].strip()

            # sécurité: cid doit être un int
            if not cid.isdigit():
                continue

            cluster_counts[int(cid)] += 1

    # maintenant: taille_cluster -> nombre_de_clusters
    dist = Counter(cluster_counts.values())
    return dist, cluster_counts

def main():
    ap = argparse.ArgumentParser(description="Compute cluster size distribution from cluster-labeled file")
    ap.add_argument("input", help="Input file (e.g. output.clstr)")
    ap.add_argument("--tsv", default="cluster_size_distribution.tsv",
                    help="Output TSV file (default: cluster_size_distribution.tsv)")
    args = ap.parse_args()

    dist, cluster_counts = cluster_size_distribution(Path(args.input))

    if not cluster_counts:
        raise SystemExit("Aucun cluster détecté: vérifie le format (tabulations, cluster id en colonne 1).")

    total_clusters = len(cluster_counts)
    total_seqs = sum(cluster_counts.values())

    print(f"Total clusters: {total_clusters}")
    print(f"Total sequences: {total_seqs}\n")

    print("cluster_size\tn_clusters\tpercent_clusters")
    for size in sorted(dist):
        n = dist[size]
        pct = 100 * n / total_clusters
        print(f"{size}\t{n}\t{pct:.2f}")

    # write TSV
    out = Path(args.tsv)
    with out.open("w", encoding="utf-8") as w:
        w.write("cluster_size\tn_clusters\n")
        for size in sorted(dist):
            w.write(f"{size}\t{dist[size]}\n")

    print(f"\n[DONE] TSV écrit: {out.resolve()}")

if __name__ == "__main__":
    main()
