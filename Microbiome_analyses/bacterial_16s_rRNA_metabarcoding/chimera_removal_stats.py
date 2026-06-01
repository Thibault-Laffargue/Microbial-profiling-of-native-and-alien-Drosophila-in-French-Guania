#!/usr/bin/env python3
import sys
import re

RE_OTU = re.compile(r"^C\d+$")

def otu_ids_from_fasta(fasta_path):
    otus = set()
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                first = line[1:].strip().split()[0]
                otus.add(first)  # chez toi c'est bien "C123"
    return otus

def read_otu_table_counts(tsv_path):
    """
    Lit une OTU table classique: 1 OTU par ligne.
    Split sur whitespace (tabs OU espaces).
    Retourne:
      counts: dict {OTU: total_reads}
      table_total: somme de toute la table (toutes OTUs)
      diagnostics: dict
    """
    counts = {}
    valid_lines = 0
    skipped_blank = 0
    skipped_bad_otu = 0
    skipped_too_short = 0

    table_total = 0

    with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        if not header:
            raise RuntimeError("OTU table vide")

        for line in f:
            if not line.strip():
                skipped_blank += 1
                continue

            parts = line.strip().split()   # tabs ou espaces
            if len(parts) < 2:
                skipped_too_short += 1
                continue

            otu = parts[0]
            if not RE_OTU.match(otu):
                skipped_bad_otu += 1
                continue

            total = 0
            for x in parts[1:]:
                try:
                    total += int(float(x))
                except ValueError:
                    # si une valeur est bizarre, on la traite comme 0
                    pass

            counts[otu] = total
            table_total += total
            valid_lines += 1

    diagnostics = {
        "valid_lines": valid_lines,
        "skipped_blank": skipped_blank,
        "skipped_bad_otu": skipped_bad_otu,
        "skipped_too_short": skipped_too_short,
    }
    return counts, table_total, diagnostics

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 chimera_removal_stats.py centroids.fasta centroids_without_chim_frogs.fasta 16s_OTU_table.tsv")
        sys.exit(1)

    centroids_fa, nochim_fa, otu_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

    all_otus = otu_ids_from_fasta(centroids_fa)
    kept_otus = otu_ids_from_fasta(nochim_fa)
    removed_otus = all_otus - kept_otus

    otu_counts, table_total, diag = read_otu_table_counts(otu_tsv)

    # Totaux calculés depuis la table (devrait être ton "nombre total de reads")
    total_reads_from_table = table_total

    # Reads retirés = somme des OTUs supprimées (si elles sont dans la table)
    removed_reads = sum(otu_counts.get(otu, 0) for otu in removed_otus)

    prop = removed_reads / total_reads_from_table if total_reads_from_table > 0 else 0.0

    # Diagnostics de correspondance FASTA vs table
    missing_from_table = [otu for otu in all_otus if otu not in otu_counts]

    print("### Résumé chimères (FROGS)")
    print(f"OTUs total (centroids.fasta)           : {len(all_otus)}")
    print(f"OTUs gardées (non chimériques)         : {len(kept_otus)}")
    print(f"OTUs retirées (chimériques)            : {len(removed_otus)}")

    print("\n### Impact en reads (depuis 16s_OTU_table.tsv)")
    print(f"Reads totaux (somme de toute la table) : {total_reads_from_table}")
    print(f"Reads retirés (OTUs chimériques)       : {removed_reads}")
    print(f"Proportion retirée                     : {prop:.6f} ({prop*100:.2f}%)")

    print("\n### Diagnostic lecture OTU table")
    print(f"Lignes OTU valides lues                : {diag['valid_lines']}")
    print(f"Lignes vides ignorées                  : {diag['skipped_blank']}")
    print(f"Lignes OTU invalides ignorées          : {diag['skipped_bad_otu']}")
    print(f"Lignes trop courtes ignorées           : {diag['skipped_too_short']}")
    print(f"OTUs dans la table (dict)              : {len(otu_counts)}")

    print("\n### Diagnostic correspondance FASTA -> table")
    print(f"OTUs de centroids absentes de la table : {len(missing_from_table)}")
    if missing_from_table:
        print("Exemples:", ", ".join(missing_from_table[:10]))

if __name__ == "__main__":
    main()
