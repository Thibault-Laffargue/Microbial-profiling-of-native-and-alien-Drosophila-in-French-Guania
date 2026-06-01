import re, csv
from pathlib import Path

tsv_in  = Path("16s_tax_table.tsv")
fasta   = Path("db/SILVA_1seqBYtax.fasta")
tsv_out = Path("16s_tax_table.tsv")  # ⚠️ écrase l'entrée

def norm_taxon(s: str) -> str:
    s = (s or "").strip()
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s.lower()

UNKNOWN_PAT = re.compile(r"^(na|n/a|unknown|unclassified|uncultured|unidentified|none|-)$", re.I)
SP_PAT = re.compile(r".*\bsp\.?\s*\d*\s*$", re.I)

def is_usable(name: str) -> bool:
    n = (name or "").strip()
    if not n:
        return False
    nn = norm_taxon(n)
    if not nn or UNKNOWN_PAT.match(nn):
        return False
    if SP_PAT.match(n):  # placeholder => on remonte
        return False
    return True

# ordre de backoff (du plus bas au plus haut)
BACKOFF = ["Species", "Genus", "Family", "Order", "Class", "Phylum"]

# maps[col][taxon_norm] = kingdom
maps = {col: {} for col in BACKOFF}

def parse_silva_tax_from_header(header: str):
    """
    header: ligne après '>' (sans le '>'), ex:
      'AB000393.1.1510 Bacteria;Pseudomonadota;...;Vibrio halioticoli'
    On récupère la partie après le 1er espace, puis split ';'.
    """
    parts = header.split(None, 1)
    if len(parts) < 2:
        return None
    tax_str = parts[1].strip()
    if not tax_str:
        return None
    taxa = [t.strip() for t in tax_str.split(";") if t.strip()]
    if len(taxa) < 2:
        return None
    return taxa

with fasta.open("r", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        taxa = parse_silva_tax_from_header(h)
        if not taxa:
            continue

        kingdom = taxa[0].strip()
        if not kingdom or UNKNOWN_PAT.match(norm_taxon(kingdom)):
            continue

        # Cas standard Bacteria/Archaea: on suppose 7 niveaux
        # Kingdom;Phylum;Class;Order;Family;Genus;Species
        if kingdom in {"Bacteria", "Archaea"} and len(taxa) >= 7:
            phylum, clazz, order, family, genus, species = taxa[1:7]
        else:
            # fallback générique : on prend les 6 derniers pour Phylum..Species
            # (utile si la taxonomie contient des niveaux intermédiaires en plus)
            if len(taxa) < 7:
                continue
            phylum, clazz, order, family, genus, species = taxa[-6:]

        rank_vals = {
            "Phylum":  phylum,
            "Class":   clazz,
            "Order":   order,
            "Family":  family,
            "Genus":   genus,
            "Species": species,
        }

        for col, val in rank_vals.items():
            vnorm = norm_taxon(val)
            if vnorm and not UNKNOWN_PAT.match(vnorm) and not SP_PAT.match(val):
                maps[col].setdefault(vnorm, kingdom)

with tsv_in.open("r", encoding="utf-8", errors="ignore", newline="") as fin:
    reader = csv.DictReader(fin, delimiter="\t")
    fields = reader.fieldnames or []

    required = ["Cluster", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    missing = [c for c in required if c not in fields]
    if missing:
        raise ValueError(f"Colonnes manquantes dans le TSV: {missing}")

    rows = []
    filled = 0

    for row in reader:
        k = (row.get("Kingdom") or "").strip()
        if k.upper() in {"", "NA"}:
            for col in BACKOFF:
                val = row.get(col) or ""
                if is_usable(val):
                    newk = maps[col].get(norm_taxon(val))
                    if newk:
                        row["Kingdom"] = newk
                        filled += 1
                    break
        rows.append(row)

with tsv_out.open("w", encoding="utf-8", newline="") as fout:
    writer = csv.DictWriter(fout, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Kingdoms filled: {filled}")
print(f"Output written to: {tsv_out}")
