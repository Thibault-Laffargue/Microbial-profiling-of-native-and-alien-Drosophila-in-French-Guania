import re, csv
from pathlib import Path

tsv_in  = Path("OL4_tax_table.tsv")
fasta   = Path("db/UNITE_1seqBYtax.fasta")
tsv_out = Path("OL4_tax_table.tsv")  # ⚠️ écrase l'entrée

def norm_taxon(s: str) -> str:
    s = (s or "").strip()
    s = s.replace("_", " ")
    s = re.sub(r"\s+", " ", s)
    return s.lower()

UNKNOWN_PAT = re.compile(r"^(na|n/a|unknown|unclassified|uncultured|unidentified|none|-)$", re.I)

# placeholder de type "... sp", "... sp.", "... sp. 1", "... sp1", "... sp.1"
SP_PAT = re.compile(r".*\bsp\.?\s*\d*\s*$", re.I)

def is_usable(name: str) -> bool:
    n = (name or "").strip()
    if not n:
        return False
    nn = norm_taxon(n)
    if not nn or UNKNOWN_PAT.match(nn):
        return False
    if SP_PAT.match(n):          # on considère "sp." comme placeholder => on remonte
        return False
    return True

# Rangs (FASTA tag -> nom de colonne TSV)
RANKS = [
    ("s", "Species"),
    ("g", "Genus"),
    ("f", "Family"),
    ("o", "Order"),
    ("c", "Class"),
    ("p", "Phylum"),
]

# Un dict par rang : maps["Species"][taxon_norm] = kingdom
maps = {col: {} for _, col in RANKS}

with fasta.open("r", encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        h = line[1:].strip()

        km = re.search(r"(?:^|[;| ])k__?([^;| ]+)", h)
        if not km:
            continue
        kingdom = km.group(1).strip()
        if not kingdom:
            continue

        for tag, col in RANKS:
            m = re.search(rf"(?:^|[;| ]){tag}__?([^;| ]+)", h)
            if not m:
                continue
            taxon_raw = m.group(1).strip()
            taxon_norm = norm_taxon(taxon_raw)
            if taxon_norm and not UNKNOWN_PAT.match(taxon_norm) and not SP_PAT.match(taxon_raw):
                maps[col].setdefault(taxon_norm, kingdom)

with tsv_in.open("r", encoding="utf-8", errors="ignore", newline="") as fin:
    reader = csv.DictReader(fin, delimiter="\t")
    fields = reader.fieldnames or []

    # On utilise les noms exacts que tu as donnés
    required = ["Cluster", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    missing = [c for c in required if c not in fields]
    if missing:
        raise ValueError(f"Colonnes manquantes dans le TSV: {missing}")

    rows = []
    filled = 0

    for row in reader:
        k = (row.get("Kingdom") or "").strip()
        if k.upper() in {"", "NA"}:
            # remonte depuis Species -> ... -> Phylum
            newk = None
            for _, col in RANKS:
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
