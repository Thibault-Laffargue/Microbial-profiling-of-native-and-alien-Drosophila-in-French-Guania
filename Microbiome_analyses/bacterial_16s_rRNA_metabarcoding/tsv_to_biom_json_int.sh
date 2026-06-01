#!/usr/bin/env python3
import json
import sys

def to_int(x: str) -> int:
    x = x.strip()
    if x == "" or x.lower() == "nan":
        return 0
    # accepte "10683", "10683.0" etc.
    return int(float(x))

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} 16s_OTU_table.tsv 16s_OTU_table.json.biom", file=sys.stderr)
        return 1

    tsv_in, biom_out = sys.argv[1], sys.argv[2]

    with open(tsv_in, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        if len(header) < 2:
            raise ValueError("TSV must have OTU column + at least one sample.")
        samples = header[1:]

        otu_ids = []
        data = []

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            otu = parts[0]
            counts = parts[1:]
            if len(counts) < len(samples):
                counts += ["0"] * (len(samples) - len(counts))

            row_idx = len(otu_ids)
            otu_ids.append(otu)

            for col_idx, v in enumerate(counts[:len(samples)]):
                iv = to_int(v)
                if iv != 0:
                    data.append([row_idx, col_idx, iv])

    biom = {
        "id": None,
        "format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "tsv_to_biom_json_int.py",
        "date": None,
        "matrix_type": "sparse",
        "matrix_element_type": "int",
        "shape": [len(otu_ids), len(samples)],
        "data": data,
        "rows": [{"id": oid, "metadata": None} for oid in otu_ids],
        "columns": [{"id": s, "metadata": None} for s in samples],
    }

    with open(biom_out, "w", encoding="utf-8") as out:
        json.dump(biom, out, ensure_ascii=False)

    print(f"[OK] wrote {biom_out} shape={len(otu_ids)}x{len(samples)} nonzero={len(data)}", file=sys.stderr)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
