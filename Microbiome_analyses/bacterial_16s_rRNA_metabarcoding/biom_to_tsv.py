#!/usr/bin/env python3
import json
import sys


def main():
    if len(sys.argv) not in (2, 3):
        print(
            f"Usage:\n  {sys.argv[0]} input.biom [output.tsv]\n\n"
            f"Example:\n  {sys.argv[0]} 16s_OTU_table_without_chim.biom 16s_OTU_table_without_chim.tsv",
            file=sys.stderr,
        )
        return 1

    biom_in = sys.argv[1]
    tsv_out = sys.argv[2] if len(sys.argv) == 3 else (biom_in.rsplit(".", 1)[0] + ".tsv")

    with open(biom_in, "r", encoding="utf-8", errors="replace") as f:
        biom = json.load(f)

    # Récupère les IDs des OTUs et des samples
    rows = biom.get("rows", [])
    cols = biom.get("columns", [])
    otu_ids = [r.get("id", "") for r in rows]
    sample_ids = [c.get("id", "") for c in cols]

    if not otu_ids or not sample_ids:
        raise ValueError("BIOM file must contain 'rows' and 'columns' with 'id' fields.")

    n_rows, n_cols = biom.get("shape", [len(otu_ids), len(sample_ids)])
    if n_rows != len(otu_ids) or n_cols != len(sample_ids):
        # On tolère, mais on garde cohérent avec rows/columns
        n_rows, n_cols = len(otu_ids), len(sample_ids)

    matrix_type = biom.get("matrix_type", "sparse")
    data = biom.get("data")

    # Initialise la matrice de comptes (int)
    table = [[0] * n_cols for _ in range(n_rows)]

    if matrix_type == "sparse":
        # Format: [ [row_idx, col_idx, value], ... ]
        if not isinstance(data, list):
            raise ValueError("Sparse BIOM must have 'data' as a list of [row, col, value].")
        for triplet in data:
            if not (isinstance(triplet, list) and len(triplet) == 3):
                continue
            r, c, v = triplet
            table[int(r)][int(c)] = int(v)

    elif matrix_type == "dense":
        # Format: [ [v00, v01, ...], [v10, v11, ...], ... ]
        if not (isinstance(data, list) and len(data) == n_rows):
            raise ValueError("Dense BIOM must have 'data' as a 2D list with the same number of rows as 'shape'.")
        for r in range(n_rows):
            row_vals = data[r]
            if not isinstance(row_vals, list):
                raise ValueError("Dense BIOM rows must be lists.")
            for c in range(min(n_cols, len(row_vals))):
                table[r][c] = int(float(row_vals[c]))  # tolère 0.0 etc.
    else:
        raise ValueError(f"Unsupported matrix_type: {matrix_type} (expected 'sparse' or 'dense')")

    # Écrit le TSV
    with open(tsv_out, "w", encoding="utf-8") as out:
        out.write("\t".join(["OTU"] + sample_ids) + "\n")
        for r, otu in enumerate(otu_ids):
            out.write(otu)
            for c in range(n_cols):
                out.write("\t" + str(table[r][c]))
            out.write("\n")

    print(f"[OK] wrote {tsv_out} shape={n_rows}x{n_cols}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
