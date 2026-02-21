#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def sample_name_from_path(path: Path) -> str:
    name = path.name
    if name.endswith('.exon.txt'):
        return name[:-len('.exon.txt')]
    return path.stem


def read_verse_counts(file_path: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    with file_path.open('r', encoding='utf-8', errors='replace') as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue

            gene_id, count_text = parts[0], parts[1]
            if gene_id.startswith('__'):
                continue

            try:
                count_value = int(float(count_text))
            except ValueError:
                continue

            counts[gene_id] = counts.get(gene_id, 0) + count_value

    return counts


def merge_counts(files: list[Path]) -> tuple[list[str], dict[str, dict[str, int]]]:
    sample_order: list[str] = []
    merged: dict[str, dict[str, int]] = {}

    for file_path in files:
        sample_name = sample_name_from_path(file_path)
        sample_order.append(sample_name)
        sample_counts = read_verse_counts(file_path)

        for gene_id, value in sample_counts.items():
            if gene_id not in merged:
                merged[gene_id] = {}
            merged[gene_id][sample_name] = value

    return sample_order, merged


def write_matrix(out_path: Path, sample_order: list[str], merged: dict[str, dict[str, int]]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open('w', newline='', encoding='utf-8') as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(['gene_id', *sample_order])

        for gene_id in sorted(merged.keys()):
            row = [gene_id]
            for sample_name in sample_order:
                row.append(merged[gene_id].get(sample_name, 0))
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description='Merge VERSE .exon.txt files into a raw count matrix CSV.')
    parser.add_argument('--out', required=True, help='Output CSV file path')
    parser.add_argument('count_files', nargs='+', help='Input VERSE count files (*.exon.txt)')
    args = parser.parse_args()

    count_paths = [Path(path) for path in args.count_files]
    sample_order, merged = merge_counts(count_paths)

    out_path = Path(args.out)
    write_matrix(out_path, sample_order, merged)


if __name__ == '__main__':
    main()
