#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path


def parse_attributes(attr_text: str) -> dict:
    attrs = {}
    for field in attr_text.strip().split(';'):
        field = field.strip()
        if not field:
            continue
        if ' ' not in field:
            continue
        key, value = field.split(' ', 1)
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def open_text_maybe_gzip(path: Path):
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return path.open('r', encoding='utf-8', errors='replace')


def build_mapping(gtf_path: Path) -> list[tuple[str, str]]:
    mapping = {}
    with open_text_maybe_gzip(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type != 'gene':
                continue

            attrs = parse_attributes(parts[8])
            gene_id = attrs.get('gene_id')
            gene_name = attrs.get('gene_name') or attrs.get('gene_symbol') or gene_id
            if not gene_id:
                continue
            if gene_id not in mapping:
                mapping[gene_id] = gene_name

    return sorted(mapping.items(), key=lambda x: x[0])


def write_csv(rows: list[tuple[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open('w', newline='', encoding='utf-8') as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(['gene_id', 'gene_symbol'])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description='Parse GTF and output gene_id to gene_symbol mapping CSV.')
    parser.add_argument('--gtf', required=True, help='Input GTF file path (.gtf or .gtf.gz)')
    parser.add_argument('--out', required=True, help='Output CSV path')
    args = parser.parse_args()

    gtf_path = Path(args.gtf)
    out_path = Path(args.out)

    rows = build_mapping(gtf_path)
    write_csv(rows, out_path)


if __name__ == '__main__':
    main()
