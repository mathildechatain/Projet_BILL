#!/usr/bin/env python3
"""Construit une table récapitulative des SNP à partir des VCF filtrés,
avec annotation ORF.

Le script lit tous les VCF SNP filtrés, extrait les informations utiles
(position, ref, alt, qual, genotype), annote chaque variant par ORF à partir
d'un GFF, puis exporte un CSV détaillé variant par variant.
"""

import argparse
import csv
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SRC_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = SRC_ROOT.parent


def resolve_project_path(path_str):
    """Résout un chemin fourni par l'utilisateur relativement à la racine du projet."""
    candidate = Path(path_str)
    if candidate.is_absolute():
        return candidate.resolve()
    return (PROJECT_ROOT / candidate).resolve()


def parse_attributes(attr):
    """Parse la colonne attributs (9e colonne) d'un GFF en dictionnaire."""
    out = {}
    for item in attr.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def parse_gene_annotations(gff_path):
    """Extrait les features de type gene du GFF pour annoter les ORF."""
    genes = []

    with gff_path.open("r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attrs = cols

            if feature != "gene":
                continue

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            attr_map = parse_attributes(attrs)

            orf = (
                attr_map.get("Name")
                or attr_map.get("locus_tag")
                or attr_map.get("ID")
                or "Unknown_ORF"
            )

            genes.append(
                {
                    "chrom": chrom,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                    "orf": orf,
                }
            )

    return genes


def find_orfs(chrom, pos, genes):
    """Retourne l'ORF (ou les ORF) recouvrant une position génomique."""
    hits = []

    for g in genes:
        if g["chrom"] == chrom and g["start"] <= pos <= g["end"]:
            hits.append(g["orf"])

    if not hits:
        return "Intergenic"

    uniq = []
    seen = set()

    for h in hits:
        if h not in seen:
            seen.add(h)
            uniq.append(h)

    return ";".join(uniq)


def parse_sample_meta(filename):
    """Extrait génération et numéro de réplication depuis le nom de fichier."""
    # ex: P25-7.trimed1000.snp_SNP_filtered.vcf
    m = re.match(r"^(P\d+)-(\d+)\.", filename)

    if not m:
        return "", ""

    return m.group(1), m.group(2)


def main():
    """Point d'entrée: parsing VCF, annotation ORF, export CSV."""

    parser = argparse.ArgumentParser(
        description="Extrait les SNP des VCF filtrés, les annote par ORF à partir du GFF, et exporte un CSV."
    )

    parser.add_argument(
        "--vcf-dir",
        default="Files/vcf_filtered_snp",
        help="Répertoire contenant les fichiers VCF SNP filtrés",
    )

    parser.add_argument(
        "--gff",
        default="Files/DQ657948.1.gff3",
        help="Fichier d'annotation GFF",
    )

    parser.add_argument(
        "--out-dir",
        default="Files/plot_files",
        help="Répertoire de sortie pour la table CSV",
    )

    args = parser.parse_args()

    vcf_dir = resolve_project_path(args.vcf_dir)
    gff_path = resolve_project_path(args.gff)
    out_dir = resolve_project_path(args.out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    genes = parse_gene_annotations(gff_path)

    vcf_files = sorted(vcf_dir.glob("*.vcf"))

    raw_rows = []

    for vcf_file in vcf_files:

        generation, replicate = parse_sample_meta(vcf_file.name)

        sample_name = vcf_file.stem

        with vcf_file.open("r") as f:

            for line in f:

                if line.startswith("##"):
                    continue

                if line.startswith("#CHROM"):

                    cols = line.rstrip("\n").split("\t")

                    if len(cols) > 9:
                        sample_name = cols[-1]

                    continue

                if line.startswith("#"):
                    continue

                cols = line.rstrip("\n").split("\t")

                if len(cols) < 8:
                    continue

                chrom = cols[0]
                pos = int(cols[1])
                ref = cols[3]
                alt = cols[4]
                qual = cols[5]

                gt = ""
                gq = ""

                if len(cols) >= 10:

                    keys = cols[8].split(":")
                    vals = cols[9].split(":")

                    format_map = dict(zip(keys, vals))

                    gt = format_map.get("GT", "")
                    gq = format_map.get("GQ", "")

                orf = find_orfs(chrom, pos, genes)

                raw_rows.append(
                    {
                        "file": vcf_file.name,
                        "sample": sample_name,
                        "generation": generation,
                        "replicate": replicate,
                        "chrom": chrom,
                        "position": pos,
                        "ref": ref,
                        "alt": alt,
                        "qual": qual,
                        "gt": gt,
                        "gq": gq,
                        "orf": orf,
                    }
                )

    raw_out = out_dir / "summary_table_snp_with_ORF.csv"

    with raw_out.open("w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "file",
                "sample",
                "generation",
                "replicate",
                "chrom",
                "position",
                "ref",
                "alt",
                "qual",
                "gt",
                "gq",
                "orf",
            ],
        )

        writer.writeheader()
        writer.writerows(raw_rows)

    print(f"CSV written: {raw_out}")
    print(f"Total SNP parsed: {len(raw_rows)}")


if __name__ == "__main__":
    main()
