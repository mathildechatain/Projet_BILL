#!/usr/bin/env python3
"""Construit une table récapitulative des variants structuraux à partir des VCF filtrés,
avec annotation ORF.
Le script lit tous les VCF SV filtrés, extrait les informations utiles
(position, type, support, VAF), annote chaque variant par ORF à partir d'un GFF,
puis exporte un CSV détaillé variant par variant. """


import argparse
import csv
import re
from pathlib import Path


def parse_info_field(info_field):
    """Parse la colonne INFO d'un VCF en dictionnaire clé/valeur."""
    info = {}
    for item in info_field.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
    return info


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
            orf = attr_map.get("Name") or attr_map.get("locus_tag") or attr_map.get("ID", "Unknown_ORF")
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
    # Si plusieurs ORF recouvrent la position, on les concatène en une chaîne unique.
    uniq = []
    seen = set()
    for h in hits:
        if h not in seen:
            seen.add(h)
            uniq.append(h)
    return ";".join(uniq)


def safe_float(x):
    """Convertit en float, sinon retourne None."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def extract_vaf(info_map, format_map):
    """Calcule la VAF selon la formule suivante VAF = DV / (DR + DV) si les champs DR et DV sont présents, sinon utilise VAF ou AF."""
    dr = safe_float(format_map.get("DR"))
    dv = safe_float(format_map.get("DV"))
    if dr is not None and dv is not None and (dr + dv) > 0:
        # Cas principal Sniffles: support variant (DV) et support référence (DR)
        return dv / (dr + dv), dr, dv

    vaf_info = safe_float(info_map.get("VAF"))
    if vaf_info is not None:
        return vaf_info, dr, dv

    af = safe_float(format_map.get("AF"))
    if af is not None:
        return af, dr, dv

    ad = format_map.get("AD")
    if ad:
        parts = ad.split(",")
        if len(parts) >= 2:
            ref = safe_float(parts[0])
            alt = safe_float(parts[1])
            if ref is not None and alt is not None and (ref + alt) > 0:
                return alt / (ref + alt), dr, dv

    return None, dr, dv


def parse_sample_meta(filename):
    """Extrait génération et numéro d'échantillon depuis le nom de fichier."""
    # ex: P25-7.trimed1000.sv_sniffles_SV_filtered.vcf
    m = re.match(r"^(P\d+)-(\d+)\.", filename)
    if not m:
        return "", ""
    return m.group(1), m.group(2)


def main():
    """Point d'entrée: parsing VCF, annotation ORF, export CSV."""
    parser = argparse.ArgumentParser(
        description="Extrait les variants structuraux des VCF filtrés, les annote par ORF à partir du GFF, et exporte un CSV."
    )
    parser.add_argument(
        "--vcf-dir",
        default="Files/vcf_filtered_sv",
        help="Répertoire contenant les fichiers VCF filtrés pour les variants structuraux",
    )
    parser.add_argument(
        "--gff",
        default="Files/DQ657948.1.gff3",
        help="Fichier d'annotation GFF pour CyHV-3",
    )
    parser.add_argument(
        "--out-dir",
        default="Files/plot_files",
        help="Répertoire de sortie pour la table CSV intermédiaire",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent # Racine du projet
    vcf_dir = (repo_root / args.vcf_dir).resolve() # Répertoire des VCF filtrés
    gff_path = (repo_root / args.gff).resolve() # Chemin vers le fichier GFF d'annotation
    out_dir = (repo_root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    genes = parse_gene_annotations(gff_path)
    # Les fichiers attendus sont de la forme:
    # P25-7.trimed1000.sv_sniffles_SV_filtered.vcf
    vcf_files = sorted(vcf_dir.glob("*.vcf"))
    
    # Table brute: une ligne par variant observé dans un échantillon.
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
                info_map = parse_info_field(cols[7])
                svtype = info_map.get("SVTYPE", "")
                svlen = info_map.get("SVLEN", "")

                format_map = {}
                if len(cols) >= 10:
                    keys = cols[8].split(":")
                    vals = cols[9].split(":")
                    format_map = dict(zip(keys, vals))

                # VAF et support de lecture, puis annotation ORF par position.
                vaf, dr, dv = extract_vaf(info_map, format_map)
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
                        "svtype": svtype,
                        "svlen": svlen,
                        "vaf": "" if vaf is None else f"{vaf:.6f}",
                        "dr": "" if dr is None else f"{dr:.0f}",
                        "dv": "" if dv is None else f"{dv:.0f}",
                        "qual": qual,
                        "orf": orf,
                    }
                )

    # Export 1: CSV complet (source principale pour les plots).
    raw_out = out_dir / "summary_table_sv_with_ORF.csv"
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
                "svtype",
                "svlen",
                "vaf",
                "dr",
                "dv",
                "qual",
                "orf",
            ],
        )
        writer.writeheader()
        writer.writerows(raw_rows)

    print(f"Intermediate CSV written: {raw_out}")
    print(f"Total variants parsed: {len(raw_rows)}")


if __name__ == "__main__":
    main()
