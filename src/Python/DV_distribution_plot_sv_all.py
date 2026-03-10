#!/usr/bin/env python3
import os
import sys
import matplotlib.pyplot as plt


def extract_all_dv(vcf_folder):
    dv_values = []
    for filename in os.listdir(vcf_folder):
        if not (filename.endswith(".vcf") and "sv_sniffles" in filename):
            continue
        filepath = os.path.join(vcf_folder, filename)
        with open(filepath, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                fmt_keys = fields[8].split(":")
                sample_vals = fields[9].split(":")
                fmt = dict(zip(fmt_keys, sample_vals))
                dv = fmt.get("DV")
                if dv is None:
                    continue
                try:
                    dv_values.append(float(dv))
                except ValueError:
                    continue
    return dv_values


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 DV_distribution_plot_sv_all.py <input_folder_vcf>")
        sys.exit(1)

    vcf_folder = sys.argv[1]
    if not os.path.isdir(vcf_folder):
        print(f"Dossier introuvable: {vcf_folder}")
        sys.exit(1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)

    list_dv_dir = os.path.join(script_dir, "List_DV")
    os.makedirs(list_dv_dir, exist_ok=True)
    plot_dir = os.path.join(repo_root, "Plots", "DV_distribution")
    os.makedirs(plot_dir, exist_ok=True)

    dv_values = extract_all_dv(vcf_folder)

    list_path = os.path.join(list_dv_dir, "List_DV_all_SV.txt")
    with open(list_path, "w") as out:
        out.write(",".join(str(v) for v in dv_values))

    if not dv_values:
        print("Aucune valeur DV trouvee dans les VCF SV.")
        print(f"Liste ecrite: {list_path}")
        sys.exit(0)

    plt.figure(figsize=(10, 6))
    plt.hist(dv_values, bins=60, color="skyblue", edgecolor="black")
    plt.title("Distribution des valeurs DV des variants structuraux (tous VCF SV)")
    plt.xlabel("DV (Depth Variant)")
    plt.ylabel("Nombre de variants")
    plt.grid(axis="y", alpha=0.75)

    plot_path = os.path.join(plot_dir, "Distribution_DV_variants_SV_all.png")
    plt.savefig(plot_path)
    plt.close()

    print(f"Nombre total de DV: {len(dv_values)}")
    print(f"Liste enregistree: {list_path}")
    print(f"Plot enregistre: {plot_path}")


if __name__ == "__main__":
    main()
