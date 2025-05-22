#!/usr/bin/env python3
import treeswift
import csv
from pathlib import Path
import subprocess
# Path to ASTRAL-Pro3 binary
ASTRAL = r"/mnt/d/College/Sophomore 2nd/CMSC 499A/ASTER-Windows/ASTER-Windows/bin/astral-pro3"
# Directory containing SimPhy output
DATA_DIR = Path("./data")
# Parameter sets
dlrates = ["0.0000000001", "0.0000000002", "0.0000000005"]
psizes = ["10000000", "50000000"]
def generate_mapping(gtrees_file, mapping_file):
    """Generates a gene-to-species mapping TSV from g_trees."""
    mapping = {}
    with open(gtrees_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            tree = treeswift.read_tree_newick(line)
            for leaf in tree.traverse_leaves():
                label = leaf.label
                species = label.split("_")[0]
                mapping[label] = species
    with open(mapping_file, "w") as outfile:
        for gene, species in mapping.items():
            outfile.write(f"{gene}\t{species}\n")
    print(f":page_facing_up: Mapping file written to: {mapping_file}")
def run():
    for rate in dlrates:
        for size in psizes:
            base = DATA_DIR / f"ntaxa-100.dlrate-{rate}.psize-{size}"
            for rep in range(1, 11):
                rep_dir = base / f"{rep:02}"
                gtrees_path = rep_dir / "g_trees"
                tagged_path = rep_dir / "astral_tagged_g_trees"
                mapping_path = rep_dir / "mapping.tsv"
                if not gtrees_path.exists():
                    print(f" Missing: {gtrees_path}, skipping.")
                    continue
                print(f"\n:file_folder: Processing: {rep_dir}")
                # Step 1: Create mapping file
                print(f" Creating mapping file...")
                generate_mapping(gtrees_path, mapping_path)
                # Step 2: Run ASTRAL-Pro tagging
                print(f" Running ASTRAL-Pro tagging...")
                with open(tagged_path, "w") as fout:
                    subprocess.run([
                        ASTRAL,
                        "-T",
                        "-a", str(mapping_path),
                        "-i", str(gtrees_path),
                    ], stdout=fout, check=True)
                print(f" Done tagging: {tagged_path.name}")
if __name__ == "__main__":
    run()