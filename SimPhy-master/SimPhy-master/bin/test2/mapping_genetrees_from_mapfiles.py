import os
from pathlib import Path
import treeswift
import csv

def read_maplg(path):
    gt_to_lt = {}
    with open(path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#') or row[0] == 'Gt_node':
                continue
            gt_node = row[0].strip("'")
            lt_node = row[1].strip("'")
            gt_to_lt[gt_node] = lt_node
    
    return gt_to_lt

def read_mapsl(path):
    lt_to_event = {}
    with open(path) as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0] == 'Lt_node':
                continue
            lt_node = row[0]
            if row[2] == 'Sp' or row[2] == 'Dup':
                lt_to_event[lt_node] = row[2]
    return lt_to_event

def tag_gene_tree(gt_path, maplg_path, mapsl_path):
    trees = treeswift.read_tree_newick(gt_path)
    tree = trees[0] if isinstance(trees, list) else trees
    gt_to_lt = read_maplg(maplg_path)
    lt_to_event = read_mapsl(mapsl_path)
    for node in tree.traverse_postorder():
        if node.is_leaf():
            continue
        label = node.label.strip("'") if node.label else None
        if label and label in gt_to_lt:
            lt_node = gt_to_lt[label]
            if lt_node in lt_to_event and lt_to_event[lt_node] == "Dup":
                node.label = "D"
            else:
                node.label = None
    return tree.newick().strip()

dlrates = ["0.0000000001", "0.0000000002", "0.0000000005"]
psizes = ["10000000", "50000000"]
base_dir = Path("data")

def main():
    for rate in dlrates:
        for size in psizes:
            parent_folder = base_dir / f"ntaxa-100.dlrate-{rate}.psize-{size}"
            for rep in range(1,11):
                rep_folder = parent_folder / f"{rep:02}"
                g_trees_path = rep_folder / "g_trees"
                tagged_g_trees_path = rep_folder / "tagged_g_trees"
                # Ensure output directory exists
                rep_folder.mkdir(parents=True, exist_ok=True)
                # Open output files
                with open(g_trees_path, "w") as g_out, open(tagged_g_trees_path, "w") as tg_out:
                    for i in range(1, 1001):
                        gene_file = rep_folder / f"g_trees{i:04}.trees"
                        maplg_file = rep_folder / f"{i:04}l1g.maplg"
                        mapsl_file = rep_folder / f"{i:04}.mapsl"
                        if not all(p.exists() for p in [gene_file, maplg_file, mapsl_file]):
                            print(f"Missing files for rate {rate}, size {size}, replicate {rep:02}, gene {i:04}, skipping.")
                            continue
                        # Append raw gene tree
                        with open(gene_file) as f:
                            g_out.write(f.read().strip() + "\n")
                        try:
                            tagged = tag_gene_tree(gene_file, maplg_file, mapsl_file)
                            tg_out.write(tagged + "\n")
                        except Exception as e:
                            print(f" Tagging failed for rate {rate}, size {size}, replicate {rep:02}, gene {i:04}: {e}")
                print(f" Finished: {rep_folder}")
if __name__ == "__main__":
    main()