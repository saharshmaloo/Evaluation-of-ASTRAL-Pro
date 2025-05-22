#!/usr/bin/env python3
import treeswift
import csv
import sys
from pathlib import Path
from typing import List, Tuple, Dict
from ete3 import Tree
def contract_nonD_edges(tree):
    for node in list(tree.traverse_postorder()):
        if node.is_root():
            continue
        parent = node.parent
        if parent and node and parent.label is None and node.label is None:
            for child in list(node.children):
                parent.add_child(child)
            parent.remove_child(node)
dlrates = ["0.0000000001", "0.0000000002", "0.0000000005"]
psizes = ["10000000", "50000000"]
base = Path("data")
def run():
    result_file = base / "comparison_results.txt"
    with open(result_file, "w") as result_out:
        result_out.write("Comparison Results\n")
        result_out.flush()
        for rate in dlrates:
            for size in psizes:
                parent = base / f"ntaxa-100.dlrate-{rate}.psize-{size}"
                for rep in range(1, 11):
                    folder = parent / f"{rep:02}"
                    ground_file = folder / "tagged_g_trees"
                    astral_file = folder / "astral_tagged_g_trees"
                    output_file = folder / "comparison.txt"
                    if not ground_file.exists() or not astral_file.exists():
                        print(f"Skipping {folder}, missing tree files.")
                        continue
                    with open(ground_file) as f1, open(astral_file) as f2, open(output_file, "w") as out:
                        num_trees = 0
                        total_rf = 0
                        total_edge_similarity = 0
                        total_reverse_edge_similarity = 0
                        for i, (gt_line, ap_line) in enumerate(zip(f1, f2), 1):
                            gt_tree = treeswift.read_tree_newick(gt_line.strip())
                            ap_tree = treeswift.read_tree_newick(ap_line.strip())
                            contract_nonD_edges(gt_tree)
                            contract_nonD_edges(ap_tree)

                            gt_newick = gt_tree.newick()
                            ap_newick = ap_tree.newick()
                            # Convert to ETE3 Tree objects
                            gt_tree_ete = Tree(gt_newick, format=1)
                            ap_tree_ete = Tree(ap_newick, format=1)
                            
                            try:
                                tree_data = ap_tree_ete.compare(gt_tree_ete, unrooted=True)
                                rf, max_rf = tree_data["rf"], tree_data["max_rf"]
                                normalized_rf = rf / max_rf if max_rf else 0
                                total_rf += normalized_rf
                                out.write(f"Tree {i:04}: RF={rf}, Normalized RF={normalized_rf:.4f}\n")
                                edge_similarity = 100 * tree_data["ref_edges_in_source"],  # how many ref edges are in source
                                reverse_edge_similarity = 100 * tree_data["source_edges_in_ref"],  # source edges in ref
                                total_edge_similarity += edge_similarity[0]
                                total_reverse_edge_similarity += reverse_edge_similarity[0]
                                out.write(f"Tree {i:04}: Edge Similarity={edge_similarity[0]:.4f}, Reverse Edge Similarity={reverse_edge_similarity[0]:.4f}\n")
                                num_trees += 1
                            except Exception as e:
                                out.write(f"Tree {i:04}: Error - {e}\n")
                                continue
                        robinson_foulds = total_rf / num_trees if num_trees > 0 else 0
                        edge_similarity = total_edge_similarity / num_trees if num_trees > 0 else 0
                        reverse_edge_similarity = total_reverse_edge_similarity / num_trees if num_trees > 0 else 0
                        out.write(f"Total Trees Evaluated: {num_trees}\n")
                        out.write(f"Average Normalized Robinson-Foulds Distance: {robinson_foulds:.4f}\n")
                        out.write(f"Average Edge Similarity: {edge_similarity:.4f}\n")
                        out.write(f"Average Reverse Edge Similarity: {reverse_edge_similarity:.4f}\n")
                        result_out.write(f"psize: {size}, rate: {rate}, rep: {folder.name}, Normalized Robinson-Foulds: {robinson_foulds:.4f}, Edge Similarity: {edge_similarity:.4f}, Reverse Edge Similarity: {reverse_edge_similarity:.4f}\n")
                        result_out.flush()
                        print(f"psize: {size}, rate: {rate}, rep: {folder.name}, Normalized Robinson-Foulds: {robinson_foulds:.4f}, Edge Similarity: {edge_similarity:.4f}, Reverse Edge Similarity: {reverse_edge_similarity:.4f}\n")

if __name__ == "__main__":
    run()