#!/usr/bin/env python
# Many thanks to ChatGPT for assisting me!

"""
AF3-pae.py - PAE-only visualization from AlphaFold3 JSON files

Generates a PAE heatmap from AlphaFold3 _full_data_0.json and 
_summary_confidences_0.json. Adds chain boundary lines if a PDB file is provided.
"""

import argparse
import json
import numpy as np
import os
import matplotlib.pyplot as plt
from Bio import PDB

# --- Argument Parsing ---
parser = argparse.ArgumentParser()
parser.add_argument("--data_file", type=str, required=True, help="AlphaFold3 _full_data_0.json file")
parser.add_argument("--pdb_file", type=str, help="Optional PDB file for chain boundaries")
args = parser.parse_args()

# --- Path Handling ---
if not args.data_file.endswith("_full_data_0.json"):
    raise ValueError("The data_file must end with '_full_data_0.json'")

prefix = args.data_file.replace("_full_data_0.json", "")
summary_file = prefix + "_summary_confidences_0.json"
output_file = os.path.splitext(args.data_file)[0]

# --- Load Data ---
with open(args.data_file, "r") as f:
    data = json.load(f)

with open(summary_file, "r") as f:
    summary_data = json.load(f)

pae_data = data.get("pae", [])
ptm = summary_data.get("ptm", 0.0)
iptm = summary_data.get("iptm", 0.0)

# --- Chain Boundary Calculation ---
cumulative_sum = []
if args.pdb_file:
    parser = PDB.PDBParser()
    structure = parser.get_structure("AFmodel", args.pdb_file)
    residue_counts = {}

    for chain in structure.get_chains():
        count = np.sum(np.fromiter((1 for res in chain.get_residues() if PDB.is_aa(res)), dtype=int))
        residue_counts[chain.get_id()] = count

    counts = list(residue_counts.values())
    cumulative_sum = [sum(counts[:i+1]) for i in range(len(counts)-1)]

# --- Plot PAE Matrix ---
matrix = np.array(pae_data)
fig, ax = plt.subplots(figsize=(7, 9))
im = ax.imshow(matrix, cmap=plt.cm.Greens_r, vmin=0, vmax=30)

# Add colorbar and labels
cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', pad=0.12)
cbar.set_label(r"Expected position error (Ångströms)"
               "\n"
               "\n"
               f"ptm={ptm:.3f}  iptm={iptm:.3f}", fontsize=18)
cbar.ax.tick_params(labelsize=18)

plt.xlabel("Scored residue", fontsize=18)
plt.ylabel("Aligned residue", fontsize=18)

for tick in cbar.ax.get_yticklabels():
    ax.tick_params(labelsize=24)

# Add chain boundary lines
for resnum in cumulative_sum:
    ax.axvline(x=resnum, color='darkred', linestyle='--', linewidth=2, alpha=0.5)
    ax.axhline(y=resnum, color='darkred', linestyle='--', linewidth=2, alpha=0.5)

fig.tight_layout()
fig.savefig(f"{output_file}_PAE.jpeg", dpi=500, bbox_inches='tight')
plt.close()

print(f"Saved: {output_file}_PAE.jpeg")
