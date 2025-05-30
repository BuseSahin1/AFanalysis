#!/usr/bin/env python
# Many thanks to ChatGPT for assisting me!


"""
AFanalysis.py - Final version with --pdb_file support

Generates PAE and pLDDT plots from AlphaFold .json files, 
with chain boundary lines if a PDB file is provided.
"""

import argparse
import json
import numpy as np
import os
import matplotlib.pyplot as plt
from Bio import PDB

# --- Argument Parsing ---
parser = argparse.ArgumentParser()
parser.add_argument("--data_file", type=str, required=True, help="AlphaFold _full_data_0.json file")
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
plddt = data.get("atom_plddts", [])
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

for resnum in cumulative_sum:
    ax.axvline(x=resnum, color='darkred', linestyle='--', linewidth=2, alpha=0.5)
    ax.axhline(y=resnum, color='darkred', linestyle='--', linewidth=2, alpha=0.5)

fig.tight_layout()
fig.savefig(f"{output_file}_PAE.jpeg", dpi=500, bbox_inches='tight')
plt.close()

# --- Plot pLDDT Graph ---
f = plt.figure(figsize=(20, 10))
ax = plt.gca()
residues = list(range(1, len(plddt) + 1))
ax.set_xlabel("Residues", fontsize=46, labelpad=10)
ax.set_ylabel("pLDDT", fontsize=46)
ax.plot(residues, plddt)
ax.tick_params(labelsize=28)
ax.get_lines()[0].set_linewidth(5)

for resnum in cumulative_sum:
    ax.axvline(x=resnum, color='darkred', linestyle='--', linewidth=2, alpha=0.5)

plt.savefig(f"{output_file}_Plddt.jpeg", dpi=500, bbox_inches='tight')
plt.close()

print(f"Saved: {output_file}_PAE.jpeg and {output_file}_Plddt.jpeg")
