#!/usr/bin/env python3
# Many thanks to ChatGPT for assisting me!

import argparse
import os
import matplotlib.pyplot as plt

def read_plddt_from_pdb(pdb_file):
    residues = []
    plddt = []
    chain_breaks = []

    prev_chain = None
    prev_resseq = None
    counter = 0

    with open(pdb_file) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            resseq = int(line[22:26])
            bfact = float(line[60:66])

            if (chain, resseq) != (prev_chain, prev_resseq):
                counter += 1
                residues.append(counter)
                plddt.append(bfact)

                if prev_chain is not None and chain != prev_chain:
                    chain_breaks.append(counter - 1)

                prev_chain = chain
                prev_resseq = resseq

    return residues, plddt, chain_breaks

def plot_plddt(residues, plddt, breaks, pdb_file):
    fig, ax = plt.subplots(figsize=(20, 10))

    # pLDDT trace
    ax.plot(residues, plddt, linewidth=5, color='blue', label='pLDDT')

    # chain‐break lines
    for b in breaks:
        ax.axvline(x=b, color='darkred', linestyle='--', linewidth=2, alpha=0.5)

    # gray dashed thresholds
    for y in (50, 70, 90):
        ax.axhline(y=y, color='gray', linestyle='--', linewidth=1.5)

    ax.set_xlabel("Residue", fontsize=46, labelpad=10)
    ax.set_ylabel("pLDDT", fontsize=46)
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 10))
    ax.tick_params(labelsize=28)

    ax.legend(loc='upper right', fontsize=24)
    plt.tight_layout()

    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    out_file = f"{pdb_name}_plddt.jpeg"
    fig.savefig(out_file, dpi=500, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved pLDDT plot to: {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Plot AF3 pLDDT directly from a PDB")
    parser.add_argument("--pdb_file", help="Path to your AF3 PDB (pLDDT in B-factor)", required=True)
    args = parser.parse_args()

    residues, plddt, chain_breaks = read_plddt_from_pdb(args.pdb_file)
    plot_plddt(residues, plddt, chain_breaks, args.pdb_file)

if __name__ == "__main__":
    main()

