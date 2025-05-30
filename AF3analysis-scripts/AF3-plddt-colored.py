#!/usr/bin/env python3
# Many thanks to ChatGPT for assisting me!

import argparse
import matplotlib.pyplot as plt

def read_plddt_from_pdb(pdb_file):
    """
    Parse a PDB file, extracting one pLDDT per residue from the B-factor column.
    Also detect chain breaks whenever the chain ID changes.
    Returns (residue_indices, plddt_list, chain_break_positions).
    """
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

            # new residue when either chain or residue number changes
            if (chain, resseq) != (prev_chain, prev_resseq):
                counter += 1
                residues.append(counter)
                plddt.append(bfact)

                if prev_chain is not None and chain != prev_chain:
                    # mark break at end of previous chain
                    chain_breaks.append(counter - 1)

                prev_chain = chain
                prev_resseq = resseq

    return residues, plddt, chain_breaks

def plot_plddt(residues, plddt, breaks, out_prefix):
    fig, ax = plt.subplots(figsize=(20, 10))

    # --- AlphaFold DB confidence bands ---
    ax.fill_between(residues, 0, 50,   color='#ff7f0e', alpha=0.2, label='Very Low (0–50)')
    ax.fill_between(residues, 50, 70,  color='#bcbd22', alpha=0.2, label='Low (50–70)')
    ax.fill_between(residues, 70, 90,  color='#17becf', alpha=0.2, label='High (70–90)')
    ax.fill_between(residues, 90, 100, color='#1f77b4', alpha=0.2, label='Very High (90–100)')
    
    # pLDDT trace
    ax.plot(residues, plddt, linewidth=5, color='black', label='pLDDT')

    # chain‐break lines
    for b in breaks:
        ax.axvline(x=b, color='darkred', linestyle='--', linewidth=2, alpha=0.5)

    # gray dashed thresholds
    for y in (50, 70, 90):
        ax.axhline(y=y, color='gray', linestyle='--', linewidth=1.5)

    # labels, limits & ticks
    ax.set_xlabel("Residue", fontsize=46, labelpad=10)
    ax.set_ylabel("pLDDT", fontsize=46)
    ax.set_ylim(0, 100)
    ax.set_yticks(range(0, 101, 10))
    ax.tick_params(labelsize=28)

    ax.legend(loc='upper right', fontsize=24)
    plt.tight_layout()

    out_file = f"{out_prefix}_plddt.jpeg"
    fig.savefig(out_file, dpi=500, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved pLDDT plot to: {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Plot AF3 pLDDT directly from a PDB")
    parser.add_argument("--pdb_file",
                        help="Path to your AF3 PDB (pLDDT in B-factor)",
                        required=True)
    parser.add_argument("--out_prefix",
                        help="Prefix for output files (default: AF3)",
                        default="AF3")
    args = parser.parse_args()

    residues, plddt, chain_breaks = read_plddt_from_pdb(args.pdb_file)
    plot_plddt(residues, plddt, chain_breaks, args.out_prefix)

if __name__ == "__main__":
    main()
