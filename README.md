# AFanalysis

Here, you can use the provided python script for the analysis of AlphaFold models, including both local predictions that produce .pkl files and ColabFold predictions that produce .json files storing quality evaluation scores. 

With this script, you can only analyze one AlphaFold model at a time.

	If only pkl/json file is provided, this script outputs:
	   - a json file containing plddt, pae, max_pae, ranking confidence, ptm and iptm scores (%s.json)
	   - draws a PAE graph showing iptm and ptm values (%s_PAE.jpeg)
	   - draws a pLDDT graph showing the plddt score for each residue (%s_Plddt.jpeg)
	   
	If only pdb file is provided, this script outputs:
	- a PyMOL-made figure showing the structure colored by chains (%s_cbc.png)
	- a PyMOL-made figure showing the structure colored by AlphaFold coloring (%s_cbaf.png)
	
	If both pkl/json and pdb files are provided, this script outputs:
	- All mentioned above and updated versions of PAE and pLDDT graphs showing the residue numbers of each monomers as lines. 

**Usage: **
```
	python AFanalysis.py --data_file <pkl file or json file>
	python AFanalysis.py --pdb_file <pdb file>
	or 
	python AFanalysis.py --data_file <pkl file or json file> --pdb_file <pdb file>
```

# AF3analysis: AlphaFold3 Model Confidence Visualizations Scripts

This repository provides tools for analyzing and visualizing model confidence scores from AlphaFold and AlphaFold3, including:

- **Predicted Aligned Error (PAE)** matrices
- **Per-residue pLDDT** confidence plots
- Optional **chain break visualization** using PDB files
- Compatible with AlphaFold3 `.json` and `.pdb` outputs

---

## Scripts

### `AF3-PAE.py`

Generates:
- A **PAE matrix** using `_full_data_0.json`
- A **pLDDT line plot** from the same JSON
- Optional vertical chain break lines (if a PDB file is provided)

**Usage:**
```
python AF3-PAE.py --data_file model_full_data_0.json --pdb_file model.pdb```

### `AF3-plddt.py`

Generates:

- A **per-residue pLDDT plot** using values from the B-factor column of an AlphaFold3-generated PDB file
- A **clean line graph**, without background confidence bands
- Optional vertical chain break lines based on changes in chain ID

**Usage:**

```
python AF3-plddt.py --pdb_file model.pdb```

### `AF3-plddt-colored.py`

Generates:

- A **per-residue pLDDT plot* using values from the B-factor column of a PDB file

- A **confidence-colored** background showing AlphaFold-style score bands:

Very Low (0–50) Orange
Low (50–70) Yellow
High (70–90) Cyan Blue
Very High (90–100) Blue

- Optional vertical chain break lines based on changes in chain ID

Usage:

```
python AF3-plddt-colored.py --pdb_file model.pdb --out_prefix model```

Recommended to use these scripts by creating a conda environment.
```
	conda create -n AFanalysis python=3.7
	conda activate AFanalysis
	conda install -c schrodinger pymol -y
	conda install matplotlib
```

