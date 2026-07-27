import os
import io
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Pristine 2D ACS-Style Molecular Structure Renderer with Subgraph PDB Mapping."
    )
    parser.add_argument("-l", "--ligand", required=True, help="Path to the input ligand reference file (.pdb).")
    parser.add_argument("-s", "--summary", required=True, help="Path to METABOLIC_FINAL_SUMMARY.xlsx file.")
    parser.add_argument("-o", "--output", default="ligand_metabolic_2D.png", help="Name of output image file.")
    parser.add_argument("-sm", "--smiles", default=None, help="SMILES string to secure exact aromaticity and bond layouts.")
    return parser.parse_args()

def main():
    args = parse_arguments()

    if not os.path.exists(args.summary):
        raise FileNotFoundError(f"Could not locate summary spreadsheet: {args.summary}")
    if not os.path.exists(args.ligand):
        raise FileNotFoundError(f"Could not locate ligand PDB structure: {args.ligand}")

    ligand_filename_base = os.path.splitext(os.path.basename(args.ligand))[0]

    # 1. Parse the metabolic summary spreadsheet
    print(f"--> Reading metabolic tracking profiles from: {args.summary}")
    df_summary = pd.read_excel(args.summary)
    active_atoms_pdb = df_summary[df_summary['Activity_Frequency_%'] > 0.0]['Atom_Name'].astype(str).str.strip().tolist()
    print(f"--> Target active metabolic centers found: {active_atoms_pdb}")

    # 2. Load ground-truth PDB to extract original atom labels
    mol_pdb = Chem.MolFromPDBFile(args.ligand, removeHs=False)
    if mol_pdb is None:
        raise ValueError(f"RDKit failed to load the structural PDB file: {args.ligand}")
    mol_pdb_heavy = Chem.RemoveHs(mol_pdb)

    # 3. Generate the pristine 2D layout molecule base
    if args.smiles:
        print("--> Constructing pristine chemical framework from canonical SMILES...")
        mol_base = Chem.MolFromSmiles(args.smiles)
        Chem.SanitizeMol(mol_base)
        
        # Cross-reference the SMILES graph with the PDB coordinates to map names
        print("--> Running Maximum Common Subgraph (MCS) isomorphism alignment layer...")
        mcs_res = rdFMCS.FindMCS([mol_base, mol_pdb_heavy], 
                                 atomCompare=rdFMCS.AtomCompare.CompareElements,
                                 bondCompare=rdFMCS.BondCompare.CompareAny)
        
        smiles_to_pdb_name = {}
        if mcs_res.numAtoms > 0:
            mcs_mol = Chem.MolFromSmarts(mcs_res.smartsString)
            match_base = mol_base.GetSubstructMatch(mcs_mol)
            match_pdb = mol_pdb_heavy.GetSubstructMatch(mcs_mol)
            
            for base_idx, pdb_idx in zip(match_base, match_pdb):
                atom_pdb = mol_pdb_heavy.GetAtomWithIdx(pdb_idx)
                p_info = atom_pdb.GetPDBResidueInfo()
                pdb_name = p_info.GetName().strip() if p_info else f"{atom_pdb.GetSymbol()}{pdb_idx+1}"
                smiles_to_pdb_name[base_idx] = pdb_name
        else:
            print("[WARNING] Isomorphism mapping failed. Falling back to index ordering.")
            for atom in mol_base.GetAtoms():
                smiles_to_pdb_name[atom.GetIdx()] = f"{atom.GetSymbol()}{atom.GetIdx()+1}"
    else:
        # Fallback if no SMILES string is provided
        print("--> No SMILES provided. Deducing bonds directly from PDB heavy atoms...")
        mol_base = Chem.DataFrame(mol_pdb_heavy)
        try: Chem.DetermineBondsOrderDetail(mol_base)
        except: pass
        Chem.SanitizeMol(mol_base)
        smiles_to_pdb_name = {}
        for atom in mol_base.GetAtoms():
            p_info = atom.GetPDBResidueInfo()
            smiles_to_pdb_name[atom.GetIdx()] = p_info.GetName().strip() if p_info else f"{atom.GetSymbol()}{atom.GetIdx()+1}"

    # 4. FIXED: Completely strip stereochemistry configurations (no wedges/dashes)
    print("--> Purging stereochemistry flags and calculating flat 2D projection...")
    Chem.RemoveStereochemistry(mol_base)
    mol_base.RemoveAllConformers()
    AllChem.Compute2DCoords(mol_base)

    # 5. Set up the Cairo graphics canvas with ACS drawing specifications
    d2d = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
    draw_opts = d2d.drawOptions()
    draw_opts.addAtomIndices = False
    draw_opts.atomLabelScale = 0.75
    draw_opts.bondLength = 30
    draw_opts.lineWidth = 2.6
    
    # Global overrides to remove text stereochemical markers
    draw_opts.addStereoAnnotation = False
    draw_opts.includeStereoComments = False

    highlight_atoms = []
    
    # 6. Apply clean labels based on the mapping dictionary
    for atom in mol_base.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        pdb_name = smiles_to_pdb_name.get(idx, f"{symbol}{idx+1}")
        
        if pdb_name in active_atoms_pdb:
            highlight_atoms.append(idx)
            # FIXED: Show ONLY the clean PDB identity label name inside the circle
            draw_opts.atomLabels[idx] = pdb_name
        elif symbol != "C":
            # Display heteroatoms cleanly (N, O, F, etc.)
            draw_opts.atomLabels[idx] = symbol
        else:
            # Leave normal skeletal carbons completely blank matching standard ACS layout rules
            draw_opts.atomLabels[idx] = ""

    # FIXED: Highlight active atoms using standard RED publication halos
    highlight_colors = {int(idx): (1.0, 0.0, 0.0) for idx in highlight_atoms}
    highlight_radii = {int(idx): 0.40 for idx in highlight_atoms}

    d2d.DrawMolecule(
        mol_base,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_colors,
        highlightAtomRadii=highlight_radii,
        highlightBonds=[]
    )
    d2d.FinishDrawing()
    
    # 7. Output to image file
    img_bytes = d2d.GetDrawingText()
    img = Image.open(io.BytesIO(img_bytes))

    fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')
    ax.imshow(img)
    ax.axis('off')
    ax.set_title(f"CYP450 2D Metabolic Structure Profile: {ligand_filename_base}", fontsize=12, fontweight='bold', pad=5)
    
    fig.savefig(args.output, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"[SUCCESS] Clean 2D ACS publication structure profile layout generated at: '{args.output}'")

if __name__ == "__main__":
    main()