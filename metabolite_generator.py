import os
import sys
import argparse
import glob
import pandas as pd
import matplotlib.pyplot as plt
import io
from PIL import Image

# Secure inline RDKit chemistry rendering dependencies
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

def parse_arguments():
    parser = argparse.ArgumentParser(description="REACCESS Phase I Automated Metabolite Generation Engine.")
    parser.add_argument("-excel", "--excel_file", required=True, help="Path to METABOLIC_FINAL_SUMMARY.xlsx or COMBINED sheet.")
    parser.add_argument("-l", "--ligand_ref", required=True, help="Path to original parent reference ligand PDB file.")
    parser.add_argument("-s", "--smiles", required=True, help="Parent reference ligand SMILES string.")
    parser.add_argument("-o", "--output_dir", default="./GENERATED_METABOLITES", help="Destination folder for generated assets.")
    return parser.parse_args()

def align_pdb_names_to_rdkit(mol_base, mol_pdb):
    """Maps custom PDB file naming formats cleanly onto RDKit index registers."""
    smiles_to_pdb_name = {}
    pdb_name_to_rdkit_idx = {}
    
    mol_pdb_heavy = Chem.RemoveHs(mol_pdb)
    mcs_res = rdFMCS.FindMCS([mol_base, mol_pdb_heavy], 
                             atomCompare=rdFMCS.AtomCompare.CompareElements, 
                             bondCompare=rdFMCS.BondCompare.CompareAny)
    
    if mcs_res.numAtoms > 0:
        mcs_mol = Chem.MolFromSmarts(mcs_res.smartsString)
        m_base = mol_base.GetSubstructMatch(mcs_mol)
        m_pdb = mol_pdb_heavy.GetSubstructMatch(mcs_mol)
        
        for b_idx, p_idx in zip(m_base, m_pdb):
            p_inf = mol_pdb_heavy.GetAtomWithIdx(p_idx).GetPDBResidueInfo()
            p_name = p_inf.GetName().strip() if p_inf else f"{mol_pdb_heavy.GetAtomWithIdx(p_idx).GetSymbol()}{p_idx+1}"
            smiles_to_pdb_name[b_idx] = p_name
            pdb_name_to_rdkit_idx[p_name] = b_idx
    else:
        for atom in mol_base.GetAtoms():
            idx = atom.GetIdx()
            p_inf = atom.GetPDBResidueInfo()
            p_name = p_inf.GetName().strip() if p_inf else f"{atom.GetSymbol()}{idx+1}"
            smiles_to_pdb_name[idx] = p_name
            pdb_name_to_rdkit_idx[p_name] = idx
            
    return smiles_to_pdb_name, pdb_name_to_rdkit_idx

def draw_molecule_acs_style(mol, filename, title_text, label_dict=None, highlight_indices=None, label_indices=None, highlight_color=(1.0, 0.0, 0.0)):
    """Generates a skeletal 2D image showing text labels and circle highlights ONLY for active/modified atoms."""
    mol_copy = Chem.Mol(mol)
    for atom in mol_copy.GetAtoms():
        atom.SetIntProp('orig_idx', atom.GetIdx())
        
    mol_drawn = Chem.RemoveHs(mol_copy)
    Chem.RemoveStereochemistry(mol_drawn)
    mol_drawn.RemoveAllConformers()
    AllChem.Compute2DCoords(mol_drawn)
    
    # Process circle highlights mapping
    atoms_to_highlight = []
    if highlight_indices is not None:
        highlight_indices_set = set(int(i) for i in highlight_indices)
        for atom in mol_drawn.GetAtoms():
            if atom.GetIntProp('orig_idx') in highlight_indices_set:
                atoms_to_highlight.append(atom.GetIdx())
                
    # Isolate authorized label targets
    label_indices_set = set(int(i) for i in label_indices) if label_indices is not None else set()
                
    d2d = rdMolDraw2D.MolDraw2DCairo(800, 800)
    draw_opts = d2d.drawOptions()
    draw_opts.addAtomIndices = False
    draw_opts.lineWidth = 2.4            
    draw_opts.fixedBondLength = 35       
    draw_opts.fillHighlights = False     
    
    # FIXED LOGIC LAYER: Render atom labels exclusively for authorized parent active centers
    for atom in mol_drawn.GetAtoms():
        idx = atom.GetIdx()
        orig_idx = atom.GetIntProp('orig_idx')
        
        if orig_idx in label_indices_set:
            if label_dict and orig_idx in label_dict and label_dict[orig_idx] != "":
                draw_opts.atomLabels[idx] = label_dict[orig_idx]
        elif atom.GetSymbol() == "C":
            # Unlabeled carbons must remain hidden skeletal vertices
            draw_opts.atomLabels[idx] = ""
        # Mapped heteroatoms (N, O, S) are intentionally left out of the text label dictionary
        # so RDKit handles their element symbols and implicit valence hydrogens natively.
            
    h_colors = {idx: highlight_color for idx in atoms_to_highlight}
    h_radii = {idx: 0.45 for idx in atoms_to_highlight}
    
    d2d.DrawMolecule(mol_drawn, highlightAtoms=atoms_to_highlight, highlightAtomColors=h_colors, highlightAtomRadii=h_radii, highlightBonds=[])
    d2d.FinishDrawing()
    
    fig, ax = plt.subplots(figsize=(6, 6), facecolor='white')
    ax.imshow(Image.open(io.BytesIO(d2d.GetDrawingText())))
    ax.axis('off')
    ax.set_title(title_text, fontsize=11, fontweight='bold', pad=8)
    fig.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig)

def main():
    args = parse_arguments()
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    print("--> [Metabolite Core] Extracting target matrices from summary logs...")
    if not os.path.exists(args.excel_file):
        print(f"[CRITICAL FAILURE] Targeted excel spreadsheet metrics log file missing: {args.excel_file}")
        sys.exit(1)
        
    df = pd.read_excel(args.excel_file)
    if 'Atom_Name' not in df.columns:
        print("[CRITICAL FAILURE] Column 'Atom_Name' missing from selected target analytical spreadsheet.")
        sys.exit(1)
        
    df_active = df[df.get('Activity_Frequency_%', pd.Series([100.0]*len(df))) > 0.0]
    active_atoms = df_active['Atom_Name'].dropna().astype(str).str.strip().unique().tolist()
    
    if not active_atoms:
        print("--> [Metabolite Core] Zero active clear SOM target entities found inside file records.")
        return
        
    print(f"--> Found {len(active_atoms)} active structural centers scheduled for metabolic expansion: {active_atoms}")

    mol_pdb = Chem.MolFromPDBFile(args.ligand_ref, removeHs=False)
    mol_base = Chem.MolFromSmiles(args.smiles)
    if mol_base is None:
        raise ValueError("Flawless structural assignment failure: Input parent SMILES parsing returned None validation.")
    mol_base = Chem.AddHs(mol_base)
    Chem.SanitizeMol(mol_base)
    
    smiles_to_pdb_name, pdb_name_to_rdkit_idx = align_pdb_names_to_rdkit(mol_base, mol_pdb)
    
    parent_highlight_indices = []
    for atom_name in active_atoms:
        if atom_name in pdb_name_to_rdkit_idx:
            parent_highlight_indices.append(pdb_name_to_rdkit_idx[atom_name])

    # Save parent reference map showing labels and circles on active indices
    parent_png_path = os.path.join(output_dir, "parent_structure_reference.png")
    draw_molecule_acs_style(mol_base, parent_png_path, 
                            "Reference Parent Structure Map: Active SOM Centers Highlighted", 
                            label_dict=smiles_to_pdb_name, 
                            highlight_indices=parent_highlight_indices,
                            label_indices=parent_highlight_indices,
                            highlight_color=(1.0, 0.0, 0.0))

    for target_atom_name in active_atoms:
        if target_atom_name not in pdb_name_to_rdkit_idx:
            print(f"    [Skipped Target] Atom {target_atom_name} could not be securely mapped to the topology index registers.")
            continue
            
        target_idx = pdb_name_to_rdkit_idx[target_atom_name]
        atom_obj = mol_base.GetAtomWithIdx(target_idx)
        element = atom_obj.GetSymbol()
        
        metabolite_name = f"parent-{target_atom_name}"
        print(f"--> Processing target asset generation for: {metabolite_name} (Element Center: {element})...")
        
        rw_mol = Chem.RWMol(mol_base)
        
        if element == 'C':
            h_neighbor_idx = None
            for neighbor in atom_obj.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    h_neighbor_idx = neighbor.GetIdx()
                    break
            
            if h_neighbor_idx is not None:
                for rw_atom in rw_mol.GetAtoms():
                    if rw_atom.GetIdx() == atom_obj.GetIdx():
                        for rw_nb in rw_atom.GetNeighbors():
                            if rw_nb.GetSymbol() == 'H':
                                rw_mol.RemoveAtom(rw_nb.GetIdx())
                                break
                        break
            
            o_idx = rw_mol.AddAtom(Chem.Atom(8))  
            rw_mol.AddBond(target_idx, o_idx, Chem.BondType.SINGLE)
            h_idx = rw_mol.AddAtom(Chem.Atom(1))  
            rw_mol.AddBond(o_idx, h_idx, Chem.BondType.SINGLE)
            
        elif element in ['N', 'S']:
            o_idx = rw_mol.AddAtom(Chem.Atom(8))  
            if element == 'N':
                rw_mol.AddBond(target_idx, o_idx, Chem.BondType.SINGLE)
                rw_mol.GetAtomWithIdx(target_idx).SetFormalCharge(atom_obj.GetFormalCharge() + 1)
                rw_mol.GetAtomWithIdx(o_idx).SetFormalCharge(-1)
            else:
                rw_mol.AddBond(target_idx, o_idx, Chem.BondType.DOUBLE)
            
        else:
            print(f"    [Skipped Type] Element {element} conversion paths are outside Phase I tracking rules boundaries.")
            continue

        metabolite_mol = rw_mol.GetMol()
        try:
            Chem.SanitizeMol(metabolite_mol)
        except Exception as san_err:
            print(f"    [Sanitization Warning] Structural convergence aborted for {metabolite_name}: {san_err}")
            continue

        mol_3d = Chem.Mol(metabolite_mol)
        AllChem.EmbedMolecule(mol_3d, randomSeed=42, maxAttempts=500)
        try:
            AllChem.MMFFOptimizeMolecule(mol_3d, maxIters=200)
        except:
            AllChem.UFFOptimizeMolecule(mol_3d, maxIters=200)

        sdf_out_path = os.path.join(output_dir, f"{metabolite_name}.sdf")
        pdb_out_path = os.path.join(output_dir, f"{metabolite_name}.pdb")
        png_out_path = os.path.join(output_dir, f"{metabolite_name}.png")
        
        sdf_writer = Chem.SDWriter(sdf_out_path)
        sdf_writer.write(mol_3d)
        sdf_writer.close()
        
        Chem.MolToPDBFile(mol_3d, pdb_out_path)
        
        # Populate the custom identification dictionary mapping exclusively for active centers
        metabolite_labels = {}
        for atom in metabolite_mol.GetAtoms():
            idx = atom.GetIdx()
            metabolite_labels[idx] = smiles_to_pdb_name.get(idx, "")

        # FIXED SETTINGS: The red highlight circle sits strictly centered on the reacting carbon atom vertex.
        # Only the parent carbon index is passed to label_indices so it retains its custom PDB text ID.
        highlight_target = [target_idx]
        metabolite_label_targets = [target_idx]

        # Draw the structure cleanly using the upgraded skeletal rendering parameters
        draw_molecule_acs_style(metabolite_mol, png_out_path, 
                                f"Generated Phase I Metabolite Structure: {metabolite_name}", 
                                label_dict=metabolite_labels, 
                                highlight_indices=highlight_target,
                                label_indices=metabolite_label_targets,
                                highlight_color=(1.0, 0.0, 0.0))

    print("\n" + "="*80)
    print(" [COMPLETE SUCCESS] ALL PHASE I METABOLITE STRUCTURES COMPILED")
    print("="*80)
    print(f" Output Location: {output_dir}/")
    print(f" Parent Diagram:  {os.path.basename(parent_png_path)}")
    print(f" Total Generated: {len(active_atoms)} (.sdf, .pdb, and clean ACS flat .png arrays)")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()