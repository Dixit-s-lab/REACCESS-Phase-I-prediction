import os
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import MDAnalysis as mda

def parse_arguments():
    """Parses command-line arguments for the script."""
    parser = argparse.ArgumentParser(
        description="Calculate distances and bond angles between Heme Fe=O (defined via Fe-S_Cys axis) and ligand heavy atom-H vectors."
    )
    parser.add_argument(
        "-t", "--topology", required=True, 
        help="Path to receptor/protein topology file (.pdb, .prmtop, etc.)"
    )
    parser.add_argument(
        "-l", "--ligand_file", required=True,
        help="Path to the ligand structure file (.pdb or .mol2) to read coordinates and auto-detect residue name."
    )
    parser.add_argument(
        "-x", "--trajectory", default=None, 
        help="Path to trajectory file (.nc, .dcd, .xtc, etc.). Optional if analyzing a single static structure."
    )
    parser.add_argument(
        "-m", "--heme", default="HEM", 
        help="Residue name of the Heme group (default: HEM)."
    )
    parser.add_argument(
        "-r", "--resolution", type=float, default=None, 
        help="X-ray resolution in Angstroms. If not provided, script attempts to extract it from the topology PDB header, falling back to 2.0 Å."
    )
    parser.add_argument(
        "-o", "--output", default="cyp_analysis_results.xlsx", 
        help="Output Excel file path (default: cyp_analysis_results.xlsx)."
    )
    return parser.parse_args()

def extract_resolution_from_pdb(pdb_path):
    """Attempts to find X-RAY RESOLUTION in a PDB file header."""
    if not pdb_path.lower().endswith('.pdb'):
        return None
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if "REMARK   2 RESOLUTION." in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            return float(parts[3])
                        except ValueError:
                            pass
    except Exception:
        pass
    return None

def auto_detect_ligand_resname(ligand_file_path):
    """Reads a ligand PDB or MOL2 file to automatically extract the residue name."""
    if not os.path.exists(ligand_file_path):
        raise FileNotFoundError(f"Ligand file not found: {ligand_file_path}")
        
    print(f"--> Reading ligand file to extract residue name: '{ligand_file_path}'")
    
    try:
        u_lig = mda.Universe(ligand_file_path)
        resnames = list(set(u_lig.atoms.resnames))
        ignore_res = ['WAT', 'HOH', 'SOL', 'TIP3']
        valid_resnames = [r for r in resnames if r not in ignore_res]
        
        if valid_resnames:
            detected_name = valid_resnames[0]
            print(f"--> Auto-detected ligand residue name: '{detected_name}'")
            return detected_name
    except Exception as e:
        print(f"Warning parsing residue name via MDAnalysis: {e}")

    if ligand_file_path.lower().endswith('.pdb'):
        with open(ligand_file_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    resname = line[17:20].strip()
                    if resname:
                        print(f"--> Auto-detected ligand residue name (fallback parser): '{resname}'")
                        return resname
                        
    raise ValueError(f"Could not automatically detect a valid residue name from '{ligand_file_path}'.")

def calculate_angle(v1, v2):
    """Calculate the angle in degrees between two vectors."""
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    if norm_v1 == 0 or norm_v2 == 0:
        return 0.0
    cos_angle = dot_product / (norm_v1 * norm_v2)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))

def add_hydrogens_rdkit(ligand_selection):
    """Converts MDAnalysis ligand selection to RDKit, adds hydrogens if missing, and returns the molecule."""
    temp_pdb = "_temp_ligand.pdb"
    ligand_selection.atoms.write(temp_pdb)
    
    mol = Chem.MolFromPDBFile(temp_pdb, removeHs=False)
    if mol is None:
        if os.path.exists(temp_pdb): os.remove(temp_pdb)
        raise ValueError("RDKit failed to read the ligand from the generated PDB snippet.")
        
    has_h = any(atom.GetSymbol() == 'H' for atom in mol.GetAtoms())
    if not has_h:
        print("   [Step 3] No hydrogens found on ligand. Adding hydrogens with RDKit and minimizing...")
        mol = Chem.AddHs(mol, addCoords=True)
        AllChem.UFFOptimizeMolecule(mol)
    else:
        print("   [Step 3] Ligand already contains hydrogens.")
        
    if os.path.exists(temp_pdb): 
        os.remove(temp_pdb)
    return mol

def process_frame(u, ligand_resname, heme_resname, r_resolution):
    """Processes a single structure frame for distance and angular alignment criteria."""
    results = []
    
    # Target Selection for Heme Iron with explicit user fallbacks
    fe_atoms = u.select_atoms(f"resname {heme_resname} and (name FE or name FE1 or name FE2 or name FE3)")
    if len(fe_atoms) == 0:
        fe_atoms = u.select_atoms("(resname HM1 or resname FE1 or resname HEM) and (name FE1 or name FE or name FE2)")
    if len(fe_atoms) == 0:
        fe_atoms = u.select_atoms("name FE or name FE1 or name FE2 or name FE3 or element FE")
        
    if len(fe_atoms) == 0:
        print(f"Error: Heme Iron (FE) could not be identified anywhere in the system.")
        return results
    
    fe_atom = fe_atoms[0]
    print(f"   [Step 1] Found Heme Iron: Atom '{fe_atom.name}' in Residue '{fe_atom.resname}' (Index: {fe_atom.index})")

    fe_pos = fe_atom.position
    
    # Multi-stage robust Cysteine Sulfur spatial selection hierarchy to prevent MDAnalysis syntax crashes
    # 1. Target CM1 residue sulfur first within coordination distance
    s_atoms = u.select_atoms(f"resname CM1 and name SG and point {fe_pos[0]} {fe_pos[1]} {fe_pos[2]} 3.5")
    
    # 2. Target standard CYS residue sulfur within coordination distance
    if len(s_atoms) == 0:
        s_atoms = u.select_atoms(f"resname CYS and name SG and point {fe_pos[0]} {fe_pos[1]} {fe_pos[2]} 3.5")
        
    # 3. Target any generic Sulfur atom located in close proximity (3.5 Angstrom window)
    if len(s_atoms) == 0:
        s_atoms = u.select_atoms(f"(name SG or name S or element S) and point {fe_pos[0]} {fe_pos[1]} {fe_pos[2]} 3.5")
        
    if len(s_atoms) == 0:
        raise ValueError(f"CRITICAL ERROR: Coordinated Cysteine Sulfur (SG/S) not discovered within 3.5 Å of Heme Iron. Cannot calculate Fe=O axis trajectory.")
    
    s_atom = s_atoms[0]
    s_pos = s_atom.position
    print(f"   [Context] Found Coordinated Cysteine Sulfur: Atom '{s_atom.name}' in Residue '{s_atom.resname}' {s_atom.resid} (Dist: {np.linalg.norm(fe_pos - s_pos):.2f} Å)")

    # Identify Ligand
    ligand_atoms = u.select_atoms(f"resname {ligand_resname}")
    if len(ligand_atoms) == 0:
        print(f"Error: Ligand with resname '{ligand_resname}' not found in the combined system.")
        return results
    
    # 2. Get heavy atoms and filter distances
    print("   [Step 2] Processing ligand heavy atoms and measuring distances from Fe...")
    rdkit_mol = add_hydrogens_rdkit(ligand_atoms)
    
    # 4. Generate Fe=O vector along the S_cys -> Fe trajectory axis
    cys_s_to_fe_axis = fe_pos - s_pos
    cys_s_to_fe_axis_normalized = cys_s_to_fe_axis / np.linalg.norm(cys_s_to_fe_axis)
    
    # Place Oxo Oxygen along the axis vector line 
    o_pos = fe_pos + (cys_s_to_fe_axis_normalized * 1.63)
    fe_o_vector = o_pos - fe_pos
    print("   [Step 4] Defined Fe=O vector strictly using the coordinated S_cys-Fe linear bond axis.")

    conf = rdkit_mol.GetConformer()
    min_dist = 6.0 - r_resolution
    max_dist = 6.0 + r_resolution
    print(f"   [Distance Filter Range]: {min_dist:.2f} Å to {max_dist:.2f} Å (6 +- {r_resolution} Å)")

    # 5, 6, 7. Vector processing loop
    print("   [Step 5 & 6] Calculating Atom-H vectors and evaluating angles...")
    for atom in rdkit_mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
            
        atom_idx = atom.GetIdx()
        atom_pos = np.array(conf.GetAtomPosition(atom_idx))
        
        dist_from_fe = np.linalg.norm(atom_pos - fe_pos)
        passes_distance = min_dist <= dist_from_fe <= max_dist
        
        neighbors = atom.GetNeighbors()
        hydrogens = [n for n in neighbors if n.GetSymbol() == 'H']
        has_hydrogen = len(hydrogens) > 0
        
        angles = []
        if has_hydrogen:
            for h_atom in hydrogens:
                h_pos = np.array(conf.GetAtomPosition(h_atom.GetIdx()))
                atom_h_vector = h_pos - atom_pos
                angle = calculate_angle(fe_o_vector, atom_h_vector)
                angles.append(angle)
        
        passes_angle = False
        if has_hydrogen:
            passes_angle = any(100.0 <= ang <= 145.0 for ang in angles)
            
        results.append({
            "Atom_Name": f"{atom.GetSymbol()}{atom.GetIdx()+1}",
            "Distance_to_Fe": round(dist_from_fe, 3),
            "Passes_Distance": passes_distance,
            "Has_H": has_hydrogen,
            "H_Count": len(hydrogens),
            "Angles": [round(a, 2) for a in angles],
            "Passes_Angle_Criteria": passes_angle,
            "Passes_Both": passes_distance and passes_angle
        })
        
    return results

def main():
    args = parse_arguments()
    
    # Auto-detect residue name from the standalone ligand file
    ligand_resname = auto_detect_ligand_resname(args.ligand_file)
    
    resolution = args.resolution
    if resolution is None:
        resolution = extract_resolution_from_pdb(args.topology)
        if resolution is not None:
            print(f"--> Extracted resolution from topology header: {resolution} Å")
        else:
            resolution = 2.0
            print(f"--> No resolution explicitly specified or found in header. Defaulting to: {resolution} Å")
    else:
        print(f"--> Using user-specified resolution: {resolution} Å")

    all_frames_data = []

    # Check operating mode (Static/Docked vs. Trajectory Loop)
    if args.trajectory is None:
        print(f"--> No trajectory file supplied. Operating in Single Structure/Docked mode.")
        print(f"--> Merging receptor structure '{args.topology}' and ligand structure '{args.ligand_file}'...")
        
        u_receptor = mda.Universe(args.topology)
        u_ligand = mda.Universe(args.ligand_file)
        u = mda.Merge(u_receptor.atoms, u_ligand.atoms)
        
        print(f"\n[Processing] Static Structure Layout...")
        frame_results = process_frame(u, ligand_resname, args.heme, resolution)
        for res in frame_results:
            res["Frame"] = 0
            all_frames_data.append(res)
    else:
        print(f"--> Trajectory file supplied: '{args.trajectory}'. Operating in Trajectory mode.")
        print(f"--> Loading complex system files into MDAnalysis...")
        u = mda.Universe(args.topology, args.trajectory, format="NCDF")
        
        for ts in u.trajectory:
            frame_idx = ts.frame
            print(f"\n[Processing] Trajectory Frame {frame_idx}...")
            frame_results = process_frame(u, ligand_resname, args.heme, resolution)
            
            for res in frame_results:
                res["Frame"] = frame_idx
                all_frames_data.append(res)

    if not all_frames_data:
        print("\n[ERROR] No structural or geometric data was compiled. Please verify your file configurations.")
        return

    print("\n[Step 8] Tabulating records and exporting data matrix...")
    df = pd.DataFrame(all_frames_data)
    cols = ["Frame", "Atom_Name", "Distance_to_Fe", "Passes_Distance", "Has_H", "H_Count", "Angles", "Passes_Angle_Criteria", "Passes_Both"]
    df = df[cols]
    
    df.to_excel(args.output, index=False)
    print(f"[SUCCESS] Output written to: {args.output}")
    
    passing_atoms = df[df["Passes_Both"] == True]
    
    print("\n" + "="*60)
    print("      METABOLIC ACTIVE ATOMS (PASS DISTANCE & ANGLE)")
    print("="*60)
    if not passing_atoms.empty:
        print(passing_atoms[["Frame", "Atom_Name", "Distance_to_Fe", "Angles"]].to_string(index=False))
    else:
        print("No atoms passed both conditions.")
    print("="*60)

if __name__ == "__main__":
    main()