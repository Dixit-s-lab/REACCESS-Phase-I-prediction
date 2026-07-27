import os
import sys
import argparse
import subprocess
import re
import shutil
import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
from PIL import Image

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

# Pre-emptively register binary NetCDF components
try:
    import scipy
    import netcdf4
except ImportError:
    pass

try:
    import MDAnalysis as mda
except ImportError:
    print("[ERROR] MDAnalysis is required. Please install it via: pip install MDAnalysis")
    sys.exit(1)

def parse_arguments():
    """Parses command line inputs matching the core engine signature."""
    parser = argparse.ArgumentParser(
        description="REACCESS Core Analytics Engine - Continuous Time Resolution Suite."
    )
    parser.add_argument("-t", "--topology", required=True, help="Path to receptor/protein topology file (.prmtop / .pdb).")
    parser.add_argument("-l", "--ligand_file", required=True, help="Path to ligand reference file (.pdb).")
    parser.add_argument("-x", "--trajectory", default=None, help="Path to trajectory file (.nc / .dcd / .xtc).")
    parser.add_argument("-m", "--heme", default="HEM", help="Residue name of Heme group.")
    parser.add_argument("-r", "--resolution", type=float, default=2.5, help="Geometric filter cut-off scaling resolution.")
    parser.add_argument("-o", "--output", default="cyp_analysis_results.xlsx", help="Output primary Excel path.")
    parser.add_argument("-s", "--smiles", default=None, help="Optional SMILES string of the ligand to guarantee correct bond orders.")
    return parser.parse_known_args()

def find_ligand_resname_from_universe(u):
    """Scans the active Universe topology to dynamically isolate the custom ligand residue code."""
    known_solvents_and_proteins = {
        'WAT', 'HOH', 'SOL', 'TIP3', 'TIP4', 'Cl-', 'Na+', 'K+', 'BR', 'CA', 'MG', 'ZN',
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'HIE', 'HID', 'HIP',
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HEM', 'HEME', 'FE1', 'CM1', 'HM1'
    }
    for resname in set(u.residues.resnames):
        if resname.strip() not in known_solvents_and_proteins and not resname.startswith('N'):
            return resname
    return "UNK"

def calculate_angle(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    if norm_v1 == 0 or norm_v2 == 0: return 0.0
    return np.degrees(np.arccos(np.clip(dot_product / (norm_v1 * norm_v2), -1.0, 1.0)))

def write_mopac_input(filename, mol, charge=0, multiplicity=1):
    """Writes a MOPAC geometry optimization input file using PM7."""
    keywords = f"PM7 charge={charge} "
    keywords += "UHF doublet PRECISE OPT" if multiplicity == 2 else "singlet PRECISE OPT"

    conf = mol.GetConformer()
    with open(filename, 'w') as f:
        f.write(f"{keywords}\n")
        f.write(f"MOPAC Automated Generation: {filename}\n\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            f.write(f"{symbol:<2} {pos.x:12.6f} 1 {pos.y:12.6f} 1 {pos.z:12.6f} 1\n")

def run_mopac(filename):
    """Invokes MOPAC and parses the final HEAT OF FORMATION (kcal/mol)."""
    base = os.path.splitext(filename)[0]
    out_file = f"{base}.out"
    if os.path.exists(out_file): os.remove(out_file)
    
    try:
        subprocess.run(["mopac", filename], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        print(f"Error running MOPAC on {filename}: {e}")
        return None
        
    if not os.path.exists(out_file): return None
        
    hof = None
    hof_re = re.compile(r"FINAL\s+HEAT\s+OF\s+FORMATION\s+=\s+(-?\d+\.\d+)\s+KCAL/MOL")
    with open(out_file, 'r') as f:
        for line in f:
            match = hof_re.search(line)
            if match:
                hof = float(match.group(1))
                break
    return hof

def main():
    args, unknown_args = parse_arguments()
    
    abs_ligand_path = os.path.abspath(args.ligand_file)
    abs_topology_path = os.path.abspath(args.topology)
    abs_trajectory_path = os.path.abspath(args.trajectory) if args.trajectory else None
    
    ligand_filename_base = os.path.splitext(os.path.basename(args.ligand_file))[0]
    output_dir = os.path.abspath(ligand_filename_base)
    os.makedirs(output_dir, exist_ok=True)
    
    local_excel_name = os.path.basename(args.output)

    print("--> [Step 1] Initializing structural universe coordinate networks...")
    if abs_trajectory_path:
        print(f"--> Trajectory file detected. Building active coordination grids...")
        fmt_tag = "NCDF" if abs_trajectory_path.lower().endswith('.nc') else None
        
        try:
            u = mda.Universe(abs_topology_path, abs_trajectory_path, format=fmt_tag)
        except ValueError as err:
            if "Supplied n_atoms" in str(err):
                print("--> [Interceptor] Mismatch caught. Trajectory has explicit waters. Rescaling memory slicing views...")
                total_traj_atoms = int(str(err).split("from ncdf (")[1].split(")")[0])
                u_topo = mda.Universe(abs_topology_path)
                num_topo_atoms = len(u_topo.atoms)
                
                u_dummy = mda.Universe.empty(total_traj_atoms, n_residues=total_traj_atoms, atom_resindex=np.arange(total_traj_atoms))
                u_dummy.load_new(abs_trajectory_path, format=fmt_tag)
                
                u = mda.Merge(u_topo.atoms)
                u.load_new(u_dummy.select_atoms(f"index 0:{num_topo_atoms-1}").positions)
                print(f"--> [Success] In-memory trajectory views sliced cleanly to {len(u.atoms)} atoms matching topology.")
            else:
                raise err
    else:
        print(f"--> Operating in Static Single Docked Pose mode. Merging molecular layers...")
        u_receptor = mda.Universe(abs_topology_path)
        u_ligand = mda.Universe(abs_ligand_path)
        u = mda.Merge(u_receptor.atoms, u_ligand.atoms)
        
    detected_resname = find_ligand_resname_from_universe(u)
    print(f"--> Dynamic Topology Core: Found ligand residue tracking code: '{detected_resname}'")

    fe_atoms = u.select_atoms(f"resname {args.heme} or resname FE1 or name FE or name FE1 or element FE")
    s_atoms = u.select_atoms(f"(resname CM1 or resname CYS or resname CYM) and (name SG or name S or element S)")
    ligand_atoms = u.select_atoms(f"resname {detected_resname}")

    if len(fe_atoms) == 0 or len(ligand_atoms) == 0:
        raise ValueError("Critical identification error: Catalyst Heme Iron or Ligand missing from coordinates.")

    fe_atom = fe_atoms[0]
    s_atom = s_atoms[0] if len(s_atoms) > 0 else None

    print("--> [Step 2] Setting up molecular templates via RDKit...")
    mol_pdb = Chem.MolFromPDBFile(abs_ligand_path, removeHs=False)
    if mol_pdb is None:
        raise ValueError(f"RDKit failed to read the ligand PDB file: {abs_ligand_path}")

    if args.smiles:
        template_mol = Chem.MolFromSmiles(args.smiles)
        if template_mol:
            template_mol = Chem.AddHs(template_mol)
            try:
                mol_receptor = Chem.AllChem.AssignBondOrdersFromTemplate(template_mol, mol_pdb)
            except Exception as e:
                print(f"--> [Warning] AssignBondOrdersFromTemplate failed: {e}. Falling back to standard PDB format.")
                mol_receptor = mol_pdb
        else:
            mol_receptor = mol_pdb
    else:
        mol_receptor = mol_pdb
        try: Chem.DetermineBondsOrderDetail(mol_receptor)
        except: pass

    conf_pdb = mol_pdb.GetConformer()
    conf_rec = mol_receptor.GetConformer()
    
    pdb_name_to_rdkit_idx = {}
    rdkit_idx_to_pdb_name = {}
    
    for r_atom in mol_receptor.GetAtoms():
        r_idx = r_atom.GetIdx()
        r_pos = conf_rec.GetAtomPosition(r_idx)
        
        best_p_atom = None
        min_dist = 999.9
        for p_atom in mol_pdb.GetAtoms():
            p_pos = conf_pdb.GetAtomPosition(p_atom.GetIdx())
            dist = (r_pos.x - p_pos.x)**2 + (r_pos.y - p_pos.y)**2 + (r_pos.z - p_pos.z)**2
            if dist < min_dist:
                min_dist = dist
                best_p_atom = p_atom
                
        if best_p_atom is not None and min_dist < 0.01:
            p_info = best_p_atom.GetPDBResidueInfo()
            name_stripped = p_info.GetName().strip() if p_info else f"{best_p_atom.GetSymbol()}{best_p_atom.GetIdx()+1}"
            
            pdb_name_to_rdkit_idx[name_stripped] = r_idx
            rdkit_idx_to_pdb_name[r_idx] = name_stripped
            if p_info:
                r_atom.SetMonomerInfo(p_info)

    geometric_records = []
    min_d, max_d = 6.0 - args.resolution, 6.0 + args.resolution
    frames_loop = u.trajectory if abs_trajectory_path else [u.trajectory[0]]
    total_trajectory_frames = len(frames_loop)
    print(f"--> [Status] Screening {total_trajectory_frames} configurations against active geometric cut-offs...")

    for f_idx, ts in enumerate(frames_loop):
        fe_pos = fe_atom.position
        cys_axis = fe_pos - s_atom.position if s_atom is not None else np.array([0.0, 0.0, 1.0])
        fe_o_vec = (cys_axis / np.linalg.norm(cys_axis)) * 1.63

        for atom in ligand_atoms:
            if atom.element in ['H', 'W'] or atom.name.strip().startswith('H'): continue
            dist_fe = np.linalg.norm(atom.position - fe_pos)
            pass_d = min_d <= dist_fe <= max_d

            bonded_h = u.select_atoms(f"(element H or name H*) and around 1.35 index {atom.index}")
            has_h = len(bonded_h) > 0
            angles = []
            
            if has_h and atom.element == 'C':
                for h_atom in bonded_h:
                    angles.append(calculate_angle(fe_o_vec, h_atom.position - atom.position))
            
            pass_a = any(100.0 <= ang <= 145.0 for ang in angles) if (has_h and atom.element == 'C') else False
            rd_idx = pdb_name_to_rdkit_idx.get(atom.name.strip(), -1)

            if rd_idx >= 0:
                valid_ang_list = [a for a in angles if 100.0 <= a <= 145.0]
                chosen_angle = valid_ang_list[0] if valid_ang_list else (angles[0] if angles else 120.0)
                pass_both = (pass_d and pass_a) if atom.element == 'C' else pass_d
                
                geometric_records.append({
                    "Frame": f_idx, 
                    "Atom_Name": atom.name.strip(), 
                    "Distance_to_Fe": round(dist_fe, 3),
                    "Passes_Distance": pass_d, 
                    "Has_H": has_h, 
                    "H_Count": len(bonded_h),
                    "Angles": [round(a, 2) for a in angles], 
                    "Parsed_Angle": chosen_angle if atom.element == 'C' else 0.0,
                    "Passes_Angle_Criteria": pass_a if atom.element == 'C' else True, 
                    "Passes_Both": pass_both,
                    "RDKit_Idx": int(rd_idx), 
                    "Element": atom.element
                })

    df_geo = pd.DataFrame(geometric_records)

    # REFACTORED LOGIC: Time series evaluation checking for continuous consecutive frame windows
    if abs_trajectory_path:
        frame_dt = getattr(u.trajectory, 'dt', 1.0)
        if frame_dt <= 0: frame_dt = 1.0
        
        min_required_frames = int(np.ceil(100.0 / frame_dt))
        print(f"--> [Occupancy Core] Trajectory storage stride detected: {frame_dt} ps/frame.")
        print(f"--> [Occupancy Core] Enforcing strict continuous duration: atoms must pass for >= {min_required_frames} CONSECUTIVE frames.")
        
        # Calculate maximum consecutive frames passing criteria for Carbons
        active_carbon_ids = []
        df_carbons_geo = df_geo[df_geo['Element'] == 'C']
        if not df_carbons_geo.empty:
            for rd_idx, group in df_carbons_geo.groupby('RDKit_Idx'):
                s = group['Passes_Both']
                max_consec = (s.groupby((~s).cumsum()).sum()).max()
                if pd.isna(max_consec): max_consec = 0
                if max_consec >= min_required_frames:
                    active_carbon_ids.append(rd_idx)

        # Calculate maximum consecutive frames passing distance criteria for Heteroatoms
        active_hetero_ids = []
        df_hetero_geo = df_geo[df_geo['Element'].isin(['N', 'S'])]
        if not df_hetero_geo.empty:
            for rd_idx, group in df_hetero_geo.groupby('RDKit_Idx'):
                s = group['Passes_Distance']
                max_consec = (s.groupby((~s).cumsum()).sum()).max()
                if pd.isna(max_consec): max_consec = 0
                if max_consec >= min_required_frames:
                    active_hetero_ids.append(rd_idx)
    else:
        print("--> [Occupancy Core] Static Pose Mode: Bypassing time thresholds. Retaining all items passing geometric constraints (Min: 1 Frame).")
        df_carbons_geo = df_geo[df_geo['Element'] == 'C']
        if not df_carbons_geo.empty:
            active_frames_c = df_carbons_geo.groupby('RDKit_Idx')['Passes_Both'].sum()
            active_carbon_ids = list(active_frames_c[active_frames_c >= 1].index)
        else:
            active_carbon_ids = []

        df_hetero_geo = df_geo[df_geo['Element'].isin(['N', 'S'])]
        if not df_hetero_geo.empty:
            active_frames_h = df_hetero_geo.groupby('RDKit_Idx')['Passes_Distance'].sum()
            active_hetero_ids = list(active_frames_h[active_frames_h >= 1].index)
        else:
            active_hetero_ids = []

    df_reactive_carbons = df_geo[(df_geo['Element'] == 'C') & (df_geo['RDKit_Idx'].isin(active_carbon_ids))]
    reactive_carbons_group = df_reactive_carbons.groupby('RDKit_Idx')
    total_unique_carbons = reactive_carbons_group.ngroups
    
    print(f"--> [Status] Isolated {total_unique_carbons} unique carbon centers passing selection criteria for quantum analysis.")
    print(f"--> [Status] Isolated {len(active_hetero_ids)} unique heteroatom centers (N/S) passing distance metrics criteria.")

    df_geo.to_excel(os.path.join(output_dir, local_excel_name), index=False)
    os.chdir(output_dir)

    bde_results = []
    summary_data = []

    # Compute baseline MOPAC values exclusively for active carbon targets
    if total_unique_carbons > 0:
        print("--> [Step 3] Computing baseline MOPAC values for radical fragment parameters...")
        h_mol = Chem.MolFromSmiles("[H]")
        h_conf = Chem.Conformer(1)
        h_conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
        h_mol.AddConformer(h_conf)
        
        write_mopac_input("H_rad.mop", h_mol, charge=0, multiplicity=2)
        h_rad_energy = run_mopac("H_rad.mop")
        if os.path.exists("H_rad.out"): os.remove("H_rad.out")
            
        write_mopac_input("parent.mop", mol_receptor, charge=0, multiplicity=1)
        parent_energy = run_mopac("parent.mop")
        if os.path.exists("parent.out"): os.remove("parent.out")

        print("--> [MOPAC Engine] Launching parallel radical optimizations loop...")
        for progress_count, (rd_idx, group) in enumerate(reactive_carbons_group, 1):
            rd_idx_int = int(rd_idx)
            if rd_idx_int < 0: continue
            atom_name_pdb = rdkit_idx_to_pdb_name.get(rd_idx_int, f"C{rd_idx_int+1}")
            print(f"    [MOPAC Optimization {progress_count}/{total_unique_carbons}] Simulating radical at: {atom_name_pdb}...")
            
            rd_atom = mol_receptor.GetAtomWithIdx(rd_idx_int)
            hydrogens = [n for n in rd_atom.GetNeighbors() if n.GetSymbol() == 'H']
            
            if not hydrogens: 
                print(f"    [Skipped] Atom {atom_name_pdb} possesses no explicit hydrogen setups.")
                continue
            h_idx_global = hydrogens[0].GetIdx()
            
            editable_mol = Chem.RWMol(mol_receptor)
            editable_mol.RemoveAtom(h_idx_global)
            radical_mol = editable_mol.GetMol()
            
            rad_base = f"radical_{atom_name_pdb}"
            rad_mop_name = f"{rad_base}.mop"
            
            write_mopac_input(rad_mop_name, radical_mol, charge=0, multiplicity=2)
            radical_energy = run_mopac(rad_mop_name)
            
            if os.path.exists(f"{rad_base}.out"): os.remove(f"{rad_base}.out")
            if radical_energy is None: 
                print(f"    [WARNING] MOPAC optimization failed to converge for atom {atom_name_pdb}.")
                continue
                
            bde = (radical_energy + h_rad_energy) - parent_energy
            
            for _, row in group.iterrows():
                bde_results.append({
                    "Frame": row['Frame'],
                    "Atom_Name": atom_name_pdb, 
                    "Distance_to_Fe": row['Distance_to_Fe'],
                    "Target_Angle": row['Parsed_Angle'],
                    "Calculated_BDE": round(bde, 3),
                    "Passes_Distance": row['Passes_Distance'],
                    "Passes_Angle_Criteria": row['Passes_Angle_Criteria'],
                    "RDKit_Idx": rd_idx_int
                })

    # Tabulate Carbon outputs
    if bde_results:
        bde_df = pd.DataFrame(bde_results)
        bde_df.to_excel(f"{local_excel_name.replace('.xlsx', '')}_BDE_summary.xlsx", index=False)
        bde_df['Passes_BDE'] = bde_df['Calculated_BDE'].between(70.0, 95.0)
        bde_df['Is_Active_Site'] = bde_df['Passes_Distance'] & bde_df['Passes_Angle_Criteria'] & bde_df['Passes_BDE']

        for idx_atom, group in bde_df.groupby('RDKit_Idx'):
            active_frames_count = group[group['Is_Active_Site'] == True]['Frame'].nunique()
            activity_percentage = (active_frames_count / total_trajectory_frames) * 100.0
            atom_name_pdb = rdkit_idx_to_pdb_name.get(int(idx_atom))
            
            summary_data.append({
                "Atom_Name": atom_name_pdb,
                "Mean_Distance": round(group['Distance_to_Fe'].mean(), 2),
                "Mean_Angle": round(group['Target_Angle'].mean(), 2),
                "Calculated_BDE": round(group['Calculated_BDE'].mean(), 2),
                "Active_Frames": active_frames_count,
                "Total_Frames": total_trajectory_frames,
                "Activity_Frequency_%": round(activity_percentage, 1)
            })
    else:
        bde_df = pd.DataFrame()
        pd.DataFrame().to_excel(f"{local_excel_name.replace('.xlsx', '')}_BDE_summary.xlsx", index=False)

    # Append active Heteroatoms directly to summary data
    for idx_atom in active_hetero_ids:
        group = df_geo[df_geo['RDKit_Idx'] == idx_atom]
        active_frames_count = group[group['Passes_Distance'] == True]['Frame'].nunique()
        activity_percentage = (active_frames_count / total_trajectory_frames) * 100.0
        atom_name_pdb = rdkit_idx_to_pdb_name.get(int(idx_atom))
        
        summary_data.append({
            "Atom_Name": atom_name_pdb,
            "Mean_Distance": round(group['Distance_to_Fe'].mean(), 2),
            "Mean_Angle": 0.0,
            "Calculated_BDE": 0.0,
            "Active_Frames": active_frames_count,
            "Total_Frames": total_trajectory_frames,
            "Activity_Frequency_%": round(activity_percentage, 1)
        })

    if not summary_data:
        print("[WARNING] Zero active sites met the tracking parameters threshold limits.")
        return

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_excel("METABOLIC_FINAL_SUMMARY.xlsx", index=False)

    # Generate all 1D and 3D charts inside the active workspace folder
    if not bde_df.empty:
        print("--> [Graphics Core] Generating 1D and 3D energy landscape profiles...")
        sns.set_theme(style="ticks")
        dist_low, dist_high = 6.0 - args.resolution, 6.0 + args.resolution
        ang_low, ang_high = 100.0, 145.0
        bde_low, bde_high = 70.0, 95.0

        unique_atoms = bde_df['Atom_Name'].unique()
        colors_palette = sns.color_palette("turbo", len(unique_atoms))
        atom_color_dict = dict(zip(unique_atoms, colors_palette))

        # Plot 1: Distance vs BDE
        fig1, ax1 = plt.subplots(figsize=(7.5, 4.5))
        sns.scatterplot(data=bde_df, x='Distance_to_Fe', y='Calculated_BDE', hue='Atom_Name', palette=atom_color_dict, s=40, edgecolor='black', alpha=0.8, ax=ax1)
        rect1 = patches.Rectangle((dist_low, bde_low), dist_high - dist_low, bde_high - bde_low, linewidth=1.5, edgecolor='red', facecolor='red', alpha=0.08, linestyle='--')
        ax1.add_patch(rect1)
        ax1.set_title('Correlation: Heme Fe Distance vs Carbon C-H BDE', fontsize=11, fontweight='bold')
        ax1.set_xlabel('Distance to Iron (Å)', fontsize=10)
        ax1.set_ylabel('C-H Bond Dissociation Energy (kcal/mol)', fontsize=10)
        ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title="Atoms")
        fig1.tight_layout()
        fig1.savefig('plot_1_distance_vs_BDE.png', dpi=300)
        plt.close(fig1)

        # Plot 2: Angle vs BDE
        fig2, ax2 = plt.subplots(figsize=(7.5, 4.5))
        sns.scatterplot(data=bde_df, x='Target_Angle', y='Calculated_BDE', hue='Atom_Name', palette=atom_color_dict, s=40, edgecolor='black', alpha=0.8, ax=ax2)
        rect2 = patches.Rectangle((ang_low, bde_low), ang_high - ang_low, bde_high - bde_low, linewidth=1.5, edgecolor='blue', facecolor='blue', alpha=0.08, linestyle='--')
        ax2.add_patch(rect2)
        ax2.set_title('Correlation: Fe=O Axis Alignment Angle vs C-H BDE', fontsize=11, fontweight='bold')
        ax2.set_xlabel('Fe=O ··· H-C Alignment Angle (Degrees)', fontsize=10)
        ax2.set_ylabel('C-H Bond Dissociation Energy (kcal/mol)', fontsize=10)
        ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title="Atoms")
        fig2.tight_layout()
        fig2.savefig('plot_2_angle_vs_BDE.png', dpi=300)
        plt.close(fig2)

        # Plot 3: 3D Landscape Plot
        fig3 = plt.figure(figsize=(9, 6.5))
        ax3 = fig3.add_subplot(111, projection='3d')
        for atom_name, group_3d in bde_df.groupby('Atom_Name'):
            ax3.scatter(group_3d['Distance_to_Fe'], group_3d['Target_Angle'], group_3d['Calculated_BDE'], label=atom_name, color=atom_color_dict[atom_name], s=50, edgecolor='k', alpha=0.75)
        corners = np.array([
            [dist_low, ang_low, bde_low], [dist_high, ang_low, bde_low], [dist_high, ang_high, bde_low], [dist_low, ang_high, bde_low],
            [dist_low, ang_low, bde_high], [dist_high, ang_low, bde_high], [dist_high, ang_high, bde_high], [dist_low, ang_high, bde_high]
        ])
        faces = [
            [corners[0], corners[1], corners[2], corners[3]], [corners[4], corners[5], corners[6], corners[7]], 
            [corners[0], corners[1], corners[5], corners[4]], [corners[2], corners[3], corners[7], corners[6]], 
            [corners[0], corners[3], corners[7], corners[4]], [corners[1], corners[2], corners[6], corners[5]]  
        ]
        cuboid = Poly3DCollection(faces, facecolors='green', linewidths=0.5, edgecolors='darkgreen', alpha=0.05)
        ax3.add_collection3d(cuboid)
        ax3.set_title('3D Geometric Profile Landscape: Distance & Angle vs C-H BDE', fontsize=12, fontweight='bold', pad=15)
        ax3.set_xlabel('Distance to Fe (Å)', labelpad=10)
        ax3.set_ylabel('Alignment Angle (°)', labelpad=10)
        ax3.set_zlabel('BDE (kcal/mol)', labelpad=10)
        ax3.legend(bbox_to_anchor=(1.12, 0.9), loc='upper left', borderaxespad=0, title="Atoms")
        fig3.tight_layout()
        fig3.savefig('plot_3_3D_landscape_BDE.png', dpi=300)
        plt.close(fig3)

    # Generating pristine 2D ACS structure map layouts
    print("--> [Graphics Core] Generating pristine 2D ACS structure map layouts...")
    df_active_atoms = summary_df[summary_df['Activity_Frequency_%'] > 0.0]
    active_atoms_list = df_active_atoms['Atom_Name'].tolist()

    mol_base = Chem.MolFromSmiles(args.smiles) if args.smiles else Chem.RemoveHs(mol_pdb)
    Chem.SanitizeMol(mol_base)
    
    smiles_to_pdb_name = {}
    if args.smiles:
        mol_pdb_heavy = Chem.RemoveHs(mol_pdb)
        mcs_res = rdFMCS.FindMCS([mol_base, mol_pdb_heavy], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareAny)
        if mcs_res.numAtoms > 0:
            mcs_mol = Chem.MolFromSmarts(mcs_res.smartsString)
            m_base = mol_base.GetSubstructMatch(mcs_mol)
            m_pdb = mol_pdb_heavy.GetSubstructMatch(mcs_mol)
            for b_idx, p_idx in zip(m_base, m_pdb):
                p_inf = mol_pdb_heavy.GetAtomWithIdx(p_idx).GetPDBResidueInfo()
                smiles_to_pdb_name[b_idx] = p_inf.GetName().strip() if p_inf else f"C{p_idx+1}"
    else:
        for atom in mol_base.GetAtoms():
            p_inf = atom.GetPDBResidueInfo()
            smiles_to_pdb_name[atom.GetIdx()] = p_inf.GetName().strip() if p_inf else f"C{atom.GetIdx()+1}"

    Chem.RemoveStereochemistry(mol_base)
    mol_base.RemoveAllConformers()
    AllChem.Compute2DCoords(mol_base)

    d2d = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
    draw_opts = d2d.drawOptions()
    draw_opts.addAtomIndices = False
    draw_opts.addStereoAnnotation = False
    draw_opts.lineWidth = 2.6
    draw_opts.bondLength = 30

    highlight_atoms = []
    for atom in mol_base.GetAtoms():
        idx = atom.GetIdx()
        pdb_name = smiles_to_pdb_name.get(idx, f"C{idx+1}")
        if pdb_name in active_atoms_list:
            highlight_atoms.append(idx)
            draw_opts.atomLabels[idx] = pdb_name
        elif atom.GetSymbol() != "C":
            draw_opts.atomLabels[idx] = atom.GetSymbol()
        else:
            draw_opts.atomLabels[idx] = ""

    h_colors = {int(idx): (1.0, 0.0, 0.0) for idx in highlight_atoms}
    h_radii = {int(idx): 0.40 for idx in highlight_atoms}
    
    d2d.DrawMolecule(mol_base, highlightAtoms=highlight_atoms, highlightAtomColors=h_colors, highlightAtomRadii=h_radii, highlightBonds=[])
    d2d.FinishDrawing()

    img = Image.open(io.BytesIO(d2d.GetDrawingText()))
    fig2d, ax2d = plt.subplots(figsize=(8, 8), facecolor='white')
    ax2d.imshow(img)
    ax2d.axis('off')
    ax2d.set_title(f"CYP450 2D Metabolic Structure Profile: {ligand_filename_base}", fontsize=12, fontweight='bold')
    fig2d.savefig(f"{ligand_filename_base}_labeled_ligand_2D.png", dpi=300, bbox_inches='tight')
    plt.close(fig2d)

    print("\n" + "="*80)
    print("               FINAL CYP450 METABOLIC ACTIVITY FREQUENCY SUMMARY")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("="*80 + "\n")

if __name__ == "__main__":
    main()