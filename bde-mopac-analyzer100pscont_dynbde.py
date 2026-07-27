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
    parser = argparse.ArgumentParser(description="REACCESS Core Engine - Dynamic BEP Thresholding Suite.")
    parser.add_argument("-t", "--topology", required=True)
    parser.add_argument("-l", "--ligand_file", required=True)
    parser.add_argument("-x", "--trajectory", default=None)
    parser.add_argument("-m", "--heme", default="HEM")
    parser.add_argument("-o", "--output", default="cyp_analysis_results.xlsx")
    parser.add_argument("-s", "--smiles", default=None)
    parser.add_argument("--min_dist", type=float, default=3.5)
    parser.add_argument("--max_dist", type=float, default=8.5)
    parser.add_argument("--min_angle", type=float, default=100.0)
    parser.add_argument("--max_angle", type=float, default=145.0)
    parser.add_argument("--alpha", type=float, default=0.5)
    return parser.parse_known_args()

def find_ligand_resname_from_universe(u):
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
    norm_v1, norm_v2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if norm_v1 == 0 or norm_v2 == 0: return 0.0
    return np.degrees(np.arccos(np.clip(dot_product / (norm_v1 * norm_v2), -1.0, 1.0)))

def write_mopac_input(filename, mol, charge=0, multiplicity=1):
    keywords = f"PM7 charge={charge} "
    keywords += "UHF doublet PRECISE OPT" if multiplicity == 2 else "singlet PRECISE OPT"
    conf = mol.GetConformer()
    with open(filename, 'w') as f:
        f.write(f"{keywords}\nMOPAC Auto-Generation\n\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol():<2} {pos.x:12.6f} 1 {pos.y:12.6f} 1 {pos.z:12.6f} 1\n")

def run_mopac(filename):
    base = os.path.splitext(filename)[0]
    out_file = f"{base}.out"
    if os.path.exists(out_file): os.remove(out_file)
    try:
        subprocess.run(["mopac", filename], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except:
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
    
    # FIXED: Explicitly mapped string variables boundary definitions
    ligand_filename_base = os.path.splitext(os.path.basename(args.ligand_file))[0]
    output_dir = os.path.abspath(ligand_filename_base)
    os.makedirs(output_dir, exist_ok=True)
    local_excel_name = os.path.basename(args.output)

    print("--> [Step 1] Initializing structural universe coordinate networks...")
    if not abs_trajectory_path and abs_topology_path.lower().endswith(('.prmtop', '.parm7', '.top', '.psf')):
        print("[CRITICAL ERROR] Standalone topology file used in static mode without trajectory matrix inputs.")
        sys.exit(1)

    if abs_trajectory_path:
        fmt_tag = "NCDF" if abs_trajectory_path.lower().endswith('.nc') else None
        try:
            u = mda.Universe(abs_topology_path, abs_trajectory_path, format=fmt_tag)
        except ValueError as err:
            if "Supplied n_atoms" in str(err):
                total_traj_atoms = int(str(err).split("from ncdf (")[1].split(")")[0])
                u_topo = mda.Universe(abs_topology_path)
                num_topo_atoms = len(u_topo.atoms)
                u_dummy = mda.Universe.empty(total_traj_atoms, n_residues=total_traj_atoms, atom_resindex=np.arange(total_traj_atoms))
                u_dummy.load_new(abs_trajectory_path, format=fmt_tag)
                u = mda.Merge(u_topo.atoms)
                u.load_new(u_dummy.select_atoms(f"index 0:{num_topo_atoms-1}").positions)
            else:
                raise err
    else:
        u = mda.Merge(mda.Universe(abs_topology_path).atoms, mda.Universe(abs_ligand_path).atoms)
        
    detected_resname = find_ligand_resname_from_universe(u)
    fe_atoms = u.select_atoms(f"resname {args.heme} or name FE or element FE")
    s_atoms = u.select_atoms(f"(resname CYS or resname CYM) and (name SG or element S)")
    ligand_atoms = u.select_atoms(f"resname {detected_resname}")

    if len(fe_atoms) == 0 or len(ligand_atoms) == 0:
        raise ValueError("Critical identification error: Catalyst Heme Iron or Ligand missing.")

    fe_atom, s_atom = fe_atoms[0], s_atoms[0] if len(s_atoms) > 0 else None

    print("--> [Step 2] Setting up molecular templates via RDKit...")
    mol_pdb = Chem.MolFromPDBFile(abs_ligand_path, removeHs=False)
    if args.smiles:
        template_mol = Chem.AddHs(Chem.MolFromSmiles(args.smiles))
        try: mol_receptor = Chem.AllChem.AssignBondOrdersFromTemplate(template_mol, mol_pdb)
        except: mol_receptor = mol_pdb
    else:
        mol_receptor = mol_pdb
        try: Chem.DetermineBondsOrderDetail(mol_receptor)
        except: pass

    pdb_name_to_rdkit_idx = {}
    rdkit_idx_to_pdb_name = {}
    conf_pdb, conf_rec = mol_pdb.GetConformer(), mol_receptor.GetConformer()
    
    for r_atom in mol_receptor.GetAtoms():
        r_idx = r_atom.GetIdx()
        r_pos = conf_rec.GetAtomPosition(r_idx)
        best_p_atom, min_dist = None, 999.9
        for p_atom in mol_pdb.GetAtoms():
            p_pos = conf_pdb.GetAtomPosition(p_atom.GetIdx())
            dist = (r_pos.x - p_pos.x)**2 + (r_pos.y - p_pos.y)**2 + (r_pos.z - p_pos.z)**2
            if dist < min_dist:
                min_dist, best_p_atom = dist, p_atom
        if best_p_atom is not None and min_dist < 0.01:
            p_info = best_p_atom.GetPDBResidueInfo()
            name_stripped = p_info.GetName().strip() if p_info else f"{best_p_atom.GetSymbol()}{best_p_atom.GetIdx()+1}"
            pdb_name_to_rdkit_idx[name_stripped] = r_idx
            rdkit_idx_to_pdb_name[r_idx] = name_stripped

    geometric_records = []
    frames_loop = u.trajectory if abs_trajectory_path else [u.trajectory[0]]
    total_trajectory_frames = len(frames_loop)
    frame_dt = getattr(u.trajectory, 'dt', 1.0) if abs_trajectory_path else 100.0
    if frame_dt <= 0: frame_dt = 1.0

    print(f"--> Screening {total_trajectory_frames} configurations using custom bounds...")
    for f_idx, ts in enumerate(frames_loop):
        fe_pos = fe_atom.position
        cys_axis = fe_pos - s_atom.position if s_atom is not None else np.array([0.0, 0.0, 1.0])
        fe_o_vec = (cys_axis / np.linalg.norm(cys_axis)) * 1.63

        for atom in ligand_atoms:
            if atom.element in ['H', 'W'] or atom.name.strip().startswith('H'): continue
            dist_fe = np.linalg.norm(atom.position - fe_pos)
            pass_d = args.min_dist <= dist_fe <= args.max_dist

            bonded_h = u.select_atoms(f"(element H or name H*) and around 1.35 index {atom.index}")
            angles = [calculate_angle(fe_o_vec, h_atom.position - atom.position) for h_atom in bonded_h]
            pass_a = any(args.min_angle <= ang <= args.max_angle for ang in angles) if (len(bonded_h) > 0 and atom.element == 'C') else False
            
            rd_idx = pdb_name_to_rdkit_idx.get(atom.name.strip(), -1)
            if rd_idx >= 0:
                valid_ang_list = [a for a in angles if args.min_angle <= a <= args.max_angle]
                chosen_angle = valid_ang_list[0] if valid_ang_list else (angles[0] if angles else 0.0)
                pass_both = (pass_d and pass_a) if atom.element == 'C' else pass_d
                
                geometric_records.append({
                    "Frame": f_idx, "Atom_Name": atom.name.strip(), "Distance_to_Fe": round(dist_fe, 3),
                    "Passes_Distance": pass_d, "Parsed_Angle": chosen_angle,
                    "Passes_Angle_Criteria": pass_a if atom.element == 'C' else True, "Passes_Both": pass_both,
                    "RDKit_Idx": int(rd_idx), "Element": atom.element
                })

    df_geo = pd.DataFrame(geometric_records)
    min_required_frames = int(np.ceil(100.0 / frame_dt)) if abs_trajectory_path else 1

    active_carbon_ids, atom_max_consec_time = [], {}
    df_carbons_geo = df_geo[df_geo['Element'] == 'C']
    if not df_carbons_geo.empty:
        for rd_idx, group in df_carbons_geo.groupby('RDKit_Idx'):
            s = group['Passes_Both']
            max_consec = (s.groupby((~s).cumsum()).sum()).max()
            max_consec = 0 if pd.isna(max_consec) else max_consec
            atom_max_consec_time[rd_idx] = max_consec * frame_dt
            if max_consec >= min_required_frames:
                active_carbon_ids.append(rd_idx)

    active_hetero_ids = []
    df_hetero_geo = df_geo[df_geo['Element'].isin(['N', 'S'])]
    if not df_hetero_geo.empty:
        for rd_idx, group in df_hetero_geo.groupby('RDKit_Idx'):
            s = group['Passes_Distance']
            max_consec = (s.groupby((~s).cumsum()).sum()).max()
            max_consec = 0 if pd.isna(max_consec) else max_consec
            if max_consec >= min_required_frames:
                active_hetero_ids.append(rd_idx)

    df_reactive_carbons = df_geo[(df_geo['Element'] == 'C') & (df_geo['RDKit_Idx'].isin(active_carbon_ids))]
    reactive_carbons_group = df_reactive_carbons.groupby('RDKit_Idx')
    
    df_geo.to_excel(os.path.join(output_dir, local_excel_name), index=False)
    os.chdir(output_dir)

    bde_results, summary_data = [], []
    if reactive_carbons_group.ngroups > 0:
        print("--> [Step 3] Launching baseline radical MOPAC optimizations loop...")
        h_mol = Chem.MolFromSmiles("[H]")
        h_conf = Chem.Conformer(1)
        h_conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
        h_mol.AddConformer(h_conf)
        write_mopac_input("H_rad.mop", h_mol, charge=0, multiplicity=2)
        h_rad_energy = run_mopac("H_rad.mop")
        
        write_mopac_input("parent.mop", mol_receptor, charge=0, multiplicity=1)
        parent_energy = run_mopac("parent.mop")

        for progress_count, (rd_idx, group) in enumerate(reactive_carbons_group, 1):
            rd_idx_int = int(rd_idx)
            atom_name_pdb = rdkit_idx_to_pdb_name.get(rd_idx_int, f"C{rd_idx_int+1}")
            rd_atom = mol_receptor.GetAtomWithIdx(rd_idx_int)
            hydrogens = [n for n in rd_atom.GetNeighbors() if n.GetSymbol() == 'H']
            if not hydrogens: continue
            
            editable_mol = Chem.RWMol(mol_receptor)
            editable_mol.RemoveAtom(hydrogens[0].GetIdx())
            
            rad_base = f"radical_{atom_name_pdb}"
            write_mopac_input(f"{rad_base}.mop", editable_mol.GetMol(), charge=0, multiplicity=2)
            radical_energy = run_mopac(f"{rad_base}.mop")
            if radical_energy is None: continue
                
            bde = (radical_energy + h_rad_energy) - parent_energy
            for _, row in group.iterrows():
                bde_results.append({
                    "Frame": row['Frame'], "Atom_Name": atom_name_pdb, "Distance_to_Fe": row['Distance_to_Fe'],
                    "Target_Angle": row['Parsed_Angle'], "Calculated_BDE": round(bde, 3),
                    "Passes_Distance": row['Passes_Distance'], "Passes_Angle_Criteria": row['Passes_Angle_Criteria'],
                    "RDKit_Idx": rd_idx_int
                })

    if bde_results:
        bde_df = pd.DataFrame(bde_results)
        bde_df.to_excel(f"{local_excel_name.replace('.xlsx', '')}_BDE_summary.xlsx", index=False)
        
        R, T = 1.987e-3, 310.15
        RT = R * T
        
        bde_df['Passes_BDE'] = False
        for i, row in bde_df.iterrows():
            bde_val = row['Calculated_BDE']
            c_time = atom_max_consec_time.get(row['RDKit_Idx'], 0.0)
            if 70.0 <= bde_val <= 95.0:
                bde_df.at[i, 'Passes_BDE'] = True
            elif bde_val > 95.0:
                tau_required = 100.0 * np.exp(args.alpha * (bde_val - 95.0) / RT)
                if c_time >= tau_required:
                    bde_df.at[i, 'Passes_BDE'] = True

        bde_df['Is_Active_Site'] = bde_df['Passes_Distance'] & bde_df['Passes_Angle_Criteria'] & bde_df['Passes_BDE']

        for idx_atom, group in bde_df.groupby('RDKit_Idx'):
            active_frames_count = group[group['Is_Active_Site'] == True]['Frame'].nunique()
            activity_percentage = (active_frames_count / total_trajectory_frames) * 100.0
            if activity_percentage > 0.0 or not abs_trajectory_path:
                summary_data.append({
                    "Atom_Name": rdkit_idx_to_pdb_name.get(int(idx_atom)), "Mean_Distance": round(group['Distance_to_Fe'].mean(), 2),
                    "Mean_Angle": round(group['Target_Angle'].mean(), 2), "Calculated_BDE": round(group['Calculated_BDE'].mean(), 2),
                    "Active_Frames": active_frames_count, "Total_Frames": total_trajectory_frames, "Activity_Frequency_%": round(activity_percentage, 1)
                })
    else:
        bde_df = pd.DataFrame()

    for idx_atom in active_hetero_ids:
        group = df_geo[df_geo['RDKit_Idx'] == idx_atom]
        active_frames_count = group[group['Passes_Distance'] == True]['Frame'].nunique()
        summary_data.append({
            "Atom_Name": rdkit_idx_to_pdb_name.get(int(idx_atom)), "Mean_Distance": round(group['Distance_to_Fe'].mean(), 2),
            "Mean_Angle": 0.0, "Calculated_BDE": 0.0, "Active_Frames": active_frames_count,
            "Total_Frames": total_trajectory_frames, "Activity_Frequency_%": round((active_frames_count / total_trajectory_frames) * 100.0, 1)
        })

    if not summary_data: return
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_excel("METABOLIC_FINAL_SUMMARY.xlsx", index=False)

    if not bde_df.empty:
        print("--> [Graphics Core] Generating 1D and 3D energy landscape profiles...")
        sns.set_theme(style="ticks")
        unique_atoms = bde_df['Atom_Name'].unique()
        atom_color_dict = dict(zip(unique_atoms, sns.color_palette("turbo", len(unique_atoms))))
        max_seen_bde = max(98.0, bde_df['Calculated_BDE'].max() + 2.0)

        fig1, ax1 = plt.subplots(figsize=(7.5, 4.5))
        sns.scatterplot(data=bde_df, x='Distance_to_Fe', y='Calculated_BDE', hue='Atom_Name', palette=atom_color_dict, ax=ax1)
        ax1.add_patch(patches.Rectangle((args.min_dist, 70.0), args.max_dist - args.min_dist, max_seen_bde - 70.0, linewidth=1.5, edgecolor='red', facecolor='red', alpha=0.05, linestyle='--'))
        fig1.tight_layout(); fig1.savefig('plot_1_distance_vs_BDE.png', dpi=300); plt.close(fig1)

        fig2, ax2 = plt.subplots(figsize=(7.5, 4.5))
        sns.scatterplot(data=bde_df, x='Target_Angle', y='Calculated_BDE', hue='Atom_Name', palette=atom_color_dict, ax=ax2)
        ax2.add_patch(patches.Rectangle((args.min_angle, 70.0), args.max_angle - args.min_angle, max_seen_bde - 70.0, linewidth=1.5, edgecolor='blue', facecolor='blue', alpha=0.05, linestyle='--'))
        fig2.tight_layout(); fig2.savefig('plot_2_angle_vs_BDE.png', dpi=300); plt.close(fig2)

        fig3 = plt.figure(figsize=(9, 6.5))
        ax3 = fig3.add_subplot(111, projection='3d')
        for atom_name, group_3d in bde_df.groupby('Atom_Name'):
            ax3.scatter(group_3d['Distance_to_Fe'], group_3d['Target_Angle'], group_3d['Calculated_BDE'], label=atom_name, color=atom_color_dict[atom_name], s=40)
        corners = np.array([
            [args.min_dist, args.min_angle, 70.0], [args.max_dist, args.min_angle, 70.0], [args.max_dist, args.max_angle, 70.0], [args.min_dist, args.max_angle, 70.0],
            [args.min_dist, args.min_angle, max_seen_bde], [args.max_dist, args.min_angle, max_seen_bde], [args.max_dist, args.max_angle, max_seen_bde], [args.min_dist, args.max_angle, max_seen_bde]
        ])
        faces = [[corners[0], corners[1], corners[2], corners[3]], [corners[4], corners[5], corners[6], corners[7]], [corners[0], corners[1], corners[5], corners[4]], [corners[2], corners[3], corners[7], corners[6]], [corners[0], corners[3], corners[7], corners[4]], [corners[1], corners[2], corners[6], corners[5]]]
        ax3.add_collection3d(Poly3DCollection(faces, facecolors='green', linewidths=0.5, edgecolors='darkgreen', alpha=0.03))
        ax3.legend(bbox_to_anchor=(1.12, 0.9), title="Atoms"); fig3.tight_layout(); fig3.savefig('plot_3_3D_landscape_BDE.png', dpi=300); plt.close(fig3)

    print("--> [Graphics Core] Generating pristine 2D ACS structure map layouts...")
    active_atoms_list = summary_df[summary_df['Activity_Frequency_%'] > 0.0]['Atom_Name'].tolist()
    mol_base = Chem.MolFromSmiles(args.smiles) if args.smiles else Chem.RemoveHs(mol_pdb)
    Chem.SanitizeMol(mol_base)
    smiles_to_pdb_name = {}
    if args.smiles:
        mol_pdb_heavy = Chem.RemoveHs(mol_pdb)
        mcs_res = rdFMCS.FindMCS([mol_base, mol_pdb_heavy], atomCompare=rdFMCS.AtomCompare.CompareElements, bondCompare=rdFMCS.BondCompare.CompareAny)
        if mcs_res.numAtoms > 0:
            m_base = mol_base.GetSubstructMatch(Chem.MolFromSmarts(mcs_res.smartsString))
            m_pdb = mol_pdb_heavy.GetSubstructMatch(Chem.MolFromSmarts(mcs_res.smartsString))
            for b_idx, p_idx in zip(m_base, m_pdb):
                p_inf = mol_pdb_heavy.GetAtomWithIdx(p_idx).GetPDBResidueInfo()
                smiles_to_pdb_name[b_idx] = p_inf.GetName().strip() if p_inf else f"C{p_idx+1}"
    else:
        for atom in mol_base.GetAtoms():
            p_inf = atom.GetPDBResidueInfo()
            smiles_to_pdb_name[atom.GetIdx()] = p_inf.GetName().strip() if p_inf else f"C{atom.GetIdx()+1}"

    Chem.RemoveStereochemistry(mol_base); mol_base.RemoveAllConformers(); AllChem.Compute2DCoords(mol_base)
    d2d = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
    d2d.drawOptions().addAtomIndices = False
    d2d.drawOptions().lineWidth = 2.6
    
    highlight_atoms = []
    for atom in mol_base.GetAtoms():
        idx = atom.GetIdx()
        pdb_name = smiles_to_pdb_name.get(idx, f"C{idx+1}")
        if pdb_name in active_atoms_list:
            highlight_atoms.append(idx)
            d2d.drawOptions().atomLabels[idx] = pdb_name
        else:
            d2d.drawOptions().atomLabels[idx] = atom.GetSymbol() if atom.GetSymbol() != "C" else ""

    d2d.DrawMolecule(mol_base, highlightAtoms=highlight_atoms, highlightAtomColors={i: (1.0, 0.0, 0.0) for i in highlight_atoms}, highlightAtomRadii={i: 0.40 for i in highlight_atoms}, highlightBonds=[])
    d2d.FinishDrawing()
    fig2d, ax2d = plt.subplots(figsize=(8, 8))
    ax2d.imshow(Image.open(io.BytesIO(d2d.GetDrawingText()))); ax2d.axis('off')
    
    # RENDER FIXED LINE: Uses the globally instantiated naming string token flawlessly
    fig2d.savefig(f"{ligand_filename_base}_labeled_ligand_2D.png", dpi=300, bbox_inches='tight'); plt.close(fig2d)

    print("\n" + "="*80)
    print("               FINAL CYP450 METABOLIC ACTIVITY FREQUENCY SUMMARY")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("="*80 + "\n")

if __name__ == "__main__":
    main()