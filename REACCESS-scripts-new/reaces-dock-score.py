import os
import sys
import argparse
import subprocess
import re
import shutil
import time
import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
from PIL import Image
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry  

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
        description="Master Automated MOPAC PM7 BDE Engine via disangle-dock-MD2.py outputs."
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

def generate_labeled_2d_depiction(
        mol,
        highlight_atoms,
        bde_labels,
        dist_labels,
        angle_labels,
        ranking_df,
        output_png,
        display_name):

    # ----------------------------
    # Remove explicit hydrogens
    # ----------------------------
    mol_2d = Chem.RemoveHs(Chem.Mol(mol))
    Chem.RemoveStereochemistry(mol_2d)
    mol_2d.RemoveAllConformers()
    AllChem.Compute2DCoords(mol_2d)

    d2d = rdMolDraw2D.MolDraw2DCairo(900,700)
    opts = d2d.drawOptions()

    opts.addAtomIndices = False
    opts.clearBackground = False
    opts.addStereoAnnotation = False

    # keep hetero atom colours
  

    opts.setAtomPalette({
        6:(0.0,0.0,0.0),      # Carbon
        7:(0.0,0.0,1.0),      # Nitrogen
        8:(1.0,0.0,0.0),      # Oxygen
        16:(1.0,0.8,0.0),     # Sulfur
        17:(0.0,0.6,0.0),     # Chlorine
    })

    try:
        opts.atomLabelFontSize = 16
    except:
        pass

    # ----------------------------
    # Remove wedge/dashed bonds
    # ----------------------------
    for bond in mol_2d.GetBonds():
        bond.SetBondDir(Chem.BondDir.NONE)

    # ----------------------------
    # Save PDB atom names
    # ----------------------------
    pdb_labels = {}

    for atom in mol_2d.GetAtoms():

        idx = atom.GetIdx()

        info = atom.GetPDBResidueInfo()

        if info:
            pdb_labels[idx] = info.GetName().strip()
        else:
            pdb_labels[idx] = atom.GetSymbol()

    # ----------------------------
    # Default atom labels
    # ----------------------------
    for atom in mol_2d.GetAtoms():

        idx = atom.GetIdx()
        symbol = atom.GetSymbol()

        # Carbon atoms: no label
        if symbol == "C":
            opts.atomLabels[idx] = ""

        # Nitrogen
        elif symbol == "N":
            opts.atomLabels[idx] = "N"

        # Sulfur
        elif symbol == "S":
            opts.atomLabels[idx] = "S"

        else:
            opts.atomLabels[idx] = symbol

    # ----------------------------
    # Highlight atoms
    # ----------------------------
    highlight_colors = {}
    highlight_radii = {}
    # ----------------------------------------------------------
    # Build rank dictionary
    # ----------------------------------------------------------
    
    rank_dict = {}
    
    for _, row in ranking_df.iterrows():
    
        rank_dict[row["Atom_Name"]] = int(row["Rank"])

    for idx in highlight_atoms:

        idx = int(idx)

        if idx >= mol_2d.GetNumAtoms():
            continue

        atom = mol_2d.GetAtomWithIdx(idx)
        

        if atom.GetSymbol() not in ["C","N","S"]:
            continue
        atom_name = pdb_labels[idx]

        rank = rank_dict.get(atom_name,999)
    
        # Rank colours
        if rank == 1:
    
            colour = (1.0,0.0,0.0)        # Red
    
        elif rank == 2:
    
            colour = (1.0,0.55,0.0)       # Orange
    
        elif rank == 3:
    
            colour = (1.0,0.90,0.0)       # Yellow
    
        else:
    
            colour = (0.10,0.70,0.20)     # Green

        highlight_colors[idx] = colour
        highlight_radii[idx] = 0.45

        # Show PDB atom ID ONLY for highlighted atoms
        opts.atomLabels[idx] = atom_name
        opts.useDefaultAtomPalette()

        opts.setAtomPalette({
            6:(0.0,0.0,0.0),     # Carbon
            7:(0.0,0.0,1.0),     # Nitrogen
            8:(1.0,0.0,0.0),     # Oxygen
            16:(0.85,0.75,0.0),  # Sulfur
            17:(0.0,0.6,0.0),    # Chlorine
        })

    # ----------------------------
    # Draw molecule
    # ----------------------------
    d2d.DrawMolecule(
        mol_2d,
        highlightAtoms=list(highlight_colors.keys()),
        highlightAtomColors=highlight_colors,
        highlightAtomRadii=highlight_radii,
        highlightBonds=[]
    )

    # ----------------------------------------------------------
    # Create ACS style figure
    # ----------------------------------------------------------
    
    d2d.FinishDrawing()
    
    img = Image.open(
        io.BytesIO(
            d2d.GetDrawingText()
        )
    )
    
    fig = plt.figure(
        figsize=(14,8),
        facecolor="white"
    )
    
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[3.3,1.3]
    )
    
    ax_img = fig.add_subplot(gs[0])
    ax_tbl = fig.add_subplot(gs[1])
    
    # -----------------------------
    # Molecule
    # -----------------------------
    
    ax_img.imshow(img)
    
    ax_img.axis("off")
    
    ax_img.set_title(
        f"{display_name}\nPredicted Metabolic Hotspots",
        fontsize=13,
        fontweight="bold"
    )
    # ----------------------------------------------------------
    # Ranking Table
    # ----------------------------------------------------------
    
    table_df = ranking_df[
        [
            "Rank",
            "Atom_Name",
            "Composite_Score"
        ]
    ].copy()
    
    table_df["Composite_Score"] = (
        table_df["Composite_Score"]
        .round(3)
    )
    
    table_df.columns = [
        "Rank",
        "Atom",
        "Score"
    ]
    
    ax_tbl.axis("off")
    
    tbl = ax_tbl.table(
    
        cellText=table_df.values,
    
        colLabels=table_df.columns,
    
        cellLoc="center",
    
        colLoc="center",
    
        loc="center"
    
    )
    
    tbl.auto_set_font_size(False)
    
    tbl.set_fontsize(10)
    
    tbl.scale(1.15,1.6)
    # ----------------------------------------------------------
    # Colour rows according to ranking
    # ----------------------------------------------------------
    
    for (row,col), cell in tbl.get_celld().items():
    
        if row == 0:
    
            cell.set_facecolor("#D9D9D9")
    
            cell.set_text_props(weight="bold")
    
        else:
    
            rank = int(table_df.iloc[row-1]["Rank"])
    
            if rank == 1:
    
                cell.set_facecolor("#FF9999")
    
            elif rank == 2:
    
                cell.set_facecolor("#FFD08A")
    
            elif rank == 3:
    
                cell.set_facecolor("#FFF59D")
    
            else:
    
                cell.set_facecolor("#C8F7C5")
    ax_tbl.set_title("Pose Ranking", fontsize=12, fontweight="bold", pad=12)
  
    plt.tight_layout()
    
    fig.savefig(
    
        output_png,
    
        dpi=300,
    
        bbox_inches="tight"
    
    )
    
    plt.close(fig)

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
    base=os.path.splitext(filename)[0]
    out_file=f"{base}.out"
    try:
        subprocess.run(["mopac", filename], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        print(f"Error running MOPAC on {filename}: {e}")
        return None
    for _ in range(60):
        if os.path.exists(out_file) and os.path.getsize(out_file)>0:
            break
        time.sleep(1)
    if not os.path.exists(out_file):
        print(f"[ERROR] Output missing: {out_file}")
        return None
    hof=None
    hof_re=re.compile(r"FINAL\s+HEAT\s+OF\s+FORMATION\s+=\s+(-?\d+\.\d+)\s+KCAL/MOL")
    with open(out_file,"r",errors="ignore") as f:
        for line in f:
            m=hof_re.search(line)
            if m:
                hof=float(m.group(1))
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

    script_name = "disangle-dock-MD2.py"
    if not os.path.exists(script_name):
        script_name = "disangle-dock-MD.py"
        
    temp_patched_script = "temp_runtime_patched_engine.py"
    with open(script_name, "r") as f:
        lines = f.readlines()
        
    new_lines = []
    for line in lines:
        if "format=\"NETCDF\"" in line or "format='NETCDF'" in line or "mda.Universe(args.topology, args.trajectory" in line:
            new_lines.append(f"        ligand_resname = '{find_ligand_resname_from_universe(mda.Universe(abs_topology_path))}'\n")
        new_lines.append(line.replace('format=\"NETCDF\"', 'format=\"NCDF\"').replace("format='NETCDF'", "format='NCDF'"))
        
    with open(temp_patched_script, "w") as f:
        f.writelines(new_lines)
        
    command = [sys.executable, temp_patched_script, "-t", abs_topology_path, "-l", abs_ligand_path]
    if abs_trajectory_path: 
        command += ["-x", abs_trajectory_path]
    if args.heme: 
        command += ["-m", args.heme]
    if args.resolution is not None: 
        command += ["-r", str(args.resolution)]
    command += ["-o", local_excel_name]
    
    print(f"--> [Step 1 Execution] Invoking Native Geometry Engine: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
    finally:
        if os.path.exists(temp_patched_script):
            os.remove(temp_patched_script)

    print("--> [Step 2] Setting up molecular templates via RDKit...")
    
    if args.smiles:
        mol_receptor = Chem.MolFromSmiles(args.smiles)
        mol_receptor = Chem.AddHs(mol_receptor)
        mol_receptor_pdb = Chem.MolFromPDBFile(abs_ligand_path, removeHs=False)
        if mol_receptor_pdb:
            try: mol_receptor = Chem.AllChem.AssignBondOrdersFromTemplate(mol_receptor, mol_receptor_pdb)
            except: mol_receptor = mol_receptor_pdb
    else:
        mol_receptor = Chem.MolFromPDBFile(abs_ligand_path, removeHs=False)
        try: Chem.DetermineBondsOrderDetail(mol_receptor)
        except: pass

    pdb_name_to_rdkit_idx = {}
    rdkit_idx_to_pdb_name = {}
    for atom in mol_receptor.GetAtoms():
        p_info = atom.GetPDBResidueInfo()
        if p_info:
            name_stripped = p_info.GetName().strip()
            pdb_name_to_rdkit_idx[name_stripped] = atom.GetIdx()
            rdkit_idx_to_pdb_name[atom.GetIdx()] = name_stripped
        else:
            pdb_name_to_rdkit_idx[f"{atom.GetSymbol()}{atom.GetIdx()+1}"] = atom.GetIdx()
            rdkit_idx_to_pdb_name[atom.GetIdx()] = f"{atom.GetSymbol()}{atom.GetIdx()+1}"

    df_geo = pd.read_excel(local_excel_name)
    df_geo['RDKit_Idx'] = df_geo['Atom_Name'].apply(lambda n: pdb_name_to_rdkit_idx.get(str(n).strip(), -1))
    df_geo['Element'] = df_geo['Atom_Name'].str.extract(r'([A-Za-z]+)')
    
    def clean_angle_val(val):
        if isinstance(val, str):
            try:
                lst = ast.literal_eval(val)
                if not lst: return 120.0
                valid_angles = [a for a in lst if 100.0 <= a <= 145.0]
                return float(valid_angles[0]) if valid_angles else float(lst[0])
            except: 
                return 120.0
        return float(val)

    df_geo['Parsed_Angle'] = df_geo['Angles'].apply(clean_angle_val)

    shutil.move(local_excel_name, os.path.join(output_dir, local_excel_name))
    os.chdir(output_dir)

    # FIXED: Identify and filter strictly for Carbon atoms that passed BOTH criteria at least once
    active_carbon_ids = df_geo[(df_geo['Element'] == 'C') & (df_geo['Passes_Both'] == True)]['RDKit_Idx'].unique()
    df_reactive_carbons = df_geo[(df_geo['Element'] == 'C') & (df_geo['RDKit_Idx'].isin(active_carbon_ids))]
    
    reactive_carbons_group = df_reactive_carbons.groupby('RDKit_Idx')
    total_unique_carbons = reactive_carbons_group.ngroups
    print(f"--> [Status] Isolated {total_unique_carbons} unique geometrically active carbon centers for quantum analysis.")

    print("--> [Step 3] Computing baseline MOPAC values for radical fragment parameters...")
    h_mol = Chem.MolFromSmiles("[H]")
    h_conf = Chem.Conformer(1)
    h_conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    h_mol.AddConformer(h_conf)
    
    write_mopac_input("H_rad.mop", h_mol, charge=0, multiplicity=2)
    h_rad_energy = run_mopac("H_rad.mop")
    # keep H_rad.out
        
    write_mopac_input("parent.mop", mol_receptor, charge=0, multiplicity=1)
    parent_energy = run_mopac("parent.mop")
    # keep parent.out

    bde_results = []
    bde_labels_map = {}

    print("--> [MOPAC Engine] Launching parallel radical optimizations loop...")
    for progress_count, (rd_idx, group) in enumerate(reactive_carbons_group, 1):
        rd_idx_int = int(rd_idx)
        if rd_idx_int < 0: continue
        atom_name_pdb = rdkit_idx_to_pdb_name.get(rd_idx_int, f"C{rd_idx_int+1}")
        print(f"    [MOPAC Optimization {progress_count}/{total_unique_carbons}] Simulating radical at: {atom_name_pdb}...")
        
        rd_atom = mol_receptor.GetAtomWithIdx(rd_idx_int)
        hydrogens = [n for n in rd_atom.GetNeighbors() if n.GetSymbol() == 'H']
        
        if not hydrogens: 
            print(f"    [Skipped] Atom {atom_name_pdb} possesses no valid explicit hydrogen configurations.")
            continue
        h_idx_global = hydrogens[0].GetIdx()
        
        editable_mol = Chem.RWMol(mol_receptor)
        editable_mol.RemoveAtom(h_idx_global)
        radical_mol = editable_mol.GetMol()
        
        rad_base = f"radical_{atom_name_pdb}"
        rad_mop_name = f"{rad_base}.mop"
        
        write_mopac_input(rad_mop_name, radical_mol, charge=0, multiplicity=2)
        radical_energy = run_mopac(rad_mop_name)
        
        # keep radical out
        if radical_energy is None: 
            print(f"    [WARNING] MOPAC calculation failed to converge for target atom {atom_name_pdb}.")
            continue
            
        bde = (radical_energy + h_rad_energy) - parent_energy
        bde_labels_map[rd_idx_int] = float(bde)
        
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

    if not bde_results:
        print("[WARNING] No active structural centers generated data points.")
        return

    bde_df = pd.DataFrame(bde_results)
    bde_df.to_excel(f"{local_excel_name.replace('.xlsx', '')}_BDE_summary.xlsx", index=False)

    bde_df['Passes_BDE'] = bde_df['Calculated_BDE'].between(50.0, 95.0)
    bde_df['Is_Active_Site'] = bde_df['Passes_Distance'] & bde_df['Passes_Angle_Criteria'] & bde_df['Passes_BDE']

    total_trajectory_frames = df_geo['Frame'].nunique()
    summary_data = []

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

    summary_df = pd.DataFrame(summary_data)
    ranking_df = summary_df.copy()
    # Carbon atoms only
    ranking_df = ranking_df[
        ranking_df["Calculated_BDE"].notna()
    ].copy()
    
    # ----------------------------------------------------------
    # Ideal values for ranking
    # ----------------------------------------------------------
    
    IDEAL_DISTANCE = ranking_df["Mean_Distance"].min()
    
    IDEAL_ANGLE = 122.5
    
    IDEAL_BDE = ranking_df["Calculated_BDE"].min()
        
    # ----------------------------------------------------------
    # Distance score
    # Minimum distance receives score = 1
    # ----------------------------------------------------------
    
    ranking_df["Distance_Score"] = (
        ranking_df["Mean_Distance"].max()
        - ranking_df["Mean_Distance"]
    )
    
    ranking_df["Distance_Score"] = (
        ranking_df["Distance_Score"]
        /
        ranking_df["Distance_Score"].max()
    )
    # ----------------------------------------------------------
    # Angle score
    # Smaller deviation from ideal angle is better
    # ----------------------------------------------------------
    
    ranking_df["Angle_Error"] = (
        ranking_df["Mean_Angle"] - 122.5
    ).abs()
    
    ranking_df["Angle_Score"] = (
        1
        - ranking_df["Angle_Error"]
          /
          ranking_df["Angle_Error"].max()
    )
    
    # ----------------------------------------------------------
    # Lower BDE receives higher score
    # ----------------------------------------------------------
    
    ranking_df["BDE_Score"] = (
        ranking_df["Calculated_BDE"].max()
        - ranking_df["Calculated_BDE"]
    )
    
    ranking_df["BDE_Score"] = (
        ranking_df["BDE_Score"]
        /
        ranking_df["BDE_Score"].max()
    )
    
    # ----------------------------------------------------------
    # Equal importance
    # ----------------------------------------------------------
    
    ranking_df["Composite_Score"] = (
          ranking_df["Distance_Score"]
        + ranking_df["Angle_Score"]
        + ranking_df["BDE_Score"]
    ) / 3.0
    
    # ----------------------------------------------------------
    # Final Rank
    # ----------------------------------------------------------
    
    ranking_df = ranking_df.sort_values(
        "Composite_Score",
        ascending=False
    ).reset_index(drop=True)
    
    ranking_df["Rank"] = np.arange(
        1,
        len(ranking_df) + 1
    )
    
    # ----------------------------------------------------------
    # Merge ranking back
    # ----------------------------------------------------------
    
    summary_df = summary_df.merge(
    
        ranking_df[
            [
                "Atom_Name",
                "Distance_Score",
                "Angle_Score",
                "BDE_Score",
                "Composite_Score",
                "Rank"
            ]
        ],
    
        on="Atom_Name",
        how="left"
    
    )
    
    summary_df.to_excel(
        "METABOLIC_FINAL_SUMMARY.xlsx",
        index=False
    )
    # ----------------------------------------------------------
    # Add hetero atoms (N, O, S) that satisfy the accessibility
    # criteria. They are not subjected to BDE calculations.
    # ----------------------------------------------------------

    hetero_df = df_geo[
        (df_geo["Element"].isin(["N", "S"])) &
        (df_geo["Passes_Distance"] == True)
    ]

    for atom_name, group in hetero_df.groupby("Atom_Name"):

        summary_df.loc[len(summary_df)] = {
            "Atom_Name": atom_name,
            "Mean_Distance": round(group["Distance_to_Fe"].mean(), 2),
            "Mean_Angle": np.nan,
            "Calculated_BDE": np.nan,
            "Active_Frames": group["Frame"].nunique(),
            "Total_Frames": total_trajectory_frames,
            "Activity_Frequency_%": round(
                group["Frame"].nunique()
                / total_trajectory_frames
                * 100,
                1
            )
        }

    summary_df.to_excel(
        "METABOLIC_FINAL_SUMMARY.xlsx",
        index=False
    )
    print("\n" + "="*80)
    print("               FINAL CYP450 METABOLIC ACTIVITY FREQUENCY SUMMARY")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("="*80 + "\n")

    df_active_atoms = summary_df[
    (
        summary_df["Calculated_BDE"].between(50.0,95.0)
    )
    |
    (
        summary_df["Calculated_BDE"].isna()
    )
    ]
    highlight_atoms_list = [pdb_name_to_rdkit_idx[name] for name in df_active_atoms['Atom_Name']]
    
    dist_labels_map = {pdb_name_to_rdkit_idx[r['Atom_Name']]: r['Mean_Distance'] for _, r in df_active_atoms.iterrows()}
    angle_labels_map = {pdb_name_to_rdkit_idx[r['Atom_Name']]: r['Mean_Angle'] for _, r in df_active_atoms.iterrows()}
    print(df_active_atoms[["Atom_Name","Calculated_BDE","Activity_Frequency_%"]])
    print("--> [Step 4] Rendering annotated structure visualization coordinates engine...")
    png_final_name = f"{ligand_filename_base}_labeled_ligand_2D.png"
    generate_labeled_2d_depiction(mol_receptor, highlight_atoms_list, bde_labels_map, dist_labels_map, angle_labels_map, ranking_df, png_final_name, ligand_filename_base)

    print("--> [Graphics Core] Generating color-coded plots with legends and boundary boxes...")
    sns.set_theme(style="ticks")
    
    res_val = args.resolution if args.resolution is not None else 2.5
    dist_low, dist_high = 6.0 - res_val, 6.0 + res_val
    ang_low, ang_high = 100.0, 145.0
    bde_low, bde_high = 50.0, 95.0

    x_dist, y_angl, z_bde = bde_df['Distance_to_Fe'].values, bde_df['Target_Angle'].values, bde_df['Calculated_BDE'].values

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

    # Plot 3: 3D Landscape Plot with Color Palette and Cuboid
    fig3 = plt.figure(figsize=(9, 6.5))
    ax3 = fig3.add_subplot(111, projection='3d')
    
    for atom_name, group in bde_df.groupby('Atom_Name'):
        ax3.scatter(group['Distance_to_Fe'], group['Target_Angle'], group['Calculated_BDE'], label=atom_name, color=atom_color_dict[atom_name], s=50, edgecolor='k', alpha=0.75)
        
    corners = np.array([
        [dist_low, ang_low, bde_low], [dist_high, ang_low, bde_low], [dist_high, ang_high, bde_low], [dist_low, ang_high, bde_low],
        [dist_low, ang_low, bde_high], [dist_high, ang_low, bde_high], [dist_high, ang_high, bde_high], [dist_low, ang_high, bde_high]
    ])
    faces = [
        [corners[0], corners[1], corners[2], corners[3]], 
        [corners[4], corners[5], corners[6], corners[7]], 
        [corners[0], corners[1], corners[5], corners[4]], 
        [corners[2], corners[3], corners[7], corners[6]], 
        [corners[0], corners[3], corners[7], corners[4]], 
        [corners[1], corners[2], corners[6], corners[5]]  
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
    
    print(f"\n[SUCCESS] Unified re-indexing pipeline completed flawlessly. Results inside: '{output_dir}/'")

if __name__ == "__main__":
    main()