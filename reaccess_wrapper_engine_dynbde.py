import os
import sys
import argparse
import glob
import subprocess
import shutil  
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
from PIL import Image
import io

# Secure inline RDKit chemistry rendering dependencies
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

def get_resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

def parse_arguments():
    parser = argparse.ArgumentParser(description="REACCESS Meta-Analysis Wrapper Orchestrator Engine.")
    parser.add_argument("-e", "--engine", required=True, choices=["MD4", "100ps", "100pscont"])
    parser.add_argument("-t", "--topology", required=True)
    parser.add_argument("-l", "--ligand_ref", required=True)
    parser.add_argument("-s", "--smiles", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-m", "--heme", default="HEM")
    parser.add_argument("-f", "--files", nargs="+", default=[])
    parser.add_argument("-d", "--directory", default=None)
    parser.add_argument("--min_dist", type=float, default=3.5)
    parser.add_argument("--max_dist", type=float, default=8.5)
    parser.add_argument("--min_angle", type=float, default=100.0)
    parser.add_argument("--max_angle", type=float, default=145.0)
    parser.add_argument("--alpha", type=float, default=0.5)
    return parser.parse_args()

def main():
    args = parse_arguments()
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    engine_scripts = {
        "MD4": "bde-mopac-analyzerMD4_dynbde.py", 
        "100ps": "bde-mopac-analyzer100ps_dynbde.py", 
        "100pscont": "bde-mopac-analyzer100pscont_dynbde.py"
    }
    
    engine_name = engine_scripts[args.engine]
    target_script = get_resource_path(engine_name)
    if not os.path.exists(target_script) and args.engine == "100pscont":
        target_script = get_resource_path("bde-mopac-analyzer100ps_dynbde.py") 

    if not os.path.exists(target_script):
        print(f"[CRITICAL ERROR] Core analytical script file missing from environment: {target_script}")
        sys.exit(1)

    files_queue = []
    if args.files:
        files_queue = [os.path.abspath(f) for f in args.files]
    elif args.directory and os.path.exists(args.directory):
        search_path = os.path.abspath(args.directory)
        files_queue = sorted(glob.glob(os.path.join(search_path, "*.nc"))) + sorted(glob.glob(os.path.join(search_path, "*.dcd"))) + sorted(glob.glob(os.path.join(search_path, "*.xtc"))) + sorted(glob.glob(os.path.join(search_path, "*.pdb")))
        files_queue = [f for f in files_queue if f != os.path.abspath(args.topology) and f != os.path.abspath(args.ligand_ref)]

    master_summary_records, master_details_records = [], []
    for idx, path_item in enumerate(files_queue, 1):
        item_name = os.path.splitext(os.path.basename(path_item))[0]
        sandbox_dir = os.path.join(output_dir, f"REACCESS_OUT_{item_name}")
        os.makedirs(sandbox_dir, exist_ok=True)
        is_trajectory = path_item.lower().endswith(('.nc', '.dcd', '.xtc'))
        
        cmd = [
            sys.executable, "-u", target_script,
            "-t", os.path.abspath(args.topology) if is_trajectory else path_item,
            "-l", os.path.abspath(args.ligand_ref) if is_trajectory else path_item,
            "-s", args.smiles, "-m", args.heme, "-o", f"cyp_results_{item_name}.xlsx",
            "--min_dist", str(args.min_dist), "--max_dist", str(args.max_dist),
            "--min_angle", str(args.min_angle), "--max_angle", str(args.max_angle), "--alpha", str(args.alpha)
        ]
        if is_trajectory: cmd += ["-x", path_item]

        print(f"--> Executing Subprocess calculations line details...")
        # FIXED: Captured stdout stream buffer explicitly
        process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        # FIXED: Flushed individual sub-script console summaries directly into wrapper log stream
        print(process.stdout, flush=True)

        if is_trajectory:
            gen_folder = os.path.abspath(os.path.splitext(os.path.basename(args.ligand_ref))[0])
        else:
            gen_folder = os.path.abspath(os.path.splitext(os.path.basename(path_item))[0])

        src_summary = os.path.join(gen_folder, "METABOLIC_FINAL_SUMMARY.xlsx")
        src_details = os.path.join(gen_folder, f"cyp_results_{item_name}_BDE_summary.xlsx")

        if os.path.exists(src_summary):
            df_s = pd.read_excel(src_summary); df_s['System_ID'] = item_name; master_summary_records.append(df_s)
            if os.path.exists(src_details):
                df_d = pd.read_excel(src_details); df_d['System_ID'] = item_name; master_details_records.append(df_d)
            for file_asset in glob.glob(os.path.join(gen_folder, "*.*")):
                shutil.move(file_asset, os.path.join(sandbox_dir, os.path.basename(file_asset)))
            try: shutil.rmtree(gen_folder)
            except: pass

    if not master_summary_records: sys.exit(1)
    
    print(f"\n" + "="*80)
    print(" COMPILING COMBINED META-ANALYSIS MASTER SUMMARIES")
    print("="*80)
    
    combined_summary_df = pd.concat(master_summary_records, ignore_index=True)
    with pd.ExcelWriter(os.path.join(output_dir, "COMBINED_METABOLIC_MASTER_SUMMARY.xlsx")) as writer:
        combined_summary_df.to_excel(writer, sheet_name="Per_System_Summary", index=False)
        if master_details_records:
            combined_details_df = pd.concat(master_details_records, ignore_index=True)
            combined_details_df.to_excel(writer, sheet_name="All_Frames_Details", index=False)
        else:
            combined_details_df = pd.DataFrame()

    global_meta_summary = combined_summary_df.groupby('Atom_Name').agg({
        'Mean_Distance': 'mean',
        'Mean_Angle': 'mean',
        'Calculated_BDE': 'mean',
        'Activity_Frequency_%': 'mean'
    }).reset_index()
    global_meta_summary.to_excel(os.path.join(output_dir, "GLOBAL_AVERAGES_METABOLIC_SUMMARY.xlsx"), index=False)

    # FIXED: Embedded 2D ACS Structure rendering engine directly inline to eliminate separate file dependency
    print("--> Generating combined master 2D ACS highlighted map layout panels...")
    master_2d_png = os.path.join(output_dir, "COMBINED_MASTER_2D_ACS_MAP.png")
    try:
        df_active_atoms = global_meta_summary[global_meta_summary['Activity_Frequency_%'] > 0.0]
        active_atoms_list = df_active_atoms['Atom_Name'].tolist()

        mol_pdb = Chem.MolFromPDBFile(os.path.abspath(args.ligand_ref), removeHs=False)
        mol_base = Chem.MolFromSmiles(args.smiles) if args.smiles else Chem.RemoveHs(mol_pdb)
        Chem.SanitizeMol(mol_base)
        
        smiles_to_pdb_name = {}
        if args.smiles and mol_pdb:
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

        Chem.RemoveStereochemistry(mol_base)
        mol_base.RemoveAllConformers()
        AllChem.Compute2DCoords(mol_base)

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

        d2d.DrawMolecule(mol_base, highlightAtoms=highlight_atoms, 
                         highlightAtomColors={i: (1.0, 0.0, 0.0) for i in highlight_atoms}, 
                         highlightAtomRadii={i: 0.40 for i in highlight_atoms}, highlightBonds=[])
        d2d.FinishDrawing()
        
        fig2d, ax2d = plt.subplots(figsize=(8, 8), facecolor='white')
        ax2d.imshow(Image.open(io.BytesIO(d2d.GetDrawingText())))
        ax2d.axis('off')
        ax2d.set_title("Aggregated Master 2D Metabolic Screen Map Profile", fontsize=12, fontweight='bold')
        fig2d.savefig(master_2d_png, dpi=300, bbox_inches='tight')
        plt.close(fig2d)
        print(f"[SUCCESS] Combined master 2D panel rendered at: {master_2d_png}")
    except Exception as img_err:
        print(f"[WARNING] Inline 2D master structural graph mapping failed: {img_err}")

    if not combined_details_df.empty:
        print("--> Generating aggregated master diagnostic plot arrays...")
        sns.set_theme(style="ticks")
        unique_atoms = combined_details_df['Atom_Name'].unique()
        atom_color_dict = dict(zip(unique_atoms, sns.color_palette("turbo", len(unique_atoms))))
        max_seen_bde = max(98.0, combined_details_df['Calculated_BDE'].max() + 2.0)

        fig1, ax1 = plt.subplots(figsize=(8, 5))
        sns.scatterplot(data=combined_details_df, x='Distance_to_Fe', y='Calculated_BDE', hue='Atom_Name', style='System_ID', palette=atom_color_dict, ax=ax1)
        ax1.add_patch(patches.Rectangle((args.min_dist, 70.0), args.max_dist - args.min_dist, max_seen_bde - 70.0, linewidth=1.5, edgecolor='red', facecolor='red', alpha=0.04, linestyle='--'))
        fig1.tight_layout(); fig1.savefig(os.path.join(output_dir, 'COMBINED_MASTER_PLOT_1_distance_vs_BDE.png'), dpi=300); plt.close(fig1)

        fig2, ax2 = plt.subplots(figsize=(8, 5))
        sns.scatterplot(data=combined_details_df, x='Target_Angle', y='Calculated_BDE', hue='Atom_Name', style='System_ID', palette=atom_color_dict, ax=ax2)
        ax2.add_patch(patches.Rectangle((args.min_angle, 70.0), args.max_angle - args.min_angle, max_seen_bde - 70.0, linewidth=1.5, edgecolor='blue', facecolor='blue', alpha=0.04, linestyle='--'))
        fig2.tight_layout(); fig2.savefig(os.path.join(output_dir, 'COMBINED_MASTER_PLOT_2_angle_vs_BDE.png'), dpi=300); plt.close(fig2)

        fig3 = plt.figure(figsize=(11, 7))
        ax3 = fig3.add_subplot(111, projection='3d')
        for (atom_name, system_id), group in combined_details_df.groupby(['Atom_Name', 'System_ID']):
            ax3.scatter(group['Distance_to_Fe'], group['Target_Angle'], group['Calculated_BDE'], color=atom_color_dict[atom_name], label=atom_name, s=40)
        corners = np.array([
            [args.min_dist, args.min_angle, 70.0], [args.max_dist, args.min_angle, 70.0], [args.max_dist, args.max_angle, 70.0], [args.min_dist, args.max_angle, 70.0],
            [args.min_dist, args.min_angle, max_seen_bde], [args.max_dist, args.min_angle, max_seen_bde], [args.max_dist, args.max_angle, max_seen_bde], [args.min_dist, args.max_angle, max_seen_bde]
        ])
        faces = [[corners[0], corners[1], corners[2], corners[3]], [corners[4], corners[5], corners[6], corners[7]], [corners[0], corners[1], corners[5], corners[4]], [corners[2], corners[3], corners[7], corners[6]], [corners[0], corners[3], corners[7], corners[4]], [corners[1], corners[2], corners[6], corners[5]]]
        ax3.add_collection3d(Poly3DCollection(faces, facecolors='green', linewidths=0.5, edgecolors='darkgreen', alpha=0.02))
        handles, labels = ax3.get_legend_handles_labels()
        unique_legend_map = dict(zip(labels, handles))
        ax3.legend(unique_legend_map.values(), unique_legend_map.keys(), bbox_to_anchor=(1.08, 0.85), title="Active Atoms")
        fig3.tight_layout(); fig3.savefig(os.path.join(output_dir, 'COMBINED_MASTER_PLOT_3_3D_landscape_BDE.png'), dpi=300); plt.close(fig3)

    print(f"\n[COMPLETE SUCCESS] Operational pipelines closed out. Master reports saved to: '{output_dir}/'")

if __name__ == "__main__":
    main()