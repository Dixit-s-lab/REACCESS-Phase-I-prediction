import os
import sys
import argparse
import glob
import subprocess
import shutil  
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
from PIL import Image
import io
import traceback

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
    parser = argparse.ArgumentParser(
        description="REACCESS Combined Docking Pose Analysis"
    )

    parser.add_argument(
        "-l", "--ligand_ref",
        required=True,
        help="Ligand PDB file"
    )

    parser.add_argument(
        "-s", "--smiles",
        required=True,
        help="Ligand SMILES"
    )

    parser.add_argument(
        "-o", "--output_dir",
        default="Combined_Analysis"
    )

    parser.add_argument(
        "--min_dist",
        type=float,
        default=3.5
    )

    parser.add_argument(
        "--max_dist",
        type=float,
        default=8.5
    )

    parser.add_argument(
        "--min_angle",
        type=float,
        default=100.0
    )

    parser.add_argument(
        "--max_angle",
        type=float,
        default=145.0
    )

    return parser.parse_args()
def add_split_legend(ax):

    handles, labels = ax.get_legend_handles_labels()

    unique = {}

    for h, l in zip(handles, labels):

        if l in ["Atom_Name", "Pose_ID"]:
            continue

        if l not in unique:
            unique[l] = h

    ax.legend(
        unique.values(),
        unique.keys(),
        title="Reactive Carbon Atoms",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        fontsize=8,
        title_fontsize=9,
        frameon=True,
        ncol=2,
        borderaxespad=0.0,
        columnspacing=1.2,
        handletextpad=0.5
    )

    
def main():
    args = parse_arguments()
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    master_summary_records = []
    master_details_records = []

    folders = sorted(
        [
            d for d in os.listdir(".")
            if os.path.isdir(d) 
        ]
    )

    print("=" * 80)
    print("READING DOCKING POSE RESULTS")
    print("=" * 80)

    for folder in folders:

        summary_file = os.path.join(
            folder,
            "METABOLIC_FINAL_SUMMARY.xlsx"
        )

        detail_files = sorted(
            glob.glob(os.path.join(folder, "*_BDE_summary.xlsx"))
        )

        if not os.path.exists(summary_file):
            print(f"Skipping {folder}")
            continue

        print(f"Reading {folder}")

        df_summary = pd.read_excel(summary_file)
        df_summary["Pose_ID"] = folder
        master_summary_records.append(df_summary)

        if detail_files:

            df_detail = pd.read_excel(detail_files[0])
            df_detail["Pose_ID"] = folder
            master_details_records.append(df_detail)

    if not master_summary_records:

        print("No docking pose results found.")
        sys.exit(1)
    
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
    # ==========================================================
    # Global Ranking (Activity Frequency)
    # ==========================================================
    
    rank_df = (
        global_meta_summary
        .sort_values(
            "Activity_Frequency_%",
            ascending=False
        )
        .reset_index(drop=True)
    )
    
    # Insert Rank column
    rank_df.insert(
        0,
        "Rank",
        np.arange(1, len(rank_df) + 1)
    )
    
    # ----------------------------------------------------------
    # Remove old Rank column if it exists
    # ----------------------------------------------------------
    
    if "Rank" in global_meta_summary.columns:
    
        global_meta_summary.drop(
            columns=["Rank"],
            inplace=True
        )
    
    # ----------------------------------------------------------
    # Merge Rank
    # ----------------------------------------------------------
    
    global_meta_summary = pd.merge(
    
        global_meta_summary,
    
        rank_df[
            [
                "Atom_Name",
                "Rank"
            ]
        ],
    
        on="Atom_Name",
    
        how="left"
    
    )
    
    # ----------------------------------------------------------
    # Move Rank after Activity_Frequency_%
    # ----------------------------------------------------------
    
    cols = [
    
        c for c in global_meta_summary.columns
    
        if c != "Rank"
    
    ]
    
    idx = cols.index("Activity_Frequency_%") + 1
    
    cols.insert(idx, "Rank")
    
    global_meta_summary = global_meta_summary[cols]
    print(global_meta_summary.columns.tolist())
    print(global_meta_summary.head())
    
    rank_df.to_excel(
    
        os.path.join(
    
            output_dir,
    
            "GLOBAL_ACTIVITY_RANKING.xlsx"
    
        ),
    
        index=False
    
    )
    global_meta_summary.to_excel(

        os.path.join(
    
            output_dir,
    
            "GLOBAL_AVERAGES_METABOLIC_SUMMARY.xlsx"
    
        ),
    
        index=False
    
    )
    
    rank_df.to_excel(
    
        os.path.join(
    
            output_dir,
    
            "GLOBAL_ACTIVITY_RANKING.xlsx"
    
        ),
    
        index=False
    
    )
    
    # FIXED: Embedded 2D ACS Structure rendering engine directly inline to eliminate separate file dependency
    print("--> Generating combined master 2D ACS highlighted map layout panels...")

    master_2d_png = os.path.join(
        output_dir,
        "COMBINED_MASTER_2D_ACS_MAP.png"
    )
    
    try:
        df_active_atoms = global_meta_summary[global_meta_summary['Activity_Frequency_%'] > 0.0]
        active_atoms_list = df_active_atoms['Atom_Name'].tolist()

        mol_pdb = Chem.MolFromPDBFile(os.path.abspath(args.ligand_ref), removeHs=False)
        mol_base = Chem.MolFromSmiles(args.smiles) if args.smiles else Chem.RemoveHs(mol_pdb)
        Chem.SanitizeMol(mol_base)
        
        # ------------------------------------------------------------------
        # Build 2D structure from SMILES and obtain PDB atom names
        # ------------------------------------------------------------------
        print("1")
        mol_base = Chem.MolFromSmiles(args.smiles)
        print("2")
        mol_base = Chem.RemoveHs(mol_base)
        print("3")

        mol_pdb = Chem.MolFromPDBFile(
            os.path.abspath(args.ligand_ref),
            removeHs=False
        )
        print("4")

        mol_pdb = Chem.RemoveHs(mol_pdb)
        print("5")

        Chem.RemoveStereochemistry(mol_base)
        print("6")
        mol_base.RemoveAllConformers()
        print("7")
        AllChem.Compute2DCoords(mol_base)
        print("8")

        mcs = rdFMCS.FindMCS(
            [mol_base, mol_pdb],
            ringMatchesRingOnly=True,
            completeRingsOnly=True
        )
        
        mcsMol = Chem.MolFromSmarts(mcs.smartsString)
        
        smiles_match = mol_base.GetSubstructMatch(mcsMol)
        pdb_match = mol_pdb.GetSubstructMatch(mcsMol)
        
        smiles_to_pdb_name = {}
        
        for smi_idx, pdb_idx in zip(smiles_match, pdb_match):
        
            atom = mol_pdb.GetAtomWithIdx(pdb_idx)
        
            info = atom.GetPDBResidueInfo()
        
            if info:
                smiles_to_pdb_name[smi_idx] = info.GetName().strip()

        Chem.RemoveStereochemistry(mol_base)
        mol_base.RemoveAllConformers()
        AllChem.Compute2DCoords(mol_base)
        
        print("9")
        d2d = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
        d2d.drawOptions().addAtomIndices = False
        #d2d.drawOptions().lineWidth = 2.6
        
        highlight_atoms = []
        print("10")

        opts = d2d.drawOptions()
        print("11")

        for atom in mol_base.GetAtoms():

            idx = atom.GetIdx()

            pdb_name = smiles_to_pdb_name.get(
                idx,
                atom.GetSymbol()
            )

            if pdb_name in active_atoms_list:

                highlight_atoms.append(idx)

                opts.atomLabels[idx] = pdb_name

            else:

                symbol = atom.GetSymbol()

                if symbol == "C":
                    opts.atomLabels[idx] = ""

                elif symbol == "N":
                    opts.atomLabels[idx] = "N"

                elif symbol == "O":
                    opts.atomLabels[idx] = "O"

                elif symbol == "S":
                    opts.atomLabels[idx] = "S"

                elif symbol == "Cl":
                    opts.atomLabels[idx] = "Cl"

                elif symbol == "F":
                    opts.atomLabels[idx] = "F"

                elif symbol == "Br":
                    opts.atomLabels[idx] = "Br"

                elif symbol == "P":
                    opts.atomLabels[idx] = "P"

                else:
                    opts.atomLabels[idx] = symbol

        # ----------------------------------------------------------
        # Highlight colours according to ranking
        # ----------------------------------------------------------
        
        rank_lookup = dict(
            zip(
                rank_df["Atom_Name"],
                rank_df["Rank"]
            )
        )
        
        highlight_colors = {}
        
        for idx in highlight_atoms:
        
            pdb_name = smiles_to_pdb_name.get(idx)
        
            rank = rank_lookup.get(pdb_name, 999)
        
            if rank == 1:
                highlight_colors[idx] = (1.00, 0.60, 0.60)      # Red
        
            elif rank == 2:
                highlight_colors[idx] = (1.00, 0.80, 0.55)      # Orange
        
            elif rank == 3:
                highlight_colors[idx] = (1.00, 1.00, 0.50)      # Yellow
        
            else:
                highlight_colors[idx] = (0.85, 0.97, 0.85)      # Green
        print("12")
        d2d.DrawMolecule(mol_base, highlightAtoms=highlight_atoms, 
                         highlightAtomColors=highlight_colors, 
                         highlightAtomRadii={i: 0.28 for i in highlight_atoms}, highlightBonds=[])
        print("13")
        d2d.FinishDrawing()
        print("14")
        
        png_bytes = d2d.GetDrawingText()

        img = Image.open(io.BytesIO(d2d.GetDrawingText()))
        fig = plt.figure(figsize=(15,9))

        gs = fig.add_gridspec(
            1,
            2,
            width_ratios=[3.5,1]
        )
        
        ax_img = fig.add_subplot(gs[0])
        ax_tbl = fig.add_subplot(gs[1])
        
        ax_img.imshow(img)
        ax_img.axis("off")
        ax_img.set_title("Potential SOMs identified using REACCESS approch", fontsize=20, fontweight="bold", pad=3)
        table_df = rank_df[
        ["Rank","Atom_Name","Activity_Frequency_%"]
    ].copy()
    
        table_df["Activity_Frequency_%"] = (
            table_df["Activity_Frequency_%"]
            .round(1)
        )
        
        table_df.columns = [
            "Rank",
            "Atom",
            "Activity (%)"
        ]
        ax_tbl.axis("off")
    
        tbl = ax_tbl.table(
            cellText=table_df.values,
            colLabels=table_df.columns,
            cellLoc="center",
            colLoc="center",
            loc="center" ,
            colWidths=[0.30, 0.30, 0.60]   # Rank | Atom | Activity
        )
        
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(15)
        tbl.scale(1.25,2.2)
        
        # ----------------------------------------------------------
        # Colour table according to ranking
        # ----------------------------------------------------------
        
        for (row, col), cell in tbl.get_celld().items():

            cell.set_edgecolor("black")
            cell.set_linewidth(0.8)
        
            if row == 0:
                cell.set_facecolor("#CFCFCF")
                cell.get_text().set_weight("normal")
                continue
        
            rank = int(table_df.iloc[row-1]["Rank"])
        
            if rank == 1:
                bg = (1.00, 0.60, 0.60)      # Red
        
            elif rank == 2:
                bg = (1.00, 0.80, 0.55)      # Orange
        
            elif rank == 3:
                bg = (1.00, 1.00, 0.50)      # Yellow
        
            else:
                bg = (0.85, 0.97, 0.85)      # Green
        
            # Colour the entire row
            cell.set_facecolor(bg)
            
        
        ax_tbl.set_title(
            "Reactive Site Ranking",
            fontsize=13,
            fontweight="bold",
            pad=3
        )
        plt.savefig(
            master_2d_png,
            dpi=400,
            bbox_inches="tight"
        )
        
        plt.close(fig)
                
        print("15")
        print(f"[SUCCESS] Combined master 2D panel rendered at: {master_2d_png}")

    except Exception:
        traceback.print_exc()

    if not combined_details_df.empty:
       print("--> Generating aggregated master diagnostic plot arrays...")
       sns.set_theme(style="ticks")
       unique_atoms = combined_details_df['Atom_Name'].unique()
       atom_color_dict = dict(zip(unique_atoms, sns.color_palette("turbo", len(unique_atoms))))
       max_seen_bde = max(98.0, combined_details_df['Calculated_BDE'].max() + 2.0)

       fig1, ax1 = plt.subplots(figsize=(9, 6))
       sns.scatterplot(data=combined_details_df, x='Distance_to_Fe', y='Calculated_BDE', hue='Atom_Name', style='Pose_ID', palette=atom_color_dict, s=35, edgecolor="black", alpha=0.85, ax=ax1)
       ax1.add_patch(patches.Rectangle((args.min_dist, 70.0), args.max_dist - args.min_dist, max_seen_bde - 70.0, linewidth=1.5, edgecolor='red', facecolor='red', alpha=0.05, linestyle='--'))
       ax1.set_title("Distance to Heme Fe vs C-H Bond Dissociation Energy", fontsize=13, fontweight='bold')
       ax1.set_xlabel("Distance to Heme Fe (Å)")
       ax1.set_ylabel("Calculated C-H BDE (kcal/mol)")

       add_split_legend(ax1)

       fig1.subplots_adjust(right=0.72)
        
       fig1.savefig(os.path.join(output_dir, "COMBINED_MASTER_PLOT_1_distance_vs_BDE.png"), dpi=300, bbox_inches="tight")
        
       plt.close(fig1)

       fig2, ax2 = plt.subplots(figsize=(9,6))

       fig2, ax2 = plt.subplots(figsize=(9,6))

       sns.scatterplot(
           data=combined_details_df,
           x='Target_Angle',
           y='Calculated_BDE',
           hue='Atom_Name',
           style='Pose_ID',
           palette=atom_color_dict,
           s=35,
           edgecolor="black",
           alpha=0.85,
           ax=ax2
       )

       ax2.add_patch(
           patches.Rectangle(
               (args.min_angle,70.0),
               args.max_angle-args.min_angle,
               max_seen_bde-70.0,
               linewidth=1.5,
               edgecolor='blue',
               facecolor='blue',
               alpha=0.05,
               linestyle='--'
           )
       )
      
       ax2.set_title("Target Angle vs C-H Bond Dissociation Energy", fontsize=13, fontweight='bold')
        
       ax2.set_xlabel("Target Angle (°)")
       ax2.set_ylabel("Calculated C-H BDE (kcal/mol)")
        
       add_split_legend(ax2)
        
       fig2.subplots_adjust(right=0.72)
        
       fig2.savefig(
            os.path.join(output_dir,
            "COMBINED_MASTER_PLOT_2_angle_vs_BDE.png"),
            dpi=300,
            bbox_inches="tight"
        )
        
       plt.close(fig2)

       fig3 = plt.figure(figsize=(16, 7))
       ax3 = fig3.add_subplot(111, projection='3d')
       # Make the plot rectangular
       try:
           ax3.set_box_aspect((1.8, 1.3, 1.0))
       except:
           pass
       for (atom_name, system_id), group in combined_details_df.groupby(['Atom_Name', 'Pose_ID']):
            ax3.scatter(group['Distance_to_Fe'], group['Target_Angle'], group['Calculated_BDE'], color=atom_color_dict[atom_name], label=atom_name, s=40)
       corners = np.array([
            [args.min_dist, args.min_angle, 70.0], [args.max_dist, args.min_angle, 70.0], [args.max_dist, args.max_angle, 70.0], [args.min_dist, args.max_angle, 70.0],
            [args.min_dist, args.min_angle, max_seen_bde], [args.max_dist, args.min_angle, max_seen_bde], [args.max_dist, args.max_angle, max_seen_bde], [args.min_dist, args.max_angle, max_seen_bde]
        ])
       faces = [[corners[0], corners[1], corners[2], corners[3]], [corners[4], corners[5], corners[6], corners[7]], [corners[0], corners[1], corners[5], corners[4]], [corners[2], corners[3], corners[7], corners[6]], [corners[0], corners[3], corners[7], corners[4]], [corners[1], corners[2], corners[6], corners[5]]]
       ax3.add_collection3d(Poly3DCollection(faces, facecolors='green', linewidths=0.5, edgecolors='darkgreen', alpha=0.02))
       handles, labels = ax3.get_legend_handles_labels()

       unique = {}
        
       for h, l in zip(handles, labels):
        
           if l not in unique:
               unique[l] = h
       ax3.set_title("Combined Metabolic Reactivity Landscape", fontsize=14, fontweight="bold")

       ax3.set_xlabel("Distance to Heme Fe (Å)", fontsize=11, labelpad=18)

       ax3.set_ylabel("Target Angle (°)", fontsize=11, labelpad=18)

       ax3.set_zlabel("Calculated C-H BDE (kcal/mol)", fontsize=11, labelpad=18)

       ax3.view_init(elev=24, azim=-60) 
       ax3.legend(unique.values(), unique.keys(), title="Reactive Carbon Atoms", loc="upper left", bbox_to_anchor=(1.20,0.55), fontsize=8, title_fontsize=9, frameon=True, ncol=2, columnspacing=1.3, handletextpad=0.4, borderaxespad=0.0)
       fig3.subplots_adjust(left=0.06, right=0.68,bottom=0.08, top=0.92)

       fig3.savefig(os.path.join(output_dir,"COMBINED_MASTER_PLOT_3_3D_landscape_BDE.png"), dpi=300, bbox_inches="tight", pad_inches=0.5)

       plt.close(fig3)

    print(
        f"\n[COMPLETE SUCCESS] Operational pipelines closed out. "
        f"Master reports saved to: '{output_dir}/'"
    )

if __name__ == "__main__":
    main()