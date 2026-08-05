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

    atom_handles = []
    atom_labels = []

    pose_handles = []
    pose_labels = []

    mode = "atoms"

    for h, l in zip(handles, labels):

        if l == "Atom_Name":
            continue

        if l == "Pose_ID":
            mode = "poses"
            continue

        if mode == "atoms":
            atom_handles.append(h)
            atom_labels.append(l)
        else:
            pose_handles.append(h)
            pose_labels.append(l)

    # Carbon atom legend
    leg_atoms = ax.legend(
        atom_handles,
        atom_labels,
        title="Carbon Atoms",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.00),
        fontsize=8,
        ncol=2,
        frameon=True,
        borderaxespad=0.0,
        columnspacing=1.0,
        handletextpad=0.4
    )

    # Pose legend placed to the RIGHT of carbon legend
    leg_pose = ax.legend(
        pose_handles,
        pose_labels,
        title="Docking Pose",
        loc="upper left",
        bbox_to_anchor=(1.27, 1.00),
        fontsize=8,
        frameon=True,
        borderaxespad=0.0
    )

    ax.add_artist(leg_atoms)
def main():
    args = parse_arguments()
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    master_summary_records = []
    master_details_records = []

    folders = sorted(
        [
            d for d in os.listdir(".")
            if os.path.isdir(d) and d.startswith("lig")
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
    # ==========================================================
    # Collect pose rankings directly from ligand summaries
    # ==========================================================
    
    pose_rank_records = []
    
    required_cols = [
        "Atom_Name",
        "Pose_ID",
        "Rank",
        "Composite_Score",
        "Mean_Distance",
        "Mean_Angle",
        "Calculated_BDE",
        "Activity_Frequency_%"
    ]
    
    for pose_name, pose_df in combined_summary_df.groupby("Pose_ID"):
    
        tmp = pose_df.copy()
    
        tmp = tmp[
            (tmp["Activity_Frequency_%"] > 0)
            &
            (tmp["Calculated_BDE"].notna())
        ].copy()
    
        if len(tmp) == 0:
            continue
    
        # keep only required columns
        tmp = tmp[required_cols]
    
        # rename for consistency
        tmp.rename(
            columns={
                "Rank":"Pose_Rank"
            },
            inplace=True
        )
    
        pose_rank_records.append(tmp)
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
    # Global Ranking based on Pose-wise Ranking
    # ==========================================================
    
    pose_rank_df = pd.concat(
        pose_rank_records,
        ignore_index=True
    )
    print(
        pose_rank_df[
            [
                "Pose_ID",
                "Atom_Name",
                "Composite_Score",
                "Pose_Rank"
            ]
        ]
    )

    # ==========================================================
    # Global ranking based only on average pose-wise rank
    # ==========================================================
    ranking_df = (
        pose_rank_df
        .groupby("Atom_Name")
        .agg(
            Average_Rank=("Pose_Rank","mean"),
            Pose_Frequency=("Pose_ID","nunique")
        )
        .reset_index()
    )
    # ---------------------------------------------
    # Sort by average rank
    # ---------------------------------------------
    
    ranking_df = ranking_df.sort_values(

        "Average_Rank",
    
        ascending=True
    
    ).reset_index(drop=True)
    # ----------------------------------------------------------
    # Keep average rank (1 decimal)
    # ----------------------------------------------------------
    
    ranking_df["Average_Rank"] = (
        ranking_df["Average_Rank"]
        .round(1)
    )
    
    # ----------------------------------------------------------
    # Rounded average rank (nearest integer)
    # ----------------------------------------------------------
    
    ranking_df["Rounded_Average_Rank"] = (
        ranking_df["Average_Rank"]
        .round(0)
        .astype(int)
    )
    
    # ----------------------------------------------------------
    # Global Rank (same as rounded average rank)
    # ----------------------------------------------------------
    
    ranking_df["Rank"] = (
        ranking_df["Rounded_Average_Rank"]
    )
       
    print(

        ranking_df[
            [
    
                "Atom_Name",
                "Average_Rank",
                "Rounded_Average_Rank",
                "Pose_Frequency",
                 "Rank"        
            ]
        ]
    
    )
       
    # ----------------------------------------
    # Merge with global summary
    # ----------------------------------------
    
    global_meta_summary = global_meta_summary.merge(

        ranking_df[
            [
                "Atom_Name",
                "Average_Rank",
                "Rounded_Average_Rank",
                "Pose_Frequency",
                "Rank"
            ]
        ],
    
        on="Atom_Name",
    
        how="left"
    
    )
    
    # ----------------------------------------
    # Save new global summary
    # ----------------------------------------
    
    global_meta_summary.to_excel(
    
        os.path.join(
    
            output_dir,
    
            "GLOBAL_AVERAGES_METABOLIC_SUMMARY.xlsx"
    
        ),
    
        index=False
    
    )
    # FIXED: Embedded 2D ACS Structure rendering engine directly inline to eliminate separate file dependency
    print("--> Generating combined master 2D ACS highlighted map layout panels...")
    master_2d_png = os.path.join(output_dir, "COMBINED_MASTER_2D_ACS_MAP.png")
    try:
        df_active_atoms = ranking_df.copy()
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
        #d2d.drawOptions().lineWidth = 2.6
        
        # ----------------------------------------------------------
        # Build rank dictionary
        # ----------------------------------------------------------

        ranking_df = ranking_df.sort_values(
            "Rank"
        ).reset_index(drop=True)
        
        rank_dict = {}
        
        for _, row in ranking_df.iterrows():
        
            rank_dict[row["Atom_Name"]] = int(row["Rank"])

        highlight_atoms = []

        highlight_colors = {}

        highlight_radii = {}

        opts = d2d.drawOptions()

        opts.setAtomPalette({

            6:(0.0,0.0,0.0),

            7:(0.0,0.0,1.0),

            8:(1.0,0.0,0.0),

            16:(0.85,0.75,0.0),

            17:(0.0,0.6,0.0),

        })


        for atom in mol_base.GetAtoms():

            idx = atom.GetIdx()

            pdb_name = smiles_to_pdb_name.get(idx,f"C{idx+1}")

            symbol = atom.GetSymbol()

            if symbol=="C":

                opts.atomLabels[idx]=""

            else:

                opts.atomLabels[idx]=symbol


            if pdb_name not in active_atoms_list:

                continue


            highlight_atoms.append(idx)

            opts.atomLabels[idx]=pdb_name

            rank = rank_dict.get(pdb_name,999)


            if rank == 1:

                colour = (1.0,0.6,0.6)

            elif rank == 2:

                colour = (1.0,0.80,0.55)

            elif rank == 3:

                colour = (1.0,1.00,0.5)

            else:

                colour = (0.85,0.97,0.85)


            highlight_colors[idx] = colour

            highlight_radii[idx] = 0.45


        d2d.DrawMolecule(

            mol_base,

            highlightAtoms=highlight_atoms,

            highlightAtomColors=highlight_colors,

            highlightAtomRadii=highlight_radii,

            highlightBonds=[]

        )
        d2d.FinishDrawing()
        
        # ----------------------------------------------------------
        # ACS Figure
        # ----------------------------------------------------------

        img = Image.open(io.BytesIO(d2d.GetDrawingText()))

        fig = plt.figure(figsize=(13,8))

        gs = fig.add_gridspec(
            1,
            2,
            width_ratios=[1.15,0.85]
        )

        axMol = fig.add_subplot(gs[0])

        axMol.imshow(img)

        axMol.axis("off")

        axMol.set_title(

            "Global Ranked Metabolic Hotspots",

            fontsize=15,

            fontweight="bold"

        )


        # ==========================================================
        # Ranking Table
        # ==========================================================

        axTbl = fig.add_subplot(gs[1])

        axTbl.axis("off")
        axTbl.set_title("Global Metabolic Ranking", fontsize=15, fontweight="bold", pad=3)


        # ==========================================================
        # Ranking Table
        # ==========================================================
        
        table_df = ranking_df[[
            "Rank",
            "Atom_Name"
        ]].copy()
        
        table_df.columns = [
            "Rank",
            "Atom"
        ]


        tbl = axTbl.table(

            cellText=table_df.values,

            colLabels=table_df.columns,

            loc="center",

            cellLoc="center"

        )


        tbl.auto_set_font_size(False)

        tbl.set_fontsize(14)

        tbl.scale(1.15,1.60)
        # Set column widths
        for (row, col), cell in tbl.get_celld().items():
        
            if col == 0:          # Rank
                cell.set_width(0.22)
        
            elif col == 1:        # Atom
                cell.set_width(0.22)
        
        # ----------------------------------------------------------
        # Colour table rows according to Rank
        # ----------------------------------------------------------

        for (row, col), cell in tbl.get_celld().items():

            if row == 0:

                cell.set_facecolor("#D9D9D9")

                cell.set_text_props(

                    weight="bold",

                    fontsize=15,

                    color="black"

                )

            else:

                rank = int(table_df.iloc[row-1]["Rank"])

                if rank == 1:

                    cell.set_facecolor("#FFBFBF")

                elif rank == 2:

                    cell.set_facecolor("#FFD699")

                elif rank == 3:

                    cell.set_facecolor("#FFFF99")

                else:

                    cell.set_facecolor("#D9F7D9")

                cell.set_text_props(

                    fontsize=14,

                    color="black"

                )

                cell.set_edgecolor("black")

                cell.set_linewidth(0.6)
        fig.tight_layout()

        fig.savefig(

            master_2d_png,

            dpi=300,

            bbox_inches="tight"

        )

        plt.close(fig)

        print(

            f"[SUCCESS] Combined ranked ACS figure saved at: {master_2d_png}"

        )
    except Exception as img_err:
        print(f"[WARNING] Inline 2D master structural graph mapping failed: {img_err}")

    if not combined_details_df.empty:
        print("--> Generating aggregated master diagnostic plot arrays...")
        sns.set_theme(style="ticks")
        unique_atoms = combined_details_df['Atom_Name'].unique()
        atom_color_dict = dict(zip(unique_atoms, sns.color_palette("turbo", len(unique_atoms))))
        bde_low = 50.0
        bde_high = 95.0
        max_seen_bde = max(bde_high, combined_details_df['Calculated_BDE'].max() + 2.0)

        fig1, ax1 = plt.subplots(figsize=(14, 5))
        ax1.set_title(
            "Distance to Heme Fe vs C-H Bond Dissociation Energy",
            fontsize=12,
            fontweight="bold"
        )
        
        ax1.set_xlabel(
            "Distance to Heme Fe (Å)",
            fontsize=11,
            fontweight="bold"
        )
        
        ax1.set_ylabel(
            "Calculated C-H BDE (kcal/mol)",
            fontsize=11,
            fontweight="bold"
        )
        scatter = sns.scatterplot(data=combined_details_df, x='Distance_to_Fe', y='Calculated_BDE', hue='Atom_Name', style='Pose_ID', palette=atom_color_dict, ax=ax1, legend="full")
        ax1.add_patch(patches.Rectangle((args.min_dist, bde_low), args.max_dist - args.min_dist, bde_high - bde_low, linewidth=1.5, edgecolor='red', facecolor='red', alpha=0.05, linestyle='--'))
        add_split_legend(ax1)

        fig1.subplots_adjust(right=0.68); fig1.savefig(os.path.join(output_dir, 'COMBINED_MASTER_PLOT_1_distance_vs_BDE.png'), dpi=300); plt.close(fig1)

        fig2, ax2 = plt.subplots(figsize=(14, 5))
        ax2.set_title(
            "Fe-O-H-C Alignment Angle vs C-H Bond Dissociation Energy",
            fontsize=12,
            fontweight="bold"
        )
        
        ax2.set_xlabel(
            "Target Angle (Degrees)",
            fontsize=11,
            fontweight="bold"
        )
        
        ax2.set_ylabel(
            "Calculated C-H BDE (kcal/mol)",
            fontsize=11,
            fontweight="bold"
        )
        scatter = sns.scatterplot(
            data=combined_details_df,
            x='Target_Angle',
            y='Calculated_BDE',
            hue='Atom_Name',
            style='Pose_ID',
            palette=atom_color_dict,
            ax=ax2,
            legend="full"
        )

        ax2.add_patch(
            patches.Rectangle(
                (args.min_angle, bde_low),
                args.max_angle - args.min_angle,
                bde_high - bde_low,
                linewidth=1.5,
                edgecolor='blue',
                facecolor='blue',
                alpha=0.05,
                linestyle='--'
            )
        )

        add_split_legend(ax2)

        fig2.subplots_adjust(right=0.68)

        fig2.savefig(
            os.path.join(
                output_dir,
                "COMBINED_MASTER_PLOT_2_angle_vs_BDE.png"
            ),
            dpi=300
        )

        plt.close(fig2)

        fig3 = plt.figure(figsize=(15, 7))
        ax3 = fig3.add_subplot(111, projection='3d')
        ax3.set_title(
            "3D Metabolic Reactivity Landscape",
            fontsize=13,
            fontweight="bold",
            pad=20
        )
        
        ax3.set_xlabel(
            "Distance to Heme Fe (Å)",
            fontsize=11,
            fontweight="bold",
            labelpad=12
        )
        
        ax3.set_ylabel(
            "Target Angle (Degrees)",
            fontsize=11,
            fontweight="bold",
            labelpad=12
        )
        
        ax3.set_zlabel(
            "Calculated C-H BDE (kcal/mol)",
            fontsize=11,
            fontweight="bold",
            labelpad=12
        )
        for (atom_name, Pose_ID), group in combined_details_df.groupby(['Atom_Name', 'Pose_ID']):
            ax3.scatter(group['Distance_to_Fe'], group['Target_Angle'], group['Calculated_BDE'], color=atom_color_dict[atom_name], label=atom_name, s=40)
        corners = np.array([
            [args.min_dist, args.min_angle, bde_low], [args.max_dist, args.min_angle, bde_low], [args.max_dist, args.max_angle, bde_low], [args.min_dist, args.max_angle, bde_low],
            [args.min_dist, args.min_angle, bde_high], [args.max_dist, args.min_angle, bde_high], [args.max_dist, args.max_angle, bde_high], [args.min_dist, args.max_angle, bde_high]
        ])
        faces = [[corners[0], corners[1], corners[2], corners[3]], [corners[4], corners[5], corners[6], corners[7]], [corners[0], corners[1], corners[5], corners[4]], [corners[2], corners[3], corners[7], corners[6]], [corners[0], corners[3], corners[7], corners[4]], [corners[1], corners[2], corners[6], corners[5]]]
        ax3.add_collection3d(Poly3DCollection(faces, facecolors='green', linewidths=0.5, edgecolors='darkgreen', alpha=0.02))
        handles, labels = ax3.get_legend_handles_labels()
        unique_legend_map = dict(zip(labels, handles))
        ax3.legend(
            unique_legend_map.values(),
            unique_legend_map.keys(),
            title="Carbon Atoms",
            loc="upper left",
            bbox_to_anchor=(1.25,1.00),
            fontsize=8,
            frameon=True,
            ncol=2,
            columnspacing=1.0,
            handletextpad=0.5
        )
        
        fig3.subplots_adjust(right=0.72); fig3.savefig(os.path.join(output_dir, 'COMBINED_MASTER_PLOT_3_3D_landscape_BDE.png'), dpi=300); plt.close(fig3)

    print(f"\n[COMPLETE SUCCESS] Operational pipelines closed out. Master reports saved to: '{output_dir}/'")

if __name__ == "__main__":
    main()