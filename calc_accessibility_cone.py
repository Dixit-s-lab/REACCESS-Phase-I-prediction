import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
try:
    import MDAnalysis as mda
    import warnings
    warnings.filterwarnings('ignore', category=UserWarning)
except ImportError:
    mda = None

def log_print(msg, end='\n', log_to_file=True):
    print(msg, end=end)
    if log_to_file and end == '\n':
        with open("accessibility_calc.log", "a", encoding="utf-8") as f:
            f.write(msg + "\n")

class AccessibilityModel:
    def __init__(self, rmin, rmax, theta_min, theta_max, grid_spacing, prot_buffer, lig_buffer):
        self.rmin = rmin
        self.rmax = rmax
        self.theta_min = np.radians(theta_min)
        self.theta_max = np.radians(theta_max)
        self.spacing = grid_spacing
        self.prot_buffer = prot_buffer
        self.lig_buffer = lig_buffer
        self.voxel_vol = self.spacing ** 3
        self.voxel_area = self.spacing ** 2

    def generate_local_grid(self):
        axis = np.arange(-self.rmax - 2.0, self.rmax + 2.0, self.spacing)
        X, Y, Z = np.meshgrid(axis, axis, axis, indexing='ij')
        coords = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
        return coords, X.shape

    def calculate_geometric_cone(self, grid_coords, fe_pos, o_pos):
        global_grid = grid_coords + fe_pos
        
        dists_to_Fe = np.linalg.norm(global_grid - fe_pos, axis=1)
        valid_dist = (dists_to_Fe >= self.rmin) & (dists_to_Fe <= self.rmax)
        
        v_O_Fe = fe_pos - o_pos
        v_O_Fe_norm = v_O_Fe / (np.linalg.norm(v_O_Fe) + 1e-10)
        
        v_O_Grid = global_grid - o_pos
        norm_O_Grid = np.linalg.norm(v_O_Grid, axis=1)
        safe_norm = np.where(norm_O_Grid == 0, 1e-10, norm_O_Grid)
        
        cos_thetas = np.dot(v_O_Grid, v_O_Fe_norm) / safe_norm
        thetas = np.arccos(np.clip(cos_thetas, -1.0, 1.0))
        valid_angle = (thetas >= self.theta_min) & (thetas <= self.theta_max)
        
        return valid_dist & valid_angle, global_grid

    def calculate_theoretical_cone(self):
        log_print(f"--> Building Theoretical Grid (Resolution: {self.spacing} Å)...")
        grid_coords, grid_shape = self.generate_local_grid()
        
        fe_pos = np.array([0.0, 0.0, 0.0])
        o_pos = np.array([0.0, 0.0, 1.62])
        
        cone_mask, _ = self.calculate_geometric_cone(grid_coords, fe_pos, o_pos)
        acc_volume = np.sum(cone_mask) * self.voxel_vol

        grid_3d = cone_mask.reshape(grid_shape)
        exposed_faces = 0
        for ax in (0, 1, 2):
            exposed_faces += np.sum(grid_3d & ~np.roll(grid_3d, 1, axis=ax))
            exposed_faces += np.sum(grid_3d & ~np.roll(grid_3d, -1, axis=ax))
        
        acc_surface_area = exposed_faces * self.voxel_area
        return acc_volume, acc_surface_area

    def process_structure_frame(self, u, fe_sel_str, lig_sel_str, prot_sel_str):
        fe_atoms = u.select_atoms(fe_sel_str)
        if len(fe_atoms) == 0:
            raise ValueError(f"No Iron (Fe) atom found using selection: '{fe_sel_str}'")
        fe_pos = fe_atoms.positions[0]

        ox_atoms = u.select_atoms(f"name O O1 O2 and around 2.5 ({fe_sel_str})")
        if len(ox_atoms) > 0:
            dists = np.linalg.norm(ox_atoms.positions - fe_pos, axis=1)
            o_pos = ox_atoms.positions[np.argmin(dists)]
        else:
            heme_n = u.select_atoms(f"(same residue as ({fe_sel_str})) and name NA NB NC ND N1 N2 N3 N4")
            cys_s = u.select_atoms(f"name SG and around 3.5 ({fe_sel_str})")
            v_FeO_norm = np.array([0,0,1])
            
            if len(heme_n) >= 3:
                centered_N = heme_n.positions - fe_pos
                _, _, vh = np.linalg.svd(centered_N)
                plane_normal = vh[2]
                if len(cys_s) > 0:
                    v_FeS = cys_s.positions[0] - fe_pos
                    if np.dot(plane_normal, v_FeS) > 0:
                        plane_normal = -plane_normal
                v_FeO_norm = plane_normal
            elif len(cys_s) > 0:
                v_FeS = cys_s.positions[0] - fe_pos
                v_FeO_norm = -v_FeS / np.linalg.norm(v_FeS)
                
            o_pos = fe_pos + 1.62 * v_FeO_norm

        grid_coords, grid_shape = self.generate_local_grid()
        raw_cone_mask, global_grid_coords = self.calculate_geometric_cone(grid_coords, fe_pos, o_pos)
        
        raw_volume = np.sum(raw_cone_mask) * self.voxel_vol

        # Raw Surface Area
        raw_grid_3d = raw_cone_mask.reshape(grid_shape)
        raw_exposed_faces = 0
        for ax in (0, 1, 2):
            raw_exposed_faces += np.sum(raw_grid_3d & ~np.roll(raw_grid_3d, 1, axis=ax))
            raw_exposed_faces += np.sum(raw_grid_3d & ~np.roll(raw_grid_3d, -1, axis=ax))
        raw_surface_area = raw_exposed_faces * self.voxel_area

        # Protein Clash Masking
        prot_atoms = u.select_atoms(prot_sel_str)
        if len(prot_atoms) > 0:
            tree = cKDTree(prot_atoms.positions)
            clash_dists, _ = tree.query(global_grid_coords[raw_cone_mask], distance_upper_bound=self.prot_buffer)
            free_mask = clash_dists > self.prot_buffer
        else:
            free_mask = np.ones(np.sum(raw_cone_mask), dtype=bool)

        acc_volume = np.sum(free_mask) * self.voxel_vol

        # Accessible Surface Area
        full_accessible_mask = np.zeros(len(grid_coords), dtype=bool)
        cone_indices = np.where(raw_cone_mask)[0]
        full_accessible_mask[cone_indices[free_mask]] = True
        acc_grid_3d = full_accessible_mask.reshape(grid_shape)
        
        exposed_faces = 0
        for ax in (0, 1, 2):
            exposed_faces += np.sum(acc_grid_3d & ~np.roll(acc_grid_3d, 1, axis=ax))
            exposed_faces += np.sum(acc_grid_3d & ~np.roll(acc_grid_3d, -1, axis=ax))
        acc_surface_area = exposed_faces * self.voxel_area

        # Ligand Metrics
        lig_atoms = u.select_atoms(lig_sel_str)
        lig_vol_in_raw = 0.0
        lig_vol_in_acc = 0.0
        lig_total_vol = 0.0
        lig_total_surf = 0.0
        lig_surf_in_raw = 0.0
        lig_surf_in_acc = 0.0
        lig_com_dist = 0.0
        lig_com_angle = 0.0
        
        if len(lig_atoms) > 0:
            lig_com = lig_atoms.center_of_mass()
            lig_com_dist = np.linalg.norm(lig_com - fe_pos)
            
            v_O_Fe_norm = (fe_pos - o_pos) / (np.linalg.norm(fe_pos - o_pos) + 1e-10)
            v_O_Lig = lig_com - o_pos
            cos_lig = np.dot(v_O_Lig, v_O_Fe_norm) / (np.linalg.norm(v_O_Lig) + 1e-10)
            lig_com_angle = np.degrees(np.arccos(np.clip(cos_lig, -1.0, 1.0)))

            lig_tree = cKDTree(lig_atoms.positions)
            
            # Volumetric overlaps
            raw_valid_coords = global_grid_coords[raw_cone_mask]
            if len(raw_valid_coords) > 0:
                raw_lig_dists, _ = lig_tree.query(raw_valid_coords, distance_upper_bound=self.lig_buffer)
                lig_vol_in_raw = np.sum(raw_lig_dists <= self.lig_buffer) * self.voxel_vol
            
            acc_valid_coords = global_grid_coords[full_accessible_mask]
            if len(acc_valid_coords) > 0:
                acc_lig_dists, _ = lig_tree.query(acc_valid_coords, distance_upper_bound=self.lig_buffer)
                lig_vol_in_acc = np.sum(acc_lig_dists <= self.lig_buffer) * self.voxel_vol

            # Total Ligand Volume & Surface Area setup
            min_b = lig_atoms.positions.min(axis=0) - self.lig_buffer - self.spacing
            max_b = lig_atoms.positions.max(axis=0) + self.lig_buffer + self.spacing
            
            x_ax = np.arange(min_b[0], max_b[0], self.spacing)
            y_ax = np.arange(min_b[1], max_b[1], self.spacing)
            z_ax = np.arange(min_b[2], max_b[2], self.spacing)
            LX, LY, LZ = np.meshgrid(x_ax, y_ax, z_ax, indexing='ij')
            lig_grid = np.column_stack((LX.ravel(), LY.ravel(), LZ.ravel()))
            
            lig_dists_all, _ = lig_tree.query(lig_grid, distance_upper_bound=self.lig_buffer)
            lig_mask = lig_dists_all <= self.lig_buffer
            lig_total_vol = np.sum(lig_mask) * self.voxel_vol
            
            # Ligand Boundary (Surface) detection
            lig_grid_3d = lig_mask.reshape(LX.shape)
            lig_exposed_mask = np.zeros_like(lig_grid_3d, dtype=bool)
            for ax in (0, 1, 2):
                lig_exposed_mask |= (lig_grid_3d & ~np.roll(lig_grid_3d, 1, axis=ax))
                lig_exposed_mask |= (lig_grid_3d & ~np.roll(lig_grid_3d, -1, axis=ax))
                
            lig_total_surf = np.sum(lig_exposed_mask) * self.voxel_area

            # Find how much of the ligand surface is inside the cone
            boundary_coords = lig_grid[lig_exposed_mask.ravel()]
            if len(boundary_coords) > 0:
                boundary_raw_mask, _ = self.calculate_geometric_cone(boundary_coords - fe_pos, fe_pos, o_pos)
                lig_surf_in_raw = np.sum(boundary_raw_mask) * self.voxel_area
                
                if len(prot_atoms) > 0:
                    clash_dists_surf, _ = tree.query(boundary_coords[boundary_raw_mask], distance_upper_bound=self.prot_buffer)
                    boundary_acc_submask = clash_dists_surf > self.prot_buffer
                    lig_surf_in_acc = np.sum(boundary_acc_submask) * self.voxel_area
                else:
                    lig_surf_in_acc = lig_surf_in_raw

        return (raw_volume, raw_surface_area, acc_volume, acc_surface_area, 
                lig_vol_in_raw, lig_vol_in_acc, lig_total_vol, lig_total_surf, 
                lig_surf_in_raw, lig_surf_in_acc, lig_com_dist, lig_com_angle)

def plot_metrics(df, out_prefix, time_col="Frame"):
    sns.set_theme(style="whitegrid")
    
    # Plot 1: Cone Accessibility Timeseries
    fig, ax1 = plt.subplots(figsize=(10, 6), dpi=300)
    sns.lineplot(data=df, x=time_col, y="Accessible_Volume", hue="System", ax=ax1, palette="viridis", lw=2)
    ax1.set_ylabel("Accessible Cone Volume (Å³)", fontweight='bold')
    ax1.set_xlabel(time_col, fontweight='bold')
    
    ax2 = ax1.twinx()
    sns.lineplot(data=df, x=time_col, y="Accessible_Surface", hue="System", ax=ax2, palette="magma", linestyle="--", lw=1.5, alpha=0.7)
    ax2.set_ylabel("Accessible Surface Area (Å²)", fontweight='bold')
    
    plt.title("Structural Accessibility Over Time", fontweight='bold', pad=15)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_Accessibility_Timeseries.png")
    plt.close()

    # Plot 2: Ligand Volume and Area Timeseries
    if "Ligand_Total_Volume" in df.columns and df["Ligand_Total_Volume"].sum() > 0:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), dpi=300, sharex=True)
        
        # Top Panel: Volumes
        melted_vol = df.melt(id_vars=["System", time_col], 
                             value_vars=["Ligand_Total_Volume", "Lig_Vol_in_Acc_Cone"],
                             var_name="Metric", value_name="Volume (Å³)")
        sns.lineplot(data=melted_vol, x=time_col, y="Volume (Å³)", hue="Metric", style="System", ax=ax1, palette="mako", lw=2)
        ax1.set_title("Ligand Volume Over Time (Total vs Accessible Cone Overlap)", fontweight='bold')
        ax1.set_ylabel("Volume (Å³)", fontweight='bold')
        
        # Bottom Panel: Surfaces
        melted_surf = df.melt(id_vars=["System", time_col], 
                              value_vars=["Ligand_Total_Surface", "Lig_Surf_in_Acc_Cone"],
                              var_name="Metric", value_name="Area (Å²)")
        sns.lineplot(data=melted_surf, x=time_col, y="Area (Å²)", hue="Metric", style="System", ax=ax2, palette="rocket", lw=2)
        ax2.set_title("Ligand Surface Area Over Time (Total vs Accessible Cone Overlap)", fontweight='bold')
        ax2.set_ylabel("Surface Area (Å²)", fontweight='bold')
        ax2.set_xlabel(time_col, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_Ligand_Metrics_Timeseries.png")
        plt.close()

    # Plot 3: Ligand Penetration % Violin Plot
    if "Pct_Lig_Vol_in_Acc" in df.columns and df["Lig_Vol_in_Acc_Cone"].sum() > 0:
        plt.figure(figsize=(9, 6), dpi=300)
        sns.violinplot(data=df, x="System", y="Pct_Lig_Vol_in_Acc", palette="Blues", inner="quartile")
        plt.title("% Ligand Volume within Accessible Cone", fontweight='bold')
        plt.ylabel("Percent (%)", fontweight='bold')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_Ligand_Pct_Penetration.png")
        plt.close()


def main():
    with open("accessibility_calc.log", "w", encoding="utf-8") as f:
        f.write("CYP450 Accessibility Cone Log\n")
        f.write("============================================================\n")

    parser = argparse.ArgumentParser(description="CYP450 Accessibility Cone Calculator (Theoretical & Structural)")
    
    parser.add_argument("--top", nargs='+', help="Topology files (PDB or PRMTOP)")
    parser.add_argument("--traj", nargs='+', help="Trajectory files (NC, DCD, XTC)")
    
    parser.add_argument("--rmin", type=float, default=3.5, help="Minimum distance from Fe (Å)")
    parser.add_argument("--rmax", type=float, default=8.5, help="Maximum distance from Fe (Å)")
    parser.add_argument("--thetamin", type=float, default=100.0, help="Minimum angle from Fe=O axis (degrees)")
    parser.add_argument("--thetamax", type=float, default=145.0, help="Maximum angle from Fe=O axis (degrees)")
    
    parser.add_argument("--grid_spacing", type=float, default=0.5, help="Resolution of the 3D grid (Å)")
    parser.add_argument("--exclude_dist", type=float, default=2.5, help="Exclusion buffer distance around protein atoms (Å)")
    parser.add_argument("--lig_vdw", type=float, default=2.0, help="Estimated vdW radius for ligand atoms (Å)")
    parser.add_argument("--dt", type=float, default=1.0, help="Timestep per frame in ps")
    
    parser.add_argument("--fe_sel", type=str, default="resname HM1 CM1 FE1 HEM and name FE Fe", help="MDAnalysis string for Iron")
    parser.add_argument("--lig_sel", type=str, default="resname LIG or resname MOL or resname UNK or resname UNL", help="MDAnalysis string for Ligand")
    parser.add_argument("--prot_sel", type=str, default="protein", help="MDAnalysis string for Protein boundaries")
    
    args = parser.parse_args()

    model = AccessibilityModel(
        rmin=args.rmin, rmax=args.rmax, 
        theta_min=args.thetamin, theta_max=args.thetamax, 
        grid_spacing=args.grid_spacing, 
        prot_buffer=args.exclude_dist, 
        lig_buffer=args.lig_vdw
    )

    if not args.top:
        log_print("\n" + "="*60)
        log_print(" 🧪 RUNNING IN PURE THEORETICAL MODEL MODE")
        log_print("="*60)
        vol, surf = model.calculate_theoretical_cone()
        log_print("\n[RESULTS]")
        log_print(f"• Theoretical Accessible Volume:       {vol:.2f} Å³")
        log_print(f"• Theoretical Accessible Surface Area: {surf:.2f} Å²\n")
        sys.exit(0)

    if mda is None:
        log_print("[ERROR] MDAnalysis is required for structural calculations.")
        sys.exit(1)

    trajs = args.traj if args.traj else [None] * len(args.top)
    if len(args.top) != len(trajs):
        log_print("[ERROR] Number of topologies must match number of trajectories.")
        sys.exit(1)

    all_data = []

    for top, traj in zip(args.top, trajs):
        sys_name = os.path.basename(top).split('.')[0]
        log_print(f"\n---> Analyzing System: {sys_name}")
        
        try:
            u = mda.Universe(top, traj) if traj else mda.Universe(top)
        except Exception as e:
            log_print(f"[ERROR] Could not load universe for {sys_name}: {e}")
            continue
            
        test_lig_atoms = u.select_atoms(args.lig_sel)
        if len(test_lig_atoms) == 0:
            log_print(f"     [WARNING] 0 atoms found for ligand selection: '{args.lig_sel}'.")
        else:
            log_print(f"     [INFO] Found {len(test_lig_atoms)} ligand atoms matching '{args.lig_sel}'.")

        n_frames = len(u.trajectory)
        
        for ts in u.trajectory:
            if ts.frame % 10 == 0:
                log_print(f"     Processing Frame {ts.frame}/{n_frames}...", end='\r', log_to_file=False)
            
            try:
                (raw_vol, raw_surf, acc_vol, acc_surf, 
                 lig_raw, lig_acc, lig_tot_vol, lig_tot_surf, 
                 lig_surf_raw, lig_surf_acc, lig_dist, lig_ang) = model.process_structure_frame(u, args.fe_sel, args.lig_sel, args.prot_sel)
                
                all_data.append({
                    "System": sys_name,
                    "Frame": ts.frame,
                    "Time_ps": ts.frame * args.dt,
                    "Lig_COM_Dist_A": lig_dist,
                    "Lig_COM_Angle_deg": lig_ang,
                    "Raw_Cone_Volume": raw_vol,
                    "Raw_Cone_Surface": raw_surf,
                    "Accessible_Volume": acc_vol,
                    "Accessible_Surface": acc_surf,
                    "Ligand_Total_Volume": lig_tot_vol,
                    "Ligand_Total_Surface": lig_tot_surf,
                    "Lig_Vol_in_Raw_Cone": lig_raw,
                    "Lig_Vol_in_Acc_Cone": lig_acc,
                    "Lig_Surf_in_Raw_Cone": lig_surf_raw,
                    "Lig_Surf_in_Acc_Cone": lig_surf_acc
                })
            except Exception as e:
                log_print(f"\n[WARNING] Skipping Frame {ts.frame}: {e}")

        log_print(f"\n     Done with {sys_name}.")

    if not all_data:
        log_print("\n[ERROR] No data generated.")
        sys.exit(1)

    df = pd.DataFrame(all_data)
    
    # Calculate % values safely avoiding division by zero
    df["Pct_Acc_Vol"] = np.where(df["Raw_Cone_Volume"] > 0, (df["Accessible_Volume"] / df["Raw_Cone_Volume"]) * 100, 0)
    df["Pct_Acc_Surf"] = np.where(df["Raw_Cone_Surface"] > 0, (df["Accessible_Surface"] / df["Raw_Cone_Surface"]) * 100, 0)
    df["Pct_Lig_Vol_in_Acc"] = np.where(df["Ligand_Total_Volume"] > 0, (df["Lig_Vol_in_Acc_Cone"] / df["Ligand_Total_Volume"]) * 100, 0)
    df["Pct_Lig_Surf_in_Acc"] = np.where(df["Ligand_Total_Surface"] > 0, (df["Lig_Surf_in_Acc_Cone"] / df["Ligand_Total_Surface"]) * 100, 0)

    df.to_csv("COMBINED_ACCESSIBILITY_METRICS.csv", index=False)
    
    log_print("\n" + "="*145)
    log_print(" 📊 STRUCTURAL ACCESSIBILITY AVERAGES (Combined)")
    log_print("="*145)
    
    summary = df.groupby("System")[
        ["Lig_COM_Dist_A", "Lig_COM_Angle_deg", 
         "Raw_Cone_Volume", "Accessible_Volume", "Pct_Acc_Vol",
         "Raw_Cone_Surface", "Accessible_Surface", "Pct_Acc_Surf",
         "Ligand_Total_Volume", "Lig_Vol_in_Acc_Cone", "Pct_Lig_Vol_in_Acc",
         "Ligand_Total_Surface", "Lig_Surf_in_Acc_Cone", "Pct_Lig_Surf_in_Acc"]
    ].mean()
    
    # Abbreviate names strictly for clean terminal printing
    summary_print = summary.rename(columns={
        "Lig_COM_Dist_A": "LigDist",
        "Lig_COM_Angle_deg": "LigAng",
        "Raw_Cone_Volume": "RawVol",
        "Accessible_Volume": "AccVol",
        "Pct_Acc_Vol": "%AccVol",
        "Raw_Cone_Surface": "RawSurf",
        "Accessible_Surface": "AccSurf",
        "Pct_Acc_Surf": "%AccSurf",
        "Ligand_Total_Volume": "LigVol",
        "Lig_Vol_in_Acc_Cone": "LigInAccV",
        "Pct_Lig_Vol_in_Acc": "%LigInV",
        "Ligand_Total_Surface": "LigSurf",
        "Lig_Surf_in_Acc_Cone": "LigInAccS",
        "Pct_Lig_Surf_in_Acc": "%LigInS"
    })
    
    log_print(summary_print.round(1).to_string())
    log_print("="*145)

    log_print("\n---> Generating intuitive plots...")
    plot_col = "Time_ps" if args.traj else "Frame"
    plot_metrics(df, "CYP450", time_col=plot_col)
    log_print("[SUCCESS] Data, logs, and plots saved to current directory.")

if __name__ == "__main__":
    main()