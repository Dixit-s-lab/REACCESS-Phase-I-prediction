import os
import sys
import glob
import shutil
import subprocess
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from PIL import Image, ImageTk

class MetaboliteGeneratorDashboard:
    def __init__(self, root):
        self.root = root
        self.root.title("REACCESS - Phase I Metabolite Explorer Dashboard")
        self.root.geometry("1340x950")
        self.root.minsize(1150, 850)
        
        self.root.rowconfigure(2, weight=1)
        self.root.columnconfigure(0, weight=1)
        
        self.active_output_dir = ""
        self.image_references = {}
        
        self.apply_scientific_theme()
        self.create_header_banner()
        self.create_control_panel()
        self.create_scrollable_viewport()
        
    def apply_scientific_theme(self):
        """Applies a uniform scientific color theme matching the primary REACCESS GUI."""
        self.root.configure(bg="#eef2f7")
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('.', background="#eef2f7", foreground="#1e293b")
        style.configure('TFrame', background="#eef2f7")
        style.configure('TLabelframe', background="#eef2f7", bordercolor="#cbd5e1", borderwidth=1)
        style.configure('TLabelframe.Label', background="#eef2f7", font=('Helvetica', 10, 'bold'), foreground="#1e40af")
        style.configure('TLabel', background="#eef2f7", font=('Helvetica', 10), foreground="#334155")
        style.configure('TEntry', fieldbackground="#ffffff", bordercolor="#cbd5e1", font=('Helvetica', 10))
        style.configure('TButton', font=('Helvetica', 9, 'bold'), foreground="#1e40af", background="#bae6fd", bordercolor="#0284c7")
        style.map('TButton', background=[('active', '#7dd3fc'), ('disabled', '#e2e8f0')])
        style.configure('Action.TButton', font=('Helvetica', 11, 'bold'), foreground="#ffffff", background="#0284c7", bordercolor="#0369a1")
        style.map('Action.TButton', background=[('active', '#0369a1'), ('disabled', '#cbd5e1')])

    def create_header_banner(self):
        """Creates a high-contrast top identity navigation banner."""
        header_frame = tk.Frame(self.root, bg="#0f172a", padx=15, pady=12)
        header_frame.grid(row=0, column=0, sticky="ew")
        tk.Label(header_frame, text="PHASE I METABOLITE GENERATOR LAYER", font=("Helvetica", 20, "bold"), bg="#0f172a", fg="#38bdf8").pack(anchor="w")
        tk.Label(header_frame, text="Automated Skeletal Expansion & High-Fidelity 3D Conformer Generation Suite", font=("Helvetica", 11, "italic"), bg="#0f172a", fg="#94a3b8").pack(anchor="w", pady=(2, 0))

    def create_control_panel(self):
        """Builds a scannable parameter input and ingestion dashboard layout."""
        ctrl_frame = ttk.LabelFrame(self.root, text=" Core Parameters & Asset Ingestion Panel ", padding=10)
        ctrl_frame.grid(row=1, column=0, sticky="ew", padx=12, pady=5)
        for i in range(5): ctrl_frame.columnconfigure(i, weight=1)

        # Row 0: Summary Spreadsheet Configurations
        ttk.Label(ctrl_frame, text="REACCESS Excel File:").grid(row=0, column=0, sticky="w", padx=5, pady=4)
        self.excel_entry = ttk.Entry(ctrl_frame, width=50)
        self.excel_entry.grid(row=0, column=1, columnspan=2, sticky="ew", padx=5, pady=4)
        ttk.Button(ctrl_frame, text="Browse...", command=self.browse_excel).grid(row=0, column=3, sticky="w", padx=2, pady=4)
        
        # High-Contrast Built-In Demo Presets Configurator Card
        preset_subframe = ttk.LabelFrame(ctrl_frame, text=" Demo Presets Engine ")
        preset_subframe.grid(row=0, column=4, rowspan=3, sticky="nsew", padx=10, pady=2)
        ttk.Button(preset_subframe, text="Load MD Meta-Summary", command=self.load_md_preset).pack(fill="x", padx=8, pady=4)
        ttk.Button(preset_subframe, text="Load Docking Summary", command=self.load_docking_preset).pack(fill="x", padx=8, pady=4)

        # Row 1: Topology Reference Inputs
        ttk.Label(ctrl_frame, text="Ref Ligand PDB (-l):").grid(row=1, column=0, sticky="w", padx=5, pady=4)
        self.lig_entry = ttk.Entry(ctrl_frame, width=50)
        self.lig_entry.grid(row=1, column=1, columnspan=2, sticky="ew", padx=5, pady=4)
        ttk.Button(ctrl_frame, text="Browse...", command=self.browse_ligand).grid(row=1, column=3, sticky="w", padx=2, pady=4)
        
        # Row 2: SMILES Verification Input Layouts
        ttk.Label(ctrl_frame, text="Ligand SMILES String:").grid(row=2, column=0, sticky="w", padx=5, pady=4)
        self.smiles_entry = ttk.Entry(ctrl_frame, width=50)
        self.smiles_entry.grid(row=2, column=1, columnspan=3, sticky="ew", padx=5, pady=4)

        # Row 3: Output Destination Paths Setup
        ttk.Label(ctrl_frame, text="Output Directory (-o):").grid(row=3, column=0, sticky="w", padx=5, pady=4)
        self.out_entry = ttk.Entry(ctrl_frame, width=50)
        self.out_entry.grid(row=3, column=1, columnspan=2, sticky="ew", padx=5, pady=4)
        ttk.Button(ctrl_frame, text="Set Dir...", command=self.browse_output).grid(row=3, column=3, sticky="w", padx=2, pady=4)
        
        self.run_btn = ttk.Button(ctrl_frame, text="LAUNCH GENERATION PIPELINE", style="Action.TButton", command=self.start_generation_thread)
        self.run_btn.grid(row=3, column=4, padx=10, pady=4, sticky="ew")

    def create_scrollable_viewport(self):
        """Builds a scrollable window grid optimized for displaying multi-system structural profile graphs."""
        self.view_container = ttk.Frame(self.root, padding=5)
        self.view_container.grid(row=2, column=0, sticky="nsew", padx=12, pady=5)
        self.view_container.rowconfigure(0, weight=1)
        self.view_container.columnconfigure(0, weight=1)
        
        self.canvas = tk.Canvas(self.view_container, bg="#ffffff", highlightthickness=1, highlightbackground="#cbd5e1")
        scrollbar = ttk.Scrollbar(self.view_container, orient="vertical", command=self.canvas.yview)
        
        self.scrollable_frame = ttk.Frame(self.canvas)
        self.scrollable_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=scrollbar.set)
        
        self.canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        
        self.placeholder_lbl = ttk.Label(self.scrollable_frame, text="Configure ingestion parameters above and launch pipeline to track generated targets...", font=("Helvetica", 11, "italic"), foreground="#64748b")
        self.placeholder_lbl.pack(pady=50, padx=50, anchor="center")

    def browse_excel(self):
        fn = filedialog.askopenfilename(filetypes=[("Excel Spreadsheet", "*.xlsx"), ("All Files", "*.*")])
        if fn: self.excel_entry.delete(0, tk.END); self.excel_entry.insert(0, fn)

    def browse_ligand(self):
        fn = filedialog.askopenfilename(filetypes=[("Ligand PDB Reference", "*.pdb"), ("All Files", "*.*")])
        if fn: self.lig_entry.delete(0, tk.END); self.lig_entry.insert(0, fn)

    def browse_output(self):
        dn = filedialog.askdirectory()
        if dn: self.out_entry.delete(0, tk.END); self.out_entry.insert(0, dn)

    def load_md_preset(self):
        self.clear_fields()
        self.excel_entry.insert(0, "./REACCESS_METASCREEN_OUT/COMBINED_METABOLIC_MASTER_SUMMARY.xlsx")
        self.lig_entry.insert(0, "cetrizine-ps4.pdb")
        self.smiles_entry.insert(0, "C1CN(CCN1CCOCC(=O)O)C(C2=CC=CC=C2)C3=CC=C(C=C3)Cl")
        self.out_entry.insert(0, "./AGGREGATED_SCREEN_METABOLITES")

    def load_docking_preset(self):
        self.clear_fields()
        self.excel_entry.insert(0, "./REACCESS_METASCREEN_OUT/REACCESS_OUT_cetrizine-ps4/METABOLIC_FINAL_SUMMARY.xlsx")
        self.lig_entry.insert(0, "cetrizine-ps4.pdb")
        self.smiles_entry.insert(0, "C1CN(CCN1CCOCC(=O)O)C(C2=CC=CC=C2)C3=CC=C(C=C3)Cl")
        self.out_entry.insert(0, "./PHASE_I_METABOLITES_POOL")

    def clear_fields(self):
        self.excel_entry.delete(0, tk.END)
        self.lig_entry.delete(0, tk.END)
        self.smiles_entry.delete(0, tk.END)
        self.out_entry.delete(0, tk.END)

    def start_generation_thread(self):
        if not self.excel_entry.get() or not self.lig_entry.get() or not self.out_entry.get():
            messagebox.showerror("Configuration Error", "All folder path selection fields must be populated.")
            return
        self.run_btn.config(state="disabled")
        threading.Thread(target=self.run_backend_subprocess, daemon=True).start()

    def run_backend_subprocess(self):
        self.active_output_dir = os.path.abspath(self.out_entry.get())
        cmd = [
            sys.executable, "metabolite_generator.py",
            "-excel", os.path.abspath(self.excel_entry.get()),
            "-l", os.path.abspath(self.lig_entry.get()),
            "-s", self.smiles_entry.get().strip(),
            "-o", self.active_output_dir
        ]
        
        try:
            process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if process.returncode == 0:
                self.root.after(0, self.render_structural_viewport_deck)
            else:
                self.root.after(0, lambda: messagebox.showerror("Subprocess Error", f"Execution failed:\n{process.stdout}"))
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Subprocess Error", f"Failed to start backend process:\n{str(e)}"))
        finally:
            self.root.after(0, lambda: self.run_btn.config(state="normal"))

    def render_structural_viewport_deck(self):
        """Sequentially loads parent layouts and flat skeletal ACS-style product files into a scrollable viewer."""
        for child in self.scrollable_frame.winfo_children():
            child.destroy()
        self.image_references.clear()
        
        parent_img = os.path.join(self.active_output_dir, "parent_structure_reference.png")
        
        # 1. Render High-Contrast Reference Parent Layout Card First
        if os.path.exists(parent_img):
            p_card = tk.Frame(self.scrollable_frame, bg="#ffffff", bd=1, relief="solid", padx=15, pady=15)
            p_card.pack(fill="x", padx=20, pady=15)
            
            tk.Label(p_card, text="Reference Parent Structure Map (Identified SOM Centers Highlighted)", font=("Helvetica", 12, "bold"), bg="#ffffff", fg="#1e40af").pack(anchor="w", pady=(0, 5))
            
            img_p = Image.open(parent_img).resize((400, 400), Image.Resampling.LANCZOS)
            ph_p = ImageTk.PhotoImage(img_p)
            self.image_references["parent"] = ph_p
            
            tk.Label(p_card, image=ph_p, bg="#ffffff").pack(pady=5)
            
            btn_bar = tk.Frame(p_card, bg="#ffffff")
            btn_bar.pack(anchor="e", pady=(5, 0))
            tk.Label(btn_bar, text="Original Parent Asset Source Files:", font=("Helvetica", 9, "italic"), bg="#ffffff", fg="#64748b").pack(side="left", padx=10)
            
            # Stylized high-contrast file loading action layouts button overrides
            pb = tk.Button(btn_bar, text="Export Source PDB File", font=("Helvetica", 9, "bold"), bg="#bae6fd", fg="#1e40af", bd=1, relief="groove", padx=8,
                           command=lambda: self.export_asset(os.path.abspath(self.lig_entry.get()), "parent_reference.pdb", [("Structure PDB", "*.pdb")]))
            pb.pack(side="left", padx=2)
        
        # 2. Render Generated Flat Skeletal Phase I Products Serially Underneath
        png_patterns = sorted(glob.glob(os.path.join(self.active_output_dir, "parent-*.png")))
        
        if not png_patterns:
            ttk.Label(self.scrollable_frame, text="No active structural metabolites were generated based on spreadsheet parameters.", font=("Helvetica", 10, "italic")).pack(pady=20)
            return
            
        for path_png in png_patterns:
            metabolite_id = os.path.splitext(os.path.basename(path_png))[0]
            
            m_card = tk.Frame(self.scrollable_frame, bg="#ffffff", bd=1, relief="solid", padx=15, pady=15)
            m_card.pack(fill="x", padx=20, pady=10)
            
            tk.Label(m_card, text=f"Metabolite Profile ID: {metabolite_id}", font=("Helvetica", 12, "bold"), bg="#ffffff", fg="#047857").grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 5))
            
            # Left Subcolumn: Flat Skeletal Mapping Image
            img_m = Image.open(path_png).resize((350, 350), Image.Resampling.LANCZOS)
            ph_m = ImageTk.PhotoImage(img_m)
            self.image_references[metabolite_id] = ph_m
            
            tk.Label(m_card, image=ph_m, bg="#ffffff").grid(row=1, column=0, padx=10, pady=5, sticky="w")
            
            # Right Subcolumn: Coordinates Export Action Controls Layout Card Panel
            actions_frame = tk.Frame(m_card, bg="#f8fafc", bd=1, relief="solid", padx=15, pady=15)
            actions_frame.grid(row=1, column=1, padx=20, pady=5, sticky="nsew")
            m_card.columnconfigure(1, weight=1)
            
            tk.Label(actions_frame, text="Export Computed 3D Coordinates", font=("Helvetica", 10, "bold"), bg="#f8fafc", fg="#334155").pack(anchor="w", pady=(0, 10))
            
            sdf_source = os.path.join(self.active_output_dir, f"{metabolite_id}.sdf")
            pdb_source = os.path.join(self.active_output_dir, f"{metabolite_id}.pdb")
            
            # Clean blue/green high-contrast download configuration setups buttons map layout controls
            sb = tk.Button(actions_frame, text="Download Metabolite SDF (.sdf)", font=("Helvetica", 9, "bold"), bg="#a7f3d0", fg="#047857", bd=1, relief="groove", width=30, pady=4,
                           command=lambda s=sdf_source, n=metabolite_id: self.export_asset(s, f"{n}.sdf", [("Structure SDF", "*.sdf")]))
            sb.pack(pady=4)
            
            tb = tk.Button(actions_frame, text="Download Metabolite PDB (.pdb)", font=("Helvetica", 9, "bold"), bg="#bfdbfe", fg="#1d4ed8", bd=1, relief="groove", width=30, pady=4,
                           command=lambda p=pdb_source, n=metabolite_id: self.export_asset(p, f"{n}.pdb", [("Structure PDB", "*.pdb")]))
            tb.pack(pady=4)
            
        self.canvas.yview_moveto(0)

    def export_asset(self, source_path, default_name, file_types):
        if not os.path.exists(source_path):
            messagebox.showerror("Export Failure", "Target computational asset file was missing or failed to converge.")
            return
        dest = filedialog.asksaveasfilename(initialfile=default_name, filetypes=file_types, defaultextension=file_types[0][1][1:])
        if dest:
            shutil.copy(source_path, dest)
            messagebox.showinfo("Export Complete", f"Asset successfully compiled at location:\n--> {os.path.basename(dest)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = MetaboliteGeneratorDashboard(root)
    root.mainloop()