import os
import sys
import subprocess
import threading
import shutil
import re
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from PIL import Image, ImageTk

class REACCESSUnifiedDashboard:
    def __init__(self, root):
        self.root = root
        self.root.title("REACCESS - Multi-Mode Unified Desktop Analytical Platform")
        self.root.geometry("1340x950")
        self.root.minsize(1150, 850)
        
        self.root.rowconfigure(2, weight=1)
        self.root.columnconfigure(0, weight=1)
        
        self.target_files_list = []
        self.active_output_dir = None
        
        self.apply_scientific_theme()
        self.create_header_banner()
        self.create_control_panel()
        self.create_workspace_notebook()
        self.create_status_bar()

    def apply_scientific_theme(self):
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
        style.configure('TNotebook', background="#cbd5e1", borderwidth=1)
        style.configure('TNotebook.Tab', font=('Helvetica', 10, 'bold'), background="#e2e8f0", foreground="#475569", padding=[14, 5])
        style.map('TNotebook.Tab', background=[('selected', '#ffffff')], foreground=[('selected', '#1e40af')])

    def create_header_banner(self):
        header_frame = tk.Frame(self.root, bg="#0f172a", padx=15, pady=12)
        header_frame.grid(row=0, column=0, sticky="ew")
        tk.Label(header_frame, text="REACCESS UNIFIED SUITE", font=("Helvetica", 20, "bold"), bg="#0f172a", fg="#38bdf8").pack(anchor="w")
        tk.Label(header_frame, text="Orchestrated Computational Interface for Phase I Site of Metabolism Screening via MD & Docking Analysis", font=("Helvetica", 11, "italic"), bg="#0f172a", fg="#94a3b8").pack(anchor="w", pady=(2, 0))

    def create_control_panel(self):
        ctrl_frame = ttk.LabelFrame(self.root, text=" Core Parameters & Asset Ingestion Panel ", padding=10)
        ctrl_frame.grid(row=1, column=0, sticky="ew", padx=12, pady=5)
        for i in range(6): ctrl_frame.columnconfigure(i, weight=1)

        ttk.Label(ctrl_frame, text="Engine Model:").grid(row=0, column=0, sticky="w", padx=5, pady=3)
        self.engine_combo = ttk.Combobox(ctrl_frame, values=["MD4 (5% Frame Occupancy)", "100ps (Cumulative Time Window)", "100pscont (Consecutive Frame Window)"], state="readonly", width=22)
        self.engine_combo.current(2)
        self.engine_combo.grid(row=0, column=1, sticky="w", padx=5)

        ttk.Label(ctrl_frame, text="Topology File (-t):").grid(row=0, column=2, sticky="w", padx=15)
        self.top_entry = ttk.Entry(ctrl_frame, width=20)
        self.top_entry.grid(row=0, column=3, sticky="ew", padx=5)
        ttk.Button(ctrl_frame, text="Browse...", command=lambda: self.browse_file(self.top_entry, [("Structure Files", "*.pdb *.prmtop *.gro"), ("All Files", "*.*")])).grid(row=0, column=4, sticky="w", padx=2)

        bounds_subframe = ttk.LabelFrame(ctrl_frame, text=" Cutoff Configuration Settings ")
        bounds_subframe.grid(row=0, column=5, rowspan=3, sticky="nsew", padx=10, pady=2)
        
        ttk.Label(bounds_subframe, text="Min Dist (Å):").grid(row=0, column=0, sticky="w", padx=3, pady=2)
        self.min_dist_entry = ttk.Entry(bounds_subframe, width=6); self.min_dist_entry.insert(0, "3.5"); self.min_dist_entry.grid(row=0, column=1, padx=3)
        
        ttk.Label(bounds_subframe, text="Max Dist (Å):").grid(row=1, column=0, sticky="w", padx=3, pady=2)
        self.max_dist_entry = ttk.Entry(bounds_subframe, width=6); self.max_dist_entry.insert(0, "8.5"); self.max_dist_entry.grid(row=1, column=1, padx=3)
        
        ttk.Label(bounds_subframe, text="Min Angle (°):").grid(row=2, column=0, sticky="w", padx=3, pady=2)
        self.min_ang_entry = ttk.Entry(bounds_subframe, width=6); self.min_ang_entry.insert(0, "100.0"); self.min_ang_entry.grid(row=2, column=1, padx=3)
        
        ttk.Label(bounds_subframe, text="Max Angle (°):").grid(row=3, column=0, sticky="w", padx=3, pady=2)
        self.max_ang_entry = ttk.Entry(bounds_subframe, width=6); self.max_ang_entry.insert(0, "145.0"); self.max_ang_entry.grid(row=3, column=1, padx=3)

        ttk.Label(bounds_subframe, text="BEP Alpha (α):").grid(row=4, column=0, sticky="w", padx=3, pady=2)
        self.alpha_entry = ttk.Entry(bounds_subframe, width=6); self.alpha_entry.insert(0, "0.5"); self.alpha_entry.grid(row=4, column=1, padx=3)

        ttk.Label(ctrl_frame, text="Ref Ligand PDB (-l):").grid(row=1, column=0, sticky="w", padx=5, pady=3)
        self.lig_entry = ttk.Entry(ctrl_frame, width=22)
        self.lig_entry.grid(row=1, column=1, sticky="ew", padx=5)
        ttk.Button(ctrl_frame, text="Browse...", command=lambda: self.browse_file(self.lig_entry, [("Ligand PDB", "*.pdb"), ("All Files", "*.*")])).grid(row=1, column=2, sticky="w", padx=2)

        ttk.Label(ctrl_frame, text="Ligand SMILES (-s):").grid(row=1, column=3, sticky="w", padx=15)
        self.smiles_entry = ttk.Entry(ctrl_frame, width=25)
        self.smiles_entry.grid(row=1, column=4, sticky="ew", padx=5)

        ttk.Label(ctrl_frame, text="Queue List:").grid(row=2, column=0, sticky="nw", padx=5, pady=8)
        list_container = ttk.Frame(ctrl_frame)
        list_container.grid(row=2, column=1, columnspan=3, sticky="nsew", padx=5, pady=3)
        list_container.rowconfigure(0, weight=1); list_container.columnconfigure(0, weight=1)
        self.files_listbox = tk.Listbox(list_container, height=4, font=("Helvetica", 9), background="#ffffff")
        self.files_listbox.grid(row=0, column=0, sticky="nsew")
        scrollbar = ttk.Scrollbar(list_container, orient="vertical", command=self.files_listbox.yview)
        scrollbar.grid(row=0, column=1, sticky="ns"); self.files_listbox.configure(yscrollcommand=scrollbar.set)
        
        list_btns = ttk.Frame(ctrl_frame)
        list_btns.grid(row=2, column=4, sticky="nw", padx=5, pady=3)
        ttk.Button(list_btns, text="Add Files...", command=self.add_individual_files_to_queue).pack(fill="x", pady=2)
        ttk.Button(list_btns, text="Add Folder...", command=self.add_entire_directory_to_queue).pack(fill="x", pady=2)
        ttk.Button(list_btns, text="Clear Queue", command=self.clear_file_listbox_queue).pack(fill="x", pady=2)

        ttk.Label(ctrl_frame, text="Destination Folder:").grid(row=3, column=0, sticky="w", padx=5, pady=3)
        self.out_entry = ttk.Entry(ctrl_frame, width=22)
        self.out_entry.grid(row=3, column=1, sticky="ew", padx=5)
        ttk.Button(ctrl_frame, text="Set Dir...", command=self.browse_output_directory).grid(row=3, column=2, sticky="w", padx=2)

        action_row = ttk.Frame(ctrl_frame)
        action_row.grid(row=3, column=3, columnspan=3, pady=5, sticky="e")
        self.launch_btn = ttk.Button(action_row, text="LAUNCH PIPELINE SCREEN", style="Action.TButton", command=self.start_pipeline_worker_thread)
        self.launch_btn.pack(side="left", padx=10)
        ttk.Button(action_row, text="LOAD PRESETS", command=self.load_built_in_demo_presets).pack(side="left", padx=5)

    def create_workspace_notebook(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.grid(row=2, column=0, sticky="nsew", padx=12, pady=5)

        self.log_tab = ttk.Frame(self.notebook); self.notebook.add(self.log_tab, text=" Execution Console Log ")
        self.log_text = tk.Text(self.log_tab, wrap="word", background="#1e293b", foreground="#f8fafc", font=("Courier", 10))
        self.log_text.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        self.log_tab.rowconfigure(0, weight=1); self.log_tab.columnconfigure(0, weight=1)

        self.table_tab = ttk.Frame(self.notebook); self.notebook.add(self.table_tab, text=" Meta Spreadsheet Summary ")
        self.table_tab.rowconfigure(0, weight=1); self.table_tab.columnconfigure(0, weight=1)
        grid_frame = ttk.Frame(self.table_tab); grid_frame.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        grid_frame.rowconfigure(0, weight=1); grid_frame.columnconfigure(0, weight=1)
        self.tree = ttk.Treeview(grid_frame, show="headings"); self.tree.grid(row=0, column=0, sticky="nsew")
        v_sc = ttk.Scrollbar(grid_frame, orient="vertical", command=self.tree.yview); v_sc.grid(row=0, column=1, sticky="ns"); self.tree.configure(yscrollcommand=v_sc.set)
        ttk.Button(self.table_tab, text="Save Excel Spreadsheet Asset (.xlsx)", command=self.download_excel_action).grid(row=1, column=0, pady=4)

        self.image_tab = ttk.Frame(self.notebook); self.notebook.add(self.image_tab, text=" Labeled 2D ACS Structure ")
        self.image_tab.rowconfigure(0, weight=1); self.image_tab.columnconfigure(0, weight=1)
        self.canvas_label = ttk.Label(self.image_tab, text="Awaiting tracking execution outputs...", font=('Helvetica', 11, 'italic'), anchor="center")
        self.canvas_label.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        ttk.Button(self.image_tab, text="Save Structural Image Panel (.png)", command=self.download_png_action).grid(row=1, column=0, pady=4)

        self.plots_tab = ttk.Frame(self.notebook); self.notebook.add(self.plots_tab, text=" Aggregated Multi-System Charts ")
        self.plots_tab.rowconfigure(0, weight=1)
        for i in range(3): self.plots_tab.columnconfigure(i, weight=1)
        self.p1_lbl = ttk.Label(self.plots_tab, anchor="center"); self.p1_lbl.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        self.p2_lbl = ttk.Label(self.plots_tab, anchor="center"); self.p2_lbl.grid(row=0, column=1, sticky="nsew", padx=2, pady=2)
        self.p3_lbl = ttk.Label(self.plots_tab, anchor="center"); self.p3_lbl.grid(row=0, column=2, sticky="nsew", padx=2, pady=2)

    def create_status_bar(self):
        status_frame = tk.Frame(self.root, bg="#cbd5e1", padx=10, pady=5)
        status_frame.grid(row=3, column=0, sticky="ew")
        self.status_lbl = tk.Label(status_frame, text="Status: Ready.", font=("Helvetica", 9), bg="#cbd5e1", fg="#334155")
        self.status_lbl.pack(side="left")
        self.progress_bar = ttk.Progressbar(status_frame, orient="horizontal", mode="determinate")
        self.progress_bar.pack(side="right", fill="x", expand=True, padx=10)

    def browse_file(self, entry_field, file_types):
        fn = filedialog.askopenfilename(filetypes=file_types)
        if fn: entry_field.delete(0, tk.END); entry_field.insert(0, fn)

    def browse_output_directory(self):
        dn = filedialog.askdirectory()
        if dn: self.out_entry.delete(0, tk.END); self.out_entry.insert(0, dn)

    def add_individual_files_to_queue(self):
        fns = filedialog.askopenfilenames(filetypes=[("Structural Run Files", "*.pdb *.nc *.dcd *.xtc"), ("All Files", "*.*")])
        for f in fns:
            if f not in self.target_files_list:
                self.target_files_list.append(f); self.files_listbox.insert(tk.END, os.path.basename(f))

    def add_entire_directory_to_queue(self):
        dn = filedialog.askdirectory()
        if dn:
            matches = sorted(glob.glob(os.path.join(dn, "*.nc"))) + sorted(glob.glob(os.path.join(dn, "*.dcd"))) + sorted(glob.glob(os.path.join(dn, "*.xtc"))) + sorted(glob.glob(os.path.join(dn, "*.pdb")))
            for f in matches:
                if f not in self.target_files_list:
                    self.target_files_list.append(f); self.files_listbox.insert(tk.END, f"[Folder Asset] -> {os.path.basename(f)}")

    def clear_file_listbox_queue(self):
        self.target_files_list.clear(); self.files_listbox.delete(0, tk.END)

    def load_built_in_demo_presets(self):
        self.top_entry.insert(0, "5TE8-A.pdb"); self.lig_entry.insert(0, "cetrizine-ps4.pdb")
        self.smiles_entry.insert(0, "C1CN(CCN1CCOCC(=O)O)C(C2=CC=CC=C2)C3=CC=C(C=C3)Cl")
        self.out_entry.insert(0, "./REACCESS_METASCREEN_OUT")

    def start_pipeline_worker_thread(self):
        if not self.target_files_list or not self.top_entry.get() or not self.lig_entry.get() or not self.out_entry.get():
            messagebox.showerror("Error", "Missing configuration fields entry items setup bounds.")
            return
        self.launch_btn.config(state="disabled"); self.log_text.delete("1.0", tk.END); self.notebook.select(0)
        threading.Thread(target=self.execute_orchestration_subprocess, daemon=True).start()

    def execute_orchestration_subprocess(self):
        # REMAPPED: Targeting the renamed routing scheduler engine
        wrapper_script = os.path.abspath("reaccess_wrapper_engine_dynbde.py")
        engine_raw = self.engine_combo.get()
        engine_mode = "MD4" if "MD4" in engine_raw else ("100ps" if "100ps (" in engine_raw else "100pscont")
        self.active_output_dir = os.path.abspath(self.out_entry.get())
        
        cmd = [
            sys.executable, "-u", wrapper_script, "-e", engine_mode,
            "-t", os.path.abspath(self.top_entry.get()), "-l", os.path.abspath(self.lig_entry.get()),
            "-s", self.smiles_entry.get().strip(), "-o", self.active_output_dir,
            "--min_dist", self.min_dist_entry.get().strip(), "--max_dist", self.max_dist_entry.get().strip(),
            "--min_angle", self.min_ang_entry.get().strip(), "--max_angle", self.max_ang_entry.get().strip(),
            "--alpha", self.alpha_entry.get().strip(), "-f"
        ] + self.target_files_list

        self.root.after(0, lambda: self.status_lbl.config(text="Status: Processing screening configurations..."))
        try:
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
            regex_sys = re.compile(r"PROCESSING TARGET SYSTEM\s*\((\d+)/(\d+)\)")
            for line in process.stdout:
                self.root.after(0, lambda l=line: self.log_text.insert(tk.END, l))
                self.root.after(0, lambda: self.log_text.see(tk.END))
                m_sys = regex_sys.search(line)
                if m_sys:
                    c, t = int(m_sys.group(1)), int(m_sys.group(2))
                    self.root.after(0, lambda v=int(((c-1)/t)*95): self.progress_bar.config(value=v))
            process.wait()
            if process.returncode == 0:
                self.root.after(0, lambda: self.progress_bar.config(value=100))
                self.root.after(0, lambda: self.status_lbl.config(text="Status: Complete."))
                self.root.after(0, self.render_completed_meta_dashboards)
            else:
                self.root.after(0, lambda: self.status_lbl.config(text="Status: Pipeline Aborted."))
                self.root.after(0, lambda: self.launch_btn.config(state="normal"))
        except Exception as err:
            self.root.after(0, lambda: self.launch_btn.config(state="normal"))

    def render_completed_meta_dashboards(self):
        self.launch_btn.config(state="normal")
        summary_excel = os.path.join(self.active_output_dir, "COMBINED_METABOLIC_MASTER_SUMMARY.xlsx")
        master_png = os.path.join(self.active_output_dir, "COMBINED_MASTER_2D_ACS_MAP.png")
        plots_names = ["COMBINED_MASTER_PLOT_1_distance_vs_BDE.png", "COMBINED_MASTER_PLOT_2_angle_vs_BDE.png", "COMBINED_MASTER_PLOT_3_3D_landscape_BDE.png"]

        if os.path.exists(summary_excel):
            df = pd.read_excel(summary_excel, sheet_name="Per_System_Summary")
            self.tree["columns"] = list(df.columns)
            for c in df.columns:
                self.tree.heading(c, text=c); self.tree.column(c, width=120, anchor="center")
            for item in self.tree.get_children(): self.tree.delete(item)
            for _, r in df.iterrows(): self.tree.insert("", tk.END, values=list(r))

        if os.path.exists(master_png):
            im2d = Image.open(master_png); im2d.thumbnail((650, 650), Image.Resampling.LANCZOS)
            ph2d = ImageTk.PhotoImage(im2d); self.canvas_label.config(image=ph2d, text=""); self.canvas_label.image = ph2d

        lbls = [self.p1_lbl, self.p2_lbl, self.p3_lbl]
        for name, lbl in zip(plots_names, lbls):
            path_p = os.path.join(self.active_output_dir, name)
            if os.path.exists(path_p):
                imp = Image.open(path_p); imp.thumbnail((340, 340), Image.Resampling.LANCZOS)
                php = ImageTk.PhotoImage(imp); lbl.config(image=php); lbl.image = php
        self.notebook.select(1)

    def download_excel_action(self):
        if self.active_output_dir: self.export_file_dialog_copy(os.path.join(self.active_output_dir, "COMBINED_METABOLIC_MASTER_SUMMARY.xlsx"), "COMBINED_METABOLIC_MASTER_SUMMARY.xlsx", [("Excel Sheet", "*.xlsx")], ".xlsx")

    def download_png_action(self):
        if self.active_output_dir: self.export_file_dialog_copy(os.path.join(self.active_output_dir, "COMBINED_MASTER_2D_ACS_MAP.png"), "COMBINED_MASTER_2D_ACS_MAP.png", [("Vector Graphic", "*.png")], ".png")

    def export_file_dialog_copy(self, source, init_file, types, ext):
        if not os.path.exists(source): return
        dest = filedialog.asksaveasfilename(initialfile=init_file, filetypes=types, defaultextension=ext)
        if dest: shutil.copy(source, dest); messagebox.showinfo("Export Success", "Asset copy generated clean.")

if __name__ == "__main__":
    root = tk.Tk(); app = REACCESSUnifiedDashboard(root); root.mainloop()