=======================================================================
               REACCESS INTUITIVE TRAJECTORY VISUALIZER
=======================================================================
Script: plot_intuitive_metrics.py

This command-line tool parses multi-trajectory Molecular Dynamics (MD) 
outputs (CSV or Excel) and generates publication-quality intuitive metrics. 
It translates complex geometric and quantum parameters into clean activity 
barcodes, continuous time-hold bar charts, and 2D spatial landscapes.

-----------------------------------------------------------------------
1. REQUIRED DEPENDENCIES
-----------------------------------------------------------------------
If not already in your REACCESS environment, install:
    conda install numpy pandas matplotlib seaborn openpyxl -y

-----------------------------------------------------------------------
2. COMPLETE COMMAND-LINE OPTIONS
-----------------------------------------------------------------------
Position 1: [FILE] 
    Path to your summary file (e.g., COMBINED_METABOLIC_MASTER_SUMMARY.csv)

Time & Hold Parameters:
    --dt            Time duration of a single frame in ps (Default: 1.0)
    --threshold     Strict continuous time hold requirement in ps (Default: 100.0)

Geometric & Quantum Cutoffs:
    --min_dist      Minimum Fe-Atom distance in Å (Default: 3.5)
    --max_dist      Maximum Fe-Atom distance in Å (Default: 8.5)
    --min_angle     Minimum Cone Angle in degrees (Default: 100.0)
    --max_angle     Maximum Cone Angle in degrees (Default: 145.0)
    --max_bde       Maximum acceptable BDE in kcal/mol (Default: 95.0)

Routing & Selection Options:
    --trajectories  Which trajectories to plot. 
                    Options: 'all', 'combined_only', or specific IDs like 'run_1 run_2'
                    (Default: all)
    --plots         Which specific graphics to generate. 
                    Options: 'all', 'barcode', 'barchart', 'dist_angle', 'dist_bde', 'angle_bde'
                    (Default: all)

-----------------------------------------------------------------------
3. EXHAUSTIVE EXAMPLES
-----------------------------------------------------------------------
A. The Standard Run (Defaults to dt=1.0, threshold=100)
   python plot_intuitive_metrics.py COMBINED_METABOLIC_MASTER_SUMMARY.csv

B. Custom Time Step & Hold Threshold (e.g., 10 ps stride, 50 ps hold)
   python plot_intuitive_metrics.py summary.csv --dt 10.0 --threshold 50.0

C. Modifying Geometric and BDE Cutoffs 
   (e.g., stricter distance 4.0 to 7.0 A, and stricter BDE of 90 kcal/mol)
   python plot_intuitive_metrics.py summary.csv --min_dist 4.0 --max_dist 7.0 --max_bde 90.0

D. "Summary Only" (Saves memory by NOT plotting individual runs)
   python plot_intuitive_metrics.py summary.csv --trajectories combined_only --plots all

E. "Specific Trajectories Only" (Only analyze Run 1 and Run 3)
   python plot_intuitive_metrics.py summary.csv --trajectories run_1 run_3

F. "Specific Plots Only" (Generate ONLY the barcode timeline and the Barchart to save space)
   python plot_intuitive_metrics.py summary.csv --plots barcode barchart

G. The "Total Customization" Command
   python plot_intuitive_metrics.py "COMBINED_METABOLIC_MASTER_SUMMARY.csv"        --dt 10.0        --threshold 100.0        --min_dist 3.0        --max_dist 8.0        --max_bde 92.5        --trajectories combined_only run_2        --plots barcode dist_bde

-----------------------------------------------------------------------
4. UNDERSTANDING THE OUTPUTS (Found in ./INTUITIVE_PLOTS_OUT/)
-----------------------------------------------------------------------
1_Metabolic_Barcode: 
    Heatmap timeline (Green = Active, Gray = Inactive). Identifies if an 
    atom is stable in the pocket or violently vibrating in and out.

2_Consecutive_Time: 
    Bar chart of the maximum unbroken frames spent active. Bars are Green 
    if they pass your --threshold, Red if they fail.

3A/3B/3C_Landscapes: 
    Grid of scatter plots per atom. The "Goldilocks Zone" is shaded green. 
    Instantly diagnoses WHY an atom failed (e.g., distance was fine but 
    angle shifted too high).
=======================================================================
