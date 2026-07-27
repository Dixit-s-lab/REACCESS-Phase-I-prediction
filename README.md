# REACCESS-Phase-I-prediction
This method uses a combination of reactivity and accessibility approach to predict potential SOM for CYP450 catalyzed metabolism of drug-like compounds.

Once all the files are present in a working directory, follow the steps given below.

1) Install the enviornment called "reaccess-env" using the command

conda create --file reaccess-env.yaml

2) Once the enviornment is installed activate it using
   conda activate reaccess-env

3) Launch the REACCESS GUI via command line using
   python3 reaccess_master_gui_dynbde.py

   Then load test files from the reaccess-test-files and you can run the calculation by "LAUNCH PIPELINE SCREEN"

4) Using command line to run REACCESS calculations

   make sure that you have all the files in the current working directory or you can run the analysis from anywhere by giving full path to the script and necessary files.
   Below is a sample command that will run the analysis on MD trajector and single docked structure (separate apo protein and docked ligand PDB files)

   python reaccess_wrapper_engine_dynbde.py -t protein.pdb -l lig1.pdb -r 2.7 -o cet-p1.xlsx -s "OC(=O)COCCN1CCN(CC1)C(C2=CC=CC=C2)C3=CC=C(Cl)C=C3"
   OR   
   python bde-mopac-analyzer100pscont-dynbde.py -t protein.pdb -l lig1.pdb -r 2.7 -o cet-p1.xlsx -s "OC(=O)COCCN1CCN(CC1)C(C2=CC=CC=C2)C3=CC=C(Cl)C=C3"
   OR
   Other bde-mopac-analyzer*.py scripts depending on the type of analyses you need.
   
   These script will call the disangle-dock-MD2.py and 2D-image-reaccess.py scripts when required and perform all the analysis and save oputput

-----------------------------------------------------------------------
COMPLETE COMMAND-LINE OPTIONS
-----------------------------------------------------------------------
Additional options include the following for which default values are used if not selected by the user via command line.
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
