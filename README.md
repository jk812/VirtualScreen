# VirtualScreen
Short Virtual Screening Workflow, for Software Dev Course TTPS

# Installation
Using the [pixi package manager](https://github.com/prefix-dev/pixi), simply type the following to install the project environment to run the workflow:
```
pixi install
```
Ensure that this command is run in the same directory as the `pixi.lock` and `pixi.toml` files.

To install the separate conda environment needed for plotting, use the following command after creating a new conda environment (tested with Python 3.14.7):
```
conda install -c conda-forge rdkit umap-learn matplotlib seaborn pandas numpy scipy mdanalysis && pip install prolif starbars statannotations
```

# Running the workflow
To run the workflow after installing the pixi environment, use the following command:
```
pixi run python main.py
```

# Config input input
This project has a configuration file input. The default configuration file, as bundled with this repository, is `workflow_config.ini`. Other names for this config file are possible, but named to be passed into the workflow script call with the `-c` flag, e.g. for a config file named `special_config.ini`:
```
pixi run python main.py -c special_config.ini
```
## Config file structure
### \[DEFAULT]
`ph` - pH specified for ligand and protein protonation
### \[workflow]
`input_protein_structure_path` - filepath to .pdb co-crystal (protein+ligand) structure. Ideally this structure should have already been cleaned to contain only the protein and ligand.
`working_directory` = filepath to directory where results will be saved (this directory will be created for the user).
`workflow_rounds` - run the workflow workflow_rounds times.
`gen_candidates_to_save` - maximum number of candidates to save for the next generation-docking round.
`cluster_distance_threshold`  - Tanimoto (ECFP4, length 1024) distance threshold for clustering generated candidates.
`docking_score_candidate_threshold` - the docking score threshold buffer used to selecting generation candidates for the next generation-docking round, i.e. only compounds that are within `docking_score_candidate_threshold` of the best docking score obtained in this round can be selected as generation candidates. This is to avoid including relatively poorly scoring molecules as generation candidates.

### \[ligand.prep]
`input_smiles` - either a SMILES or a .smi file containing as many SMILES as the user desires.
`ph_deviation` - pH deviation of +/- `ph_deviation` log units for ligand protonation state enumeration.
`enumerate_only_unassigned` - only enumerate undefined stereocentres of True, if False enumerate all stereocentres.
`max_opt_iter` - when generating 3D molecules, this is the maximum number of optimisation steps for the 3D geometry optimisation.

### \[docking]
`exhaustiveness` - number of Monte Carlo samplings performed by Autodock Vina.
`n_poses` - only save the best `n_poses`poses.
`grid_buffer` - buffer in Angstroms to add to the receptor grid dimensions around the co-crystal input structure cognate ligand.



# Repository structure
- The files required to run the workflow are `main.py`, `unit_functions.py` and `workflow_functions.py`. `main.py` is the script used to run the workflow. `workflow_functions.py` contain the workflow-level functions used by `main.py` and `unit_functions.py` contains the functions used by `workflow_functions.py`.
- `workflow_config.ini` is a configuration file read in by `main.py`. See above for its file structure.
- `smiles.smi` is a SMILES file that currently contains one SMILES string. This is read into `main.py` by the instruction of `workflow_config.ini` with the `input_smiles` argument.
- `workflow_io.py` contains the functions for reading and writing the `workflow_config.ini` file into `main.py`.
