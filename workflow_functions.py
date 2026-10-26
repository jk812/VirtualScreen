import os
import shutil

from unit_functions import *
from workflow_io import read_config_file


## PROTEIN PREPARATION
def prepare_protein(input_protein_path:str, working_directory:str, pH=7.4):
    """Convenience function which fixes a protein structure with PDBFixer if any issues are present, then protonates the protein with
    pdb2pqr (propka under the hood) at a given pH, and converts the protonated structure file (.pqr used in this case) to a .pdbqt file
    ready for docking with Autodock Vina.

    Args:
        input_protein_path (str): Input protein PDB filepath.
        working_directory (str): The directory where the prepared protein files will be saved.
        pH (float, optional): The pH at which to assign protonated states. Defaults to 7.4 (physiological pH).

    Returns:
    str: The path where the prepared protein .pdbqt file has been saved.
    """
    # create protein structure output directory within the working directory
    output_directory_path = os.path.join(working_directory, "protein")
    os.makedirs(output_directory_path)

    # run protein preparation
    fixed_protein_path = fix_protein(input_protein_path, output_directory_path)
    protonated_protein_pqr_path = protonate_protein(fixed_protein_path, output_directory_path, pH=pH)
    pdbqt_path = pdb_to_pdbqt(protonated_protein_pqr_path, output_directory_path)

    return pdbqt_path

### LIGAND PREPARATION
def prepare_ligands(smiles_list:list[str], generation_round_dir:str, pH=7.4, pH_deviation=2, enumerate_only_unassigned=True, max_iter=1000) -> list[str]:
    """Convenience function which generates the protonation states and tautomers of a list of SMILES, then enumerates
    stereoisomers, generating 3D structures for each stereoisomer and optimising the structures, and saving .pdbqt files
    for docking with Autodock Vina.

    Args:
        smiles_list (list[str]): List of SMILES strings.
        generation_round_dir (str): Directory for the current generation round.
        pH (float, optional): pH to protonate ligands. Defaults to 7.4.
        pH_deviation (int, optional): Deviation from the specified pH to consider for ligand protonation. Defaults to 2.
        enumerate_only_unassigned (bool, optional): If set (default), stereocenters which have specified stereochemistry will not be perturbed unless they are part of a relative stereo group. Defaults to True.
        max_iter (int, optional): Maximum number of iterations to optimise a 3D structure. Defaults to 1000.
    
    Returns:
        tuple[list[str], str]: List of absolute filepaths for each generated ligand .pdbqt file and the path to the enumerated ligands CSV file.
    """
    ligand_directory = os.path.join(generation_round_dir, "ligand_prep")
    os.makedirs(ligand_directory)

    # For below: progress_apply is from tqdm.pandas() and applies a function whilst showing a progress bar
    ligand_df = pd.DataFrame({"annalog": smiles_list})

    print("Protonating molecules...")
    ligand_df["protonated"] = ligand_df["annalog"].progress_apply(protonate_ligand, pH=pH, pH_deviation=pH_deviation)
    ligand_df = ligand_df.explode("protonated", ignore_index=True)

    print("Enumerating stereoisomers...")
    ligand_df["isomers"] = ligand_df["protonated"].progress_apply(enumerate_stereoisomers, only_unassigned=enumerate_only_unassigned)
    ligand_df = ligand_df.explode("isomers", ignore_index=True)
    invalid_smiles = ligand_df.isnull().sum().sum() # when generating stereoisomers, NaN can be returned if a SMILES cannot be converted to 3D, thus it is considered nonsensical.
    print(f"{invalid_smiles} invalid SMILES detected.")
    if(invalid_smiles > 0):
        ligand_df = ligand_df.dropna()
        print("Invalid SMILES removed.")

    print("Deduplicating isomers...")
    len_before_deduplication = ligand_df.shape[0]
    ligand_df = ligand_df.drop_duplicates("isomers", ignore_index=True)
    len_after_deduplication = ligand_df.shape[0]
    print(f"{len_before_deduplication} molecules identified after enumeration, {len_after_deduplication} after deduplication.")
    print(f"{len_before_deduplication-len_after_deduplication} duplicates removed.")

    print("Generating 3D structures per SMILES...") # this is done more to check the structures later on, if the user wishes
    mols_3d = ligand_df["isomers"].progress_apply(generate_3d_molecule, max_iter=max_iter)
    if(mols_3d.isnull().any()): # it is possible for 3D generation to fail despite the check during stereo enumeration.
        mols_3d, ligand_df = clean_3d_molecule_generation(mols_3d, ligand_df, ligand_directory)
    mols_3d = mols_3d.to_list() # convert to list for later use in writing .pdbqt and .sdf files

    sdf_path = os.path.join(ligand_directory, "prepped_ligands.sdf")
    print(f"Saving 3D structures to {sdf_path}")
    write_3d_mols_to_sdf(mols_3d, ligand_df["isomers"].to_list(), sdf_path)

    ligand_info_path = os.path.join(ligand_directory, "prepped_ligands.csv")
    print(f"Saving enumerated ligand smiles to {ligand_info_path}")
    ligand_df.insert(0, "sdf_number", ligand_df.index+1) # Add a column for the sdf number, which is the index of the molecule in the .sdf file (1-indexed)
    ligand_df.to_csv(ligand_info_path, index=False)
    
    print("Writing .pdbqt files of each ligand...")
    ligand_pdbqt_filepath_list = write_multiple_ligand_pdbqt_files(mols_3d, ligand_directory)
    print(f"{len(ligand_pdbqt_filepath_list)} molecules have been prepped for docking.")
    
    return ligand_pdbqt_filepath_list, ligand_info_path


### DOCKING FUNCTIONS
def vina_setup_receptor(vina_object:Vina, protein_with_ligand_struct_path:str, protein_pdbqt_path:str, grid_buffer:float=5) -> Vina:
    """Initialise and setup a Vina object with a receptor structure and receptor grid.

    Args:
        vina_object (Vina): An instance of the Vina class from the vina package.
        protein_with_ligand_struct_path (str): Filepath to a protein-ligand co-crystal structure (PDB format) to setup receptor grid.
        protein_pdbqt_path (str): Filepath to the prepared protein .pdbqt file.
        grid_buffer (float, optional): Buffer size in Angstroms to add to the receptor grid dimensions around cognate ligand. Defaults to 5.

    Returns:
        Vina: The initialised Vina object with the receptor structure and receptor grid set up.
    """
    box_centre, box_dimensions = get_receptor_grid(protein_with_ligand_struct_path, buffer=grid_buffer)
    vina_object.set_receptor(protein_pdbqt_path)
    vina_object.compute_vina_maps(center=box_centre, box_size=box_dimensions)
    return vina_object

def docking(ligand_pdbqt_filepath_list:list[str], prepared_vina_object:Vina, n_poses:int, exhaustiveness:int, smiles_df_path:str, generation_round_dir:str):
    """Convenience function to perform docking and save docking results.

    Args:
        ligand_pdbqt_filepath_list (list[str]): List of prepared ligand PDBQT filepaths.
        prepared_vina_object (Vina): A set up Vina object that already has a receptor and receptor grid set up, as done with vina_setup_receptor.
        n_poses (int): Vina parameter - Number of poses to generate for docking. This controls how many poses are returned after docking.
        exhaustiveness (int): Vina parameter - Exhaustiveness of the global search. Higher values lead to more thorough searches but take longer.
        smiles_df_path (str): Filepath to file containing enumerated isomer SMILES made earlier in the workflow.
        generation_round_dir (str): Directory for the current generation round.
    
    Returns:
        str: Filepath to the saved .csv file containing SMILES and docking energies.
    """
    docking_results_directory = os.path.join(generation_round_dir, "docking")
    os.makedirs(docking_results_directory)

    docking_energy_list, docked_pose_pdbqt_list = dock_and_score(ligand_pdbqt_filepath_list, prepared_vina_object, n_poses, exhaustiveness)
    docking_results_output_path = save_docking_results(smiles_df_path, docking_results_directory, docking_energy_list, docked_pose_pdbqt_list)

    return docking_results_output_path

### CLUSTERING AND CANDIDATE SELECTION FUNCTIONS
def cluster_and_select_candidates(results_filepath:str, generation_round_dir:str, clustering_threshold=0.4, selection_method="docking_score", threshold_buffer=2, candidates_to_save=5):
    """Cluster the docking results and select candidates based on the specified selection method. Clustered results are saved to a CSV file,
    along with the selected candidates for the next generation-docking round. The selection method can be based on docking scores or random
    selection.

    Args:
        results_filepath (str): Path to the CSV file containing docking results.
        generation_round_dir (str): Directory for the current generation round.
        clustering_threshold (float, optional): Tanimoto distance threshold for Butina clustering. Defaults to 0.4.
        selection_method (str, optional): Method for selecting candidates. Options are "docking_score" or "random". Defaults to "docking_score".
        threshold_buffer (int, optional): Buffer for docking score thresholding when using "docking_score" selection method (final threshold is set to the best docking score plus this buffer). Defaults to 2.
        candidates_to_save (int, optional): Number of candidates to save total, where each candidate comes from a different cluster. Defaults to 5.

    Returns:
        str: Path to the CSV file containing the selected candidates.
    """
    clustering_results_directory = os.path.join(generation_round_dir, "clustering")
    os.makedirs(clustering_results_directory)

    results_df_with_clusters = cluster(results_filepath, clustering_results_directory, clustering_threshold)

    if selection_method == "docking_score":
        selected_candidates = select_candidates_docking_score(results_df_with_clusters, threshold_buffer, candidates_to_save)
    elif selection_method == "random":
        selected_candidates = select_candidates_random(results_df_with_clusters, candidates_to_save)
    else:
        raise ValueError("Invalid selection method. Choose 'docking_score' or 'random'.")

    selected_candidates_filepath = os.path.join(clustering_results_directory, "generation_candidates.csv")
    selected_candidates.to_csv(selected_candidates_filepath, index=False)
    return selected_candidates_filepath

def get_generation_candidates_smi(generation_candidates_filepath:str, generation_round_dir:str) -> str:
    """Convert the previous generation-docking round's generation candidates file into a .smi file for the current
    generation-docking round. The filepath of the generation candidates .smi file is returned to be passed into
    ANNalog.

    Args:
        generation_candidates_filepath (str): Path to the CSV file containing generation candidates acquired at the end of the previous round.
        generation_round_dir (str): Directory for the current generation round.
    Returns:
        str: Path to the generated .smi file containing generation candidates.
    """

    gen_candidates_list = pd.read_csv(generation_candidates_filepath)["smiles"].to_numpy()
    gen_candidates_smi_filepath = os.path.join(generation_round_dir, "generation_candidates.smi")
    np.savetxt(gen_candidates_smi_filepath, gen_candidates_list, fmt="%s") # save the generation candidates to a text file for use in the next round of generation and docking
    return gen_candidates_smi_filepath

### RUN FULL WORKFLOW
def run_workflow(config_file_path:str='workflow_config.ini', selection_method:str="docking_score"):
    """Runs the entire generation-docking workflow based on the settings in the input config file.

    Args:
        config_file_path (str): Path to the configuration file. Defaults to 'workflow_config.ini'.
        selection_method (str, optional): Method for selecting candidates. Options are "docking_score" or "random". Defaults to "docking_score".
    """
    config = read_config_file(config_file_path) # reading in the workflow configuration file and storing it in a WorkflowInput object
    input_smiles = config.input_smiles # initial input SMILES (SMILES string or .smi filepath) for the first round of generation and docking

    # protein preparation
    protein_pdbqt_path = prepare_protein(config.input_protein_structure_path, config.workflow_results_dir, pH=config.pH)
    v=Vina(sf_name="vina", verbosity=0) # using the default vina scoring function, and setting the verbosity to 0 to silence Vina output
    v=vina_setup_receptor(v, config.input_protein_structure_path, protein_pdbqt_path, grid_buffer=config.grid_buffer)

    # set up directory filepath for the workflow
    run_directory = os.path.join(config.workflow_results_dir, config.gen_algo) # naming the directory to save results to after the generation algorithm used
    os.makedirs(run_directory)
    shutil.copyfile(config.config_file_path, os.path.join(run_directory, os.path.basename(config.config_file_path))) # move the workflow configuration file to the workflow run directory for record keeping

    # ligand generation and docking loop
    selected_candidates_filepath = None # initialise the selected candidates filepath to None, which will be updated after each round of generation and docking.
    for round_num in range(1, config.workflow_rounds+1): # 1-indexing to make printing nicer
        print(f"Starting round {round_num} of {config.workflow_rounds}...")
        print("-"*60)

        generation_round_dir = os.path.join(run_directory, f"round_{round_num}")
        os.makedirs(generation_round_dir)

        if(round_num > 1): # for subsequent rounds, the input SMILES is the generation candidates from the previous round
            input_smiles = get_generation_candidates_smi(selected_candidates_filepath, generation_round_dir)

        print("Reading in molecules...")
        smiles_to_dock = generate_and_deduplicate(config.num_of_mols_to_gen, input_smiles, generation_round_dir, method=config.gen_algo, internal_multiplier=config.internal_gen_multiplier)

        print("Preparing ligands for docking...")
        ligand_pdbqt_filepaths, enumerated_smiles_path = prepare_ligands(smiles_to_dock,
                                                                        generation_round_dir,
                                                                        pH=config.pH,
                                                                        pH_deviation=config.pH_deviation,
                                                                        enumerate_only_unassigned=config.enumerate_only_unassigned,
                                                                        max_iter=config.max_opt_iter)

        print("Docking...")
        docking_results_filepath = docking(ligand_pdbqt_filepaths, v, config.n_poses, config.exhaustiveness, enumerated_smiles_path, generation_round_dir)

        selected_candidates_filepath = cluster_and_select_candidates(docking_results_filepath, generation_round_dir, selection_method=selection_method, candidates_to_save=config.gen_candidates_to_save, clustering_threshold=config.cluster_distance_threshold, threshold_buffer=config.docking_score_candidate_threshold)
    print("-"*60+"\n"+"-"*60)
    print("Workflow completed.")