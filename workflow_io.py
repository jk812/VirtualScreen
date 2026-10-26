import configparser
import os


class WorkflowInput:
    """Class to hold the workflow input parameters."""
    def __init__(self, config: configparser.ConfigParser):
        # workflow settings
        workflow_settings = config["workflow"]
        self.workflow_rounds = workflow_settings.getint("workflow_rounds") # run the workflow workflow_rounds times
        self.gen_candidates_to_save = workflow_settings.getint("gen_candidates_to_save") # maximum number of candidates to save for the next generation-docking round
        self.cluster_distance_threshold = workflow_settings.getfloat("cluster_distance_threshold") # distance threshold for clustering generated candidates
        self.docking_score_candidate_threshold = workflow_settings.getfloat("docking_score_candidate_threshold") # threshold buffer for the docking score of a candidate to be considered for the next generation
        # input and output directories for the workflow
        self.workflow_results_dir = workflow_settings["working_directory"] # working directory for the workflow results, which will contain subdirectories for each round of generation and docking.
        self.input_protein_structure_path = workflow_settings["input_protein_structure_path"] # protein structure to dock to

        ## preparation settings (protein+ligand)
        prep_settings = config["ligand.prep"]
        self.input_smiles = prep_settings["input_smiles"] # PDB ligand ID: 9JU, bound to PDB ID 5ZTY
        self.pH = prep_settings.getfloat("pH") # physiological pH - specified for ligand and protein protonation
        self.pH_deviation = prep_settings.getfloat("pH_deviation") # deviation of +/- 2 log units for ligand protonation
        self.enumerate_only_unassigned = prep_settings.getboolean("enumerate_only_unassigned") # only enumerate undefined stereocentres.
        self.max_opt_iter = prep_settings.getint("max_opt_iter") # when generating 3D molecules, this is the maximum number of optimisation steps for the 3D geometry optimisation

        ## docking settings
        docking_settings = config["docking"]
        self.exhaustiveness = docking_settings.getint("exhaustiveness") # increased from the default of 8
        self.n_poses = docking_settings.getint("n_poses") # only save the best pose per ligand.
        self.grid_buffer = docking_settings.getfloat("grid_buffer") # buffer in Angstroms to add to the receptor grid dimensions around the co-crystal input structure cognate ligand.

        ## generation settings
        generation_settings = config["generation"]
        self.gen_algo = generation_settings["gen_algo"] # generation method for ANNalog. Options are: "beam" (default), "BF-beam", and "sampling". See ANNalog documentation for more details.
        self.num_of_mols_to_gen = generation_settings.getint("num_of_mols_to_gen")
        self.internal_gen_multiplier = generation_settings.getfloat("internal_gen_multiplier") # internally, generate 2 times the number of molecules to account for duplicates, which will be removed later

def generate_default_workflow_config():
    """Generates a default configuration file named 'workflow_config.ini' with predefined settings for the workflow,
    ligand preparation, docking, and generation parameters.
    """
    config = configparser.ConfigParser()
    config["DEFAULT"] = {"pH": 7.4}
    config['workflow'] = {"input_protein_structure_path": os.path.join("Structures", "refined_cnr2_human_5ZTY_inactive.pdb"),
                        "working_directory": os.path.join("results", "workflow"),
                        'workflow_rounds': '5',
                        'gen_candidates_to_save': '5',
                        'cluster_distance_threshold': '0.45',
                        'docking_score_candidate_threshold': '2'}
    config['ligand.prep'] = {"input_smiles": "smiles.smi",
                            "pH_deviation": "2",
                            "enumerate_only_unassigned": "True",
                            "max_opt_iter": "1000"}
    config['docking'] = {"exhaustiveness": "24",
                        "n_poses": "1",
                        "grid_buffer": "5"}

    with open('workflow_config.ini', 'w') as configfile:
        config.write(configfile)

    return config

def read_config_file(config_file_path:str='workflow_config.ini') -> WorkflowInput:
    """Reads the configuration file and returns a WorkflowInput object.

    Args:
        config_file_path (str): Path to the configuration file.

    Returns:
        WorkflowInput: WorkflowInput object containing the typed configuration settings.
    """
    
    config = configparser.ConfigParser()
    config.read(config_file_path)
    workflow_input = WorkflowInput(config)  # Create a WorkflowInput object to assign types and store settings
    workflow_input.config_file_path = config_file_path  # Store the path to the configuration file in the WorkflowInput object
    return workflow_input