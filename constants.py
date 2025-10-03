import yaml, os


def load_config():
    # Load default template
    with open("config.yaml") as f:
        config = yaml.safe_load(f)

    # If user has a local config, override defaults
    if os.path.exists("config.local.yaml"):
        with open("config.local.yaml") as f:
            local_config = yaml.safe_load(f)
        config.update(local_config)

    return config

CONFIG = load_config()
BASE_PATH = CONFIG["base_path"]
PYTHON_ENV = CONFIG["python_env"]

STARTING_DATA_PATH = f"{BASE_PATH}/starting_data"
DEMOGRAPHICS_CSV = f"{STARTING_DATA_PATH}/Pangenomes Populations - VCF.csv"

DEFAULT_SUBJECTS_FILE = f"{STARTING_DATA_PATH}/subjects_files/subjects.txt"
DEFAULT_CAPACITY_FILE = f"{STARTING_DATA_PATH}/capacity_files/capacities_log.txt"

DEAD_WEIGHT = 10000000000
VERBOSE = True

PANGENOME_NPY = "pangenome.npy"
PANGENOME_POSITIONS_NPY = "pangenome_positions.npy"
PANGENOME_SUBJECTS_NPY = "pangenome_subjects.npy"

AFS_FREQ_NPY = "allele_frequencies.npy"

ONEK_NPY = "1000g_phased.npy"
ONEK_MASKED_NPY = "1000g_phased_masked.npy"
ONEK_SUBJECTS_NPY = "1000g_subjects.npy"
ONEK_POSITIONS_NPY = "1000g_positions.npy"

EXPERIMENT_PATH = f"{BASE_PATH}/experiments"
PLOT_OUT_PATH = f"{BASE_PATH}/plots"

AF_PATH = f"{BASE_PATH}/starting_data/chr21/population_af_freqs.pickle"

PANGENOME_HAPLOTYPES = f"{BASE_PATH}/starting_data/chr21/pangenome.npy"
PANGENOME_SUBJECTS = f"{BASE_PATH}/starting_data/chr21/pangenome_subjects.npy"
PANGENOME_SUPPORT = f"{BASE_PATH}/starting_data/chr21/pangenome_support.pickle"

POPULATION_AF_FREQS = f"{BASE_PATH}/starting_data/chr21/population_af_freqs.pickle"

NUM_1000G_SUBJECTS = 3202

ALL_POPULATIONS = ["East Asian Ancestry", "American Ancestry", "South Asian Ancestry", "African Ancestry"]

READ_SUBJECTS = ["HG00138", "HG00635", "HG01112", "HG01600", "HG02698"]
TYPES_OF_RECONSTRUCTION = ["new_haplotypes"]

READ_POPULATIONS_FILE = f"{BASE_PATH}/starting_data/read_populations.csv"
READ_SUBJECTS_FILE = f"{BASE_PATH}/starting_data/all_samples.txt"

STATS_FILE_SCHEMA = [
    ("total_alignments", 0),
    ("total_primary", 0),
    ("total_secondary", 0),
    ("total_aligned", 0),
    ("total_perfect", 0),
    ("total_gapless", 0),
    ("total_paired", 0),
    ("total_properly_paired", 0),
    ("alignment_score", 1),
    ("mapping_quality", 1),
    ("junk", 2),
    ("junk", 2),
    ("junk", 2),
    ("junk", 2),
    ("junk", 2),
    ("junk", 2),
]

PANGENIE_COLUMNS = ["quality","allele_frequency","unique_kmers","missing_alleles","total_baseline","total_baseline_biallelic","total_baseline_nonref","total_intersection","correct_all",
                    "wrong_all","not_typed_all","correct_biallelic","wrong_biallelic","not_typed_biallelic","correct_non-ref","wrong_non-ref","not_typed_non-ref",
                    "nr_correct_all","nr_wrong_all","nr_not_typed_all","nr_not_in_callset_all","nr_correct_biallelic","nr_wrong_biallelic","nr_not_typed_biallelic",
                    "nr_not_in_callset_biallelic","nr_correct_non-ref","nr_wrong_non-ref","nr_not_typed_non-ref","nr_not_in_callset_non-ref"]

INTERESTED_FILES_GAPSCORE = [
    "no_obfuscation",
    "left_inference",
    "center_inference",
    "right_inference"
]

TYPES_OF_SVS = ["complex", "insertion", "deletion", "snp"]
SIZE_OF_SVS = ["large", "midsize", "small"]

AGGREGATION_DICTIONARY = {
    "subject": "index",
    "capacity": "index",
    "region": "index",
    "population": "index",
    "population_code": "index",

    "max_utility_loss": "sum",
    "max_pmi_gain": "sum",
    "utility_loss": "sum",
    "pmi_gain": "sum",
    "utility_loss_normalized": "mean",
    "pmi_gain_normalized": "mean",
    "number_of_moves": "sum",
    "percent_shared_snps": "mean",

    "maf_wd": "mean",
    "maf_kl": "mean",
    "ld_euclidean": "mean",

    "af_loss": "mean",
    "af_loss_snp_only": "mean",
    "af_complex": "mean",
    "af_loss_snps_0_01": "mean",
    "af_loss_snps_0_05": "mean",
    "af_loss_snps_0_1": "mean",
    "af_loss_snps_0_5": "mean",

    "wd_af": "mean",
    "wd_af_snp_only": "mean",
    "wd_af_complex": "mean",
    "wd_af_snps_0_01": "mean",
    "wd_af_snps_0_05": "mean",
    "wd_af_snps_0_1": "mean",
    "wd_af_snps_0_5": "mean",

    "changed_alleles_percentage": "mean",
    "mean_difference": "mean",
    "number_of_changed_alleles": "mean",
    "number_of_alleles_touched": "mean",
    "number_of_alleles_touched_percentage": "mean",

    "ld": "mean",
}

TOTAL_UTILITY_JSON = f"{BASE_PATH}/starting_data/utility_loss.json"

SUBPOPS = ["All", "East Asian Ancestry", "American Ancestry", "African Ancestry", "South Asian Ancestry"]

SUBPOPS_TEST = ["All", "East Asian Ancestry", "American Ancestry"]

MAX_LD_DISTANCE_KBP=1000

READS_DIR=f"{BASE_PATH}/reads_fastq"