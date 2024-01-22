import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_FOLDER = os.path.join(BASE_DIR, "input")

OUTPUT_FOLDER = os.path.join(BASE_DIR, "output")

LOGS_FOLDER = os.path.join(BASE_DIR, "logs")

STDERR_FOLDER = os.path.join(LOGS_FOLDER, "e")

STDOUT_FOLDER = os.path.join(LOGS_FOLDER, "o")

DATA_FOLDER = os.path.join(BASE_DIR, "data")

MOLS_FOLDER = os.path.join(DATA_FOLDER, "bound_mols")

INTERS_FOLDER = os.path.join(DATA_FOLDER, "bound_mol_inters")

EXP_FOLDER = os.path.join(DATA_FOLDER, "exp")

MATS_FOLDER = os.path.join(DATA_FOLDER, "supp_mats")

SEGMENT_FOLDER = os.path.join(DATA_FOLDER, "segments")

STRUCTURE_FOLDER = os.path.join(DATA_FOLDER, "structures")

ASYM_FOLDER = os.path.join(STRUCTURE_FOLDER, "asymmetric")

ASSEMBLY_FOLDER = os.path.join(STRUCTURE_FOLDER, "assembly")