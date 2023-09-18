import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_FOLDER = os.path.join(BASE_DIR, "input")

OUTPUT_FOLDER = os.path.join(BASE_DIR, "output")

DATA_FOLDER = os.path.join(BASE_DIR, "data")

MOLS_FOLDER = os.path.join(DATA_FOLDER, "bound_mols")

INTERS_FOLDER = os.path.join(DATA_FOLDER, "bound_mol_inters")

EXP_FOLDER = os.path.join(DATA_FOLDER, "exp")

MATS_FOLDER = os.path.join(DATA_FOLDER, "supp_mats")

SEGMENT_FOLDER = os.path.join(DATA_FOLDER, "segment_data")