
### IMPORTS

import os
import sys
import copy
import gzip
import math
import time
import scipy
import pickle
import shutil
import logging
import argparse
import Bio.SeqIO
import numpy as np
import pandas as pd
import configparser
from Bio import PDB
from Bio import AlignIO
import scipy.stats as stats

from prointvar.pdbx import PDBXreader, PDBXwriter
from prointvar.sifts import SIFTSreader
from prointvar.dssp import DSSPrunner, DSSPreader

import varalign.align_variants
import varalign.alignments

from urllib.error import HTTPError
from urllib.error import URLError

from config import BASE_DIR, INPUT_FOLDER, OUTPUT_FOLDER, DATA_FOLDER, MOLS_FOLDER, INTERS_FOLDER, EXP_FOLDER, MATS_FOLDER, SEGMENT_FOLDER

### SETTING UP LOGGER

logging.basicConfig(filename = "ligysis.log", format = '%(asctime)s %(name)s [%(levelname)-8s] - %(message)s', level = logging.INFO)

log = logging.getLogger("LIGYSIS")

## READ CONFIG FILE

config = configparser.ConfigParser()
config_path = os.path.join(BASE_DIR, "ligysis_config.txt")
config.read(config_path) # assuming this program is being executed one level above where this script and its config file are located

clean_pdb_bin = config["binaries"].get("clean_pdb_bin")                     # location of clean_pdb.py.
clean_pdb_python_bin = config["binaries"].get("clean_pdb_python_bin")       # location of python binary to run clean_pdb.py.
dssp_bin = config["binaries"].get("dssp_bin")                               # location of DSSP binary.
biolip_data = config["dbs"].get("biolip_data")                              # location of dictionary containing information about ligands in Biolip.
ensembl_sqlite_path = config["dbs"].get("ensembl_sqlite")                   # location of a local copy of ENSEMBL mappings from UniProt Accession to genome (sqlite)
gnomad_vcf = config["dbs"].get("gnomad_vcf")                                # location of gnomAD VCF. This database is not updated.
pdb_db_path = config["dbs"].get("pdb_db_path")                              # location of local copy of PDB. This should be updated, but it might be missing some files.
sifts_db_path = config["dbs"].get("sifts_db_path")                          # location of a local copy of SIFTS. This database might not be updated, current version is Feb 2023.
sifts_db_path2 = config["dbs"].get("sifts_db_path2")                        # location of a local backup copy of SIFTS. This database might not be updated, current version is Jul 2023.
swissprot = config["dbs"].get("swissprot")                                  # location of MY local SwissProt copy. This database is not updated, current version is Nov 2021.
MSA_fmt = config["formats"].get("MSA_fmt")                                  # MSA format used for all calculations. Currently only working with Stockholm format.
struc_fmt = config["formats"].get("struc_fmt")                              # structure format used for all calculations. Currently, only working format is PDB.
experimental_methods = config["other"].get("experimental_methods")          # experimental method used to determine structures. Default is X-ray.
experimental_methods_list = experimental_methods.split(",")                 # splits string into list of strings.
lig_clust_method = config["other"].get("lig_clust_method")                  # linkage method used for the hierarchical clustering.
lig_clust_metric = config["other"].get("lig_clust_metric")                  # similarity metric used defined for ligand binding site definition (later transformed into distance metric).
cons_t_h = float(config["thresholds"].get("cons_t_h"))                      # conservation score upper threshold to consider position highly divergent. Currently only working with Shenkin divergence score.
cons_t_l = float(config["thresholds"].get("cons_t_l"))                      # conservation score lower threshold to consider position highly conserved. Currenyly only working with Shenkin divergence score.
cons_ts = [cons_t_l, cons_t_h]                                              # groups both thresholds into a list.
jackhmmer_n_it = int(config["thresholds"].get("jackhmmer_n_it"))            # number of iterations to perform on remote-homologue search by jackHMMR
lig_clust_dist = float(config["thresholds"].get("lig_clust_dist"))          # distance threshold where hierarchical tree is cut for ligand binding site definition.
MES_t = float(config["thresholds"].get("MES_t"))                            # Missense Enrichment Score threshold to consider a position missense-depleted, or enriched.
max_retry = int(config["thresholds"].get("max_retry"))                      # number of maximum attempts to make to retrieve a certain piece of data from PDBe API.
resolution = float(config["thresholds"].get("resolution"))                  # resolution threshold to consider a structure high-resolution.
sleep_time = float(config["thresholds"].get("sleep_time"))                  # time to sleep between queries to the PDBe API.

### LISTS

bbone_atoms = ["C", "CA", "N", "O"]

pdb_resnames = [
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"
]

aas_1l= [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
    "-"
]

chimera_cmd_args = [
    "ksdssp", "background solid white", "~dis",
    "sel ~@/color=white", "dis sel", "namesel lois",
    "~sel"
]

exp_data_cols = [
    "resolution", "resolution_low", "resolution_high",
    "r_factor", "r_free", "r_work", "experimental_method", # columns of interest from experimental data table
    "structure_determination_method"
]

dssp_cols = [
    "PDB_ResNum", "SS", "ACC", "KAPPA",     # selects subset of columns
    "ALPHA", "PHI", "PSI", "RSA"
]

consvar_class_colours = [
    "royalblue", "green", "grey", "firebrick", "orange"
]

### FUNCTIONS

## NEW SECTION

def get_status_code_data(status_code_file):
    """
    Reads status code file and returns dictionaries
    for accession and segment status codes.
    """
    accs_status_dict = {}
    segs_status_dict = {}
    with open(status_code_file, "r") as f:
        for line in f:
            d = line.strip().split("\t")
            if "_" in d[0]:
                segs_status_dict[d[0]] = d[1]
            else:
                accs_status_dict[d[0]] = d[1]
    return accs_status_dict, segs_status_dict

def cp_sqlite(wd, og_path = ensembl_sqlite_path):
    """
    Copies ensembl_cache.sqlite to execution directory.
    """
    hidden_var_dir = os.path.join(wd, ".varalign")
    sqlite_name = os.path.basename(og_path)
    if not os.path.isdir(hidden_var_dir):
        os.mkdir(hidden_var_dir)
    else:
        pass
    cp_path = os.path.join(hidden_var_dir, sqlite_name)
    shutil.copy(og_path, cp_path)
    return cp_path

def rm_sqlite(cp_path):
    """
    Removes ensembl_cache.sqlite from execution directory.
    """
    hidden_var_dir = os.path.dirname(cp_path)
    os.remove(cp_path)
    os.rmdir(hidden_var_dir)

def dump_pickle(data, f_out):
    """
    Dumps data to pickle.
    """
    with open(f_out, "wb") as f:
        pickle.dump(data, f)

def load_pickle(f_in):
    """
    Loads data from pickle.
    """
    with open(f_in, "rb") as f:
        data = pickle.load(f)
    return data

## CLEANING PDB FILES

def run_clean_pdb(pdb_path):
    """
    Runs pdb_clean.py.
    """
    args = [
        clean_pdb_python_bin, clean_pdb_bin, pdb_path
    ]
    exit_code = os.system(" ".join(args))
    if exit_code != 0:
        log.critical("{} was not cleaned with command: {}".format(pdb_path, " ".join(args)))

## TRANSFORMING PDB FILES

def fmt_mat_in(mat_in):
    """
    Formats a transformation matrix in a way that can be
    used by the transformation function.
    """
    mat_out = [
        [mat_in[0][0], mat_in[1][0], mat_in[2][0]],
        [mat_in[0][1], mat_in[1][1], mat_in[2][1]],
        [mat_in[0][2], mat_in[1][2], mat_in[2][2]],
        [mat_in[0][-1], mat_in[1][-1], mat_in[2][-1]]
    ]
    matrix_dict = {
        "rotation": [mat_out[0], mat_out[1], mat_out[2]],
        "translation": mat_out[3]
    }
    return matrix_dict

def get_segments_dict(supp_data, acc):
    """
    Given a superimposition table for a protein, creates
    a dictionary with information about the different
    structure coverage segments.
    """
    segment_dict = {}
    for idx, row in supp_data.iterrows():
        segment_dict[idx] = {}
        clusters = row[acc]["clusters"]
        for cluster in clusters:
            for member in cluster:
                pdb_id = member["pdb_id"]
                chain_id = member["auth_asym_id"]
                unit_id = "{}_{}".format(pdb_id, chain_id)
                segment_dict[idx][unit_id] = member
    return segment_dict

def get_segment_membership(supp_data, acc):
    """
    Given a superimposition table for a protein, creates
    a dictionary that indicates which pdb chains are included
    within each structural coverage segment.
    """
    segment_membership = {}
    for idx, row in supp_data.iterrows():
        segment_membership[idx] = []
        clusters = row[acc]["clusters"]
        for cluster in clusters:
            for member in cluster:
                pdb_id = member["pdb_id"]
                chain_id = member["auth_asym_id"]
                unit_id = "{}_{}".format(pdb_id, chain_id)
                segment_membership[idx].append(unit_id)
    return segment_membership

def parse_pdb_file(pdb_path, fmt):
    """
    :param pdb_path:
    :return: biopython's structure object
    """
    if fmt == "pdb":
        parser = PDB.PDBParser()
    elif fmt == "cif":
        parser = PDB.MMCIFParser()
        
    strct_name = os.path.basename(pdb_path).split(".")[0][3:]
    structure = parser.get_structure(strct_name, pdb_path)
    return structure

def apply_transformation(structure, matrix, output_path, chain_id, fmt):
    """
    Transforms structure based on the transformation matrix
    :param structure: biopython's structure object
    :param matrix: transformation matrix dict
    :return: transformed structure
    """
    rotation = matrix["rotation"]
    translation = matrix["translation"]
    
    if fmt == "pdb":
        io = PDB.PDBIO()
    elif fmt == "cif":
        io = PDB.MMCIFIO()
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    for atom in residue:
                        atom.transform(rotation, translation)
    
                io.set_structure(chain)
                io.save(output_path)

def pdb_transform(pdb_path, output_path, matrix_raw, chain_id, fmt = struc_fmt):
    """
    Applies transformation matrix to the PDB file and writes the new PDB to file
    :param pdb_path: path to the input PDB file
    :param matrix_path: path to the transformation matrix file
    :param output_path: path to the output PDB file
    """
    matrix_rf = fmt_mat_in(matrix_raw)
    
    structure = parse_pdb_file(pdb_path, fmt)

    structure = apply_transformation(structure, matrix_rf, output_path, chain_id, fmt)

def transform_all_files(pdb_files, matrices, chains, raw_dir, clean_dir, trans_dir):
    """
    Given a set of pdb files, matrices, and chains, uncompresses, cleans and transforms
    the coordinates according to a transformation matrix
    """
    for i, pdb_in in enumerate(pdb_files):
        pdb_root, _ = os.path.splitext(os.path.basename(pdb_in))
        pdb_out = os.path.join(raw_dir, pdb_root[3:] + ".pdb")
        #pdb_out = os.path.join(raw_dir, os.path.basename(pdb_in)[3:].replace(".ent.gz", ".pdb"))
        pdb_id = os.path.basename(pdb_in)[3:7]
        if os.path.isfile(pdb_out):
            pass
        else:
            with gzip.open(pdb_in, "rb") as f_in:
                with open(pdb_out, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        file_clean_from = os.path.join(raw_dir, pdb_root[3:] + ".clean.pdb")
        file_clean_to = os.path.join(clean_dir, pdb_root[3:] + ".clean.pdb")
        #file_clean_from = os.path.join(raw_dir, os.path.basename(pdb_out).replace(".pdb", ".clean.pdb"))
        #file_clean_to = os.path.join(clean_dir, os.path.basename(pdb_out).replace(".pdb", ".clean.pdb"))
        if os.path.isfile(file_clean_to):
            pass
        else:
            run_clean_pdb(pdb_out)
            os.remove(os.path.join(raw_dir, pdb_root[3:] + ".pdb" + ".breaks"))
            os.remove(os.path.join(raw_dir, pdb_root[3:] + ".pdb" + ".break_residues"))
            shutil.move(file_clean_from, file_clean_to)
            log.info("{} cleaned".format(pdb_id))
        transformed_out = os.path.join(trans_dir, pdb_root[3:] + "_{}_trans.pdb".format(chains[i]))
        #transformed_out = os.path.join(trans_dir, os.path.basename(pdb_out).replace(".pdb", "_{}_trans.pdb".format(chains[i])))
        if os.path.isfile(transformed_out):
            pass
        else:
            pdb_transform(file_clean_to, transformed_out, matrices[i], chains[i])
            log.info("{}_{} transformed".format(pdb_id, chains[i]))

## EXPERIMENTAL DATA AND VALIDATION

def get_experimental_data(pdb_ids, exp_data_dir, out):
    """
    Retrieves experimental data from PDBe rest API
    and saves subset of columns to pandas dataframe, and pickle.
    """
    exp_data_dfs = []
    for pdb_id in pdb_ids:
        pdb_expd_out = os.path.join(exp_data_dir, "{}_exp_data.json".format(pdb_id))
        if os.path.isfile(pdb_expd_out):
            exp_data_df = pd.read_json(pdb_expd_out, convert_axes = False, dtype = False)
            exp_data_dfs.append(exp_data_df)
        else:
            try:
                exp_data = pd.read_json("https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/{}".format(pdb_id), convert_axes = False, dtype = False)
                exp_data_dict = exp_data.loc[0, pdb_id]
                exp_data_df = pd.DataFrame({k: v for k, v in exp_data_dict.items() if k in exp_data_cols}, index = [0])
                exp_data_df["pdb_id"] = pdb_id
                exp_data_df.to_json(pdb_expd_out)
                exp_data_dfs.append(exp_data_df)
                log.debug("Experimental data retrieved for {}".format(pdb_id))
                time.sleep(sleep_time)
            except HTTPError as e:
                log.error("Experimental data not retrieved for {}: {}".format(pdb_id, e))
                continue
            except URLError as u:
                log.error("Experimental data not retrieved for {}: {}".format(pdb_id, u))
                continue
    master_exp_data_df = pd.concat(exp_data_dfs)
    master_exp_data_df.reset_index(drop = True, inplace = True)
    master_exp_data_df.to_pickle(out)
    return master_exp_data_df

## SIMPLIFYING PDB FILES

def get_simple_pdbs(trans_dir, simple_dir):
    """
    This function simplifies a group of PDB files that have been
    transformed and are imposed in space. It will only keep the
    ATOM records of the first file. For the rest, it will only
    save the HETATM records, corresponding to ligands.
    """
    pdb_files = [os.path.join(trans_dir, f) for f in os.listdir(trans_dir) if f.endswith("_trans.pdb")]
    first_simple = os.path.join(simple_dir, os.path.basename(pdb_files[0]))
    if os.path.isfile(first_simple):
        pass
    else:
        shutil.copy(pdb_files[0], first_simple)
    for pdb_in in pdb_files[1:]:
        pdb_id = os.path.basename(pdb_in)[:6]
        pdb_out = os.path.join(simple_dir, os.path.basename(pdb_in))
        if os.path.isfile(pdb_out):
            continue
        pdb_df = PDBXreader(inputfile = pdb_in).atoms(format_type = "pdb", excluded=())
        hetatm_df = pdb_df.query('group_PDB == "HETATM"') #[pdb_df.group_PDB == "HETATM"]
        if len(hetatm_df) == 0:
            log.warning("No HETATM records in {}".format(pdb_id))
            continue
        hetatm_df = hetatm_df.replace({"label_alt_id": ""}, " ")
        w = PDBXwriter(outputfile = pdb_out)
        w.run(hetatm_df, format_type = "pdb", category = "auth")
        log.debug("{} simplified".format(pdb_id))

## RELATIVE INTERSECTION, AND METRIC FUNCTIONS

def get_intersect_rel_matrix(binding_ress):
    """
    Given a set of ligand binding residues, calcualtes a
    similarity matrix between all the different sets of ligand
    binding residues.
    """
    inters = {i: {} for i in range(len(binding_ress))}
    for i in range(len(binding_ress)):
        inters[i][i] = intersection_rel(binding_ress[i], binding_ress[i])
        for j in range(i+1, len(binding_ress)):
            inters[i][j] = intersection_rel(binding_ress[i], binding_ress[j])
            inters[j][i] = inters[i][j]
    return inters # implement different distances, e.g., jaccard_sim

def intersection_rel(l1, l2):
    """
    Calculates relative intersection.
    """
    len1 = len(l1)
    len2 = len(l2)
    I_max = min([len1, len2])
    I = len(list(set(l1).intersection(l2)))
    return I/I_max

def jaccard_sim(l1, l2):
    """
    Calculates Jaccard Similiarity index.
    """
    I = len(list(set(l1).intersection(l2)))
    U = (len(set(l1)) + len(set(l2))) - I
    return float(I) / U

## PDBe-KB REST API CALLS PARSING FUNCTIONS

def get_lig_pdbs(acc):
    """
    For a given UniProt protein accession, it returns two lists:
    The first one contains all PDB IDs containing a protein-ligand interaction.
    The second one contains all the unique ligand identifiers.
    """
    df = pd.read_json("https://www.ebi.ac.uk/pdbe/graph-api/uniprot/ligand_sites/{}".format(acc), convert_axes = False, dtype = False)
    d = df.loc["data", acc]
    dd = pd.DataFrame(d)
    ligs = dd.accession.unique().tolist() # accession here is ligand name
    all_lig_pdbs = []
    for lig in ligs:
        dd_lig = pd.DataFrame(dd.query('accession == @lig').residues.tolist()[0])#[dd.accession == lig].residues.tolist()[0])
        all_lig_pdbs.extend(dd_lig.allPDBEntries.tolist()[0])
    all_un_lig_pdbs = list(set(all_lig_pdbs))
    return sorted(all_un_lig_pdbs)#, sorted(ligs) #return ligs here and then log.warning (no ligs in acc are relevant)

def get_lig_lab_dict(df_mols):
    """
    Creates a dictionary with bm_ids as keys,
    and ligand names as values given a dataframe
    indicating bound molecules in a PDB.
    """
    lig_lab_dict = {}
    for idx, row in df_mols.iterrows():
        lig_lab = "{}_{}_{}".format(row["chem_comp_id"], str(row["author_residue_number"]), row["chain_id"])
        lig_lab_dict[idx] = lig_lab
    return lig_lab_dict

def get_bound_mols(pdb_id, bound_mols_dir):
    """
    This function processes the json output of an API call
    and returns a formatted dataframe of the bound molecules
    of a given pdb structure.
    """
    bound_mols_out = os.path.join(bound_mols_dir, "{}_bound_mols.json".format(pdb_id))
    if os.path.isfile(bound_mols_out):
        bound_mols_df = pd.read_json(bound_mols_out, convert_axes = False, dtype = False)
        bound_mols_df.index = bound_mols_df.index.astype(int) # not sure I need this
    else:
        try:
            bound_mols_df = pd.read_json("https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecules/{}".format(pdb_id), convert_axes = False, dtype = False)
        except:
            log.error("Bound molecules data were not recovered for {}".format(pdb_id))
            return pd.DataFrame()
        idx = bound_mols_df.index.tolist()
    df_rows = []
    for i in idx:
        bm_id = bound_mols_df.loc[i, pdb_id]["bm_id"]
        bound_mol_i_df = pd.DataFrame(bound_mols_df.loc[i, pdb_id])
        bound_mol_i_dff = pd.DataFrame(bound_mol_i_df.loc["ligands", "composition"])
        bound_mol_i_dff["bmid"] = bm_id
        df_rows.append(bound_mol_i_dff)
    df_mols = pd.concat(df_rows)
    df_mols.set_index("bmid", inplace = True)
    return df_mols

def get_bound_mol_inters(pdb_id, bm_id, bound_mol_inters_dir):
    """
    This function returns a table containing the protein-ligand interactions
    for a given PDB ID and a biomolecule ID.
    """
    bound_mol_inters_out = os.path.join(bound_mol_inters_dir, "{}_{}_inters.json".format(pdb_id, bm_id))
    if os.path.isfile(bound_mol_inters_out):
        bound_inters_df = pd.read_json(bound_mol_inters_out, convert_axes = False, dtype = False)
        bound_inters_df.index = bound_inters_df.index.astype(int)
    else:
        try:
            bound_inters_df = pd.read_json("https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecule_interactions/{}/{}".format(pdb_id, bm_id), convert_axes = False, dtype = False) # I THINK WE SHOULD SAVE THESE
        except:
            log.critical("Bound molecule interactions data were not recovered for {}_{}".format(pdb_id, bm_id))
            return pd.DataFrame()
        bound_inters_df.to_json(bound_mol_inters_out)
    bound_inters_dff = pd.DataFrame(bound_inters_df.loc[0, pdb_id]["interactions"])
    begin_rws, end_rws, inter_rws = [[], [], []]
    for i in range(len(bound_inters_dff)):
        begin_rws.append(pd.DataFrame.from_dict(bound_inters_dff.loc[i,"begin"], orient = "index").T) # DO I EVEN NEED BEGIN?
        end_rws.append(pd.DataFrame.from_dict(bound_inters_dff.loc[i,"end"], orient = "index").T)
    begin_cols_df = pd.concat(begin_rws) # DO I EVEN NEED BEGIN?
    end_cols_df = pd.concat(end_rws)
    begin_cols_df.columns = [coli + "_begin" for coli in begin_cols_df.columns.tolist()] # DO I EVEN NEED BEGIN?
    end_cols_df.columns = [coli + "_end" for coli in end_cols_df.columns.tolist()]
    inters_df = pd.concat([begin_cols_df, end_cols_df], axis = 1)
    inters_df["lig_lab"] = inters_df.chem_comp_id_begin + "_" + inters_df.author_residue_number_begin.astype(str) + "_" + inters_df.chain_id_begin

    inters_df = inters_df.query('lig_lab == @lig_lab_dict[@bm_id]').copy() # filtering out interactions where ligand is not same bm_id. Don't know why they are even here

    return inters_df

def get_fingerprints(pdb_id, bound_mols_dir, bound_mol_inters_dir, out = None):
    """
    This function returns a dictionary with a ligand ID as key,
    and a list containing all the protein residues that interact with it
    as a value.
    """
    df_mols = get_bound_mols(pdb_id, bound_mols_dir)
    if df_mols.empty:
        return {}
    lig_lab_dict = get_lig_lab_dict(df_mols)
    bm_ids = df_mols.index.tolist()
    bm_id_dfs = []
    trying_counter = {bm_id: 1 for bm_id in bm_ids}
    bm_ids_copy = copy.deepcopy(bm_ids)
    while bm_ids_copy:
        bm_id = bm_ids_copy[0]
        try:
            if trying_counter[bm_id] > max_retry:
                log.warning("Maximum number of attempts: {} made to obtain fingerprints for {} in {}".format(max_retry, bm_id, pdb_id))
                bm_ids_copy.remove(bm_id) # trying while loop
            else:
                log.debug("Attempt #{} to get interactions for {} in {}".format(trying_counter[bm_id], bm_id, pdb_id))
                bm_id_df = get_bound_mol_inters(pdb_id, bm_id, lig_lab_dict, bound_mol_inters_dir)
                bm_id_df["bm_id"] = bm_id
                bm_id_dfs.append(bm_id_df)
                bm_ids_copy.remove(bm_id) # trying while loop
        except HTTPError as e:
            log.debug("Attempt #{} to retrieve interactions for {} in {} failed".format(str(trying_counter[bm_id]), bm_id, pdb_id))
            trying_counter[bm_id] += 1
    
    if bm_id_dfs == []: # Interactions retrieval failed for ALL bm_ids of this pdb
        log.error("None of the interactions were retrieved from GRAPH-API for {}".format(pdb_id))
        return {}
    else: # data has been retrieved about bms
        pdb_id_inters_df = pd.concat(bm_id_dfs)
        
        lig_lab_inters = {}
        lig_labs = pdb_id_inters_df.lig_lab.unique().tolist()
        for lig_lab in lig_labs:
            lig_lab_rows = pdb_id_inters_df.query('lig_lab == @lig_lab') #[pdb_id_inters_df.lig_lab == lig_lab]
            resnums = lig_lab_rows.author_residue_number_end.astype(str).tolist() # changing to string so SIFTS mapping works
            resnames = lig_lab_rows.chem_comp_id_end.tolist()
            chainids = lig_lab_rows.chain_id_end.tolist()
            lig_lab_inters[lig_lab] = list(zip(resnums, resnames, chainids)) #KEEPS INFORMATION OF NUMBER, NAME, AND CHAIN OF LIGAND-INTERACTING RESIDUES

        if out != None:
            with open(out, "wb") as f:
                pickle.dump(lig_lab_inters, f)

        log.debug("Fingerprints obtained for {}".format(pdb_id))

        return lig_lab_inters

def filter_fingerprints(pdb_id, lig_lab_inters, acc, segment_i_chains):
    """
    Given the protein-ligand fingerpints dictionary
    for a PDB, filters out those interactions with
    other proteins that are not the target protein.
    Chain filtering needs to happen even if it is
    only with one chain that the ligand is interacting.
    """
    lig_lab_inters_filt = copy.deepcopy(lig_lab_inters)
    for lig_name, fingerprint in lig_lab_inters_filt.items():
        inter_chains = sorted(list(set([inter[2] for inter in fingerprint])))
        n_chains = len(inter_chains)
        if n_chains > 1:
            log.debug("{} interacts with {} chains in {}".format(lig_name, n_chains, pdb_id))
        t_chains = [ch.split("_")[1] for ch in segment_i_chains if pdb_id in ch]
        lig_lab_inters_filt[lig_name] = [inter for inter in lig_lab_inters_filt[lig_name] if inter[2].split("_")[0] in t_chains]
        if len(lig_lab_inters_filt[lig_name]) == 0:
            log.debug("{} does not interact with {} in {}".format(lig_name, acc, pdb_id))
    return lig_lab_inters_filt

def get_fingerprints_dict(acc, fps_dir, out, all_ligs_pdbs, segment_i_chains, bound_mols_dir, bound_mol_inters_dir):
    """
    Returns a dictionary for all protein ligand interactions
    across all structures for a given UniProt accession and
    saves dictionary to pickle format. Can filter results using
    the pdb_set argument, containing the list of pdbs that infor-
    mation is required about.
    """

    all_fingerprints = {}
    all_ligs_pdbs_copy = copy.deepcopy(all_ligs_pdbs)
    while all_ligs_pdbs_copy: # trying while loop
        pdb_id = all_ligs_pdbs_copy[0]
        fp_out = os.path.join(fps_dir, "{}_fps.pkl".format(pdb_id))
        try:
            if os.path.isfile(fp_out):
                with open(fp_out, "rb") as f:
                    lig_lab_inters = pickle.load(f)
                    log.debug("Fingerprints read for {}".format(pdb_id))
            else:
                lig_lab_inters = get_fingerprints(pdb_id, bound_mols_dir, bound_mol_inters_dir, fp_out)
                time.sleep(sleep_time)
            if lig_lab_inters == {}:
                all_fingerprints[pdb_id] = lig_lab_inters
                pass
            else:
                lig_lab_inters_filt = filter_fingerprints(pdb_id, lig_lab_inters, acc, segment_i_chains)
                all_fingerprints[pdb_id] = lig_lab_inters_filt

        except HTTPError as e:
            all_fingerprints[pdb_id] = {}
            log.warning("No fingerprints were retrieved for {}".format(pdb_id))
        all_ligs_pdbs_copy.remove(pdb_id) 
        
    with open(out, "wb") as f:
        pickle.dump(all_fingerprints, f)

    return all_fingerprints

def get_labs(fingerprints_dict):
    """
    Returns all ligand labels from fingerprints dict.
    """
    return [k1 + "_" + k2 for k1, v1 in fingerprints_dict.items() for k2 in v1.keys()]

def get_inters(fingerprints_dict):
    """
    Returns all ligand fingerprints from fingerprints dict.
    """
    return [v2 for v1 in fingerprints_dict.values() for v2 in v1.values()]

def filter_non_relevant_ligs(fingerprints_dict, relevant_ligs):
    """
    Removes highly frequent ligands in the PDBe that might
    not be biologically relevant, but crystallisation artifacts.
    """
    fingerprints_dict_filt = {}
    for k1, v1 in fingerprints_dict.items(): #k1 is pdb_id
        fingerprints_dict_filt[k1] = {}
        for k2, v2 in v1.items(): #k2 is lig_name
            lig_name = k2.split("_")[0]
            if lig_name in relevant_ligs:
                fingerprints_dict_filt[k1][k2] = v2
            else:
                continue
        if fingerprints_dict_filt[k1] == {}:
            log.debug("No relevant ligands found in {}".format(k1))
    return fingerprints_dict_filt

def filter_non_protein_inters(fingerprints_dict, acc, aa_resnames):
    """
    Removes non-protein-ligand interactions from
    the fingerprints dictionary. Only relevant ligands-
    protein interactions remain now.
    """
    fingerprints_dict_filt = {}
    #fingerprints_dict_filt_v2 = {}
    for k1, v1 in fingerprints_dict.items():
        fingerprints_dict_filt[k1] = {}
        for k2, v2 in v1.items():
            inter_ress = [inter_res for inter_res in v2 if inter_res[1] in aa_resnames]
            if len(inter_ress) == 0:
                log.debug("{} does not interact with any protein atoms of {} in {}".format(k2, acc, k1))
                continue
            else:
                fingerprints_dict_filt[k1][k2] = inter_ress
                #fingerprints_dict_filt_v2[k1 + "_" + k2] = inter_ress
    return fingerprints_dict_filt#, fingerprints_dict_filt_v2 # possibly do not need to return second dict

def get_pdb_lig_freq(lig_id):
    """
    Returns the absolute frequency of a ligand in PDBe structures,
    i.e., how many structures present this ligand.
    """
    return len(pd.read_json("https://www.ebi.ac.uk/pdbe/graph-api/compound/in_pdb/{}".format(lig_id)), convert_axes = False, dtype = False)

## GETTING PDB - UNIPROT MAPPINGS WITH SIFTS

def get_mapping_from_sifts(pdb_id):
    """
    Given a PDB ID, returns a dictionary with the PDB-UniProt
    residue mapping for each chain. It does this by parsing
    a SIFTs .xml.gz file via ProIntVar.
    """
    input_sifts = os.path.join(sifts_db_path, "{}/{}.xml.gz".format(pdb_id[1:3], pdb_id))

    if not os.path.isfile(input_sifts):
        input_sifts = os.path.join(sifts_db_path2, "{}.xml".format(pdb_id))  # this one is not g-zipped
        if not os.path.isfile(input_sifts):
            log.error("SIFTS data could not be found for {}".format(pdb_id)) # this is now the local download made by JSU
            return {}
    try:
        sifts_df = SIFTSreader(inputfile = input_sifts).read()
        try:
            pdb2up = sifts_df[["PDB_dbResNum", "UniProt_dbResNum", "PDB_dbChainId"]]
        except KeyError as e:
            log.error("PDB-UniProt mapping columns were not found in SIFTS table for {}".format(pdb_id))
            return {}
        pdb2up = pdb2up[pdb2up.PDB_dbResNum != "null"].dropna()
        chains = pdb2up.PDB_dbChainId.unique().tolist()
        mapping_dict = {}
        for chain in chains:
            chain_df = pdb2up[pdb2up.PDB_dbChainId == chain]
            mapping_dict[chain] = {}
            mapping_dict[chain] = {chain_df.loc[i, "PDB_dbResNum"]: int(chain_df.loc[i, "UniProt_dbResNum"]) for i in chain_df.index} #removed int for cases such as '63A'
        return mapping_dict
    except IOError as e:
        log.error("SIFTS data could not be read for {}".format(pdb_id))
        return {}

def get_up_mapping_from_prointvar(fps_dict, mappings_dict):
    """
    Given a protein-ligand fingerprints dictionary, and a sifts 
    PDB-UniProt dictionary, returns the fingerprints dictionary
    with the UniProt resnums instead of PDB ones.
    """
    fps_sifted_dict = {}
    fps_sifted_dict_v2 = {}
    for pdb_id, lig_inters in fps_dict.items():
        if len(lig_inters) == 0:
            continue
        fps_sifted_dict[pdb_id] = {}
        mappings_pdb = mappings_dict[pdb_id]
        for lig_name, inters in lig_inters.items():
            fps_sifted_dict[pdb_id][lig_name] = []
            
            v2_k = pdb_id + "_" + lig_name
            fps_sifted_dict_v2[v2_k] = []
            lig_finger = fps_dict[pdb_id][lig_name]

            fps_sifted_dict[pdb_id][lig_name] = []
            fps_sifted_dict_v2[v2_k] = []

            for i in lig_finger:
                try:
                    finger_chain = i[2]
                    finger_chain_id = finger_chain.split("_")[0]
                    mappings_chain = mappings_pdb[finger_chain_id] # assumes ligand chain is the same as protein (might not always be true)
                    fps_sifted_dict[pdb_id][lig_name].append(mappings_chain[i[0]])
                    fps_sifted_dict_v2[v2_k].append(mappings_chain[i[0]])
                except KeyError as e:
                    log.error("Residue {} in chain {} of {} has no SIFTS mapping".format(i, finger_chain_id, pdb_id))
                    continue


    return fps_sifted_dict, fps_sifted_dict_v2

def get_mol_type_dict(pdb_id):
    """
    For a given PDB, returns a dictionary
    where keys are chain ids, and values
    is a list containing the types of the
    molecules within that chain.
    """
    d = pd.read_json("https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{}".format(pdb_id), convert_axes = False, dtype = False)
    mol_type_dict = {}
    all_chains = []
    for i in d.index:
        d_i = d.loc[i, pdb_id]
        ent_id = d_i["entity_id"]
        in_chains = d_i["in_chains"]
        mol_type = d_i["molecule_type"]
        mol_type_dict[ent_id] = {}
        for chain in in_chains:
            all_chains.append(chain)
            mol_type_dict[ent_id][chain] = mol_type
    all_chains = sorted(list(set(all_chains)))
    all_chains_types = {ch: [] for ch in all_chains}
    for v1 in mol_type_dict.values():
        for k2, v2 in v1.items():
            all_chains_types[k2].append(v2)
    all_chains_types = {k: sorted(list(set(v))) for k, v in all_chains_types.items()}
    return all_chains_types

## CHIMERA COLOURING FUNCTIONS, AND VISUALISATION

def get_chimera_data(cluster_id_dict):
    """
    Gets chimera atom specs, binding site ids, and paths
    to pdb files to generate the attribute files later, and
    eventually colour models.
    """
    chimera_atom_specs, bs_ids, pdb_paths = [[], [], []]
    for k, v in cluster_id_dict.items():
        ld = k.split("_") # stands for lig data
        pdb_id, lig_resname, lig_resnum, lig_chain_id = [ld[0], ld[1], ld[2], ld[3]] # not sure why sometimes chain ID is A_1, and sometimes just A. Q9UKK9, 5qjj.
            
        bs_id = str(v)
        pdb_path = "{}_{}_trans.pdb".format(pdb_id, lig_chain_id)
        chimera_atom_spec = (
            ":"+ str(lig_resnum) +
            "."+ lig_chain_id +
            "&#/name==" + pdb_path
        )
        chimera_atom_specs.append(chimera_atom_spec)
    return chimera_atom_specs

def write_chimera_attr(attr_out, chimera_atom_specs, cluster_ids):
    """
    Writes chimera attribute file, which indicates in which cluster
    each ligand is found, and assigns it a colour.
    """
    with open(attr_out, "w") as out:
        out.write("attribute: binding_site\n")
        out.write("match mode: 1-to-1\n")
        out.write("recipient: residues\n")
        for i in range(len(cluster_ids)):
            out.write("\t" + chimera_atom_specs[i] + "\t" + str(cluster_ids[i]) + "\n")
    return

def write_chimera_command(chimera_script_out, cmds, cluster_ids, sample_colors):
    """
    Writes chimera command file, which colours all models according
    to the attribute file.
    """
    with open(chimera_script_out, "w") as out:
        out.write("# neutral colour for everything not assigned a cluster\n")
        out.write("colour white\n")

        out.write("# colour each binding site\n")
        for i in range(0, len(sorted(list(set(cluster_ids))))):
            out.write("colour {} :/binding_site=={}\n".format(",".join(list(map(str, list(sample_colors[i])))), i))
        out.write("### SOME FORMATTING ###\n")
        out.write("\n".join(cmds))
    return

## CLUSTER ANALYSIS UTILS

def get_cluster_membership(cluster_id_dict):
    """
    Creates a dictionary indicating to which cluster
    each ligand binds to.
    """
    membership_dict = {}
    for k, v in cluster_id_dict.items():
        if v not in membership_dict:
            membership_dict[v] = []
        membership_dict[v].append(k)
    return membership_dict

def get_all_cluster_ress(membership_dict, fingerprints_dict):
    """
    Given a membership dict and a fingerprint dictionary,
    returns a dictionary that indicates the protein residues
    forming each binding site.
    """
    binding_site_res_dict = {}
    for k, v in membership_dict.items():
        if k not in binding_site_res_dict:
            binding_site_res_dict[k] = []
        for v1 in v:
            binding_site_res_dict[k].extend(fingerprints_dict[v1])
    binding_site_res_dict = {k: sorted(list(set(v))) for k, v in binding_site_res_dict.items()}
    return binding_site_res_dict

def get_residue_bs_membership(cluster_ress):
    """
    Returns a dictionary indicating to which ligand binding
    site each ligand binding residue is found in. A residue
    might contribute to more than one adjacent binding site.
    """
    all_bs_ress = []
    for v in cluster_ress.values():
        all_bs_ress.extend(v)
    all_bs_ress = sorted(list(set(all_bs_ress)))
    
    bs_ress_membership_dict = {}
    for bs_res in all_bs_ress:
        bs_ress_membership_dict[bs_res] = []
        for k, v in cluster_ress.items():
            if bs_res in v:
                bs_ress_membership_dict[bs_res].append(k) # which binding site each residue belongs to
    return bs_ress_membership_dict

## DSSP

def run_dssp(pdb_path, trans_dir, dssp_dir):
    """
    Runs DSSP unbound, and saves formatted resulting output dataframe.
    pdb_path is a relative path to input pdb file, which absoulute path
    results from joining it with trans_dir.
    """
    
    struc_in = os.path.join(trans_dir, pdb_path)
    pdb_root, _ = os.path.splitext(pdb_path)
    dssp_out = os.path.join(dssp_dir, pdb_root + ".dssp")
    dssp_pickle = os.path.join(dssp_dir, pdb_root + ".pkl") # output pickle filepath
    #dssp_out = os.path.join(dssp_dir, pdb_path.replace("pdb", "dssp"))
    #dssp_pickle = os.path.join(dssp_dir, pdb_path.replace("pdb", "pkl")) # output pickle filepath
    if os.path.isfile(dssp_out):
        pass
    else:
        DSSPrunner(inputfile = struc_in, outputfile = dssp_out).write(run_unbound = True)            # runs DSSP

    if os.path.isfile(dssp_pickle):
        return
    else:
        dssp_data = DSSPreader(inputfile = dssp_out).read()            # reads DSSP output
        dssp_data = dssp_data.rename(index = str, columns = {"RES": "PDB_ResNum"})
        dssp_data.PDB_ResNum = dssp_data.PDB_ResNum.astype(str)
        dssp_data.to_pickle(dssp_pickle)#, index = False)

def get_dssp_data(trans_dir, dssp_dir, sifts_mapping_dict, out):
    """
    Given a dir with transformed files, output dssp dir, and
    sifts mapping dictionary, runs DSSP for all structures and
    returns a dataframe with DSSP data from all structures.
    """
    trans_files = [f for f in os.listdir(trans_dir) if f.endswith("_trans.pdb")]
    all_dssp_dfs = []
    for trans_file in trans_files:
        # skip those files which do not include target protein, therefore, have no SIFTS mapping?
        pdb_id, chain_id, _ = trans_file.split("_")
        try:
            mol_type_dict = get_mol_type_dict(pdb_id)
        except:
            log.error("Molecule type data could not be retrieved for {}".format(pdb_id))
            continue
        chain_mol_type = mol_type_dict[chain_id]
        if "polypeptide(L)" not in chain_mol_type: # NOT RUNNING DSSP FOR ONLY NON-PROTEIN CHAINS
            log.warning("Chain {} in {} includes {}. Skipping DSSP.".format(chain_id, pdb_id, chain_mol_type))
            continue
        else:
            try:
                run_dssp(trans_file, trans_dir, dssp_dir)
            except:
                pdb_in = os.path.join(trans_dir, trans_file)
                pdb_df = PDBXreader(inputfile = pdb_in).atoms(format_type = "pdb", excluded=())
                pdb_df_atom = pdb_df.query('group_PDB == "ATOM"')
                un_atoms = pdb_df_atom.label_atom_id.unique().tolist()
                if un_atoms == ["CA"]:
                    log.warning("Only CA atoms found in {}. Skipping DSSP".format(pdb_in)) # DSSP crashes, for example 3jcr_C. Only CA atoms.
                elif set(un_atoms).issubset(set(bbone_atoms)):
                    log.warning("Only backbone atoms found in {}. Skipping DSSP".format(pdb_in)) # incomplete residues, CA + rest of backbone, but no side chains, e.g., 6yw7
                else:
                    log.error("Unknown DSSP Error for {}".format(pdb_in)) 
                continue
            trans_root, _ = os.path.splitext(trans_file)
            dssp_df = pd.read_pickle(os.path.join(dssp_dir, trans_root + ".pkl"))
            #dssp_df = pd.read_pickle(os.path.join(dssp_dir, trans_file.replace("pdb", "pkl")))
            dssp_df.PDB_ResNum = dssp_df.PDB_ResNum.astype(str)
            try:
                dssp_df["UniProt_ResNum"] = dssp_df["PDB_ResNum"].map(sifts_mapping_dict[pdb_id][chain_id])
            except KeyError as e:
                dssp_df["UniProt_ResNum"] = np.nan
                log.error("SIFTS mapping could not be found for chain {} in {}".format(chain_id, pdb_id))
            dssp_df["PDB_ID"] = pdb_id
            all_dssp_dfs.append(dssp_df)
    if all_dssp_dfs == []: # list of dssp dfs is empty, e.g., all structures have CA atoms only (Q12840)
        return pd.DataFrame()
    else:
        master_dssp_df = pd.concat(all_dssp_dfs)
        master_dssp_df.reset_index(drop = True, inplace = True)
        master_dssp_df.to_pickle(out)#, index = False)
        return master_dssp_df

## MSA UTILS

def get_best_from_segment_data(segment_data):
    """
    Gets best pdb_id and chain_id for a given
    segment of structural coverage of a protein
    from the superimposition data from PDBe-KB.
    """
    for v in segment_data.values():
        if v["is_representative"] == True:
            #best_pdb_id = v["pdb_id"]
            #best_chain_id = v["auth_asym_id"]
            return v
        else:
            continue

def get_best_seq_SOLR(best_pdb_id, best_chain_id):
    """
    Obtains the sequence for the best PDB,
    and best chain for a given accession. The
    pdb_id and chain_id passed must already be
    the best representations for the protein.
    """
    slr_query_df = pd.read_json("https://www.ebi.ac.uk/pdbe/search/pdb/select?q=pdb_id:{}&fl=chain_id+molecule_sequence&wt=json".format(best_pdb_id), convert_axes = False, dtype = False)
    seq_dict = {j: i["molecule_sequence"] for i in slr_query_df.loc["docs", "response"] for j in i["chain_id"]} # now it takes into account multiple chains with the same sequence
    best_seq = seq_dict[best_chain_id]
    return best_seq

def get_best_struct_seq(acc, segment, out, best = None):
    """
    Using pdbe APIs, gets the sequence from the best
    structure for a given protein and prints it in
    fasta format.
    """
    if best != None:
        best_pdb_id = best["pdb_id"]
        best_entity_id = best["entity_id"]
        best_chain_id = best["auth_asym_id"]
        try:
            d = pd.read_json("https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages/domains/{}/{}".format(best_pdb_id, best_entity_id), convert_axes = False, dtype = False)
            best_seq = d.loc["sequence", best_pdb_id]
            log.debug("Getting sequence from SEGMENT REPRESENTATIVE structure through DOMAINS GRAPH-API for Segment {} of {}".format(str(segment), acc))
        except HTTPError as e:
            log.warning("Could not retrieve sequence from DOMAINS GRAPH-API for entity {}, chain {} in {} for Segment {} of {}".format(str(best_entity_id), best_chain_id, best_pdb_id, str(segment), acc))
            try:
                best_seq = get_best_seq_SOLR(best_pdb_id, best_chain_id)
                log.debug("Getting sequence from SEGMENT REPRESENTATIVE structure from Solr-based query system for Segment {} of {}".format(str(segment), acc))
            except:
                log.warning("Could not retrieve sequence from SOLR-query API for entity {}, chain {} in {} for Segment {} of {}".format(str(best_entity_id), best_chain_id, best_pdb_id, str(segment), acc))
                try:
                    df = pd.read_json("https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{}".format(best_pdb_id), convert_axes = False, dtype = False)
                    mols_df = pd.DataFrame(df[best_pdb_id].tolist())
                    best_seq = mols_df.query('entity_id == @best_entity_id').sequence.tolist()[0]
                    log.debug("Getting sequence from SEGMENT REPRESENTATIVE structure from molecules GRAPH-API for Segment {} of {}".format(str(segment), acc))
                except:
                    log.error("Could not get sequence for Segment {} of {}".format(str(segment), acc))

    else:
        best_structure = pd.read_json("https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/{}".format(acc), convert_axes = False, dtype = False).loc[0, acc]
        best_pdb_id = best_structure["pdb_id"]
        best_chain_id = best_structure["chain_id"]   
        best_seq = get_best_seq_SOLR(best_pdb_id, best_chain_id)
        log.debug("Getting sequence from BEST_STRUCTURES through GRAPH-API and Solr-based query system for Segment {} of {}".format(str(segment), acc))

    best_seq_id = "{}_{}_{}".format(acc, best_pdb_id, best_chain_id)
    with open (out, "w") as f:
        f.write(">{}\n{}\n".format(best_seq_id, best_seq))
    return best_seq_id

def jackhmmer(seq, hits_out, hits_aln, n_it = jackhmmer_n_it, seqdb = swissprot):
    """
    Runs jackhmmer on an input seq for a number of iterations and returns exit code, should be 0 if all is ok.
    """
    args = ["jackhmmer", "--acc", "-N", str(n_it), "-o", hits_out, "-A", hits_aln, seq, seqdb]
    cmd = " ".join(args)
    exit_code = os.system(cmd)
    if exit_code != 0:
        log.critical("Jackmmer did not work with command: {}".format(cmd))

def add_acc2msa(aln_in, aln_out, query_id, fmt_in = MSA_fmt):
    """
    Modifies AC field of jackhmmer alignment in stockholm format.
    
    :param aln_in: path of input alignment
    :type aln_in: str, required
    :param aln_out: path of output alignment
    :type aln_in: str, required
    :param fmt_in: input and output MSA format
    :type aln_in: str, defaults to stockholm
    """
    aln = Bio.SeqIO.parse(aln_in, fmt_in)
    recs = []
    for rec in aln:
        if rec.id == query_id:
            continue
        else:
            rec.annotations["accession"] = rec.id.split("|")[1]
            recs.append(rec)
    Bio.SeqIO.write(recs, aln_out, fmt_in)

def get_target_prot_cols(msa_in, query_id, msa_fmt = MSA_fmt): 
    """
    Returns list of MSA col idx that are popualted on the protein target.
    """
    seqs = [str(rec.seq) for rec in Bio.SeqIO.parse(msa_in, msa_fmt) if rec.id == query_id]
    occupied_cols = [i+1 for seq in seqs for i, el in enumerate(seq) if el != "-"]
    return sorted(list(set(occupied_cols)))

def get_human_subset_msa(aln_in, human_msa_out, fmt_in = MSA_fmt):
    """
    Creates a subset MSA containing only human sequences.
    """
    msa = Bio.AlignIO.read(aln_in, fmt_in)
    human_recs = []
    for rec in msa:
        if "HUMAN" in rec.name:
            human_recs.append(rec)
    Bio.SeqIO.write(human_recs, human_msa_out, fmt_in)

def generate_subset_aln(aln_in, aln_fmt, df, aln_out = None):
    """
    Creates a subset MSA containing only human sequences that present
    missense variants and returns the path of such MSA.
    """
    seqs_ids = df.source_id.unique().tolist()
    aln = Bio.SeqIO.parse(aln_in, aln_fmt)
    variant_seqs = [rec for rec in aln if rec.id in seqs_ids]
    n_variant_seqs = len(variant_seqs)
    if n_variant_seqs == 0:
        return ""
    else:
        log.info("There are {} sequences with variants for {}".format(str(n_variant_seqs), aln_in))
    if aln_out == None:
        pref, fmt = aln_in.split(".")
        aln_out =  pref + "_variant_seqs." + fmt
    Bio.SeqIO.write(variant_seqs, aln_out, aln_fmt)
    return aln_out

## SHENKIN CALCULATION FUNCTIONS 

def get_freqs(i_col, col):
    """
    Calculates amino acid frequences for a given MSA column.
    """
    abs_freqs = {
        "A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "Q": 0, "E": 0, "G": 0, "H": 0, "I": 0,
        "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0, "-": 0
    }
    non_standard_aas = {}
    for aa in col:
        aa = aa.upper()
        if col.count("-") == len(col):
            abs_freqs["-"] = 1
            return abs_freqs
        if aa in aas_1l:
            abs_freqs[aa] += 1
        else:
            if aa not in non_standard_aas:
                non_standard_aas[aa] = 0
            non_standard_aas[aa] += 1
    all_ns_aas = sum(non_standard_aas.values())
    if all_ns_aas != 0:
        log.warning("Column {} presents non-standard AAs: {}".format(str(i_col), non_standard_aas))
    rel_freqs = {k: v/(len(col) - all_ns_aas) for k, v in abs_freqs.items()}
    return rel_freqs

def get_entropy(freqs):
    """
    Calculates Shannon's entropy from a set of aa frequencies.
    """
    S = 0
    for f in freqs.values():
        if f != 0:
            S += f*math.log2(f)
    return -S

def get_shenkin(i_col, col):
    """
    Calculates Shenkin score for an MSA column.
    """
    S = get_entropy(get_freqs(i_col, col))
    return round((2**S)*6,2)

def in_columns(aln_in, infmt):
    """
    Returns dictionary in which column idx are the key
    and a list containing all aas aligned to that column
    is the value.
    """
    aln = Bio.AlignIO.read(aln_in, infmt)
    n_cols = len(aln[0])
    cols = {}
    for col in range(1,n_cols+1):
        cols[col] = []
    for row in aln:
        seq = str(row.seq)
        for i in range(0,len(seq)):
            cols[i+1].append(seq[i])
    return cols

def get_stats(col):
    """
    For a given MSA column, calculates some basic statistics
    such as column residue occupancy ang gaps.
    """
    n_seqs = len(col)
    gaps = col.count("-")
    occ = n_seqs - gaps
    occ_pct = round(100*(occ/n_seqs), 2)
    gaps_pct = round(100-occ_pct, 2)
    return occ, gaps, occ_pct, gaps_pct

def calculate_shenkin(aln_in, aln_fmt, out = None):
    """
    Given an MSA, calculates Shenkin ans occupancy, gap
    percentage for all columns.
    """
    cols = in_columns(aln_in, aln_fmt)
    scores = []
    occ = []
    gaps = []
    occ_pct = []
    gaps_pct = []
    for k, v in cols.items():
        scores.append(get_shenkin(k, v))
        stats = (get_stats(v))
        occ.append(stats[0])
        gaps.append(stats[1])
        occ_pct.append(stats[2])
        gaps_pct.append(stats[3])
    df = pd.DataFrame(list(zip(list(range(1,len(scores)+1)),scores, occ,gaps, occ_pct, gaps_pct)), columns = ["col", "shenkin", "occ", "gaps", "occ_pct", "gaps_pct"])
    if out != None:
        df.to_pickle(out)
    return df

def format_shenkin(shenkin, prot_cols, out = None):
    """
    Formats conservation dataframe and also
    calculates two normalised versions of it.
    """
    shenkin_filt = shenkin[shenkin.col.isin(prot_cols)].copy()
    shenkin_filt.index = range(1, len(shenkin_filt) + 1) # CONTAINS SHENKIN SCORE, OCCUPANCY/GAP PROPORTION OF CONSENSUS COLUMNS
    min_shenkin = min(shenkin_filt.shenkin)
    max_shenkin = max(shenkin_filt.shenkin)
    shenkin_filt.loc[:, "rel_norm_shenkin"] = round(100*(shenkin_filt.shenkin - min_shenkin)/(max_shenkin - min_shenkin), 2) # ADDING NEW COLUMNS WITH DIFFERENT NORMALISED SCORES
    shenkin_filt.loc[:, "abs_norm_shenkin"] = round(100*(shenkin_filt.shenkin - 6)/(120 - 6), 2)
    if out != None:
        shenkin_filt.to_pickle(out)#, index = False)
    return shenkin_filt

## VARIANT DATA

def format_variant_table(df, col_mask, vep_mask = ["missense_variant"], tab_format = "gnomad"):
    """
    Formats variant table, by gettint rid of empty rows that are not human sequences,
    changning column names and only keeping those variants that are of interest and
    are present in columns of interest.
    """
    df_filt = df.copy(deep = True)
    df_filt.reset_index(inplace = True)
    if tab_format == "gnomad":
        df_filt.columns = [" ".join(col).strip() for col in df_filt.columns.tolist()]
    df_filt.columns = [col.lower().replace(" ", "_") for col in df_filt.columns.tolist()]
    df_filt = df_filt[df_filt.source_id.str.contains("HUMAN")]
    df_filt = df_filt.dropna(subset = ["vep_consequence"])
    df_filt = df_filt[df_filt.vep_consequence.isin(vep_mask)]
    df_filt = df_filt[df_filt.alignment_column.isin(col_mask)]
    return df_filt

def get_missense_df(aln_in, variants_df, shenkin_aln, prot_cols, aln_out, aln_fmt = MSA_fmt, get_or = True):
    """
    Generates a dataframe for the subset of human sequences with variants
    mapping to them. Calculates shenkin, and occupancy data, and then
    enrichment in variants.
    """
    variants_aln = generate_subset_aln(aln_in, aln_fmt, variants_df, aln_out)
    if variants_aln == "":
        return pd.DataFrame()
    variants_aln_info = calculate_shenkin(variants_aln, aln_fmt)
    variants_aln_info = variants_aln_info[variants_aln_info.col.isin(prot_cols)]
    vars_df = pd.DataFrame(variants_df.alignment_column.value_counts().reindex(prot_cols, fill_value = 0).sort_index()).reset_index()
    vars_df.index = range(1, len(prot_cols) + 1)
    vars_df.columns = ["col", "variants"]
    merged = pd.merge(variants_aln_info, vars_df, on = "col", how = "left")
    merged.index = range(1, len(vars_df) + 1)
    merged["shenkin"] = shenkin_aln["shenkin"]
    merged["rel_norm_shenkin"] = shenkin_aln["rel_norm_shenkin"] 
    merged["abs_norm_shenkin"] = shenkin_aln["abs_norm_shenkin"]
    if get_or == True:
        merged_or = get_OR(merged)
        return merged_or
    else:
        return merged

def add_miss_class(df, miss_df_out = None, cons_col = "shenkin", MES_t = MES_t, cons_ts = cons_ts, colours = consvar_class_colours):
    """
    Adds two columns to missense dataframe. These columns will put columns
    into classes according to their divergence and missense enrichment score.
    It also adds a column that will help colour MSA columns according to their
    classifications.
    """
    for i in df.index:
        if df.loc[i, cons_col] <= cons_ts[0] and df.loc[i, "oddsratio"] < MES_t:
            df.loc[i, "miss_class"] = "CMD"
        elif df.loc[i, cons_col] <= cons_ts[0] and df.loc[i, "oddsratio"] > MES_t:
            df.loc[i, "miss_class"] = "CME"
        elif df.loc[i, cons_col] >= cons_ts[1] and df.loc[i, "oddsratio"] < MES_t:
            df.loc[i, "miss_class"] = "UMD"
        elif df.loc[i, cons_col] >= cons_ts[1] and df.loc[i, "oddsratio"] > MES_t:
            df.loc[i, "miss_class"] = "UME"
        else:
            df.loc[i, "miss_class"] = "None"
    coloring = {
        "CMD": colours[0],
        "CME": colours[1],
        "UMD": colours[3],
        "UME": colours[4],
        "None": colours[2]
    }
    df["miss_color"] =  df.miss_class.map(coloring)
    
    if miss_df_out != None:
            df.to_pickle(miss_df_out)
    return df

def merge_shenkin_df_and_mapping(shenkin_df, mapping_df, aln_ids):
    """
    Merges conservation, and variation table with MSA-UniProt
    mapping table, so conservation and variation data
    are mapped to UniProt residues.
    """
    shenkin_df = shenkin_df.rename(index = str, columns = {"col": "alignment_column"}) # renaming columns to be consistent with other StruVarPi dataframes
    prot_mapping = mapping_df.copy(deep = True).loc[aln_ids]
    prot_mapping.columns = prot_mapping.columns.droplevel(1)
    prot_mapping.reset_index(inplace = True)
    prot_mapping = prot_mapping.rename_axis(None, axis = "columns")
    prot_mapping.rename(index = None, columns = {"Alignment": "MSA_column", "Protein_position": "UniProt_ResNum"}, inplace = True)
        
    # Merging the VarAlign data to the Pfam alignment to gain conservation and variation data for the whole family...
    mapped_data = pd.merge(
        prot_mapping[["MSA_column", "UniProt_ResNum"]], shenkin_df,
        left_on = "MSA_column", right_on = "alignment_column"
    ).drop("MSA_column", axis = 1)
    return mapped_data

def get_OR(df, variant_col = "variants"):
    """
    Calculates OR, ln(OR) and associated p-value and CI,
    given a missense dataframe with variants and occupancy.
    """
    tot_occ = sum(df.occ)
    tot_vars = sum(df[variant_col])
    idx = df.index.tolist()
    for i in idx:
        i_occ = df.loc[i,"occ"]
        i_vars = df.loc[i,variant_col]
        rest_occ = tot_occ - i_occ
        rest_vars = tot_vars - i_vars
        if i_occ == 0:
            oddsr = np.nan 
            pval = np.nan
            se_or = np.nan
            log.debug("0 occupancy. Returning np.nan")
        else:
            if i_vars == 0:
                i_occ += 0.5
                i_vars += 0.5
                rest_occ += 0.5
                rest_vars += 0.5
                log.debug("0 variants. Adding 0.5 to each cell")
            oddsr, pval = stats.fisher_exact([[i_vars, rest_vars], [i_occ, rest_occ]])
            vals = [i_vars, rest_vars, i_occ, rest_occ]
            se_or = 1.96*(math.sqrt(sum(list(map((lambda x: 1/x), vals)))))
        df.loc[i, "oddsratio"] = round(oddsr, 2)
        df.loc[i, "pvalue"] = round(pval, 2)
        df.loc[i, "se_OR"] = round(se_or, 2)
    return df

### COLORS ###

# This code I did not write, I grabbed it from the following URL:

# https://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors

import colorsys
import itertools
from fractions import Fraction
from typing import Iterable, Tuple

def zenos_dichotomy() -> Iterable[Fraction]:
    """
    http://en.wikipedia.org/wiki/1/2_%2B_1/4_%2B_1/8_%2B_1/16_%2B_%C2%B7_%C2%B7_%C2%B7
    """
    for k in itertools.count():
        yield Fraction(1,2**k)

def fracs() -> Iterable[Fraction]:
    """
    [Fraction(0, 1), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4), Fraction(1, 8), Fraction(3, 8), Fraction(5, 8), Fraction(7, 8), Fraction(1, 16), Fraction(3, 16), ...]
    [0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 0.0625, 0.1875, ...]
    """
    yield Fraction(0)
    for k in zenos_dichotomy():
        i = k.denominator # [1,2,4,8,16,...]
        for j in range(1,i,2):
            yield Fraction(j,i)

# can be used for the v in hsv to map linear values 0..1 to something that looks equidistant
# bias = lambda x: (math.sqrt(x/3)/Fraction(2,3)+Fraction(1,3))/Fraction(6,5)

HSVTuple = Tuple[Fraction, Fraction, Fraction]
RGBTuple = Tuple[float, float, float]

def hue_to_tones(h: Fraction) -> Iterable[HSVTuple]:
    for s in [Fraction(6,10)]: # optionally use range
        for v in [Fraction(8,10),Fraction(5,10)]: # could use range too
            yield (h, s, v) # use bias for v here if you use range

def hsv_to_rgb(x: HSVTuple) -> RGBTuple:
    return colorsys.hsv_to_rgb(*map(float, x))

flatten = itertools.chain.from_iterable

def hsvs() -> Iterable[HSVTuple]:
    return flatten(map(hue_to_tones, fracs()))

def rgbs() -> Iterable[RGBTuple]:
    return map(hsv_to_rgb, hsvs())

def rgb_to_css(x: RGBTuple) -> str:
    uint8tuple = map(lambda y: int(y*255), x)
    return "rgb({},{},{})".format(*uint8tuple)

def css_colors() -> Iterable[str]:
    return map(rgb_to_css, rgbs())

### EXECUTING CODE

if __name__ == '__main__': ### command to run form command line: python3.6 fragsys_pdbe.py acc

    ### INITIATING LOGGER

    log.info("Logging initiated")

    ### PARSING COMMAND LINE ARGUMENTS

    parser = argparse.ArgumentParser(description = "Clusters ligands and defines binding sites.")
    parser.add_argument("up_acc", type = str, help = "UniProt accession number of the protein of interest.")
    parser.add_argument("--override", help = "Override any previously generated files.", action = "store_true")
    parser.add_argument("--override_variants", help = "Override any previously generated files (ONLY VARIANTS SECTION).", action = "store_true")
    #parser.add_argument("--no_transform", help = "Cleans, transforms, and simplifies structures.", action = "store_true")
    #parser.add_argument("--no_experimental", help = "Downloads, and processes experimental data for all structures.", action = "store_true")
    #parser.add_argument("--no_variants", help = "Retrieves Human variants form MSA and generates tables.", action = "store_true")

    args = parser.parse_args()

    acc = args.up_acc
    override = args.override
    override_variants = args.override_variants
    #run_transform = args.transform
    #run_experimental = args.experimental
    #run_variants = args.variants

    for arg, value in sorted(vars(args).items()):
        log.info("Argument %s: %r", arg, value)

    ### RETRIEVES ALL SUPERPOSITION MATRICES FOR PDB IDS IN acc, EXCEPT IF THERE ARE NOT ANY SOLVED STRUCTURES

    supp_mat_out = os.path.join(MATS_FOLDER, "{}_supp_mat.json".format(acc)) # had to change to json because of pickle issue
    if os.path.isfile(supp_mat_out):
        matrices_df = pd.read_json(supp_mat_out, convert_axes = False, dtype = False)
        log.info("Matrix table was read for {} and contains {} chains".format(acc, str(len(matrices_df))))
    else:
        try:
            matrices_df = pd.read_json("http://ftp.ebi.ac.uk/pub/databases/pdbe-kb/superposition/{}/{}/{}.json".format(acc[0], acc, acc), convert_axes = False, dtype = False).T
            matrices_df.to_json(supp_mat_out)
            log.info("Matrix table contains data from {} chains for {}".format(str(len(matrices_df)), acc))
        except HTTPError as e:
            log.warning("No structures solved for {}. Exiting programme".format(acc))
            print("{}\t{}".format(acc, str(15)), flush = True)
            sys.exit(0)

    ### RETRIEVES ALL LIGAND-BINDING PDB IDS FOR acc, EXCEPT IF THERE ARE NOT ANY LIGAND-BINDING STRUCTURES

    biolip_dict = load_pickle(biolip_data)
    try:
        all_ligs_pdbs = list(biolip_dict[acc].keys()) # NOW DONE WITH BIOLIP
        n_all_ligs_pdbs = len(all_ligs_pdbs)
        log.info("There are {} ligand-binding structures for {}".format(str(n_all_ligs_pdbs), acc))
    except KeyError as e:
        log.warning("No ligand-binding structures for {}. Exiting programme".format(acc))
        print("{}\t{}".format(acc, str(16)), flush = True)
        sys.exit(0)

    ### CREATES WORKING DIRECTORY FOR acc

    wd = os.path.join(OUTPUT_FOLDER, acc)

    if not os.path.isdir(wd):
        os.mkdir(wd)

    ### READING SUPERPOSITION DATA FROM GRAPH-API. CONTAINS INFO ABOUT SEGMENTS.

    segment_data_out = os.path.join(SEGMENT_FOLDER, "{}_segments.json".format(acc)) # had to change to json because of pickle issue
    if os.path.isfile(segment_data_out):
        supp_data = pd.read_json(segment_data_out, convert_axes = False, dtype = False)
        log.info("Segment data is being read from json file")
    else:
        try:
            supp_data = pd.read_json("https://www.ebi.ac.uk/pdbe/graph-api/uniprot/superposition/{}".format(acc), convert_axes = False, dtype = False)
            log.info("Segment data is being read from API")
            supp_data.to_json(segment_data_out)
        except HTTPError as e:
            log.warning("Superposition data could not be obtained from GRAPH-API for {}. Exiting programme".format(acc))
            print("{}\t{}".format(acc, str(17)), flush = True)
            sys.exit()

    supp_data.index = range(1, len(supp_data)+1)
    segments = supp_data.index.tolist()
    n_segments = len(supp_data)
    log.info("{} presents {} different structure coverage segments".format(acc, n_segments))

    segment_data = get_segments_dict(supp_data, acc)
    segment_chains = get_segment_membership(supp_data, acc)
    segment_pdbs = {k: list(set([vv.split("_")[0] for vv in v])) for k, v in segment_chains.items()}

    ### STARTING LOOPING THROUGH ALL PROTEIN SEGMENTS

    for segment in segments:

        seg_id = "{}_{}".format(acc, str(segment))

        try:

            ### CREATES SEGMENT DIRECTORY, AND SUBDIRECTORIES

            segment_dir = os.path.join(wd, str(segment))
            raw_dir = os.path.join(segment_dir, "raw")
            clean_dir = os.path.join(segment_dir, "clean")
            trans_dir = os.path.join(segment_dir, "trans")
            simple_dir = os.path.join(segment_dir, "simple")
            fps_dir = os.path.join(segment_dir, "fingerprints")
            dssp_dir = os.path.join(segment_dir, "dssp")
            variants_dir = os.path.join(segment_dir, "variants")
            results_dir = os.path.join(segment_dir, "results")

            dirs = [
                segment_dir, raw_dir, clean_dir, trans_dir,
                simple_dir, fps_dir, dssp_dir,
                variants_dir, results_dir,
            ]
            
            for dirr in dirs:
                if not os.path.isdir(dirr):
                    os.mkdir(dirr)

            ### CHECKS IF FINAL RESULTS TABLE EXISTS, AND IF SO, SKIPS TO NEXT SEGMENT

            final_table_out = os.path.join(results_dir, "{}_{}_{}_{}_results_table.pkl".format(acc, str(segment), experimental_methods, str(resolution)))
            if os.path.isfile(final_table_out) and not override:
                print("{}\t{}".format(seg_id, str(0)), flush = True)
                log.info("Results available for Segment {} of {}".format(str(segment), acc))
                continue

            log.info("Starting to process Segment {} of {}".format(str(segment), acc))

            segment_df = matrices_df.query('index in @segment_chains[@segment]') # subsets matrices_df to select segment rows

            if len(segment_df) == 0: # happened for 8au0 of O94901. Brand new of 19/07/2023.
                log.warning("Segment {} of {} presents no chains in supp data".format(str(segment), acc))
                print("{}\t{}".format(seg_id, str(18)), flush = True)
                continue

            log.info("Segment {} of {} presents {} chains".format(str(segment), acc, str(len(segment_df))))

            pdb_ids = segment_df.pdb_id.tolist()

            pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id in all_ligs_pdbs] # filters out pdb_ids that do not present BioLiP-defined LOIs

            pdb_files = [os.path.join(pdb_db_path, pdb_id[1:3], "pdb{}.ent.gz".format(pdb_id)) for pdb_id in pdb_ids]

            #if run_experimental:

            ### GETTING EXPERIMENTAL DATA FROM ALL STRUCTURES

            experimental_out = os.path.join(results_dir, "{}_{}_strs_exp.pkl".format(acc, str(segment)))

            if override or not os.path.isfile(experimental_out):
                exp_data_df = get_experimental_data(pdb_ids, EXP_FOLDER, experimental_out)
                log.info("Obtained experimental data")
            else:
                
                exp_data_df = pd.read_pickle(experimental_out)
                log.debug("Loaded experimental data")
                pass
            log.info("Experimental data processed for Segment {} of {}".format(str(segment), acc))

            ### FILTERS OUT PDB IDS, AND FILES THAT ARE NOT FOUND IN LOCAL DATABASE

            files2remove, ids2remove = [[], []]
            for i, pdb_file in enumerate(pdb_files):
                try:
                    assert os.path.isfile(pdb_file)
                except AssertionError as e:
                    log.error("{} was not found in local database".format(pdb_ids[i]))
                    files2remove.append(pdb_file) # saving pdbs not in database so they are removed later. removing whilst for loop is not a good idea
                    ids2remove.append(pdb_ids[i])
                    
            pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id not in ids2remove]
            pdb_files = [pdb_file for pdb_file in pdb_files if pdb_file not in files2remove]

            if len(pdb_files) == 0:
                log.error("None of the structures for Segment {} of {} are present in local database".format(str(segment), acc)) #this is actually a segment EC
                print("{}\t{}".format(seg_id, str(2)), flush = True)
                continue

            ### NEW FROM 07/2023 FILTERS OUT PDB IDS THAT DO NOT MEET CRITERION: @experimental_methods AND @resolution

            log.info("{} experimental methods are being used".format(experimental_methods_list))
            log.info("Structures with resolution < {} are being used".format(str(resolution)))

            if experimental_methods == "ALL" and math.isinf(resolution):
                pass
            else:
                if "experimental_method" not in exp_data_df or "resolution" not in exp_data_df:
                    log.warning("Medhod or resolution missing. No quality structures returned.")
                    print("{}\t{}".format(seg_id, str(3)), flush = True)
                    continue

                if experimental_methods == "ALL":
                    pass
                else:
                    exp_data_df = exp_data_df.query('experimental_method in @experimental_methods_list')

                if math.isinf(resolution):
                    pass
                else:
                    exp_data_df = exp_data_df.query('resolution < @resolution')

                good_pdbs = exp_data_df.pdb_id.tolist()

                if good_pdbs == []:
                    log.warning("None of the structures meet quality threshold for Segment {} of {}".format(str(segment), acc))
                    print("{}\t{}".format(seg_id, str(4)), flush = True)
                    continue

                files2remove, ids2remove = [[], []] # now because of exp method and resolution
                for i, pdb_id in enumerate(pdb_ids):
                    if pdb_id not in good_pdbs:
                        log.warning("{} did not meet quality standards".format(pdb_id))
                        files2remove.append(pdb_files[i]) # saving pdbs not meeting QC so they are removed later. removing whilst for loop is not a good idea
                        ids2remove.append(pdb_id)

                pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id not in ids2remove] # now only desired experimental method and resolution are kept
                pdb_files = [pdb_file for pdb_file in pdb_files if pdb_file not in files2remove] # now only desired experimental method and resolution are kept

            log.info("Segment {} of {} presents {} high quality chains".format(str(segment), acc, str(len(pdb_files))))
            
            ### CHECKING AMOUNT OF PDB FILES AND PDB IDS ARE THE SAME

            try:
                assert len(pdb_ids) == len(pdb_files)
            except AssertionError as e:
                log.critical("Number of pdb ids ({}) not equal to pdb_files ({}) for Segment {} of {}".format(str(len(pdb_ids)), str(len(pdb_files)), str(segment), acc))

            segment_df = segment_df.query('pdb_id in @pdb_ids') # filtering segment dataframe, so it only includes transformation data of those tructures present in local copy of PDB
            matrices = segment_df.matrix.tolist()
            chains = segment_df.auth_asym_id.tolist()

            ### CHECKING PDB IDS AGREE WITH SEGMENT DF PDB IDS

            try:
                assert sorted(pdb_ids) == sorted(segment_df.pdb_id.tolist())
            except AssertionError as e:
                log.critical("Filtered PDBs do not agree with those from Segment dataframe for Segment {} of {}".format(str(segment), acc))
                continue

            try:
                assert segment_df.pdb_id.tolist() == [pdb_file.split("/")[-1][3:7] for pdb_file in pdb_files]
            except AssertionError as e:
                log.critical("PDB file IDs do not agree with those from Segment dataframe for Segment {} of {}".format(str(segment), acc))
                continue 

            ### OBTAINING PROTEIN-LIGAND FINGERPRINTS
            
            fps_out = os.path.join(results_dir, "{}_{}_ligs_fingerprints.pkl".format(acc, str(segment))) #fps: will stand for fingerprints. update with main_dir and so on.
            pdb_set = segment_pdbs[segment]
            log.info("There are {} unique PDBs for Segment {} of {}".format(str(len(pdb_set)), str(segment), acc))
            all_ligs_pdbs_segment = [pdb for pdb in all_ligs_pdbs if pdb in pdb_set] # filtering pdbs so only data about segment is retrieved
            log.info("There are {} unique ligand-binding PDBs for Segment {} of {}".format(str(len(all_ligs_pdbs_segment)), str(segment), acc))

            if override or not os.path.isfile(fps_out):
                if all_ligs_pdbs_segment == []:
                    print("{}\t{}".format(seg_id, str(5)), flush = True)
                    log.warning("Segment {} of {} does not present any ligand-binding structures".format(str(segment), acc))
                    continue
                else:
                    lig_fps = get_fingerprints_dict(acc, fps_dir, fps_out, all_ligs_pdbs_segment, segment_chains[segment], MOLS_FOLDER, INTERS_FOLDER) # GETS ALL LIGAND FINGERPRINTS FROM LIG-BOUND CONTAINING PDBS
            else:
                with open(fps_out, "rb") as f:
                    lig_fps = pickle.load(f) # stands for ligand fingerprints

            ### CHECKING THAT THERE ARE FINGERPRINTS. THERE SHOULD ALWAYS BE AT THIS POINT.

            if lig_fps == {}: # if there are not any segment fingerprints (no ligands bound)
                print("{}\t{}".format(seg_id, str(6)), flush = True)
                log.warning("ACHTUNG! No fingerprints found for Segment {} of {}".format(str(segment), acc))
                continue

            ### FILTER OUT NON-SEGMENT LIGAND INTERACTIONS

            n_acc_chain_inters = sum([len(inters) for ligs in lig_fps.values() for inters in ligs.values()]) # number of ligands interactions with target UniProt Accession chains for all structures in a segment
            
            if n_acc_chain_inters == 0: #there are ligands bound, but not to the target protein chains, across all ligands bound in all structures for a segment
                print("{}\t{}".format(seg_id, str(7)), flush = True)
                log.warning("None of the ligands in Segment {} of {} interact with any of the target chains".format(str(segment), acc))
                continue

            ### FILTER OUT NON-LOIS 

            #lig_fps_filt1 = filter_non_relevant_ligs(lig_fps, relevant_ligs) # FILTERS OUT NON-LOIs INTERACTIONS. THIS IS WHERE WE WOULD CHECK FOR LIGAND SIZE AS WELL.

            lig_fps_filt1 = lig_fps # avoiding filtering, as we are already using BioLiP LOIs

            n_relevant_ligs = sum([len(rel_ligs) for rel_ligs in lig_fps_filt1.values()])

            if n_relevant_ligs == 0: #there are ligands bound, but they are all non-relevant, across all ligands bound in all structures for a segment
                print("{}\t{}".format(seg_id, str(8)), flush = True)
                log.warning("All ligands in Segment {} of {} are not relevant".format(str(segment), acc))
                continue

            ### FILTER OUT NON-PROTEIN INTERACTIONS

            lig_fps_filt2 = filter_non_protein_inters(lig_fps_filt1, acc, pdb_resnames) # FILTERS OUT NON-PROTEIN INTERACTIONS

            n_protein_inters = sum([len(inters) for rel_ligs in lig_fps_filt2.values() for inters in rel_ligs.values()]) # number of protein-RELEVANT ligands interactions for all structures in a segment

            if n_protein_inters == 0: #there are ligands bound, but they are all non-relevant, across all ligands bound in all structures for a segment
                print("{}\t{}".format(seg_id, str(9)), flush = True)
                log.warning("None of the ligand interactions in Segment {} of {} involve target protein atoms".format(str(segment), acc))
                continue

            ### TODO: IMPLEMENT ONLY-SIDECHAIN INTERACTIONS ###
            ###                                             ###
            ###                                             ###
            ### TODO: IMPLEMENT ONLY-SIDECHAIN INTERACTIONS ###

            log.info("Ligand fingerprints obtained for Segment {} of {}".format(str(segment), acc))

            ### GETTING PDB-UNIPROT SIFTS MAPPING

            sifts_out = os.path.join(results_dir, "{}_{}_strs_sifts.pkl".format(acc, str(segment)))

            if override or not os.path.isfile(sifts_out):
                sifts_mapping = {}
                for pdb_id in all_ligs_pdbs_segment: # before it was just pdb_ids, so no information would be retrieved of those with no structure in local PDB copy
                    sifts_mapping[pdb_id] = get_mapping_from_sifts(pdb_id)
                with open(sifts_out, "wb") as f:
                    pickle.dump(sifts_mapping, f)
                log.info("Obtained SIFTS mappings table")
            else:
                with open(sifts_out, "rb") as f:
                    sifts_mapping = pickle.load(f)
                log.info("Loaded SIFTS mappings table")

            lig_fps_filt2 = {k: v for k, v in lig_fps_filt2.items() if sifts_mapping[k] != {}} #filter out those fingerprints form structures where SIFTS mapping was not retrieved

            lig_fps_filt2_sifted, lig_fps_filt2_sifted_v2 = get_up_mapping_from_prointvar(lig_fps_filt2, sifts_mapping)

            ### FILTERING OUT FINGERPRINTS THAT ARE EMPTY DUE TO LACK OF SIFTS MAPPING ###

            for k1, v1 in lig_fps_filt2_sifted_v2.items():
                if v1 == []:
                    log.warning("{} resulted in an empty fingerprint".format(k1))
                else:
                    pass
            lig_fps_filt2_sifted_v2 = {k: v for k, v in lig_fps_filt2_sifted_v2.items() if v != []} # removes empty fingerprints from dict

            lig_fps_filt2_sifted = {k1: {k2: v2 for k2, v2 in v1.items() if v2 != []} for k1, v1 in lig_fps_filt2_sifted.items()} # removes empty fingerprints from dict
            lig_fps_filt2_sifted = {k: v for k, v in lig_fps_filt2_sifted.items() if v != {}} # removes pdb entries from dict if no fingerprints remain in dict
            
            # it could be that none of the fingerprints have SIFTS mappings, and therefore dictionaries are empty

            if lig_fps_filt2_sifted == {}:
                print("{}\t{}".format(seg_id, str(10)), flush = True)
                log.warning("No relevant fingerprints present SIFTS mapping for Segment {} of {}".format(str(segment), acc))
                continue

            ### DONE ###

            log.info("PDB-UniProt mappings performed for Segment {} of {}".format(str(segment), acc))

            ### REMOVING REDUNDANT LIGANDS: LIGANDS WITH SAME LIGAND ID AND SAME FINGERPRINT. possibly not needed.

            ### TODO ###

            ### CLUSTERING LIGANDS INTO BINDING SITES

            lig_sifted_inters = get_inters(lig_fps_filt2_sifted)
            lig_sifted_inters = [sorted(list(set(i))) for i in lig_sifted_inters] # making sure each residue is present only once (O00214 problematic with saccharides)

            lig_labs = get_labs(lig_fps_filt2_sifted)
            n_ligs = len(lig_labs)
            log.info("There are {} relevant ligands for Segment {} of {}".format(str(n_ligs), str(segment), acc))
            irel_mat_out = os.path.join(results_dir, "{}_{}_{}_{}_irel_matrix.pkl".format(acc, str(segment), experimental_methods, str(resolution)))

            if override or not os.path.isfile(irel_mat_out):
                irel_matrix = get_intersect_rel_matrix(lig_sifted_inters) # this is a measure of similarity, probs want to save this
                dump_pickle(irel_matrix, irel_mat_out)
                log.info("Calcualted intersection matrix")
            else:
                load_pickle(irel_mat_out)
                log.info("Loaded intersection matrix")
            if n_ligs == 1:
                cluster_ids = [0]
            else:
                irel_df = pd.DataFrame(irel_matrix)
                dist_df = 1 - irel_df # distance matrix in pd.Dataframe() format
                condensed_dist_mat = scipy.spatial.distance.squareform(dist_df) # condensed distance matrix to be used for clustering
                linkage = scipy.cluster.hierarchy.linkage(condensed_dist_mat, method = lig_clust_method, optimal_ordering = True)
                cut_tree = scipy.cluster.hierarchy.cut_tree(linkage, height = lig_clust_dist)
                cluster_ids = [int(cut) for cut in cut_tree]
            cluster_id_dict = {lig_labs[i]: cluster_ids[i] for i in range(len(lig_labs))} #dictionary indicating membership for each lig
            
            log.info("Ligand clustering realised for Segment {} of {}".format(str(segment), acc))

            ### TREE VISUALISATION AND COLOURING SECTION ###

            bs_colors = list(itertools.islice(rgbs(), 300)) # new_colours

            ### TODO ###

            #if run_transform:

            ### TRANSFORMATION OF PROTEIN-LIGAND INTERACTIONS CONTAINING PDBs

            transform_all_files(pdb_files, matrices, chains, raw_dir, clean_dir, trans_dir) # this is wrong. pdb_files and chains do not match

            log.info("Structures cleaned and transformed for Segment {} of {}".format(str(segment), acc))

            ### SIMPLIFYING PDB FILES (1 MODEL PROTEIN COORDINATES + HETATM FOR THE REST OF THEM)

            get_simple_pdbs(trans_dir, simple_dir) # right now, does not print ligands if they are actual amino acids. Could fix passing ligand data for each structure. fingerprints dict

            log.info("Structures simplified for Segment {} of {}".format(str(segment), acc))

            ### CHIMERA COLOURING SCRIPT AND ATTRIBUTE WRITING

            chimera_atom_specs = get_chimera_data(cluster_id_dict)

            attr_out = os.path.join(results_dir, "{}_{}_{}_{}_pdbe_kb_scipy_{}_{}.attr".format(acc, str(segment), experimental_methods, str(resolution), lig_clust_method, lig_clust_dist))
            
            if override or not os.path.isfile(attr_out):
                
                write_chimera_attr(attr_out, chimera_atom_specs, cluster_ids)

            chimera_script_out = os.path.join(results_dir, "{}_{}_{}_{}_pdbe_kb_scipy_{}_{}.com".format(acc, str(segment), experimental_methods, str(resolution), lig_clust_method, lig_clust_dist))

            ### IMPLEMENT CHIMERA OPENING SCRIPT: opens only those PDBs that are actually binding ligands. could be less than 50% of total chains

            if override or not os.path.isfile(chimera_script_out):

                write_chimera_command(chimera_script_out, chimera_cmd_args, cluster_ids, bs_colors)

            log.info("Chimera attributes and script generated for Segment {} of {}".format(str(segment), acc))

            ### BINDING SITE MEMBERSHIP PROCESSING

            membership_out = os.path.join(results_dir, "{}_{}_{}_{}_bss_membership.pkl".format(acc, str(segment), experimental_methods, str(resolution)))
            cluster_ress_out = os.path.join(results_dir, "{}_{}_{}_{}_bss_ress.pkl".format(acc, str(segment), experimental_methods, str(resolution)))
            bs_mm_dict_out = os.path.join(results_dir, "{}_{}_{}_{}_ress_bs_membership.pkl".format(acc, str(segment), experimental_methods, str(resolution)))
                
            if override or not os.path.isfile(membership_out):
                membership = get_cluster_membership(cluster_id_dict) # which LBS ligands belong to
                dump_pickle(membership, membership_out)
                log.info("Calculated binding site membership")
            else:
                membership = load_pickle(membership_out)
                log.debug("Loaded binding site membership")

            if override or not os.path.isfile(cluster_ress_out):
                cluster_ress = get_all_cluster_ress(membership, lig_fps_filt2_sifted_v2) # residues that form each LBS 
                dump_pickle(cluster_ress, cluster_ress_out) 
                log.info("Calculated binding site composition") 
            else:
                cluster_ress = load_pickle(cluster_ress_out)
                log.info("Loaded binding site composition")

            if override or not os.path.isfile(bs_mm_dict_out):
                bs_ress_membership_dict = get_residue_bs_membership(cluster_ress)
                log.info("Calcualted residue membership")
                dump_pickle(bs_ress_membership_dict, bs_mm_dict_out)  
            else:
                bs_ress_membership_dict = load_pickle(bs_mm_dict_out)
                log.debug("Loaded residue membership")

            ### RUNNING DSSP FOR ALL STRUCTURES

            master_dssp_out = os.path.join(results_dir, "{}_{}_{}_{}_strs_dssp.pkl".format(acc, str(segment), experimental_methods, str(resolution)))
            if override or not os.path.isfile(master_dssp_out):
                dssp_data = get_dssp_data(trans_dir, dssp_dir, sifts_mapping, master_dssp_out)
                log.info("Obtained DSSP data")
            else:
                dssp_data = pd.read_pickle(master_dssp_out)
                log.info("Loaded DSSP data")
            if dssp_data.empty:
                log.warning("There is no DSSP data for Segment {} of {}".format(str(segment), acc))
            else:
                dsspd_filt = dssp_data.query('UniProt_ResNum == UniProt_ResNum and AA != "X" and RSA == RSA').copy()
                dsspd_filt.SS = dsspd_filt.SS.fillna("C")
                dsspd_filt.SS = dsspd_filt.SS.replace("", "C")
                dsspd_filt.UniProt_ResNum = dsspd_filt.UniProt_ResNum.astype(int)

                AA_dict_out = os.path.join(results_dir, "{}_{}_{}_{}_ress_AA.pkl".format(acc, str(segment), experimental_methods, str(resolution)))
                RSA_dict_out = os.path.join(results_dir, "{}_{}_{}_{}_ress_RSA.pkl".format(acc, str(segment), experimental_methods, str(resolution)))
                SS_dict_out = os.path.join(results_dir, "{}_{}_{}_{}_ress_SS.pkl".format(acc, str(segment), experimental_methods, str(resolution)))
                rsa_profs_out = os.path.join(results_dir, "{}_{}_{}_{}_bss_RSA_profiles.pkl".format(acc, str(segment), experimental_methods, str(resolution)))

                if override or not os.path.isfile(AA_dict_out):
                    ress_AA_dict = {
                        up_resnum: dsspd_filt.query('UniProt_ResNum == @up_resnum').AA.mode()[0] # gets dict per UP residue and meore frequent AA.
                        for up_resnum in dsspd_filt.UniProt_ResNum.unique().tolist()
                    }   
                    dump_pickle(ress_AA_dict, AA_dict_out)
                    log.info("Calculated residue AA dictionary")
                else:
                    ress_AA_dict = load_pickle(AA_dict_out)
                    log.debug("Loaded residue AA dictionary")

                if override or not os.path.isfile(RSA_dict_out):
                    ress_RSA_dict = {
                        up_resnum: round(dsspd_filt.query('UniProt_ResNum == @up_resnum').RSA.mean(), 2) # gets dict per UP residue and mean RSA.
                        for up_resnum in dsspd_filt.UniProt_ResNum.unique().tolist()
                    }
                    dump_pickle(ress_RSA_dict, RSA_dict_out)  
                    log.info("Calculated residue RSA dictionary")
                else:
                    ress_RSA_dict = load_pickle(RSA_dict_out)
                    log.debug("Loaded residue RSA dictionary")
                
                if override or not os.path.isfile(SS_dict_out):
                    ress_SS_dict = {
                        up_resnum: dsspd_filt.query('UniProt_ResNum == @up_resnum').SS.mode()[0] # gets dict per UP residue and more frequent SS.
                        for up_resnum in dsspd_filt.UniProt_ResNum.unique().tolist()
                    }
                    dump_pickle(ress_SS_dict, SS_dict_out)
                    log.info("Calculated residue SS dictionary")
                else:
                    ress_SS_dict = load_pickle(SS_dict_out)
                    log.debug("Loaded residue SS dictionary")

                if override or not os.path.isfile(rsa_profs_out):
                    rsa_profiles = {}
                    for k, v in cluster_ress.items():
                        rsa_profiles[k] = []
                        for v2 in v:
                            if v2 in ress_RSA_dict:
                                rsa_profiles[k].append(ress_RSA_dict[v2])
                            else:
                                log.warning("Cannot find RSA data for UP residue {} in Segment {} of {}".format(str(v2), str(segment), acc))
                    dump_pickle(rsa_profiles, rsa_profs_out)
                    log.info("Calculated bigand binding site RSA profiles")
                else:
                    rsa_profiles = load_pickle(rsa_profs_out)
                    log.debug("Loaded ligand binding site RSA profiles")

                log.info("DSSP data processed for Segment {} of {}".format(str(segment), acc))

            ### GENERATE ALIGNMENT                                                 

            #if run_variants:

            seq_out = os.path.join(variants_dir, "{}_{}.fasta".format(acc, str(segment)))

            if override_variants or not os.path.isfile(seq_out):
                best = get_best_from_segment_data(segment_data[segment])
                best_seq_id = get_best_struct_seq(acc, segment, seq_out, best) # change how we get seq. must be representative of the segment
                log.info("Generated sequence file")
            else:
                best_seq_id = [rec.id for rec in Bio.SeqIO.parse(seq_out, "fasta")][0]
                log.debug("Sequence file already existed")
                pass

            hits_out = os.path.join(variants_dir, "{}_{}.out".format(acc, str(segment)))
            hits_aln = os.path.join(variants_dir, "{}_{}.sto".format(acc, str(segment)))

            if override_variants or not os.path.isfile(hits_out) or not os.path.isfile(hits_aln):
                jackhmmer(seq_out, hits_out, hits_aln, n_it = jackhmmer_n_it, seqdb = swissprot)
                log.info("Generated MSA")
            else:
                log.debug("MSA already existed")
                pass
            
            n_seqs = len(AlignIO.read(hits_aln, MSA_fmt))

            if n_seqs == 1: # need RESULTS TABLE even if there is no MSA data
                print("{}\t{}".format(seg_id, str(11)), flush = True)
                log.critical("No sequences were found by jackHMMER for Segment {} of {}. Finishing here".format(str(segment), acc))
                continue
            else:
                log.info("{} sequences were found by jackHMMER for Segment {} of {}".format(str(n_seqs), str(segment), acc))
            hits_aln_rf = os.path.join(variants_dir, "{}_{}_rf.sto".format(acc, str(segment)))

            if override_variants or not os.path.isfile(hits_aln_rf):
                add_acc2msa(hits_aln, hits_aln_rf, best_seq_id)
                log.info("Formatted MSA")
            else:
                log.debug("MSA already formatted")
                pass

            log.info("MSA realised for Segment {} of {}".format(str(segment), acc))

            ### CONSERVATION ANALYSIS

            prot_cols = prot_cols = get_target_prot_cols(hits_aln, best_seq_id)
            shenkin_out = os.path.join(variants_dir, "{}_{}_rf_shenkin.pkl".format(acc, str(segment)))
            if override_variants or not os.path.isfile(shenkin_out):
                shenkin = calculate_shenkin(hits_aln_rf, "stockholm", shenkin_out)
                log.info("Calculated conservation data")
            else:
                shenkin = pd.read_pickle(shenkin_out)
                log.debug("Loaded conservation data")
            
            shenkin_filt_out = os.path.join(variants_dir, "{}_{}_rf_shenkin_filt.pkl".format(acc, str(segment)))
            if override_variants or not os.path.isfile(shenkin_filt_out):
                shenkin_filt = format_shenkin(shenkin, prot_cols, shenkin_filt_out)
                log.info("Filtered conservation data")
            else:
                shenkin_filt = pd.read_pickle(shenkin_filt_out)
                log.debug("Conservation data already filtered")

            log.info("Conservation scores calculated for Segment {} of {}".format(str(segment), acc))

            ### VARIATION, POSSIBLY NEED TO IMPLEMENT CLINVAR AS WELL

            aln_obj = Bio.AlignIO.read(hits_aln_rf, "stockholm") #crashes if target protein is not human!
            aln_info_path = os.path.join(variants_dir, "{}_{}_rf_info_table.p.gz".format(acc, str(segment)))
            if override_variants or not os.path.isfile(aln_info_path):
                aln_info = varalign.alignments.alignment_info_table(aln_obj)
                aln_info.to_pickle(aln_info_path)
                log.info("Generated MSA info table")
            else:
                aln_info = pd.read_pickle(aln_info_path)
                log.debug("Loaded MSA info table")
            
            log.info("There are {} sequences in MSA for Segment {}".format(len(aln_info), str(segment)))

            indexed_mapping_path = os.path.join(variants_dir, "{}_{}_rf_mappings.p.gz".format(acc, str(segment)))
            if override_variants or not os.path.isfile(indexed_mapping_path):
                indexed_mapping_table = varalign.align_variants._mapping_table(aln_info) # now contains all species
                indexed_mapping_table.to_pickle(indexed_mapping_path) # important for merging later on
                log.info("Generated MSA mapping table")
            else:
                indexed_mapping_table = pd.read_pickle(indexed_mapping_path)
                log.debug("Loaded MSA mapping table")    

            aln_info_human = aln_info[aln_info.species == "HUMAN"]

            if len(aln_info_human) > 0:
                log.info("There are {} HUMAN sequences in the MSA for Segment {} of {}".format(len(aln_info_human), str(segment), acc))
            
                human_hits_msa = os.path.join(variants_dir, "{}_{}_rf_human.sto".format(acc, str(segment)))
                
                if override_variants or not os.path.isfile(human_hits_msa):
                    get_human_subset_msa(hits_aln_rf, human_hits_msa)
                else:
                    pass

                ### copy ensemble SQLite to directory where this is being executed
                cp_path = cp_sqlite(wd)
                log.debug("ENSEMBL_CACHE SQLite copied correctly")

                variant_table_path = os.path.join(variants_dir, "{}_{}_rf_human_variants.p.gz".format(acc, str(segment)))
                if override_variants or not os.path.isfile(variant_table_path):
                    #vcf_out_path = os.path.join(main_dir , "results", "{}_alignment_variants.vcf".format(seg_id))
                    try:
                        variants_table = varalign.align_variants.align_variants(aln_info_human, path_to_vcf = gnomad_vcf,  include_other_info = False, write_vcf_out = False)     
                    except ValueError as e:
                        print("{}\t{}".format(seg_id, str(12)), flush = True)
                        variants_table = pd.DataFrame()
                        log.warning("No variants were retrieved for Segment {} of {}".format(str(segment), acc))

                    variants_table.to_pickle(variant_table_path)

                else:
                    variants_table = pd.read_pickle(variant_table_path)

                ### remove ensembl SQLite from directory where this is being executed
                rm_sqlite(cp_path)
                log.debug("ENSEMBL_CACHE SQLite removed correctly")

                if variants_table.empty: # variant table is empty. E.g., P03915. Only 3 human sequences. They are all mitochondrial (not in gnomAD)
                    pass

                else:
                    # in order to be able to read the vcf and parse the DB, the ensemble.cache.sqlite file must be in the ./.varalign directory

                    human_miss_vars = format_variant_table(variants_table, prot_cols) # GET ONLY MISSENSE VARIANTS ROWS
                    human_miss_vars_msa_out = os.path.join(variants_dir, "{}_{}_rf_human_missense_variants_seqs.sto".format(acc, str(segment)))

                    miss_df_out = os.path.join(results_dir, "{}_{}_missense_df.pkl".format(acc, str(segment)))
                    
                    if override or not os.path.isfile(miss_df_out): # we leave it as override and not override_variants to fix the wrong pseudocounts
                        missense_variants_df = get_missense_df(
                            hits_aln_rf, human_miss_vars,
                            shenkin_filt, prot_cols, human_miss_vars_msa_out
                        )

                        if missense_variants_df.empty:
                            print("{}\t{}".format(seg_id, str(13)), flush = True)
                            log.warning("No missense variants found for MSA of Segment {} of {}".format(str(segment), acc))
                            pass

                        else:
                            missense_variants_df = add_miss_class(
                                missense_variants_df, miss_df_out,
                                cons_col = "abs_norm_shenkin",
                            )
                            log.info("Calculated missense dataframe")
                    else:
                        missense_variants_df = pd.read_pickle(miss_df_out)
                        log.debug("Loaded missense dataframe")

                    if missense_variants_df.empty:
                        pass
                        #log.info("No missense variants found for MSA of Segment {} of {}".format(str(segment), acc))
                            
                    else:
                        # ADDS COLUMNS FROM MISSENSE DF TO SHENKIN FILT DF, CONSERVATION AND VARIATION DATA ABOUT HUMAN VARIANT SUB MSA
                        shenkin_filt.loc[:, "human_shenkin"] = missense_variants_df.shenkin
                        shenkin_filt.loc[:, "human_occ"] = missense_variants_df.occ
                        shenkin_filt.loc[:, "human_gaps"] = missense_variants_df.gaps
                        shenkin_filt.loc[:, "human_occ_pct"] = missense_variants_df.occ_pct
                        shenkin_filt.loc[:, "human_gaps_pct"] = missense_variants_df.gaps_pct
                        shenkin_filt.loc[:, "variants"] = missense_variants_df.variants
                        shenkin_filt.loc[:, "oddsratio"] = missense_variants_df.oddsratio
                        shenkin_filt.loc[:, "pvalue"] = missense_variants_df.pvalue
                        shenkin_filt.loc[:, "se_OR"] = missense_variants_df.se_OR

            else:
                print("{}\t{}".format(seg_id, str(14)), flush = True)
                log.warning("No human sequences for Segment {} of {}".format(str(segment), acc))
                pass

            shenkin_mapped_out = os.path.join(results_dir, "{}_{}_ress_consvar.pkl".format(acc, str(segment)))
            if override or not os.path.isfile(shenkin_mapped_out): # we leave it as override and not override_variants to fix the wrong pseudocounts
                aln_ids = list(set([seqid[0] for seqid in indexed_mapping_table.index.tolist() if acc in seqid[0]])) # THIS IS EMPTY IF QUERY SEQUENCE IS NOT FOUND
                n_aln_ids = len(aln_ids)
                if n_aln_ids != 1:
                    log.warning("There are {} sequences matching accession for Segment {} in {}".format(str(n_aln_ids), str(segment), acc))
                mapped_data = merge_shenkin_df_and_mapping(shenkin_filt, indexed_mapping_table, aln_ids)
                mapped_data.to_pickle(shenkin_mapped_out)
            else:
                mapped_data = pd.read_pickle(shenkin_mapped_out)
            log.info("Conservation + variant data obtained for Segment {} of {}".format(str(segment), acc))

            ### GENERATE SUMMARY TABLES

            if not dssp_data.empty:
                mapped_data["AA"] = mapped_data.UniProt_ResNum.map(ress_AA_dict)
                mapped_data["RSA"] = mapped_data.UniProt_ResNum.map(ress_RSA_dict)
                mapped_data["SS"] = mapped_data.UniProt_ResNum.map(ress_SS_dict)
            else:
                log.warning("Results table will not contain DSSP columns for Segment {} of {}".format(str(segment), acc))
                pass

            mapped_data["binding_sites"] = mapped_data.UniProt_ResNum.map(bs_ress_membership_dict)
            mapped_data.to_pickle(final_table_out)

            log.info("Results available for Segment {} of {}".format(str(segment), acc))
            #else:
            #    log.info("Not running variants for Segment {} of {}".format(str(segment), acc))

            print("{}\t{}".format(seg_id, str(0)), flush = True)

        except Exception as e:
            print("{}\t{}".format(seg_id, str(1)), flush = True)
            log.error("Segment {} of {} failed due to {}".format(str(segment), acc, e))
            raise

    log.info("THE END")



# /cluster/gjb_lab/2394007/LIGYSIS_PDB

# python3.6 ./../../ligysis.py --transform --experimental --variants --override O55234

