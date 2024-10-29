# Installation

These are the complete full instructions to install **LIGYSIS** from scratch.

## Installation of DSSP

DSSP is incompatible with all other environments, and so must go on its environment of its own. You can install locally as well, all we will need is the location of the executable. The version of the libboost library must be specified to be this one, otherwise, dssp will not run.

```
conda create -n DSSP salilab::dssp=3.0.0 libboost=1.73.0
```

## Installation of HMMER

Installing HMMER on its own environment. Just easier...

```
conda create -n HMMER hmmer=3.4
```

## Installation of STAMP

The following instructions are to install STAMP. For more information refer to the [STAMP installation instructions](https://www.compbio.dundee.ac.uk/downloads/stamp/INSTALL).

```
# download STAMP
curl -O https://www.compbio.dundee.ac.uk/downloads/stamp/stamp.4.4.2.tar.gz

# decompress STAMP
tar -xzvf stamp.4.4.2.tar.gz

# change directory
cd stamp.4.4.2
```
To install STAMP, run the BUILD script in this directory using:
```
# building STAMP
./BUILD <system-type>
```
where \<system-type\> is one of:

- linux
- osx 
- dec
- sgi
- sun

The executables will be installed in bin/\<system-type\>/.

For more information refer to the [STAMP manual](https://www.compbio.dundee.ac.uk/manuals/stamp.4.4/stamp.html)

## Installation of LIGYSIS

The first step to install **LIGYSIS** is to Git Clone the repository.

```
# git clone LIGYSIS from repository
git clone -b revamped https://github.com/JavierSanchez-Utges/ligysis_custom.git
```

### Installation of pdbe-arpeggio

The following commands will crate an environment and install pdbe-arpeggio.

```
# install pdbe-arpeggio environment
conda create -n ARPEGGIO python=3.9 gemmi openbabel biopython -c conda-forge

# activating pdbe-arpeggio environment
conda activate ARPEGGIO

# install pdbe-arpeggio
pip install pdbe-arpeggio

# test pdbe-arpeggio with help function
pdbe-arpeggio -h
```

### Installation of DEEP_LEARNING environment

```
# change directory to environments directory
cd ligysis_custom/ENVS

# install DEEP_LEARNING environment
conda env create -f DEEP_LEARNING.yml
```

### Installation of CLEAN_PDB environment

```
# install CLEAN_PDB environment
conda create -n CLEAN_PDB python=2.7.15 biopython=1.74
```

## Installation of VarAlign

The following instructions are to install VarAlign. Fore more information refer to the [VarAlign repository](https://github.com/bartongroup/SM_VarAlign/tree/JSU_branch).

```
# install LIGYSIS environment
conda env create -f LIGYSIS.yml

# VarAlign installation (Grabbed from URL)

# change directory to main working directory
cd ../..

# git clone specific branch of VarAlign from repository
git clone -b JSU_branch https://github.com/bartongroup/SM_VarAlign.git

# change directory to VarAlign directory
cd SM_VarAlign

# activate varalign_env environment
conda activate LIGYSIS

# install VarAlign
pip install .
```

## Installation of ProIntVar

The following instructions are to install ProIntVar. Fore more information refer to the [ProIntVar repository](https://github.com/bartongroup/ProIntVar/tree/JSU_branch).

```
# ProIntVar installation (Grabbed from URL)

# change directory to main working directory
cd ..

# git clone specific branch of ProIntVar from repository
git clone -b JSU_branch https://github.com/bartongroup/ProIntVar.git

### change directory to ProIntVar directory
cd ProIntVar

# pip install ProIntVar dependencies
pip install -r requirements.txt

#then
python setup.py install
```

## Configuration of ProIntVar

ProIntVar needs to be configured in order to run DSSP. The path to the DSSP binary must be added.

```
# change directory to VarAlign directory
cd ../SM_VarAlign

# set up ProIntVar configuration
ProIntVar-config-setup new_config.ini

# edit the following values in new_config.ini
dssp_bin = /path/to/anaconda/envs/bin/mkdssp

# Update the settings according to user preferences and push them
ProIntVar-config-load new_config.ini
```

## Configuration of LIGYSIS

The configuration file can be found [here](ligysis_config.txt). It includes the necessary paths to databases and binaries needed to run the pipeline. It looks like this:

```
### LIGYSIS CONFIG FILE ###

[paths]

## BINARIES

stamp_bin = /path/to/stamp.4.4.2/bin/linux/stamp
transform_bin = /path/to/stamp.4.4.2/bin/linux/transform
clean_pdb_python_bin = /path/to/miniconda/envs/CLEAN_PDB/bin/python
clean_pdb_bin = ./clean_pdb.py
arpeggio_python_bin =/path/to/miniconda/envs/ARPEGGIO/bin/python
arpeggio_bin = /path/to/miniconda/envs/ARPEGGIO/bin/pdbe-arpeggio

## DATABASES

ensembl_sqlite =/path/to/.varalign/ensembl_cache.sqlite              ### WHAT DO WE DO ABOUT THIS?
gnomad_vcf =/path/to/gnomad/gnomad.exomes.r2.0.1.sites.vcf.gz
swissprot = /path/to/swissprot.fasta

## DIRECTORIES

stampdir = /path/to/stamp.4.4.2/defs/

### END OF CONFIG FILE ###
```

These are mock paths and need to be replaced with the right ones for the pipeline to run successfully.

## Patching pdbe-arpeggio

It is necessary to apply patch to pdbe-arpeggio. pdbe-arpeggio crashes when CIF files do not present the `_chem_comp.` fields. See GitHub issue [here](https://github.com/PDBeurope/arpeggio/issues/20). This will be the case for all the preferred biological assemblies in the PDBe. [This patch](OTHER/JSU_patch.diff) is a quick solve for this until the issue is appropriately approached. It simply comments the lines that extract the `_chem_comp.` information, which is not needed in LIGYSIS.

This is how the patch is applied:

```sh
# change directory to where pdbe-arpeggio is installed
cd /path/to/environments/ARPEGGIO/lib/python3.9/site-packages/

# apply patch
git apply /path/to/ligysis_custom/OTHER/JSU_patch.diff
```

This should have commented these 7 lines of code in the `arpeggio/core/interactions.py` file.

```
+        # self.component_types = protein_reader.get_component_types(filename)

+        # result_entry['bgn']['label_comp_type'] = self.component_types[utils.get_residue_name(contact.bgn_atom)]

+        # result_entry['end']['label_comp_type'] = self.component_types[utils.get_residue_name(contact.end_atom)]

+        # result_entry['bgn']['label_comp_type'] = self.component_types[utils.get_residue_name(contact.bgn_res)]

+        # result_entry['end']['label_comp_type'] = self.component_types[utils.get_residue_name(contact.end_res)]

+        # result_entry['bgn']['label_comp_type'] = self.component_types[utils.get_residue_name(contact.bgn_atom)]

+        # result_entry['end']['label_comp_type'] = self.component_types[utils.get_residue_name(contact.end_res)]
```
