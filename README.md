# LIGYSIS

This repository contains the original version of our original ligand site analysis **LIGYSIS** pipeline, which was employed to analyse all biologically relevant protien-ligand interactions on the PDBe, which results are served in the **LIGYSIS** [web server](https://www.compbio.dundee.ac.uk/ligysis/). The code for the web server can be found [here](https://github.com/JavierSanchez-Utges/LIGYSIS-web).

**LIGYSIS** makes extensive use of the [PDBe-KB](https://www.ebi.ac.uk/pdbe/pdbe-kb/) and [PDBe](https://www.ebi.ac.uk/pdbe/) [APIs](https://www.ebi.ac.uk/pdbe/pdbe-rest-api) to retrieve transformation matrices, all PDB IDs mapping to a given UniProt ID and more. LIGYSIS analyses protein-ligand contacts from BioLiP defined biologically relevant ligands from [PISA](https://www.ebi.ac.uk/pdbe/pisa/)-defined preferred biological assemblies. These protein-ligand interactions are used to cluster ligands together into binding sites, which are then characterised by using evolutionary divergence, human genetic variation and structural features such as [Relative Solvent Accessibility](https://en.wikipedia.org/wiki/Relative_accessible_surface_area) (RSA) or secondary structure.

## Pipeline methodology

The pipeline can be summarised in the following steps:

1. Transformation matrices are retrieved for each chain mapping to a [UniProt](https://www.uniprot.org/) ID from the [PDBe-KB FTP site](http://ftp.ebi.ac.uk/pub/databases/pdbe-kb/superposition/P/P00517/P00517.json).
2. Segment superposition data, i.e., segment-clustered protein chains and representative chains, are retrieved from the PDBe GRAPH API [superposition endpoint](https://www.ebi.ac.uk/pdbe/graph-api/uniprot/superposition/P00517)

   **Note:** The following steps are repeated in a loop now for each Segment within the protein.
3. Retrieve experimental data for each structure from the PDBe REST API [experiment endpoint](https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/1svh)
4. Filter out structures according to resolution (Current dataset considers <i>ALL</i> structures).
5. Download of biological assemblies from PDBe via [ProIntVar](https://github.com/bartongroup/ProIntVar).
6. Protein-ligand interactions calculation with [pdbe-rpeggio](https://github.com/PDBeurope/arpeggio).
7. Mapping of PDB residues to UniProt through [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/), encoded in the [<i>mmCIF</i> format](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdbx-mmcif).
8. Ligand clustering into binding sites using [SciPy](https://scipy.org/).
9. Transformation of all chains employing PDBe-KB matrices with [BioPython](https://biopython.org/) (`atom.transform(rotation, translation)`).
10. <i>Simplification</i> of superposed chains. This consists in keeping protein atoms only for one of the superposed structures, and heteroatoms for the rest. This is done to generate a lower-weight superposition (all ligands to a single protein scaffold).
11. Generation of [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) superposition visualisation script.
12. [Relative solvent accessibility](https://en.wikipedia.org/wiki/Relative_accessible_surface_area) (RSA) and secondary structure elements calculation with [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) via [ProIntVar](https://github.com/bartongroup/prointvar).
13. Multiple sequence alignment with [jackhmmer](http://hmmer.org/).
14. Shenkin amino acid divergence score calculation [[1](https://doi.org/10.1002/prot.340110408), [2](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335)].
15. Missense enrichment score calculation with [VarAlign](https://github.com/bartongroup/SM_varalign) [[3](https://www.biorxiv.org/content/10.1101/127050v2), [4](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3783), [5](https://www.nature.com/articles/s42003-024-06117-5)].
16. RSA-based clustering label and functional score calculation [[6](https://www.nature.com/articles/s42003-024-05970-8)].

## Dependencies

Third party dependencies for these notebooks include:
- [pdbe-arpeggio](https://github.com/PDBeurope/arpeggio) [(GNU GPL v3.0 License)](https://github.com/harryjubb/arpeggio/blob/master/LICENSE)
- [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) [(Boost Software License)](https://swift.cmbi.umcn.nl/gv/dssp/)
- [Hmmer](http://hmmer.org/) [(BSD-3 Clause License)](http://eddylab.org/software/hmmer/Userguide.pdf)
- [ProIntVar](https://github.com/bartongroup/prointvar) [(MIT License)](https://github.com/bartongroup/ProIntVar/blob/master/LICENSE.md)
- [ProteoFAV](https://github.com/bartongroup/ProteoFAV) [(MIT License)](https://github.com/bartongroup/ProteoFAV/blob/master/LICENSE.md)
- [VarAlign](https://github.com/bartongroup/SM_varalign) [(MIT License)](https://github.com/bartongroup/SM_VarAlign/blob/master/LICENSE)

Other standard python libraries:
- [Biopython](https://biopython.org/) [(BSD 3-Clause License)](https://github.com/biopython/biopython/blob/master/LICENSE.rst)
- [Keras](https://keras.io/) [(Apache v2.0 License)](https://github.com/keras-team/keras/blob/master/LICENSE)
- [Numpy](https://numpy.org/) [(BSD 3-Clause License)](https://github.com/numpy/numpy/blob/main/LICENSE.txt)
- [Pandas](https://pandas.pydata.org/) [(BSD 3-Clause License)](https://github.com/pandas-dev/pandas/blob/main/LICENSE)
- [Scipy](https://scipy.org/) [(BSD 3-Clause License)](https://github.com/scipy/scipy/blob/main/LICENSE.txt)
- [Scikit-learn](https://scikit-learn.org/stable/) [(BSD 3-Clause License)](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING)
- [Tensorflow](https://www.tensorflow.org/) [(Apache v2.0 License)](https://github.com/tensorflow/tensorflow/blob/master/LICENSE)

## Installation

For complete installation instructions refer [here](INSTALL.md).

### Downloading SwissProt

This is the database used for our analysis, but can be changed according to the user purposes, e.g. TrEMBL. What is important is to add the correct path in the [configuration file](ligysis_config.txt). To download SwissProt, follow the next steps.

```sh
# download SwissProt in fasta format (88MB)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# decrompress the file
gzip -d uniprot_sprot.fasta.gz
```

### Downloading gnomAD v2.1

This is the database used for our analysis, but can be changed according to the user purposes, e.g. v > 2.1. What is important is to add the correct path in the [configuration file](ligysis_config.txt). To download gnomAD v2.1, follow the next steps.
```
# download gnomAD Exomves vcf (large file 58GB)
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
```

For more information, refer to [gnomAD](https://gnomad.broadinstitute.org/).

After downloading gnomAD, it is required to run VEP on it, as VarAlign uses its annotations. Information on how to do so [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html).

## Execution

XXX


## Citation

If you use **LIGYSIS** pipeline, please cite:

**Utgés JS**, MacGowan SA, Ives CM, Barton GJ. Classification of likely functional class for ligand binding sites identified from fragment screening. Commun Biol. 2024 Mar 13;7(1):320. doi: [10.1038/s42003-024-05970-8](https://www.nature.com/articles/s42003-024-05970-8). PMID: 38480979; PMCID: PMC10937669.

**Utgés JS** & Barton GJ. Comparative evaluation of methods for the prediction of protein-ligand binding sites, 08 August 2024, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-4849153/v1](https://doi.org/10.21203/rs.3.rs-4849153/v1)

## References
