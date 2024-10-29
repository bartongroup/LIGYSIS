# LIGYSIS

This repository contains the original version of our original ligand site analysis **LIGYSIS** pipeline, which was employed to analyse all biologically relevant protien-ligand interactions on the [PDBe](https://www.ebi.ac.uk/pdbe/) [[1](https://europepmc.org/article/MED/31691821)], which results are served in the **LIGYSIS** [web server](https://www.compbio.dundee.ac.uk/ligysis/). The code for the web server can be found [here](https://github.com/JavierSanchez-Utges/LIGYSIS-web).

**LIGYSIS** makes extensive use of the [PDBe-KB](https://www.ebi.ac.uk/pdbe/pdbe-kb/) [[2](https://academic.oup.com/nar/article/48/D1/D344/5580911)] and PDBe [APIs](https://www.ebi.ac.uk/pdbe/pdbe-rest-api) [[3](https://academic.oup.com/bioinformatics/article/37/21/3950/6291664)] to retrieve transformation matrices, all PDB IDs mapping to a given [UniProt](https://www.uniprot.org/) [[4](https://academic.oup.com/nar/article/49/D1/D480/6006196?login=true)] ID and more. LIGYSIS analyses protein-ligand contacts from BioLiP defined biologically relevant ligands from [PISA](https://www.ebi.ac.uk/pdbe/pisa/)-defined preferred biological assemblies [[5](https://www.sciencedirect.com/science/article/pii/S0022283607006420)]. These protein-ligand interactions are used to cluster ligands together into binding sites, which are then characterised by using evolutionary divergence, human genetic variation and structural features such as [Relative Solvent Accessibility](https://en.wikipedia.org/wiki/Relative_accessible_surface_area) (RSA) [[6](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635#:~:text=It%20is%20defined%20as%20a,1%5D%E2%80%93%5B5%5D.)] or secondary structure.

## Pipeline methodology

The pipeline can be summarised in the following steps:

1. Transformation matrices are retrieved for each chain mapping to a UniProt ID from the [PDBe-KB FTP site](http://ftp.ebi.ac.uk/pub/databases/pdbe-kb/superposition/P/P00517/P00517.json).
2. Segment superposition data, i.e., segment-clustered protein chains and representative chains, are retrieved from the PDBe GRAPH API [superposition endpoint](https://www.ebi.ac.uk/pdbe/graph-api/uniprot/superposition/P00517) [[7](https://pubs.aip.org/aca/sdy/article/11/3/034701/3294234)].

   **Note:** The following steps are repeated in a loop now for each Segment within the protein.
3. Retrieve experimental data for each structure from the PDBe REST API [experiment endpoint](https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/1svh)
4. Filter out structures according to resolution (Current dataset considers <i>ALL</i> structures).
5. Download of biological assemblies from PDBe via [ProIntVar](https://github.com/bartongroup/ProIntVar) [[8](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3783)].
6. Protein-ligand interactions calculation with [pdbe-rpeggio](https://github.com/PDBeurope/arpeggio) [[9](https://www.sciencedirect.com/science/article/pii/S0022283616305332?via%3Dihub)].
7. Mapping of PDB residues to UniProt through [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) [[10](https://academic.oup.com/nar/article/47/D1/D482/5184711)], encoded in the [<i>mmCIF</i> format](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner%E2%80%99s-guide-to-pdbx-mmcif).
8. Ligand clustering into binding sites using [SciPy](https://scipy.org/).
9. Transformation of all chains employing PDBe-KB matrices with [BioPython](https://biopython.org/) (`atom.transform(rotation, translation)`).
10. <i>Simplification</i> of superposed chains. This consists in keeping protein atoms only for one of the superposed structures, and heteroatoms for the rest. This is done to generate a lower-weight superposition (all ligands to a single protein scaffold).
11. Generation of [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) [[11](https://onlinelibrary.wiley.com/doi/10.1002/pro.3943)] superposition visualisation script.
12. Relative solvent accessibility (RSA) and secondary structure elements calculation with [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) [[12](https://onlinelibrary.wiley.com/doi/10.1002/bip.360221211)] via ProIntVar.
13. Multiple sequence alignment with [jackhmmer](http://hmmer.org/) [[13](https://academic.oup.com/bioinformatics/article/14/9/755/259550)].
14. Shenkin amino acid divergence score calculation [[14](https://doi.org/10.1002/prot.340110408), [15](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335)].
15. Missense enrichment score calculation with [VarAlign](https://github.com/bartongroup/SM_varalign) [[16](https://www.biorxiv.org/content/10.1101/127050v2), [17](https://www.nature.com/articles/s42003-024-06117-5)].
16. RSA-based clustering label and functional score calculation [[18](https://www.nature.com/articles/s42003-024-05970-8)].

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

1. Armstrong DR, Berrisford JM, Conroy MJ, Gutmanas A, Anyango S, Choudhary P, Clark AR, Dana JM, Deshpande M, Dunlop R, Gane P, Gáborová R, Gupta D, Haslam P, Koča J, Mak L, Mir S, Mukhopadhyay A, Nadzirin N, Nair S, Paysan-Lafosse T, Pravda L, Sehnal D, Salih O, Smart O, Tolchard J, Varadi M, Svobodova-Vařeková R, Zaki H, Kleywegt GJ, Velankar S. PDBe: improved findability of macromolecular structure data in the PDB. Nucleic Acids Res. 2020 Jan 8;48(D1):D335-D343. doi: [10.1093/nar/gkz990](https://europepmc.org/article/MED/31691821). PMID: 31691821; PMCID: PMC7145656.

2. PDBe-KB consortium. PDBe-KB: a community-driven resource for structural and functional annotations. Nucleic Acids Res. 2020 Jan 8;48(D1):D344-D353. doi: [10.1093/nar/gkz853](https://academic.oup.com/nar/article/48/D1/D344/5580911). PMID: 31584092; PMCID: PMC6943075.

3. Nair S, Váradi M, Nadzirin N, Pravda L, Anyango S, Mir S, Berrisford J, Armstrong D, Gutmanas A, Velankar S. PDBe aggregated API: programmatic access to an integrative knowledge graph of molecular structure data. Bioinformatics. 2021 Nov 5;37(21):3950-3952. doi: [10.1093/bioinformatics/btab424](https://academic.oup.com/bioinformatics/article/37/21/3950/6291664). PMID: 34081107; PMCID: PMC8570819.

4. UniProt Consortium. UniProt: the universal protein knowledgebase in 2021. Nucleic Acids Res. 2021 Jan 8;49(D1):D480-D489. doi: [10.1093/nar/gkaa1100](https://academic.oup.com/nar/article/49/D1/D480/6006196?login=true). PMID: 33237286; PMCID: PMC7778908.

5. Krissinel E, Henrick K. Inference of macromolecular assemblies from crystalline state. J Mol Biol. 2007 Sep 21;372(3):774-97. doi: [10.1016/j.jmb.2007.05.022](https://www.sciencedirect.com/science/article/pii/S0022283607006420). Epub 2007 May 13. PMID: 17681537.

6. Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum allowed solvent accessibilites of residues in proteins. PLoS One. 2013 Nov 21;8(11):e80635. doi: [10.1371/journal.pone.0080635](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635#:~:text=It%20is%20defined%20as%20a,1%5D%E2%80%93%5B5%5D.). PMID: 24278298; PMCID: PMC3836772.

7. Ellaway JIJ, Anyango S, Nair S, Zaki HA, Nadzirin N, Powell HR, Gutmanas A, Varadi M, Velankar S. Identifying protein conformational states in the Protein Data Bank: Toward unlocking the potential of integrative dynamics studies. Struct Dyn. 2024 May 17;11(3):034701. doi: [10.1063/4.0000251](https://pubs.aip.org/aca/sdy/article/11/3/034701/3294234). PMID: 38774441; PMCID: PMC11106648.

8. MacGowan SA, Madeira F, Britto-Borges T, Warowny M, Drozdetskiy A, Procter JB, Barton GJ. The Dundee Resource for Sequence Analysis and Structure Prediction. Protein Sci. 2020 Jan;29(1):277-297. doi: [10.1002/pro.3783](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.3783). Epub 2019 Nov 28. PMID: 31710725; PMCID: PMC6933851.

9. Jubb HC, Higueruelo AP, Ochoa-Montaño B, Pitt WR, Ascher DB, Blundell TL. Arpeggio: A Web Server for Calculating and Visualising Interatomic Interactions in Protein Structures. J Mol Biol. 2017 Feb 3;429(3):365-371. doi: [10.1016/j.jmb.2016.12.004](https://www.sciencedirect.com/science/article/pii/S0022283616305332?via%3Dihub). Epub 2016 Dec 10. PMID: 27964945; PMCID: PMC5282402.

10. Velankar S, Dana JM, Jacobsen J, van Ginkel G, Gane PJ, Luo J, Oldfield TJ, O'Donovan C, Martin MJ, Kleywegt GJ. SIFTS: Structure Integration with Function, Taxonomy and Sequences resource. Nucleic Acids Res. 2013 Jan;41(Database issue):D483-9. doi: [10.1093/nar/gks1258](https://academic.oup.com/nar/article/47/D1/D482/5184711). Epub 2012 Nov 29. PMID: 23203869; PMCID: PMC3531078.

11. Pettersen EF, Goddard TD, Huang CC, Meng EC, Couch GS, Croll TI, Morris JH, Ferrin TE. UCSF ChimeraX: Structure visualization for researchers, educators, and developers. Protein Sci. 2021 Jan;30(1):70-82. doi: [10.1002/pro.3943](https://onlinelibrary.wiley.com/doi/10.1002/pro.3943). Epub 2020 Oct 22. PMID: 32881101; PMCID: PMC7737788.

12. Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers. 1983 Dec;22(12):2577-637. doi: [10.1002/bip.360221211](https://onlinelibrary.wiley.com/doi/10.1002/bip.360221211). PMID: 6667333.

13. Eddy SR. Profile hidden Markov models. Bioinformatics. 1998;14(9):755-63. doi: [10.1093/bioinformatics/14.9.755](https://academic.oup.com/bioinformatics/article/14/9/755/259550). PMID: 9918945.
   
14. Shenkin PS, Erman B, Mastrandrea LD. Information-theoretical entropy as a measure of sequence variability.
Proteins. 1991; 11(4):297–313. Epub 1991/01/01. [https://doi.org/10.1002/prot.340110408](https://doi.org/10.1002/prot.340110408)
PMID: 1758884.

15. **Utgés JS**, Tsenkov MI, Dietrich NJM, MacGowan SA, Barton GJ. Ankyrin repeats in context with human population variation. PLoS Comput Biol. 2021 Aug 24;17(8):e1009335. doi: [10.1371/journal.pcbi.1009335](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009335). PMID: 34428215; PMCID: PMC8415598.
 
16. MacGowan, SA, Madeira, F, Britto-Borges, T, Schmittner, MS, Cole, C, & Barton, GJ (2017). Human missense variation is constrained by domain structure and highlights functional and pathogenic residues. bioRxiv, 127050. [https://doi.org/10.1101/127050](https://www.biorxiv.org/content/10.1101/127050v2).

17. MacGowan SA, Madeira F, Britto-Borges T, Barton GJ. A unified analysis of evolutionary and population constraint in protein domains highlights structural features and pathogenic sites. Commun Biol. 2024 Apr 11;7(1):447. doi: [10.1038/s42003-024-06117-5](https://www.nature.com/articles/s42003-024-06117-5). PMID: 38605212; PMCID: PMC11009406.
   
18. **Utgés JS**, MacGowan SA, Ives CM, Barton GJ. Classification of likely functional class for ligand binding sites identified from fragment screening. Commun Biol. 2024 Mar 13;7(1):320. doi: [10.1038/s42003-024-05970-8](https://www.nature.com/articles/s42003-024-05970-8). PMID: 38480979; PMCID: PMC10937669.
