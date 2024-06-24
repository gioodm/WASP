# WASP: Protein Functional Annotation using AlphaFold structures
  
Welcome to the official repository for the paper *WASP: A pipeline for functional annotation of proteins using AlphaFold structural models*. **WASP**, **<ins>W</ins>hole-proteome <ins>A</ins>nnotation through <ins>S</ins>tructural-homology <ins>P</ins>ipeline**, is a python-based software designed for comprehensive organism annotation at the whole-proteome level based on structural homology.

WASP is a user-friendly command-line tool that only requires the NCBI taxonomy ID of the organism of interest as an input. Using the computational speed of Foldseek [[1](https://doi.org/10.1038/s41587-023-01773-0)], WASP generates a graphical representation of reciprocal hits between the organism protein query and the AlphaFold database [[2](https://doi.org/10.1038/s41586-021-03819-2), [3](https://doi.org/10.1093/nar/gkab1061)], enabling downstream robust functional enrichment and statistical testing. WASP annotates uncharacterised proteins using multiple functional descriptors, including GO terms, Pfam domains, PANTHER family classification and CATH superfamilies, Rhea IDs and EC numbers.\\
Additionally, WASP provides a module to map native proteins to orphan reactions in genome-scale models based on structural homology.

<p align="center">
<img src="fig.png" alt="drawing" width="600"/>
</p>


<!---
If you find WASP helpful in your research, please cite us:

    @article{,
      author = {},
      title = {},
      journal = {},
      volume = {},
      number = {},
      pages = {},
      year = {},
      doi = {},
      URL = {}
    }  
-->


## Table of Contents

1. [Install](#1-install)
    - [Requirements](#11-requirements)
    - [Quickstart](#12-quickstart)
2. [Run](#2-run)
    - [Whole-proteome annotation](#21-whole-proteome-annotation)
    - [GEM gap-filling module](#22-gem-gap-filling-module)
3. [References](#3-references)


## 1. Install

### 1.1 Requirements

- Python >= 3.9
- Foldseek >= 8
- gsutil

The manuscript results were obtained using Python 3.10.14 and Foldseek 8-ef4e960.

### 1.2 Quickstart

Start by cloning the repository, create the WASP environment and install all python requirements:

```sh
mkdir WASP
cd WASP/
git clone https://github.com/gioodm/WASP.git

conda create -n WASP python==3.10 -y
conda activate WASP
pip install -r requirements.txt
```

Then, install Foldseek in the `/bin` of the project root. Follow the installation instructions at [Foldseek GitHub](https://github.com/steineggerlab/foldseek).

Example for a Linux AVX2 build:


```sh
# Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

```

Install gsutil and initialise the gcloud CLI, following instructions for your machine at [Google Cloud Storage Documentation](https://cloud.google.com/storage/docs/gsutil_install).


## 2. Run

### 2.1 Whole-proteome annotation

#### Usage

```sh
./run.sh [-h] [-e evalue_threshold] [-b bitscore_threshold] [-n max_neighbours] [-N max_neighbours_NaN] taxid
```

First, identify the NCBI taxonomy ID of your organism of interest. You can find the ID at [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy). You can use AlphaFold DB to check how many structures are linked to that ID by searching the same ID on the AFDB search bar. This requires gsutil installed.

Once you have identified the taxid, run the pipeline as follows (e.g., for organism S. cerevisiae S288c, taxid: 559292):

```sh
    chmod +x run.sh

    ./run.sh 559292

```

This requires `gsutil` installed.

Additional parameters can be customised, including:

    `-e` evalue_threshold: set the evalue threshold (default: 10e-10)
    `-b` bitscore threshold: set the bitscore threshold (default: 50)
    `-n` max_neighbours: set the max number of neighbours (default: 10)
    `-s` step: set step to add to max neighbours (n) in additional iterations (default: 10)
    `-i` iterations: set number of iterations to perform (default: 3)

Usage examples:

```sh
    ./run.sh -e 1e-50 -b 200 -n 5 -i 5 559292
    ./run.sh -s 5 559292
```

<!---
If you want to use a custom dataset to predict function (for example, a newly sequenced genome or a set of proteins from different organisms which cannot downloaded in bulk using the above AlphaFold Bucket gsutil command) you can do so by creating a custom tarred folder containing the protein structures (in either .cif.gz or .pdb.gz format) and placing it in a folder called proteomes/ within the WASP folder. Then you can run WASP with:

```
    chmod +x run.sh

    ./run.sh custom_folder.tar

```


### 2.2 GEM gap-filling module

#### Pre-processing

Potential modifications to the python scripts might be needed depending on the GEM format.

Some pre-processing steps are required to obtain a standardised input file for the WASP pipeline.

find_orphans.py : identifies orphan reactions in the GEM (accepted extensions: .xml, .sbml, .json and .mat) using the Python3 cobrapy module. Annotation in different formats (accepted: BiGG, Rhea, EC number, KEGG, PubMed, MetaNetX) present in the model is retrieved;

rxn2code.py : each reaction annotation is mapped to the corresponding Rhea reaction ID and/or EC number, when available. Output file includes: 1. reaction id, 2. reaction extended name, 3. Rhea IDs, 4. EC numbers

If MetaNetX codes are present, the reac_xref.tsv file (retrieved from: https://www.metanetx.org/mnxdoc/mnxref.html) is needed.
The final gaps_file must contain all the 4 columns; if the reaction extended name is not present in the model, an empty column should be present. If no Rhea/EC IDs are identified, empty columns should be present instead.

#### Usage

```
./run_GEM.sh [-h] [-e evalue_threshold] [-b bitscore_threshold] [-t tmscore] taxid gaps_file
```

*gaps_file has to be generated manually, check pre-processing
-->





## References

[1] van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, et al. Fast and accurate protein structure
search with Foldseek. Nature Biotechnology. 2023. doi: https://doi.org/10.1038/s41587-023-01773-0

[2] Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, et al. Highly accurate protein structure prediction
with AlphaFold. Nature. 2021;596(7873):583-9. doi: https://doi.org/10.1038/s41586-021-03819-2

[3] Varadi M, Bertoni D, Magana P, Paramval U, Pidruchna I, Radhakrishnan M, et al. AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. Nucleic Acids Research. 2023;52(D1):D368–D375. doi: http://dx.doi.org/10.1093/nar/gkad1011
