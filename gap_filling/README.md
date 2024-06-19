# WASP: Protein Functional Annotation using AlphaFold structures
  
This is the official repository for the paper [*WASP: A pipeline for functional annotation of proteins using AlphaFold structural models*](link). **WASP**, **<ins>W</ins>hole-proteome <ins>A</ins>nnotation through <ins>S</ins>tructural-homology <ins>P</ins>ipeline**, is a python-based command-line software designed for comprehensive organism annotation at the whole-proteome level based on structural homology.

[To use CLEAN to inference the EC number for any amino acid sequence, we included the pretrained weights for both the 70% and 100% identity clustering split of SwissProt (expertly reviewed portion of the UniProt, total ~220k training data). User can follow the instruction below on how to install and inference with CLEAN. We also provide full training scripts.

<p  align="center">

<img  src="CLEAN.png"  alt="drawing"  width="600"/>

</p>


If you find WASP helpful in your research, please cite us:

    @article{doi:10.1126/science.adf2465,
      author = {Tianhao Yu  and Haiyang Cui  and Jianan Canal Li  and Yunan Luo  and Guangde Jiang  and Huimin Zhao },
      title = {Enzyme function prediction using contrastive learning},
      journal = {Science},
      volume = {379},
      number = {6639},
      pages = {1358-1363},
      year = {2023},
      doi = {10.1126/science.adf2465},
      URL = {https://www.science.org/doi/abs/10.1126/science.adf2465}
    }  
]: #

## 1. Install

### 1.1 Requirements

Python >= 3.6; PyTorch >= 1.11.0; CUDA >= 10.1
Manuscript result was obtained using Python 3.10.4; PyTorch 1.11.0; CUDA 11.3; fair-esm 1.0.2

  

### 1.2 Quickstart


```

cd CLEAN/app/
conda create -n clean python==3.10.4 -y
conda activate clean
pip install -r requirements.txt

```


