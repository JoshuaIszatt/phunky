# Phunky
[![DOI](https://zenodo.org/badge/805134858.svg)](https://doi.org/10.5281/zenodo.13981019)

Long read assembly for phages:
![Phage pipeline](pipeline.png)
**Figure 1:** Rough phage assembly pipeline, dotted line indicates where processing begins.

[![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)](https://pypi.org/project/phunky/)

## Install dependencies
1. Download the yml environment file:
```
wget https://anaconda.org/JoshIszatt/phunky/2024.05.24.145220/download/phunky.yml
```

2. Create the environment
```
conda env create --file phunky.yml
```

3. Activate the environment:
```
conda activate phunky
```

4. Install phunky pipeline:
```
pip install phunky
```

5. Optional: Setup the checkv database 
```
checkv download_database /path/to/checkv-db
```

```
export CHECKVDB=/path/to/checkv-db
```

## Usage
Open a python terminal and enter:
```py
import phunky
dir(phunky)
```

Quick phage assembly:

```py
import phunky

phunky.assembly_pipeline('example_reads.bam', 'output_directory/')
```

Batch phage assembly:
```py
import phunky
phunky.batch_assembly_pipeline('reads_directory/', 'output_directory/')
```

## Citation:
There is a CITATION file in this repository (.CFF), however Phunky may also be cited like so:
```
Iszatt, Joshua James. (2024). PHUNKY long read assembly pipeline (Version v0.3) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.13981020
```

## Dependencies:
  - python>=3
  - checkv==1.0.3
  - biopython==1.83
  - bbmap==39.06
  - pandas==2.2.1
  - matplotlib==3.8.4
  - flye==2.9.3
  - porechop_abi==0.5.0
  - nanoplot==1.42.0
  - filtlong==0.2.1


## General todo
* Create conda distribution
* Add hash key function
