# Phunky
Long read assembly for phages:
![Phage pipeline](pipeline.png)
**Figure 1:** Rough phage assembly pipeline, dotted line indicates where processing begins.

[![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)](https://pypi.org/project/phunky/)

## Usage
Install:
```
pip install phunky
```

Open a python terminal and enter:
```py
import phunky
dir(phunky)
```

Quick phage assembly:
```py
import phunky
phunky.phage_assembly_pipeline('example_reads.bam', 'output_directory/')
```

Batch phage assembly:
```py
import phunky
phunky.batch_phage_assembly_pipeline('reads_directory/', 'output_directory/')
```


### dependencies:
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


## Todo
* Create conda distribution
* Add logging function
* Add hash key function
* Add multiprocessing
