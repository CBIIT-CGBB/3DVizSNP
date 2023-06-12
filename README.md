***

<p align="center">
  <img width="666" alt="Logo" src="https://github.com/CBIIT-CGBB/3DVizSNP/assets/56087985/83be1453-ce64-42dc-accb-d06a0c95bcf3">
</p>


***

</br>

## Table of Contents
- [Introduction](#Introduction)
- [Installation](#Installation)
- [Usage](#Usage)
- [Dependencies](#Dependencies)
- [Acknowledgments](#Acknowledgments)
- [License](#License)
- [Citation](#Citation)

## Introduction

- **Single nucleotide polymorphisms (SNPs)** are characterized by a change in the single-nucleotide of the DNA sequence occurring at a specific position. Rapid advancements in next-generation sequencing (NGS) technologies have eased the identification of these point mutations. They account for approximately 90% of genetic variations in humans and are found to be associated with several diseases. (Collins, et al. [1998](https://genome.cshlp.org/content/8/12/1229.short)) Tools such as [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) and [OpenCRAVAT](https://opencravat.org) have been developed to predict the deleterious effects of these variants. However, the consequences of SNPs on protein structure and function are poorly understood. Hence, there is a great need to integrate this data with the available structural data. Mapping of SNPs to the protein structure and its visualization can significantly help us in studying various disease-causing mutations.

- **[iCn3D](https://www.ncbi.nlm.nih.gov/Structure/icn3d/)** is a powerful, web-based 3D viewer used for representing biomolecular structures. It facilitates visualization in 1D, 2D, and complex 3D while allowing synchronization of the selection across various structural displays. One of the important features of iCn3D is its ability to generate a shareable link that includes all the custom filters and user-provided labels. Anyone with this link can reproduce the same display of the structure. Annotations can either be directly extracted and mapped to the protein structure from various NCBI databases (such as dbSNP, ClinVar, Conserved Domain Database, and others), or users can submit their own annotations via custom tracks. Finally, it facilitates the display of both sequence-structure and structure-structure alignments with corresponding superposition.

- **3DVizSNP**: 3DVizSNP extracts SNPs from VCF files, submits them to the [Ensembl VEP API](https://rest.ensembl.org/#VEP), generates iCn3D links to load the mutation in either the PDB structure or the AlphaFold structure, along with a table combining the various IDs and SIFT and PolyPhen scores from VEP. It can be run locally as a python script (the [pysam](https://pysam.readthedocs.io/en/latest/installation.html) library is required) or can be accessed as a webserver at [https://analysistools.cancer.gov/3dvizsnp](https://analysistools.cancer.gov/3dvizsnp).

## Installation
3DVizSNP is developed in Python >= 3 (see the Dependencies section for more details). We have provided a conda environment file in this repository that you can use to install the required dependencies. Most other standard libraries should already be available on your system.

- Step 1: Obtain 3DVizSNP on your local machine by cloning the repository.

  ```
  git clone https://github.com/CBIIT-CGBB/3DVizSNP.git 
  ```

- Step 2: Create and activate the conda environment using our provided file. After this step, you should be ready to use the 3DVizSNP tool on the command-line.

  ```
  cd 3DVizSNP
  conda env create -f environment_3DVizSNP.yml
  conda activate environment_3DVizSNP
  ```

## Usage

- To run the script locally you need a VCF file that is bgzipped (.gz extension) and the associated tabix index file (.tbi extension) from which the variants are to be extracted. 

- You can test the script on a sample VCF file provided in the `example` directory.

  ```
  python3 3DVizSNP.py -v example/test.vcf.gz 
  ```

- The output HTML and .csv files will be saved in the `example` directory.

- Use the -h flag to see the available command-line options.

## Dependencies
- Python version >= 3
- Required modules:
  * pysam
  * pandas, other standard python packages

## Acknowledgments
- This project began at the [ISMB 2022 Hackathon](https://github.com/hackathonismb/VizSNP-St) sponsored by the **International Society for Computational Biology** and the **National Center for Biotechnology Information (NCBI)**.

## License
Licensed under MIT License - Copyright (c) 2022 hackathonismb (Refer LICENSE file for more details)

## Citation
If you use 3DVizSNP in your research, please cite as: Sierk M, Ratnayake S, Wagle MM, Chen B, Park B, Wang J, Youkharibache P & Meerzaman D (2023). "3DVizSNP: A Tool for Rapidly Visualizing Missense Mutations Identified in High Throughput Experiments in iCn3D." <i>BMC Bioinformatics</i>, 24, 244.
