***

<img width="666" alt="Logo" src="https://github.com/CBIIT-CGBB/3DVizSNP/assets/56087985/83be1453-ce64-42dc-accb-d06a0c95bcf3">


***

</br>

## Table of Contents
- [Introduction](#Introduction)
- [Usage](#Usage)
- [Dependencies](#Dependencies)
- [Acknowledgments](#Acknowledgments)
- [License](#License)

## Introduction

- **Single nucleotide polymorphisms (SNPs)** are characterized by a change in the single-nucleotide of the DNA sequence occurring at a specific position. Rapid advancements in next-generation sequencing (NGS) technologies have eased the identification of these point mutations. They account for approximately 90% of genetic variations in humans and are found to be associated with several diseases. (Collins, et al. [1998](https://genome.cshlp.org/content/8/12/1229.short)) Tools such as [VEP](https://useast.ensembl.org/info/docs/tools/vep/index.html) and [OpenCRAVAT](https://opencravat.org) have been developed to predict the deleterious effects of these variants. However, the consequences of SNPs on protein structure and function are poorly understood. Hence, there is a great need to integrate this data with the available structural data. Mapping of SNPs to the protein structure and its visualization can significantly help us in studying various disease-causing mutations.

- **[iCn3D](https://www.ncbi.nlm.nih.gov/Structure/icn3d/)** is a powerful, web-based 3D viewer used for representing biomolecular structures. It facilitates visualization in 1D, 2D, and complex 3D while allowing synchronization of the selection across various structural displays. One of the important features of iCn3D is its ability to generate a shareable link that includes all the custom filters and user-provided labels. Anyone with this link can reproduce the same display of the structure. Annotations can either be directly extracted and mapped to the protein structure from various NCBI databases (such as dbSNP, ClinVar, Conserved Domain Database, and others), or users can submit their own annotations via custom tracks. Finally, it facilitates the display of both sequence-structure and structure-structure alignments with corresponding superposition.

- **3DVizSNP**: 3DVizSNP extracts SNPs from VCF files, submits them to the [Ensembl VEP API](https://rest.ensembl.org/#VEP), generates iCn3D links to load the mutation in either the PDB structure or the AlphaFold structure, along with a table combining the various IDs and SIFT and PolyPhen scores from VEP. It can be run locally as a python script (the [pysam](https://pysam.readthedocs.io/en/latest/installation.html) library is required) or can be accessed as a webserver at [https://analysistools.cancer.gov/3dvizsnp](https://analysistools.cancer.gov/3dvizsnp).

## Usage

- To run the script locally you need a VCF file that is bgzipped (.gz extension) and the associated tabix index file (.tbi extension) from which the variants are to be extracted. 

- You can test the script on a sample VCF file provided in the `example` directory.

  ```
  python3 VizSNPSt.py -v test.vcf.gz 
  ```

- The output HTML and .csv files are available in the `example` directory.

- Use the -h flag to see the available command-line options.

## Dependencies
- Python version >= 3
- Required modules:
  * pysam (install with conda: `conda install -c bioconda pysam`)
  * pandas, other common python packages

## Acknowledgments
- This project began at the [ISMB 2022 Hackathon](https://github.com/hackathonismb/VizSNP-St) sponsored by the **International Society for Computational Biology** and the **National Center for Biotechnology Information (NCBI)**.

## License
Licensed under MIT License - Copyright (c) 2022 hackathonismb (Refer LICENSE file for more details)
