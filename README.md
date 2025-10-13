# TYK2 Deep Mutational Scanning by Octant and Bristol Meyers Squibb

Welcome! This repository contains analyses, code, and summary statistics related to and reproducing the manuscript [Deep mutational scanning to characterize the mechanistic role of TYK2 in immune signaling and disease](https://www.biorxiv.org/content/10.1101/2025.10.11.681520v1).

<center> <img src="./aux/header.jpeg" width="500"/> </center>

This repository contains the following:

```
├── aux
├── barcode_maps
├── docker
├── fastq
├── mapped_counts
├── paper
├── sra
├── src
└── sumstats
```

Code to conduct statistical modeling via negative binomial generalized linear mixed modeling in R is located in [src](./src/), and resulting summary statistics for the manuscript are located in [sumstats](./sumstats/). The [paper](./paper) directory contains several Quarto markdown files, each describing one major figure in the manuscript, as well as auxiliary plotting files.

To build an environment with all utilities installed, execute the following to create a docker image containing a nix development environment:

```
cd docker
docker build -t tyk2-dms-image .
```

Then, run an instance of the image, enter the initiated container, and open the nix environment:

```
docker run --rm -dit \
  --name tyk2-dms-container \
  -v ~/bms-dms:/bms-dms \
  tyk2-dms-image:latest

docker exec -it tyk2-dms-container bash

nix develop
```

### Other Supporting Data
  -  [Biorxiv DOI 2025.10.11.681520](https://doi.org/10.1101/2025.10.11.681520): Manuscript pre-print
  -  [SRA BioProject PRJNA1291213](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1291213): Raw FASTQ files
  -  [Zenodo 15347448](https://zenodo.org/records/15347448): Oligo-barcode maps and raw barcode counts [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15347448.svg)](https://doi.org/10.5281/zenodo.15347448)
