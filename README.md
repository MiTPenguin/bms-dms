# TYK2 Deep Mutational Scanning by Octant and Bristol Meyers Squibb

Welcome! This repository contains analyses, code, and summary statistics related to and reproducing the manuscript `Deep mutational scanning to characterize the mechanistic role of TYK2 in immune signaling and disease`. This repository contains the following:

```
├── data
├── docker
├── paper
├── src
├── sumstats
└── validation
```

Code to conduct statistical modeling via negative binomial generalized linear mixed modeling in R is located in [src](./src/), and resulting summary statistics for the manuscript are located in [sumstats](./sumstats/). The [paper](./paper) directory contains several Quarto markdown files, each describing one major figure in the manuscript, as well as auxiliary plotting files.

To build an environment with all utilities installed, execute the following to create a docker image containing a nix development environment:

```
cd docker
docker build -t analysis-image .
```

Then, run an instance of the image, enter the initiated container, and open the nix environment:

```
docker run --rm -dit \
  --name analysis-container \
  -v ~/bms-dms:/bms-dms \
  analysis-image:latest

docker exec -it analysis-container bash

nix develop
```
