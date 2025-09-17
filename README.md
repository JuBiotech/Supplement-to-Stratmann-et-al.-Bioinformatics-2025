# Online Supplementary Information for 13CFLUX -- Third generation high-performance engine for isotopically (non)stationary <sup>13</sup>C metabolic flux analysis

With this script you will be able to reproduce figures and data from the Supplementary information. 

First, get the docker image by
```shell
docker pull jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux
```
and clone this repository.

Then run the examples by
```shell
docker run -v .:/task jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux:latest /task/run.sh
```
while being in the root directory of this repository.

The output will be written to the `out` directory. 

On an AMD EPYC 9334 execution takes about 1 hour


## Sampling (S.6.2)

Because Bayesian inference is much more computationally costly than the other tasks,
we provide the raw samples on [zenodo](https://doi.org/10.5281/zenodo.17100887).
After downloading extract them to the `data` subdirectory
Additionally, the samples can be technically be reproduced by calling
```shell
docker run -v .:/task jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux:latest /task/run-sampling.sh
```
For our run, we used 4 nodes with 128 cores to run 512 chains in parallel. For this we used slurm and the slurm scripts can be reproduced by scripts/S6\_0-create-SLURM-scripts.py.
We advise using the slurm scripts and a hpc-system.
