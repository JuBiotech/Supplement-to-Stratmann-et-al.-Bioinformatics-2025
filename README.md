# Online Supplementary Information for: 13CFLUX -- Third generation high-performance engine for isotopically (non)stationary <sup>13</sup>C metabolic flux analysis

With this repository you will be able to reproduce figures and data from the supplementary information (SI). 

First, get the `13CFlUX` docker image by

```shell
docker pull jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux
```
and clone this repository.

Then generate the results by tunning
```shell
docker run -v .:/task jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux:latest /task/run.sh
```
while being in the root directory of this repository.


**Output:*** The script generates figures S.1-6 and S.8-9 as `png` and `svg`. The files will be written to the `out` directory.
The command line output will look [like this](console_output.txt)

**Run time** On an AMD EPYC 9334 execution takes about 1 hour


## Results of SI Section S.6.2

In section S.6 of the SI, Bayesian INST <sup>13</sup>C-MFA we do Bayesian inference, which is computationally much more costly than the other tasks.
Therefore, we provide the raw samples on [zenodo](https://doi.org/10.5281/zenodo.17100887).
Download and extract them to the `data` subdirectory, in order to be able to reproduce Figure S.10


Alternatively, the samples can be technically be reproduced by calling
```shell
docker run -v .:/task jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux /task/run-sampling.sh
```
For our run, we used 4 nodes with 128 cores to run 512 chains in parallel. For this we used slurm and the slurm scripts can be reproduced by `scripts/S6\_0-create-SLURM-scripts.py`.
We advise using the slurm scripts and a hpc-system.
