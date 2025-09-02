# Online Supplementary Information for 13CFLUX -- Third generation high-performance simulator for <sup>13</sup>C metabolic flux analysis

With this script you will be able to reproduce most figures and data from the Supplementary information. 

First, get the docker image by
```shell
docker pull jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux
```
and clone this repository.

To reproduce the results from S.6.2 you either need to download the results from zenodo and unpack them in the `data` directory, or generate the samples as described in section S.6.2

Then run the examples by
```shell
docker run -v .:/task jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux:latest /task/run.sh
```
while being in the root directory of this repository.

The output will be written to the `out` directory. 

On an AMD EPYC 9334 execution takes about 1 hour


