# Supplement-to-Stratmann-et-al.-Bioinformatics-2025


How to run it:

First, get the docker image by
```shell
docker pull jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux
```

Then run the examples by
```shell
docker run -v .:/task jugit-registry.fz-juelich.de/ibg-1/modsim/fluxomics/13cflux:latest /task/run.sh
```
while being in the root directory of this repository.
