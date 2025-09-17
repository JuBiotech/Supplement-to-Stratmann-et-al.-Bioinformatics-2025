#!/bin/bash

PYTHON_EXEC="/opt/python/cp313-cp313/bin/python"

cd /task
${PYTHON_EXEC} -u scripts/S6_1-INST-mcmc.py  --rng 6 --fml models/Syn.fml --configs bayes --path output --n_samples 24000 --iteration 0 --thinning 100 --n_chains 1 --n_temps 128 --draws_per_exchange_attempt 200 --step_size 1
${PYTHON_EXEC} -u scripts/S6_1-INST-mcmc.py  --rng 7 --fml models/Syn.fml --configs bayes --path output --n_samples 24000 --iteration 0 --thinning 100 --n_chains 1 --n_temps 128 --draws_per_exchange_attempt 200 --step_size 1
${PYTHON_EXEC} -u scripts/S6_1-INST-mcmc.py  --rng 8 --fml models/Syn.fml --configs bayes --path output --n_samples 24000 --iteration 0 --thinning 100 --n_chains 1 --n_temps 128 --draws_per_exchange_attempt 200 --step_size 1
${PYTHON_EXEC} -u scripts/S6_1-INST-mcmc.py  --rng 9 --fml models/Syn.fml --configs bayes --path output --n_samples 24000 --iteration 0 --thinning 100 --n_chains 1 --n_temps 128 --draws_per_exchange_attempt 200 --step_size 1
