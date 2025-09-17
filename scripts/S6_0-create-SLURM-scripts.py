if __name__ == "__main__":

    template = """#!/bin/bash -x
#SBATCH --account=TODO
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=jureca$RUN.%j
#SBATCH --error=jureca$RUN.%j
#SBATCH --time=23:59:59
#SBATCH --partition=dc-cpu

apptainer run /p/project/loki/jadebeck1/fluxtopus.sif \
python scripts/S6_1-INST-mcmc.py --fml models/Syn.fml --configs bayes --path measured_poolsizes_run$RUN --n_samples 12000 \
--iteration 0 --rng $RNG --n_chains 1 --n_temps $N_TEMPS --configs freeflux_synthetic_data_measured_poolsizes \
--mcmc_algo GaussianHitAndRun --thinning 100 --step_size $STEP_SIZE --draws_per_exchange_attempt 100
"""

    step_sizes = [1]
    n_temps = [128]

    call_script = ''

    for i in range(len(n_temps)):
        for j in range(len(step_sizes)):
            script = template
            script = script.replace('$N_TEMPS', str(n_temps[i]))
            script = script.replace('$STEP_SIZE', str(step_sizes[j]))
            script = script.replace('$RUN', str(i*len(step_sizes)+j))

            with open(f'slurm_{i*len(step_sizes)+j}.sh', 'w') as f:
                f.write(script)

            call_script += f'sbatch ./slurm_{i*len(step_sizes)+j}.sh\n'

        with open(f'run-slurm-sampling.sh', 'w') as f:
            f.write(call_script)
