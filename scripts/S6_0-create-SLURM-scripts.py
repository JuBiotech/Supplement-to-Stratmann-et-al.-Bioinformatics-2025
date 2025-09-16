if __name__ == "__main__":

    template = """#!/bin/bash -x
#SBATCH --account=TODO
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=jureca$RUN.%j
#SBATCH --error=jureca$RUN.%j
#SBATCH --time=23:59:59
#SBATCH --partition=dc-cpu

apptainer run x3cflux.sif \
python scripts/S6_1-INST-mcmc.py --fml models/S6\ -\ MCMC/syn_perturbed.fml --path measured_poolsizes_run$RUN --n_samples 10000 \
--iteration 0 --rng 50251 --n_chains 4 --n_temps $N_TEMPS --configs freeflux_synthetic_data_measured_poolsizes \
--mcmc_algo GaussianHitAndRun --thinning 100 --step_size $STEP_SIZE --draws_per_exchange_attempt 100
"""

    step_sizes = [2, 1, 0.5]
    n_temps = [8, 32]

    call_script = ''

    for i in range(len(n_temps)):
        for j in range(len(step_sizes)):
            script = template
            script = script.replace('$N_TEMPS', str(n_temps[i]))
            script = script.replace('$STEP_SIZE', str(step_sizes[j]))
            script = script.replace('$RUN', str(i*len(step_sizes)+j))

            with open(f'jureca{i*len(step_sizes)+j}.sh', 'w') as f:
                f.write(script)

            call_script += f'sbatch ./jureca{i*len(step_sizes)+j}.sh\n'

        with open(f'run.sh', 'w') as f:
            f.write(call_script)
