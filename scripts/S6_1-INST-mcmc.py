import os

os.environ['OPENBLAS_NUM_THREADS'] = '10'
from datetime import datetime
import hopsy
import numpy as np
import time as timer
import click
import x3cflux
import logging
from pathlib import Path
import pickle
from multiprocessing import set_start_method


class Wrapper(x3cflux.HopsyModel):
    def __init__(self, *args, **kwargs):
        self.x3cflux_time = 0
        self.n_eval = 0
        super().__init__(*args, **kwargs)

    def log_density(self, x):
        start = timer.perf_counter()
        ld = super().log_density(x)
        self.x3cflux_time += timer.perf_counter() - start
        self.n_eval += 1
        return ld

    def reset_statistics(self):
        self.x3cflux_time = 0
        self.n_eval = 0

def setup(**kwargs):
    configs = kwargs.pop('configs').split(',')
    rng_seed = kwargs.pop('rng')
    n_chains = kwargs.pop('n_chains')
    n_temps = kwargs.pop('n_temps')
    draws_per_exchange_attempt = kwargs.pop('draws_per_exchange_attempt')
    step_size = kwargs.pop('step_size')
    warm_up = kwargs.pop('warm_up')
    mcmc_algo = kwargs.pop('mcmc_algo')
    fml = kwargs.pop('fml')
    if mcmc_algo.upper() == "ADAPTIVEMETROPOLIS":
        proposal_type = hopsy.AdaptiveMetropolisProposal
    elif mcmc_algo.upper() == "GAUSSIANHITANDRUN":
        proposal_type = hopsy.GaussianHitAndRunProposal
    else:
        raise RuntimeError(f"unkown algo type {mcmc_algo}")
    logging.info(f'using proposal {proposal_type}')

    logging.info(f'configs are {str(configs)}')
    simulator = x3cflux.create_simulator_from_fml(fml, config_names=configs)
    ineq_sys = simulator.parameter_space.inequality_system
    if hopsy.is_polytope_empty(ineq_sys.matrix, ineq_sys.bound):
        message = "flux space is empty/infeasible."
        logging.error(message)
        raise RuntimeError(message)
    problem = hopsy.Problem(ineq_sys.matrix, ineq_sys.bound, Wrapper(simulator))
    problem = hopsy.add_box_constraints(problem, lower_bound=-30, upper_bound=250)
    problem.starting_point = x3cflux.get_parameters(simulator.parameter_space,
                                                         simulator.configurations[0].parameter_entries)
    logging.info(f"starting rounding using {hopsy.LP().settings.backend}")
    problem = hopsy.round(problem)
    logging.info("finished rounding")
    if proposal_type == hopsy.AdaptiveMetropolisProposal:
        proposal = proposal_type(problem, warm_up=warm_up)
    else:
        proposal = proposal_type(problem)

    mcs = [
        hopsy.MarkovChain(
            proposal=proposal,
            problem=problem,
        )
        for _ in range(n_chains)
    ]
    sync_rngs = [hopsy.RandomNumberGenerator(seed=rng_seed + 1, stream=r) for r in range(n_chains)]
    if (n_temps != 0):
        temperature_ladder = [1.0 - float(n) / (n_temps - 1) for n in range(n_temps)]
        # Creates one parallel tempering ensemble for each replicate.
        # Each ensemble will have len(temperature_ladder) chains.
        mcs = hopsy.create_py_parallel_tempering_ensembles(
            markov_chains=mcs,
            temperature_ladder=temperature_ladder,
            sync_rngs=sync_rngs,
            draws_per_exchange_attempt=draws_per_exchange_attempt,
        )

    for mc in mcs:
        if hasattr(mc.proposal, "stepsize"):
            coldness = 1
            if hasattr(mc, "coldness"):
                coldness = mc.coldness
            max_step_size = 10
            mc.proposal.stepsize = max_step_size + coldness * (step_size - max_step_size)
            logging.info(f'set stepsize of mc with coldness {coldness} to {mc.proposal.stepsize}')

    rngs = [hopsy.RandomNumberGenerator(rng_seed, stream=r) for r, _ in enumerate(mcs)]

    return mcs, rngs, sync_rngs


def load(path,
         iteration_to_load,
         n_temps,
         draws_per_exchange_attempt,
         n_chains,
         ):
    with open(Path(path) / Path(f"chains_{iteration_to_load}.plk"), 'rb') as f:
        mcs = pickle.load(f)
    with open(Path(path) / Path(f"rngs_{iteration_to_load}.plk"), 'rb') as f:
        rngs = pickle.load(f)
    with open(Path(path) / Path(f"sync_rngs_{iteration_to_load}.plk"), 'rb') as f:
        sync_rngs = pickle.load(f)

    if (n_temps != 0):
        temperature_ladder = [1.0 - float(n) / (n_temps - 1) for n in range(n_temps)]
        # Creates one parallel tempering ensemble for each replicate.
        # Each ensemble will have len(temperature_ladder) chains.
        dummy_mcs = [hopsy.MarkovChain(mcs[0].problem, starting_point=mcs[0].proposal.state) for _ in range(n_chains)]
        _mcs = hopsy.create_py_parallel_tempering_ensembles(
            markov_chains=dummy_mcs,
            temperature_ladder=temperature_ladder,
            sync_rngs=sync_rngs,
            draws_per_exchange_attempt=draws_per_exchange_attempt,
        )
        for i in range(len(_mcs)):
            _mcs[i].markov_chain = mcs[i]
        mcs = _mcs

    return mcs, rngs, sync_rngs


@click.command()
@click.option("--fml", type=str, required=True, help="path to fml")
@click.option("--configs", type=str, default='default',
              help="fml config. Give comma separated list for parallel experiments")
@click.option("--path", type=str, required=True, help="where to write results and checkpoints")
@click.option("--n_samples", type=int, required=True, help="how many samples for this iteration")
@click.option("--thinning", type=int, default=1, help="thinning constant")
@click.option("--iteration", type=int, required=True, help="what iteration is currently being sampled")
@click.option("--rng", type=int, required=True, help="make sure to use a different seed on every run")
@click.option("--n_chains", type=int, required=True, help="how many chains to use")
@click.option("--n_temps", type=int, default=0,
              help="If positive, parallel tempering will be used with this many temps for each chain")
@click.option("--draws_per_exchange_attempt", type=int, default=1000, help="for parallel tempering")
@click.option("--step_size", type=float, default=1, help="step size")
@click.option("--progress_bar", type=bool, default=False, help="whether to show hopsy progressbar during sampling")
@click.option("--warm_up", type=int, default=10000, help="param for warm up of adaptive metropolis")
@click.option("--mcmc_algo", type=str, default="GaussianHitAndRun", help="mcmc algorithm")
def command(
        fml: str,
        configs: str,
        path: str,
        n_samples: int,
        thinning: int,
        iteration: int,
        rng: int,
        n_chains: int,
        n_temps: int,
        draws_per_exchange_attempt,
        step_size: float,
        progress_bar: bool,
        warm_up: int,
        mcmc_algo: str,
):
    set_start_method("spawn")
    logging.basicConfig(encoding='utf-8', level=logging.INFO)
    Path(path).mkdir(exist_ok=True)
    with open(Path(path) / Path('run_params.yml'), 'a+') as f:
        lines = [
            f'start_date: {str(datetime.today())}\n',
            f'fml: {str(fml)}\n',
            f'configs: {str(configs)}\n',
            f'path: {str(path)}\n',
            f'n_samples: {str(n_samples)}\n',
            f'thinning: {str(thinning)}\n',
            f'iteration: {str(iteration)}\n',
            f'rng: {str(rng)}\n',
            f'n_chains: {str(n_chains)}\n',
            f'n_temps: {str(n_temps)}\n',
            f'draws_per_exchange_attempt: {str(draws_per_exchange_attempt)}\n',
            f'step_size: {str(step_size)}\n',
            f'progress_bar: {str(progress_bar)}\n',
            f'warm_up: {str(warm_up)}\n',
            f'mcmc_algo: {str(mcmc_algo)}\n',
            f'\n'
        ]

        f.writelines(lines)
    if iteration == 0:
        # set up problem
        mcs, rngs, sync_rngs = setup(
            fml=fml,
            configs=configs,
            rng=rng,
            n_chains=n_chains,
            n_temps=n_temps,
            draws_per_exchange_attempt=draws_per_exchange_attempt,
            step_size=step_size,
            mcmc_algo=mcmc_algo,
            warm_up=warm_up,
        )
        logging.info('setup markov chains')
    else:
        # load previous iteration
        mcs, rngs, sync_rngs = load(path,
                                    iteration - 1,
                                    n_temps=n_temps,
                                    draws_per_exchange_attempt=draws_per_exchange_attempt,
                                    n_chains=n_chains
                                    )
        logging.info('loaded previous iteration')

    for mc in mcs:
        # resets statistics on model calls, so that we only count calls made during sampling
        logging.info('resetting statistics')
        mc.problem.model.reset_statistics()
    logging.info('start time', timer.ctime())
    start = timer.perf_counter()

    acc, samples = hopsy.sample(mcs,
                                rngs,
                                n_samples=n_samples,
                                thinning=thinning,
                                n_procs=len(mcs),
                                progress_bar=progress_bar
                                )
    end = timer.perf_counter()

    runtime = end - start
    logging.info(f'runtime: {runtime}s')
    time_per_sample = runtime / (n_samples * thinning)
    logging.info(f"time per sample {time_per_sample}")

    if n_temps != 0:
        mcs = [mc.markov_chain for mc in mcs]
        for mc in mcs:
            print('n_eval', mc.problem.model.n_eval, 'time', mc.problem.model.x3cflux_time)
    else:
        for mc in mcs:
            print('n_eval', mc.problem.model.n_eval, 'time', mc.problem.model.x3cflux_time)

    with open(Path(path) / Path(f"chains_{iteration}.plk"), 'wb') as f:
        pickle.dump(mcs, f)

    with open(Path(path) / Path(f"rngs_{iteration}.plk"), 'wb') as f:
        pickle.dump(rngs, f)

    with open(Path(path) / Path(f"sync_rngs_{iteration}.plk"), 'wb') as f:
        pickle.dump(sync_rngs, f)

    with open(Path(path) / Path(f'runtime_{iteration}.csv'), 'w') as f:
        f.write("runtime: " + str(runtime) + " seconds\n")
        f.write("time per iteration: " + str(time_per_sample) + " seconds / iteration \n")
        for i, mc in enumerate(mcs):
            f.write(f"chain {i} (coldness={mc.coldness}) 13cflux time: " + str(mc.problem.model.x3cflux_time) + " seconds \n")
        for i, mc in enumerate(mcs):
            f.write(f"chain {i} (coldness={mc.coldness}) 13cflux it/s: " + str(mc.problem.model.n_eval/mc.problem.model.x3cflux_time) + "\n")

    np.savez_compressed(
        str(Path(path) / Path(f'sampling_results_{iteration}')),
        samples=samples, acceptance_rates=acc
    )


if __name__ == '__main__':
    command()
