#!/usr/bin/env python
# coding: utf-8

import x3cflux
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.stats as stats
import os
import util


def benchmark_ed(
    network_data,
    configuration,
    batch_sizes,
    substr_mix_collections,
    parallel_batch_size,
):
    parallel_times = np.zeros(len(batch_sizes))
    sequential_times = np.zeros(len(batch_sizes))
    for i, batch_size in enumerate(batch_sizes):
        print(f"Computing ED for {batch_size} mixtures")
        loc_times = []

        substrates_mixtures, _ = x3cflux.compute_mixture_samples(
            substr_mix_collections[:1], batch_size
        )
        fixed_substrates = []
        for substrate in configuration.substrates:
            if substrate.metabolite_name not in [
                mix.metabolite_name for mix in substr_mix_collections[:1]
            ]:
                fixed_substrates.append(substrate)

        size = batch_size if batch_size < parallel_batch_size else parallel_batch_size
        print(f"\tbatch-wise")
        for j in range(int(batch_size / (parallel_batch_size + 1)) + 1):
            simulator = x3cflux.create_simulator_from_inputs(
                network_data,
                configuration,
                substrates_mixtures[j * size : min(batch_size, (j + 1) * size)],
                fixed_substrates,
            )
            last_params = x3cflux.get_parameters(
                simulator.parameter_space, configuration.parameter_entries
            )
            start = time.perf_counter()
            simulator.compute_multi_jacobians(last_params)
            loc_times.append(time.perf_counter() - start)
        parallel_times[i] = np.sum(loc_times)

        loc_times = np.zeros(batch_size)
        print(f"\tsequential")
        for j in range(batch_size):
            mix_meas_config = x3cflux.MeasurementConfiguration(
                configuration.name if configuration.name else "",
                configuration.comment,
                configuration.stationary,
                substrates_mixtures[j] + fixed_substrates,
                configuration.measurements,
                configuration.net_flux_constraints,
                configuration.exchange_flux_constraints,
                configuration.pool_size_constraints,
                configuration.parameter_entries,
            )

            simulator = x3cflux.create_simulator_from_data(
                network_data, mix_meas_config
            )
            last_params = x3cflux.get_parameters(
                simulator.parameter_space, configuration.parameter_entries
            )
            start = time.perf_counter()
            simulator.compute_multi_jacobians(last_params)
            loc_times[j] += time.perf_counter() - start
        sequential_times[i] = np.sum(loc_times)

        print(
            f"\t{parallel_times[i]:0.2f}s (batch-wise) vs {sequential_times[i]:0.2f}s (sequential), speedup = {sequential_times[i] / parallel_times[i]:0.2f}"
        )

    return parallel_times, sequential_times


def plot_ed_benchmark(parallel_times, sequential_times, batch_sizes, postfix):
    slope, intercept, _, _, _ = stats.linregress(
        np.log10(batch_sizes), np.log10(parallel_times)
    )
    x = np.arange(0, 4.3, step=0.1)
    plt.plot(
        np.power(10, x),
        np.power(10, slope * x + intercept),
        linestyle="--",
        color="tab:blue",
    )
    plt.errorbar(
        batch_sizes,
        parallel_times,
        label=f"batch-wise ({slope:0.1f} * x - {-intercept:0.1f})",
        marker="o",
        ls="",
        color="tab:blue",
    )
    slope, intercept, _, _, _ = stats.linregress(
        np.log10(batch_sizes), np.log10(sequential_times)
    )
    x = np.arange(0, 4.3, step=0.1)
    plt.plot(
        np.power(10, x),
        np.power(10, slope * x + intercept),
        linestyle="--",
        color="tab:orange",
    )
    plt.errorbar(
        batch_sizes,
        sequential_times,
        label=f"sequential ({slope:0.1f} * x - {-intercept:0.1f})",
        marker="o",
        ls="",
        color="tab:orange",
    )
    plt.xlabel("number of designs")
    plt.ylabel("$\\leftarrow$ runtime [s]")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim((1e1, 2e4))
    plt.legend()
    plt.tight_layout()
    filename = "../out/figure_s06"
    print(f"saving figure to {filename}")
    plt.savefig(filename + ".png", dpi=150)
    plt.savefig(filename + ".svg")


util.print_box(f"Executing {os.path.basename(__file__)}")
plt.rcParams["font.size"] = 16
x3cflux.logging.level = 0

data = x3cflux.FluxMLParser().parse("../models/EC.fml")
substr_mix_collections = x3cflux.parse_tracer_mixtures(
    "../models/mixture.fml"
)
batch_sizes_ed = (
    [10, 25, 50, 75, 100]
    + [200 * i for i in range(1, 6)]
    + [2000 * i for i in range(1, 6)]
)

config_names = [i.name for i in data.configurations]
config_index = config_names.index("a_STAT")

parallel_times, sequential_times = benchmark_ed(
    data.network_data,
    data.configurations[config_index],
    batch_sizes_ed,
    substr_mix_collections,
    75,
)

plot_ed_benchmark(parallel_times[1:], sequential_times[1:], batch_sizes_ed[1:], "stat")
