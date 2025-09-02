#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import arviz
import x3cflux
import util
import os


def plot(simulator, samples, n_times=100, n_per_row=4, until=1):
    meas_names = simulator.measurement_names[0]
    fig, axs = plt.subplots(
        len(meas_names) // n_per_row + int(len(meas_names) % n_per_row != 0),
        n_per_row,
        figsize=(22, 11),
        sharex=True,
        sharey=True,
    )
    time_grid = np.linspace(0, until, n_times)
    for sample in samples:
        label_meas, _ = simulator.compute_measurements(sample, time_stamps=time_grid)
        for k, meas_name in enumerate(meas_names):
            ax = axs[k // n_per_row, k % n_per_row]
            for datum in simulator.configurations[0].measurements:
                if datum.name == meas_name:
                    meas_datum = datum
            meas_size = len(meas_datum.data.values[0])

            group_meas = np.array(
                label_meas[k * len(time_grid) : (k + 1) * len(time_grid)]
            )
            for j in range(meas_size):
                ax.plot(
                    time_grid,
                    group_meas[:, j],
                    linewidth=1,
                    alpha=0.01,
                    color=plt.rcParams["axes.prop_cycle"].by_key()["color"][j],
                )

    for k, meas_name in enumerate(meas_names):
        ax = axs[k // n_per_row, k % n_per_row]
        curr_ts = simulator.measurement_time_stamps[k]
        for datum in simulator.configurations[0].measurements:
            if datum.name == meas_name:
                meas_datum = datum
        meas_size = len(meas_datum.data.values[0])
        for j in range(meas_size):
            n = meas_datum.specification.weights[j]
            ax.errorbar(
                curr_ts,
                np.array(meas_datum.data.values)[:, j],
                np.array(meas_datum.data.standard_deviations)[:, j],
                label=f"M+{n}",
                linestyle=" ",
                linewidth=3,
                marker=".",
                capsize=5,
                color=plt.rcParams["axes.prop_cycle"].by_key()["color"][j],
            )
        if k == 0:
            ax.legend(ncol=8, bbox_to_anchor=(5, 1.25), fontsize=15)
        ax.set_title(meas_datum.metabolite_name, fontsize=20)
        xticks = [0, 25, 50, 75, 100, 125, 150]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks, fontsize=15)
    fig.text(0.5, 0.05, "time [h]", ha="center", va="center", fontsize=18)
    fig.text(
        0.09,
        0.5,
        "fractional enrichment",
        ha="center",
        va="center",
        rotation="vertical",
        fontsize=18,
    )
    plt.xlim((0, until))
    plt.ylim((0, 1))


util.print_box(f"Executing {os.path.basename(__file__)}")
x3cflux.logging.level = 0
simulator = x3cflux.create_simulator_from_fml("../models/Syn.fml", "bayes")
mle = x3cflux.get_parameters(
    simulator.parameter_space, simulator.configurations[0].parameter_entries
)

samples = np.zeros((4, 24000, 60))

for i in range(4):
    for j in range(2):
        fn = f"../data/results_{i}_{j}.npz"
        print(f"Reading data from {fn}")
        samples[i, j * 12000 : (j + 1) * 12000] = np.load(fn)["samples"][0]

idx = [2, 5, 1, 26, 54]
arviz.plot_pair(
    arviz.from_dict(
        {
            simulator.parameter_space.free_parameter_names[i]: samples[
                :, :, i
            ].flatten()
            for i in idx
        }
    ),
    kind="kde",
    textsize=25,
    marginals=True,  # marginal_kwargs={"quantiles": [0.05, 0.95]},
    var_names=[simulator.parameter_space.free_parameter_names[i] for i in idx],
    reference_values={
        simulator.parameter_space.free_parameter_names[i]: mle[i] for i in idx
    },
    reference_values_kwargs=dict(
        marker="X",
        markersize=25,
        color="darkgrey",
        markerfacecolor="tab:red",
        markeredgecolor="grey",
    ),
)
plt.tight_layout()
filename = "../out/figure_s10"
print(f"saving figure to {filename}")
plt.savefig(filename + ".png", dpi=300)
plt.savefig(filename + ".svg")

print("computing posterior predictives")
plot(
    simulator,
    samples[0, np.random.randint(0, 24_000, 500)],
    n_times=100,
    until=150,
    n_per_row=5,
)
filename = "../out/figure_s11"
print(f"saving figure to {filename}")
plt.savefig(filename + ".png", dpi=300)
plt.savefig(filename + ".svg")
