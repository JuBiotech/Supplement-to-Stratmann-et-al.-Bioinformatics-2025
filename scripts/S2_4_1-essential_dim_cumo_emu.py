#!/usr/bin/env python
# coding: utf-8

import x3cflux
import hopsy
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.stats as stats
import util
import os


def time_measurement_simulation(measurements, network_data, exp_config, samples, method="emu"):
    new_exp_config = x3cflux.MeasurementConfiguration(
        exp_config.name,
        exp_config.comment,
        exp_config.stationary,
        exp_config.substrates,
        measurements,
        exp_config.net_flux_constraints,
        exp_config.exchange_flux_constraints,
        exp_config.pool_size_constraints,
        exp_config.parameter_entries
    )

    new_simulator = x3cflux.create_simulator_from_data(network_data, new_exp_config, sim_method=method)
    num_unknowns = 0
    cascaded_system = new_simulator.builder.build(new_simulator.parameter_space.compute_parameters(samples[0]))
    num_levels = len(cascaded_system.solve())
    for i in range(num_levels):
        if exp_config.stationary:
            num_unknowns += cascaded_system.get_level_system(i).rhs.shape[0]
        else:
            num_unknowns += cascaded_system.get_level_system(i).initial_value.shape[0]

    loc_times = np.zeros(samples.shape[1])
    for k in range(samples.shape[1]):
        start = time.perf_counter()
        new_simulator.compute_loss(samples[k])
        loc_times[k] = (time.perf_counter() - start) * 1000.

    return num_unknowns, loc_times.mean(), loc_times.std()


def time_increasing_measurement_simulations(network_data, config, samples):
    meas_list = list(sorted(config.measurements, key=lambda x: x.name))
    n_param_meas = np.sum([1 if meas.name.startswith("f") else 0 for meas in meas_list])

    times = []
    for type in ["Tandem MS", "MS"]:
        print(f"\tFor {type} measurements")
        cumo_times = np.zeros((len(config.measurements) - n_param_meas + 1, 3))
        emu_times = np.zeros((len(config.measurements) - n_param_meas + 1, 3))
        cumo_times[0] = time_measurement_simulation(meas_list[:n_param_meas], network_data, config, samples, "cumomer")
        emu_times[0] = time_measurement_simulation(meas_list[:n_param_meas], network_data, config, samples, "emu")

        incr_meas_list = meas_list[:n_param_meas]
        for i in range(n_param_meas, len(config.measurements)):
            curr_meas = meas_list[i]
            print(f"\t\tAdding measurement {i}/{len(config.measurements)}: {curr_meas.name}")

            if type == "Tandem MS":
                n = curr_meas.num_atoms
                m = n - 1
                spec = x3cflux.MSMSSpecification("1" * n, "0" + "1" * m,
                                                 sum([m * [i] for i in range(n + 1)], []),
                                                 sum([list(range(m)) for _ in range(n + 1)], []))
                ts = curr_meas.data.time_stamps
                data = x3cflux.MeasurementDataSet(ts,
                                                  len(ts) * [np.zeros(n * m)],
                                                  len(ts) * [np.zeros(n * m)], [])
                incr_meas_list.append(x3cflux.MSMSMeasurement(curr_meas.name, curr_meas.auto_scalable,
                                                              curr_meas.metabolite_name, curr_meas.num_atoms, spec,
                                                              data))
            else:
                incr_meas_list.append(curr_meas)

            cumo_times[i - n_param_meas + 1] = time_measurement_simulation(incr_meas_list, network_data, config, samples, "cumomer")
            emu_times[i - n_param_meas + 1] = time_measurement_simulation(incr_meas_list, network_data, config, samples, "emu")
        times.append((cumo_times, emu_times))

    return times


def plot_solution_times(msms_times_stat, ms_times_stat, msms_times_inst, ms_times_inst):
    figsize = plt.rcParams["figure.figsize"]
    fig, axs = plt.subplots(2, 2, figsize=(figsize[0]*2, figsize[1]*2), sharex="col", sharey="row")

    # IST
    axs[0, 0].set_title("MS")
    axs[0, 0].errorbar(ms_times_stat[0][1:, 0], ms_times_stat[0][1:, 1], yerr=ms_times_stat[0][1:, 2], label="cumomer",
                    marker="o", ls="", color="tab:blue")
    result = stats.linregress(ms_times_stat[0][1:, 0], ms_times_stat[0][1:, 1])
    print(f"cumomer: r-squared = {result.rvalue}")
    x = np.linspace(ms_times_stat[0][1:, 0].min()-200, ms_times_stat[0][1:, 0].max()+100, 100)
    axs[0, 0].plot(x, result.slope * x + result.intercept, ls="--", color="tab:blue")
    axs[0, 0].errorbar(ms_times_stat[1][1:, 0], ms_times_stat[1][1:, 1], yerr=ms_times_stat[1][1:, 2], label="EMU",
                    marker="o", ls="", color="tab:orange")
    result = stats.linregress(ms_times_stat[1][1:,0], ms_times_stat[1][1:, 1])
    print(f"EMU: r-squared = {result.rvalue}")
    x = np.linspace(ms_times_stat[1][1:, 0].min()-50, ms_times_stat[1][1:, 0].max()+100, 100)
    axs[0, 0].plot(x, result.slope * x + result.intercept, ls="--", color="tab:orange")
    axs[0, 0].set_ylabel("$\\leftarrow$ IST runtime [ms]")
    axs[0, 0].set_xlim((200, 1500))
    #axs[0, 0].set_ylim((0, 3))
    axs[0, 0].legend()

    axs[0, 1].set_title("MS/MS")
    axs[0, 1].errorbar(msms_times_stat[0][1:, 0], msms_times_stat[0][1:, 1], yerr=msms_times_stat[0][1:, 2], label="cumomer",
                    marker="o", ls="", color="tab:blue")
    result = stats.linregress(msms_times_stat[0][1:, 0], msms_times_stat[0][1:, 1])
    print(f"cumomer: r-squared = {result.rvalue}")
    x = np.linspace(msms_times_stat[0][1:, 0].min()-100, msms_times_stat[0][1:, 0].max()+100, 100)
    axs[0, 1].plot(x, result.slope * x + result.intercept, ls="--", color="tab:blue")
    axs[0, 1].errorbar(msms_times_stat[1][1:, 0], msms_times_stat[1][1:, 1], yerr=msms_times_stat[1][1:, 2], label="EMU",
                    marker="o", ls="", color="tab:orange")
    result = stats.linregress(msms_times_stat[1][1:, 0], msms_times_stat[1][1:, 1])
    print(f"EMU: r-squared = {result.rvalue}")
    x = np.linspace(msms_times_stat[1][1:, 0].min()-100, msms_times_stat[1][1:, 0].max()+100, 100)
    axs[0, 1].plot(x, result.slope * x + result.intercept, ls="--", color="tab:orange")
    axs[0, 1].set_xlim((200, 1950))
    #axs[0, 1].set_ylim((0, 3))

    # INST
    axs[1, 0].errorbar(ms_times_inst[0][1:, 0], ms_times_inst[0][1:, 1], yerr=ms_times_inst[0][1:,2], label="cumomer",
                    marker="o", ls="", color="tab:blue")
    result = stats.linregress(ms_times_inst[0][1:, 0], ms_times_inst[0][1:, 1])
    print(f"cumomer: r-squared = {result.rvalue}")
    x = np.linspace(ms_times_inst[0][1:, 0].min()-200, ms_times_inst[0][1:, 0].max()+100, 100)
    axs[1, 0].plot(x, result.slope * x + result.intercept, ls="--", color="tab:blue")
    axs[1, 0].errorbar(ms_times_inst[1][1:, 0], ms_times_inst[1][1:, 1], yerr=ms_times_inst[1][1:,2], label="EMU",
                    marker="o", ls="", color="tab:orange")
    result = stats.linregress(ms_times_inst[1][1:, 0], ms_times_inst[1][1:, 1])
    print(f"EMU: r-squared = {result.rvalue}")
    x = np.linspace(ms_times_inst[1][1:, 0].min()-50, ms_times_inst[1][1:, 0].max()+100, 100)
    axs[1, 0].plot(x, result.slope * x + result.intercept, ls="--", color="tab:orange")
    axs[1, 0].set_xlabel("essential dimension")
    axs[1, 0].set_ylabel("$\\leftarrow$ INST runtime [ms]")
    #axs[1, 0].set_ylim((0, 200))
    axs[1, 0].set_xlabel("essential dimension")

    axs[1, 1].errorbar(msms_times_inst[0][1:, 0], msms_times_inst[0][1:, 1], yerr=msms_times_inst[0][1:, 2], label="cumomer",
                    marker="o", ls="", color="tab:blue")
    result = stats.linregress(msms_times_inst[0][1:, 0], msms_times_inst[0][1:, 1])
    print(f"cumomer: r-squared = {result.rvalue}")
    x = np.linspace(msms_times_inst[0][1:, 0].min()-100, msms_times_inst[0][1:, 0].max()+100, 100)
    axs[1, 1].plot(x, result.slope * x + result.intercept, ls="--", color="tab:blue")
    axs[1, 1].errorbar(msms_times_inst[1][1:, 0], msms_times_inst[1][1:, 1], yerr=msms_times_inst[1][1:,2], label="EMU",
                    marker="o", ls="", color="tab:orange")
    result = stats.linregress(msms_times_inst[1][1:, 0], msms_times_inst[1][1:, 1])
    print(f"EMU: r-squared = {result.rvalue}")
    x = np.linspace(msms_times_inst[1][1:, 0].min()-100, msms_times_inst[1][1:,0].max()+100, 100)
    axs[1, 1].plot(x, result.slope * x + result.intercept, ls="--", color="tab:orange")
    #axs[1, 1].set_ylim((0, 200))
    axs[1, 1].set_xlabel("essential dimension")

    #axs[1].bar([1, 4], [msms_times[0][:, 0].max(), ms_times[0][:, 0].max()], width=1, label="cumomer")
    #axs[1].bar([2, 5], [msms_times[1][:, 0].max(), ms_times[1][:, 0].max()], width=1, label="emu")
    #axs[1].set_xticks([1.5, 4.5], ["Tandem-MS", "MS"])
    #lower = min(ms_times[0][1:, 0].min(), ms_times[1][1:, 0].min()) - 10
    #upper = max(ms_times[0][1:, 0].max(), ms_times[1][1:, 0].max()) + 10
    #axs[1].set_ylabel("essential dimension")

    plt.tight_layout()
    filename = "../out/figure_s05"
    print(f"saving to {filename}")
    plt.savefig(filename + ".png", dpi=150)
    plt.savefig(filename + ".svg")


def plot_total_essential_dimension(msms_times_stat, ms_times_stat, msms_times_inst, ms_times_inst):
    figsize = plt.rcParams["figure.figsize"]
    fig, axs = plt.subplots(2, 2, figsize=(figsize[0]*2, figsize[1]*2), sharex="col", sharey="row")

    axs[0, 0].set_title("MS")
    for i in range(1, ms_times_stat[0].shape[0]):
        axs[0, 0].bar(3*i, ms_times_stat[0][i, 1], width=1, color="tab:blue", label="cumomer")
        axs[0, 0].bar(3*i + 1, ms_times_stat[1][i, 1], width=1, color="tab:orange", label="EMU")
        if i == 1:
            axs[0, 0].legend()
    axs[0, 0].errorbar([3*i for i in range(1, ms_times_stat[0].shape[0])], ms_times_stat[0][1:, 1], yerr=ms_times_stat[0][1:, 2], c="black", ls="")
    axs[0, 0].errorbar([3*i+1 for i in range(1, ms_times_stat[1].shape[0])], ms_times_stat[1][1:, 1], yerr=ms_times_stat[1][1:, 2], c="black", ls="")
    axs[0, 0].set_xticks([])
    axs[0, 0].set_ylabel("$\\leftarrow$ IST runtime [ms]")

    axs[0, 1].set_title("MS/MS")
    for i in range(1, msms_times_stat[0].shape[0]):
        axs[0, 1].bar(3*i, msms_times_stat[0][i, 1], width=1, color="tab:blue", label="cumomer")
        axs[0, 1].bar(3*i + 1, msms_times_stat[1][i, 1], width=1, color="tab:orange", label="EMU")
    axs[0, 1].errorbar([3*i for i in range(1, msms_times_stat[0].shape[0])], msms_times_stat[0][1:, 1], yerr=msms_times_stat[0][1:, 2], c="black", ls="")
    axs[0, 1].errorbar([3*i+1 for i in range(1, msms_times_stat[1].shape[0])], msms_times_stat[1][1:, 1], yerr=msms_times_stat[1][1:, 2], c="black", ls="")
    axs[0, 1].set_xticks([])

    for i in range(1, ms_times_inst[0].shape[0]):
        axs[1, 0].bar(3*i, ms_times_inst[0][i, 1], width=1, color="tab:blue", label="cumomer")
        axs[1, 0].bar(3*i + 1, ms_times_inst[1][i, 1], width=1, color="tab:orange", label="EMU")
    axs[1, 0].errorbar([3*i for i in range(1, ms_times_inst[0].shape[0])], ms_times_inst[0][1:, 1], yerr=ms_times_inst[0][1:, 2], c="black", ls="")
    axs[1, 0].errorbar([3*i+1 for i in range(1, ms_times_inst[1].shape[0])], ms_times_inst[1][1:, 1], yerr=ms_times_inst[1][1:, 2], c="black", ls="")
    axs[1, 0].set_xticks([])
    axs[1, 0].set_ylabel("$\\leftarrow$ INST runtime [ms]")

    for i in range(1, msms_times_inst[0].shape[0]):
        axs[1, 1].bar(3*i, msms_times_inst[0][i, 1], width=1, color="tab:blue", label="cumomer")
        axs[1, 1].bar(3*i + 1, msms_times_inst[1][i, 1], width=1, color="tab:orange", label="EMU")
    axs[1, 1].errorbar([3*i for i in range(1, msms_times_inst[0].shape[0])], msms_times_inst[0][1:, 1], yerr=msms_times_inst[0][1:, 2], c="black", ls="")
    axs[1, 1].errorbar([3*i+1 for i in range(1, msms_times_inst[1].shape[0])], msms_times_inst[1][1:, 1], yerr=msms_times_inst[1][1:, 2], c="black", ls="")
    axs[1, 1].set_xticks([])

    fig.supxlabel("measurement configurations", y=0.03)
    plt.tight_layout()
    filename = "../out/figure_s04"
    print(f"saving to {filename}")
    plt.savefig(filename + ".png", dpi=150)
    plt.savefig(filename + ".svg")


util.print_box(f"Executing {os.path.basename(__file__)}")
plt.rcParams["font.size"] = 16
x3cflux.logging.level = 0

data = x3cflux.FluxMLParser().parse("../models/EC.fml")
config_names = [i.name for i in data.configurations]
config_index_stat = config_names.index("a_STAT")
config_index_inst = config_names.index("a_INST")

n_samples = 1000

print("Evaluating IST Problem")
simulator = x3cflux.create_simulator_from_data(data.network_data, data.configurations[config_index_stat])
ineq_sys = simulator.parameter_space.inequality_system
problem = hopsy.Problem(ineq_sys.matrix, ineq_sys.bound)
problem = hopsy.add_box_constraints(problem, 12 * [-200.] + 23 * [1e-6], 12 * [200.] + 23 * [1000])
print(f"\tDrawing {n_samples} with hopsy")
samples = hopsy.sample(hopsy.MarkovChain(problem, hopsy.UniformCoordinateHitAndRunProposal),
                       hopsy.RandomNumberGenerator(42),
                       n_samples=n_samples, thinning=int(60 ** 2 // 6))[1].reshape((n_samples, -1))

msms_times_stat, ms_times_stat = time_increasing_measurement_simulations(data.network_data, data.configurations[config_index_stat], samples)

print("Evaluating INST Problem")
simulator = x3cflux.create_simulator_from_data(data.network_data, data.configurations[config_index_inst])
ineq_sys = simulator.parameter_space.inequality_system
problem = hopsy.Problem(ineq_sys.matrix, ineq_sys.bound)
problem = hopsy.add_box_constraints(problem,
                                    12 * [-200.] + 23 * [1e-6] + 56 * [1e-6],
                                    12 * [200.] + 23 * [1000] + 56 * [1e-1])
print(f"\tDrawing {n_samples} with hopsy")
samples = hopsy.sample(hopsy.MarkovChain(problem, hopsy.UniformCoordinateHitAndRunProposal),
                       hopsy.RandomNumberGenerator(42),
                       n_samples=n_samples, thinning=int(60 ** 2 // 6))[1].reshape((n_samples, -1))

msms_times_inst, ms_times_inst = time_increasing_measurement_simulations(data.network_data, data.configurations[config_index_inst], samples)

plot_total_essential_dimension(msms_times_stat, ms_times_stat, msms_times_inst, ms_times_inst)

plot_solution_times(msms_times_stat, ms_times_stat, msms_times_inst, ms_times_inst)
