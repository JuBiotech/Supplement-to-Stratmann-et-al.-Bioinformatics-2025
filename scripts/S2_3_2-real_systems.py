#!/usr/bin/env python
# coding: utf-8
import x3cflux
import hopsy
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import time

plt.rcParams["font.size"] = 16

def supremum_norm(x, y):
    errors = np.zeros(len(x))
    for i in range(len(x)):
        errors[i] = np.linalg.norm(x[i] - y[i])
    return errors.max()


def run_benchmark(simulator, reference_solutions, time_grid, samples):
    n_samples = len(samples)
    errors = np.zeros(n_samples)
    times = np.zeros(n_samples)
    for i in range(n_samples):
        start = time.perf_counter()
        label_meas, _ = simulator.compute_measurements(samples[i], time_stamps=time_grid)
        times[i] = time.perf_counter() - start

        n_times = len(time_grid)
        reference_trajectory = []
        simulator_trajectory = []
        for j in range(len(time_grid)):
            x = reference_solutions[i][0][j]
            y = label_meas[j]
            for k in range(len(simulator.measurement_names[0])):
                x = np.concatenate((x, reference_solutions[i][0][j + k * n_times]))
                y = np.concatenate((y, label_meas[j + k * n_times]))
            reference_trajectory.append(x)
            simulator_trajectory.append(y)
        errors[i] = supremum_norm(np.array(reference_trajectory), np.array(simulator_trajectory))
    return errors, times


simulator = x3cflux.create_simulator_from_fml("../models/Syn.fml", "a")
ineq_sys = simulator.parameter_space.inequality_system
simulator.builder.solver.num_max_steps = 1_000_000

n_samples = 10
problem = hopsy.Problem(ineq_sys.matrix, ineq_sys.bound)
samples = x3cflux.run_uniform_sampling(simulator, n_samples)[0]

print("Compute high-accuracy solution (this might take a while)")
simulator.builder.set_solver("sdirk")
simulator.builder.solver.relative_tolerance = 1e-12
simulator.builder.solver.absolute_tolerance = 1e-15
time_grid = np.linspace(0., 610., num=1_000)
reference_solutions = [simulator.compute_measurements(samples[i], time_stamps=time_grid)
                           for i in range(n_samples)]

bdf_errors, bdf_times = [], []
sdirk_errors, sdirk_times = [], []
max_tol = 12
print("Compute solutions with increasing strict tolerances")
for rel_tol in range(3, max_tol):
    print(f"\trel tol = {np.power(10., -rel_tol):.0e}")
    for solver_name in ["bdf", "sdirk"]:
        simulator.builder.set_solver(solver_name)
        simulator.builder.solver.relative_tolerance = np.power(10., -rel_tol)
        simulator.builder.solver.absolute_tolerance = np.power(10., -(rel_tol + 3))
        loc_errors, times = run_benchmark(simulator, reference_solutions, time_grid, samples)
        if solver_name=="bdf":
            bdf_errors.append((loc_errors.mean(), loc_errors.std()))
            bdf_times.append(np.mean(times))
        else:
            sdirk_errors.append((loc_errors.mean(), loc_errors.std()))
            sdirk_times.append(np.mean(times))

fig, ax = plt.subplots(1, 1)
slope, intercept, _, _, _ = stats.linregress(-np.arange(3, max_tol, step=1), np.log10(np.array(bdf_errors)[:, 0]))
x = -np.arange(2, max_tol+1, step=0.1)
ax.plot(np.power(10, x), np.power(10, slope * x + intercept),
        linestyle="--", color="tab:blue")
ax.errorbar(np.power(10., -np.arange(3, max_tol, step=1)), np.array(bdf_errors)[:, 0], np.array(bdf_errors)[:, 1], label=f"BDF ({slope:0.1f} * x - {-intercept:0.1f})",
          marker="o", linestyle="", color="tab:blue")
slope, intercept, _, _, _ = stats.linregress(-np.arange(3, max_tol, step=1), np.log10(np.array(sdirk_errors)[:, 0]))
x = -np.arange(2, max_tol+1, step=0.1)
ax.plot(np.power(10, x), np.power(10, slope * x + intercept),
        linestyle="--", color="tab:orange")
ax.errorbar(np.power(10., -np.arange(3, max_tol, step=1)), np.array(sdirk_errors)[:, 0], np.array(sdirk_errors)[:, 1], label=f"SDIRK ({slope:0.1f} * x + {intercept:0.1f})",
          marker="o", linestyle="", color="tab:orange")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim((6e-12, 2e-3))
ax.set_ylim((6e-11, 2e-2))
ax.set_xticks(np.power(10., -np.arange(3, max_tol, step=2)))
ax.set_yticks(np.power(10., -np.arange(2, 11, step=2)))
ax.set_xlabel("local error tolerance")
ax.set_ylabel("estimated global error")
ax.invert_xaxis()
plt.legend()
plt.tight_layout()
plt.savefig("../data/figure_s2.png", dpi=150)
