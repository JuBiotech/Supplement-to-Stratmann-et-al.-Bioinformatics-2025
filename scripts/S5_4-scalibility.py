#!/usr/bin/env python
# coding: utf-8

import x3cflux
import numpy as np
import pandas as pd
from time import perf_counter
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

C_X3 = '#023d6b'

data = x3cflux.FluxMLParser().parse("../models/EC.fml")
num_runs = 100
tols = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
times = pd.DataFrame(index=range(num_runs), columns=tols)
config_names = [i.name for i in data.configurations]
config_index = config_names.index("a_INST")

simulator_base = x3cflux.create_simulator_from_data(data.network_data, data.configurations[config_index])
params = x3cflux.get_parameters(simulator_base.parameter_space, data.configurations[config_index].parameter_entries)

for tol in tols:
    simulator_base.builder.solver.relative_tolerance = tol
    simulator_base.builder.solver.absolute_tolerance = tol*1e-3

    for i in range(num_runs):
        before = perf_counter()
        r = simulator_base.compute_measurements(params)
        after = perf_counter()
        time = after-before
        times.loc[i, tol] = time

xs = -np.log10(times.columns.to_numpy(dtype=float))
means = np.array(times.mean(axis=0).values, dtype=float)
ys = np.log10(means)
model = LinearRegression()
model.fit(xs.reshape(-1, 1),  ys)
x3_y=model.predict([[1],[14]])

times.to_csv("../out/figure_s9_raw.csv")

plt.errorbar(xs, times.mean(axis=0).values, yerr=times.std(axis=0).values, fmt='o', color=C_X3, label="13CFLUX(v3)")
plt.yscale('log')
plt.legend()
locs = [2,4,6,8,10,12]
plt.xticks(locs,[f"$10^{{-{i}}}$" for i in locs])
xlims=plt.xlim()
ylims=plt.ylim()
plt.plot([1,14],np.power(10, x3_y), '--', color=C_X3)
plt.xlim(xlims)
plt.ylim(ylims)
plt.xlabel(r'$tol_{rel}$')
plt.ylabel('← runtime [ms]')
plt.tight_layout()
plt.savefig('../out/figure_s9.png', dpi=150)
plt.savefig('../out/figure_s9.svg')
