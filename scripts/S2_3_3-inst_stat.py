import x3cflux
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import sys
import os
import util


def mask_to_num(mask: str):
    changes = []
    for i in np.arange(1,len(mask)):
        if mask[i]!=mask[i-1]:
            changes.append(i)
    state=False
    last_on = None
    out=''
    if mask[0]=='1':
        if not changes:
            return out
        state=True
        last_on=0
        out='1'
    for change in changes:
        if state:
            if last_on < i:
                out += f"-{change}"
            last_on = i
            state = False
        else:
            if not last_on is None:
                out += ","
            out += f"{change+1}"
            state = True
    if state:
        out += f"-{len(mask)}"
    return f"[{out}]"
    
util.print_box(f"Executing {os.path.basename(__file__)}")
plt.rcParams['font.size'] = 16
C_X3 = '#023d6b'
x3cflux.logging.level = 0

data = x3cflux.FluxMLParser().parse("../models/EC.fml")
config_names = [i.name for i in data.configurations]

# INST
config = data.configurations[config_names.index("a_INST")]
simulator_inst = x3cflux.create_simulator_from_data(data.network_data, config)
simulator_inst.builder.solver.relative_tolerance = 1e-6
simulator_inst.builder.solver.absolute_tolerance = 1e-9
params = x3cflux.get_parameters(simulator_inst.parameter_space, config.parameter_entries)
param_names = simulator_inst.parameter_space.parameter_names

# STAT
config_stat = data.configurations[config_names.index("a_STAT")]
simulator_stat = x3cflux.create_simulator_from_data(data.network_data, config_stat)
params_stat = x3cflux.get_parameters(simulator_stat.parameter_space, config_stat.parameter_entries)

samples = pd.read_csv("../data/poolsize_samples.csv", header=None).to_numpy()

total = len(samples)
print(f"simulating {total} INST samples")
r = [None]*total
for i,s in enumerate(samples):
    params[35:91] = s
    temp = simulator_inst.compute_measurements(params = params, time_stamps=[10_000_000])
    r[i] = temp
    # build & print bar
    util.progress_bar(i+1, total)
print("")

inst_names = {x.name: f"{x.metabolite_name}{mask_to_num(x.specification.mask)}" for x in simulator_inst.configurations[0].measurements if x.name in simulator_inst.measurement_names[0]}
stat_names = {x.name: f"{x.metabolite_name}{mask_to_num(x.specification.mask)}" for x in simulator_stat.configurations[0].measurements if x.name in simulator_stat.measurement_names[0]}
stat_meas_order = np.array([stat_names[x] for x in simulator_stat.measurement_names[0]])
inst_meas_order = np.array([inst_names[x] for x in simulator_inst.measurement_names[0]])
stat_order = np.argsort(stat_meas_order)
inst_order = np.argsort(inst_meas_order)

r_stat = simulator_stat.compute_measurements(params=params_stat)

diffs = []
for sample_i in range(len(samples)):
    if r[sample_i] == None:
        continue
    m = 0
    for inst_i, stat_i in zip (inst_order, stat_order):
        m=np.max([m, np.max(np.abs(r[sample_i][0][inst_i] - r_stat[0][stat_i]))])
    diffs.append(m)
print(f"Maximum Difference between any INST and the corresponding STAT measurement: {np.max(diffs):.2f}")

plt.figure(figsize=(6.4,4.8))
plt.hist(diffs, bins=np.logspace(-13,-5,25), weights=np.ones(len(diffs))/ len(diffs), color=C_X3)
plt.xscale('log')
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.minorticks_off()
plt.xlim([0.5e-13, 1.5e-6])
plt.xlabel('maximal difference to stationary solution ←')
ylims=plt.ylim()
plt.vlines(1e-6, ymin=ylims[0], ymax=ylims[1], color='k', linestyle='dashed', label='simulator\ntolerance')
plt.ylim(ylims)
plt.tight_layout()
filename = "../out/figure_s03"
print(f"saving figure to {filename}")
plt.savefig(filename + ".png", dpi=150)
plt.savefig(filename + ".svg")
