import x3cflux
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from time import perf_counter
import sys

plt.rcParams['font.size'] = 16
C_X3 = '#023d6b'

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

data = x3cflux.FluxMLParser().parse("../models/EC.fml")
config_names = [i.name for i in data.configurations]

# INST
config = data.configurations[config_names.index("a_INST")]
simulator_inst = x3cflux.create_simulator_from_data(data.network_data, config)
simulator_inst.builder.solver.relative_tolerance = 1e-6
simulator_inst.builder.solver.absolute_tolerance = 1e-9
params = x3cflux.get_parameters(simulator_inst.parameter_space, config.parameter_entries)
param_names = simulator_inst.parameter_space.parameter_names

samples = pd.read_csv("../data/poolsize_samples.csv", header=None).to_numpy()

total = len(samples)
print(f"simulating {total} INST samples")
bar_len = 80
time = np.zeros(total,)
num_repeats = 10
for i,s in enumerate(samples):
    params[35:91] = s
    before = perf_counter()
    for j in range(num_repeats):
        temp = simulator_inst.compute_measurements(params = params, time_stamps=[10_000_000])
    after = perf_counter()
    time[i] = (after-before)/num_repeats
    # build & print bar
    idx = i + 1
    filled = int(bar_len * idx / total)
    bar = '#' * filled + '-' * (bar_len - filled)
    pct = idx / total * 100
    sys.stdout.write(f'\r[{bar}] {pct:5.1f}% ({idx}/{total})')
    sys.stdout.flush()
print("")

plt.hist(time, bins=40, color=C_X3)
plt.xlabel('← runtime [ms]')
plt.tight_layout()
plt.savefig('../out/figure_s8.png', dpi=150)
plt.savefig('../out/figure_s8.svg')
