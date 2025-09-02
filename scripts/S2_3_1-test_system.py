#!/usr/bin/env python
# coding: utf-8
import x3cflux
import numpy as np
import matplotlib.pyplot as plt
from sympy import Function, dsolve, Derivative, Symbol
from sympy.abc import x
import os
import util


def global_error(jac, rhs, solution):
    a = Function("a")
    b = Function("b")
    system = [Derivative(b(x), x) - jac[0, 0] * b(x) - jac[0, 1] * a(x),
              Derivative(a(x), x) - jac[1, 0] * b(x) - jac[1, 1] * a(x) - rhs]
    result = dsolve(system, ics={a(0.): 0.01109, b(0.): 0.01109})
    errors = np.linalg.norm([np.array([float(result[0].rhs.evalf(subs={x: t})),
                                       float(result[1].rhs.evalf(subs={x: t}))])
                             - solution(t)
                             for t in np.linspace(0., 1., 100)],
                            axis=1)
    return np.max(errors)

util.print_box(f"Executing {os.path.basename(__file__)}")
plt.rcParams['font.size'] = 16
plt.rcParams['figure.figsize'][0] = 10
x3cflux.logging.level = 0
simulator = x3cflux.create_simulator_from_fml("../models/linea.fml")
simulator.parameter_space.free_parameter_names

tau = Symbol("tau")
A = Symbol("A")
a = Function("a")
b = Function("b")
print("Solving system manually")
system = [Derivative(b(x), x) + (tau + 1) * b(x) - (tau +1) * a(x),
          A * Derivative(a(x), x) - tau * b(x) + (tau + 1) * a(x) - 1]
result = dsolve(system, ics={a(0.): 0.01109, b(0.): 0.01109})

data = {}
for name in ["bdf", "sdirk"]:
    simulator.builder.set_solver(name)
    data.update({name: {}})
    print(f"Solving system with x3cflux using {name} solver")
    for j in [3, 6, 9]:
        simulator.builder.solver.relative_tolerance = np.power(10., -j)
        simulator.builder.solver.absolute_tolerance = np.power(10., -(j+3))
        cond, errors = [], []
        print(f"\tsolving for solver tolerance 1e{-j}")
        for i in np.linspace(1, 3, 25):
            val = 10 ** i
            new_params = [1., val, 1., 1 / val]
            casc_sys = simulator.builder.build(simulator.parameter_space.compute_parameters(new_params))
            level_sys = casc_sys.get_level_system(0)
            jac = level_sys.jacobian.todense()
            cond.append(np.linalg.cond(jac))
            solution = casc_sys.solve()
            errors.append(global_error(jac, level_sys.evaluate_inhomogeneity(0.)[1], solution[0]))
            data[name].update({j: np.array(list(zip(cond, errors)))})

for j in [3, 6, 9]:
    plt.axhline(np.power(10., -j), color="black", linestyle="-" if j < 4 else "--" if j < 7 else ":")
for name in data:
    for j in data[name]:
        cond, errors = data[name][j].T
        display_name = "BDF" if name == "bdf" else "SDIRK"
        plt.loglog(cond, errors, label=f"{display_name}: " + "$tol_{rel} = 10^{-" + str(j) + "}$",
                   color= "tab:blue" if name == "bdf" else "tab:orange",
                   linestyle="-" if j < 4 else "--" if j < 7 else ":")
plt.legend()
plt.xlabel("condition number")
plt.ylabel("global error")
plt.xlim((2e2, 2e6))
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
plt.tight_layout()
filename = "../out/figure_s01"
print(f"saving figure to {filename}")
plt.savefig(filename + ".png", dpi=150)
plt.savefig(filename + ".svg")
