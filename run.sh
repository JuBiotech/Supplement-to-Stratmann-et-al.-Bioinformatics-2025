#!/bin/bash

PYTHON_EXEC="/opt/python/cp313-cp313/bin/python"

mkdir -p /task/out

cd /task/scripts
${PYTHON_EXEC} -u S2_3_1-test_system.py
${PYTHON_EXEC} -u S2_3_2-real_systems.py
${PYTHON_EXEC} -u S2_3_3-inst_stat.py
${PYTHON_EXEC} -u S2_4_1-essential_dim_cumo_emu.py
${PYTHON_EXEC} -u S2_4_2-mile_ed.py
${PYTHON_EXEC} -u S6_2-mcmc_postprocessing.py
