"""
This example runs the Electrolyzer on its own and generates polarization curves.
"""
import numpy as np
import matplotlib.pyplot as plt

from electrolyzer import Stack, electrolyzer_model


n_cells = 100  # number of cells in stack
cell_area = 1000  # cell area, cm^2
temperature = 60  # temperature
max_current = 2000
fact_fit_act = 100
fact_fit_ohm = 100
fit_current_interval = .1
fit_temp_interval = 5
fit_temp_min = 50
fit_temp_max = 70


stack_dict = {
    "n_cells": n_cells,
    "cell_area": cell_area,
    "temperature": temperature,
    "max_current": max_current,
    "dt": 1,
    "fact_fit_act": fact_fit_act,
    "fact_fit_ohm": fact_fit_ohm,
    "fit_current_interval": fit_current_interval,
    "fit_temp_interval": fit_temp_interval,
    "fit_temp_min": fit_temp_min,
    "fit_temp_max": fit_temp_max,
}

elec = Stack.from_dict(stack_dict)

cur = np.linspace(0, 2500, 100)
p_fit = elec.calc_stack_power(
    electrolyzer_model((elec.calc_stack_power(cur),
                        temperature), *elec.fit_params)
)
p_actual = elec.calc_stack_power(cur)
voltage = elec.cell.calc_cell_voltage(cur, temperature)
fit_error = p_actual - p_fit

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.plot(cur / (1000), fit_error / elec.stack_rating_kW * 100)
ax1.grid(True)
ax1.set_ylabel("error (% rated power)")
ax1.set_title("Power error (actual - fit)")
ax2.plot(cur / (1000), p_actual / elec.stack_rating_kW)
ax2.plot(cur / (1000), p_fit / elec.stack_rating_kW)
ax2.grid(True)
ax2.set_title("power normalized by stack rating")
ax2.set_ylabel("normalized power")
ax2.set_xlabel("Current density (A/cm^2)")
ax2.legend(["actual", "fit"])

plt.figure()
plt.plot(cur / 1000, voltage)
plt.grid(True)
plt.xlabel("Current density (A/cm^2)")
plt.ylabel("Voltage (V)")
plt.title("Polarization Curve")
plt.show()
