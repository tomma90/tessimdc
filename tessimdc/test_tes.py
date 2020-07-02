import numpy as np
import matplotlib.pyplot as plt
import tessimdc as tes

# Initiate TES - This create an instance of a DC-biased TES with 
# some default values to bias the TES at 0.5xRn
tes_1 = tes.TesDcModel()
# You can see the default values using the print_info() method:
tes_1.print_info()
# Default values can be modified. Example to modify bias current:
# tes_1.modify_bias_current(100.e-6)

# Define time length [s] and sampling rate [Hz]
time = 10.
sampling_hz = 10000
step = int(time * sampling_hz + 1)

# Basic Constant Input
t = np.linspace(0., time, step)
Ib = np.ones_like(t)*tes_1.bias_current
#P = np.ones_like(t)*tes_1.optical_loading_power
Tb = np.ones_like(t)*tes_1.temperature_focal_plane

# Basic Input with 1% sinusoidal power fluctuation
# t = np.linspace(0., time, step)
# Ib = np.ones_like(t)*tes_1.bias_current
P = np.ones_like(t)*tes_1.optical_loading_power + 0.01*tes_1.optical_loading_power*np.sin(2.*np.pi*0.25*t)
# Tb = np.ones_like(t)*tes_1.temperature_focal_plane

# Solve
I, T = tes.TesRungeKuttaSolver(t, Ib, P, Tb, tes_1)
# TES resistance during simulation
R = tes_1.resistance_vs_temperature(T)
# Alpha during simulation
alpha = tes_1.alpha_calculation(T)
# Loop gain during simulation
loop_g = tes_1.loop_gain_calculation(T, I)
# Time constant during simulation
tau = tes_1.time_constant_calculation(T, I)
# Current responsivity during simulation (the method is approximated for slowly varying input)
S_I = tes_1.current_responsivity_calculation(T, I)

# Plot
plt.figure()
plt.plot(t, T)
plt.figure()
plt.plot(t, I)
plt.figure()
plt.plot(t, R)
plt.figure()
plt.plot(t, alpha)
plt.figure()
plt.plot(t, loop_g)
plt.figure()
plt.plot(t, tau)
plt.figure()
plt.plot(t, S_I)
plt.show()
