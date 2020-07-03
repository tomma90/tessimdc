import numpy as np
import matplotlib.pyplot as plt
import tessimdc as tes

# Initiate TES - This create an instance of a AC-biased TES with 
# some default values to bias the TES at 0.5xRn
tes_ac = tes.TesAcModel()
# You can see the default values using the print_info() method:
tes_ac.print_info()
# Default values can be modified. Example to modify bias current:
#tes_ac.modify_bias_current_amplitude(0.)

# Define time length [s] and sampling rate [Hz]. 
# The AC case requires much higher sampling rate than DC case because Mux system works at ~ few MHz!!!
time = 0.1
sampling_hz = 20.e6
step = int(time * sampling_hz + 1)

# Basic Constant Input
t = np.linspace(0., time, step)
P = np.ones_like(t)*tes_ac.optical_loading_power
Tb = np.ones_like(t)*tes_ac.temperature_focal_plane

# Basic Input with 1% sinusoidal power fluctuation
# t = np.linspace(0., time, step)
#P = np.ones_like(t)*tes_ac.optical_loading_power + 0.01*tes_ac.optical_loading_power*np.sin(2.*np.pi*10.*t)
# Tb = np.ones_like(t)*tes_ac.temperature_focal_plane

# Solve
I, T = tes.TesAcRungeKuttaSolver(t, P, Tb, tes_ac)
# TES resistance during simulation
R = tes_ac.resistance_vs_temperature(T)

# These functions are not optimized yet for AC case. They are simply imported from the DC case. Use with caution!!!
# Alpha during simulation
alpha = tes_ac.alpha_calculation(T)
# Loop gain during simulation
loop_g = tes_ac.loop_gain_calculation(T, I)
# Time constant during simulation
tau = tes_ac.time_constant_calculation(T, I)
# Current responsivity during simulation (the method is approximated for slowly varying input)
S_I = tes_ac.current_responsivity_calculation(T, I)

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
