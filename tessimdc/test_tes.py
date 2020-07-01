import numpy as np
import matplotlib.pyplot as plt
import tessimdc as tes

# Initiate TES
tes_1 = tes.DcBiasedTesParameters()

# Define time length [s] and sampling rate [Hz]
time = 10.
sampling_hz = 10000
step = int(time * sampling_hz + 1)

# Basic Input with 1% sinusoidal power fluctuation
t = np.linspace(0., time, step)
Ib = np.ones_like(t)*tes_1.bias_current
P = np.ones_like(t)*tes_1.optical_loading_power + 0.01*tes_1.optical_loading_power*np.sin(2.*np.pi*0.25*t)
Tb = np.ones_like(t)*tes_1.temperature_focal_plane

# Solve
I, T = tes.TesRungeKuttaSolver(t, Ib, P, Tb)

# Plot
plt.figure()
plt.plot(t, T)
plt.figure()
plt.plot(t, I)
plt.show()
