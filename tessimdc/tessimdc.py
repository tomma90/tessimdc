# TES bolo model (see K.Irwin, G.C.Hilton, Transition Edge Sensors, 2005)

import numpy as np

PI = np.pi

class TesDcModel:

	def __init__(self):
		self.squid_input_inductor = 65.e-6 # SQUID input inductor [henry]
		self.shunt_resistor = 0.02 # shunt resistor [ohm]
		self.temperature_focal_plane = 0.1 # focal plane temperature [kelvin]
		self.tes_normal_resistance = 1. # TES noraml resistance [ohm]
		self.tes_log_sensitivity_alpha = 100. # TES alpha = dlogR/dlogT
		self.tes_leg_thermal_carrier_exponent = 4. # phonons (for electrons change to 2)
		self.tes_normal_time_constant = 33.e-3 # thermal time-constant in normal state [seconds]
		self.optical_loading_power = 0.5e-12 # optical loading power [watts]
		self.tes_saturation_power = 2.5*0.5e-12 # tes saturation power [watts]
		self.tes_transition_temperature = 1.71*0.1 # tes transition temperature [kelvin]
		self.tes_leg_thermal_conductivity = 2.5*0.5e-12/((1.71*0.1)**4.-0.1**4.)*4.*(1.71*0.1)**(4.-1.) # tes thermal conductivity [watts/kelvin]
		self.tes_heat_capacity = 33.e-3*2.5*0.5e-12/((1.71*0.1)**4.-0.1**4.)*4.*(1.71*0.1)**(4.-1.) # tes heat capacity [joule/kelvin]
		self.bias_current = 33.e-6 # bias current [ampere]	
	
	def print_info(self):
		print('\n')
		print('TES Model from Irwin and Hilton, 2005.')
		print('Default values for the DC biased TES:')
		print('\n')
		print('Normal resistance:%e ohm'%self.tes_normal_resistance)
		print('Shunt resistor:%e ohm'%self.shunt_resistor)
		print('SQUID input inductor:%e H'%self.squid_input_inductor)
		print('TES Log Sensitivity at 0.5 ohm (for arctan approx):%e'%self.tes_log_sensitivity_alpha)
		print('Focal plane temperature:%e K'%self.temperature_focal_plane)
		print('TES critical temperature:%e K (1.71 x bath temperature)'%self.tes_transition_temperature)
		print('TES leg thermal carrier:%e (n=b+1, b = 3 for phonons or 1 for electrons)'%self.tes_leg_thermal_carrier_exponent)
		print('TES normal time-constant:%e s'%self.tes_normal_time_constant)
		print('Total optical loading power expected:%e W'%self.optical_loading_power)
		print('Saturation power expected:%e W (Psat = 2.5 x optical expected)'%self.tes_saturation_power)
		print('TES leg thermal conductivity:%e W/K (P_sat*n*Tc^(n-1)/(Tc^n-Tb^n))'%self.tes_leg_thermal_conductivity)
		print('TES heat capacity:%e J/K (tau*G)'%self.tes_heat_capacity)
		print('Bias current:%e A'%self.bias_current)
		print('\n')
		print('To change one or more parameter use the modify methods.')
		print('Example: to modify bias_current use modify_bias_current. Same sintax for other parameters.')
		print('\n')

	def modify_squid_input_inductor(self, squid_input_inductor):
		self.squid_input_inductor = squid_input_inductor  

	def modify_shunt_resistor(self, shunt_resistor):
		self.shunt_resistor = shunt_resistor

	def modify_temperature_focal_plane(self, temperature_focal_plane):
		self.temperature_focal_plane = temperature_focal_plane  

	def modify_tes_normal_resistance(self, tes_normal_resistance):
		self.tes_normal_resistance = tes_normal_resistance

	def modify_tes_log_sensitivity_alpha(self, tes_log_sensitivity_alpha):
		self.tes_log_sensitivity_alpha = tes_log_sensitivity_alpha

	def modify_tes_leg_thermal_carrier_exponent(self, tes_leg_thermal_carrier_exponent):
		self.tes_leg_thermal_carrier_exponent = tes_leg_thermal_carrier_exponent

	def modify_tes_normal_time_constant(self, tes_normal_time_constant):
		self.tes_normal_time_constant = tes_normal_time_constant

	def modify_optical_loading_power(self, optical_loading_power):
		self.optical_loading_power = optical_loading_power

	def modify_tes_saturation_power(self, tes_saturation_power):
		self.tes_saturation_power = tes_saturation_power

	def modify_tes_transition_temperature(self, tes_transition_temperature):
		self.tes_transition_temperature = tes_transition_temperature

	def modify_tes_leg_thermal_conductivity(self, tes_leg_thermal_conductivity):
		self.tes_leg_thermal_conductivity = tes_leg_thermal_conductivity

	def modify_tes_heat_capacity(self, tes_heat_capacity):
		self.tes_heat_capacity = tes_heat_capacity

	def modify_bias_current(self, bias_current):
		self.bias_current = bias_current


	def resistance_vs_temperature(self, temperature):

		alpha = self.tes_log_sensitivity_alpha
		Tc = self.tes_transition_temperature
		Rn = self.tes_normal_resistance

		steepness_rt = alpha * PI * (Tc * Tc + 1.) * Rn / 2. / Tc
	
		return Rn * (np.arctan((temperature - Tc) * steepness_rt) + PI / 2.) / PI


	def differential_equations(self, current_temperature, bias_current, loading_power, bath_temperature):

		I = current_temperature[0]
		T = current_temperature[1]
		R = self.resistance_vs_temperature(T)
		Rs = self.shunt_resistor
		G = self.tes_leg_thermal_conductivity
		n = self.tes_leg_thermal_carrier_exponent
		L = self.squid_input_inductor
		C = self.tes_heat_capacity

		V = bias_current * R * Rs / (R + Rs)
		Pb = G * (T**n - bath_temperature**n) / (n * T**(n - 1.))
		Pr = V * V / R

		return [(V - I * Rs - I * R) / L, (-Pb + Pr + loading_power) / C]


	def dR_dT(self, temperature):
		
		h_t = 1.e-6
		
		return (self.resistance_vs_temperature(temperature + h_t) - self.resistance_vs_temperature(temperature - h_t)) / 2. / h_t 


	def alpha_calculation(self, temperature):
		
		return temperature * self.dR_dT(temperature) / self.resistance_vs_temperature(temperature)


	def loop_gain_calculation(self, temperature, current):
		
		return self.alpha_calculation(temperature) * current * current * self.resistance_vs_temperature(temperature) / self.tes_leg_thermal_conductivity / temperature

	
	def time_constant_calculation(self, temperature, current):
		
		return self.tes_normal_time_constant / (self.loop_gain_calculation(temperature, current) + 1.)


	def current_responsivity_calculation(self, temperature, current):
		
		return - 1. / current / self.resistance_vs_temperature(temperature) * self.loop_gain_calculation(temperature, current) / (self.loop_gain_calculation(temperature, current) + 1.)


class TesAcModel:

	def __init__(self):
		self.squid_input_inductor = 10.e-9 # SQUID input inductor [henry]
		self.shunt_resistor = 0.02 # shunt resistor [ohm]
		self.temperature_focal_plane = 0.1 # focal plane temperature [kelvin]
		self.tes_normal_resistance = 1. # TES noraml resistance [ohm]
		self.tes_log_sensitivity_alpha = 100. # TES alpha = dlogR/dlogT
		self.tes_leg_thermal_carrier_exponent = 4. # phonons (for electrons change to 2)
		self.tes_normal_time_constant = 33.e-3 # thermal time-constant in normal state [seconds]
		self.optical_loading_power = 0.5e-12 # optical loading power [watts]
		self.tes_saturation_power = 2.5*0.5e-12 # tes saturation power [watts]
		self.tes_transition_temperature = 1.71*0.1 # tes transition temperature [kelvin]
		self.tes_leg_thermal_conductivity = 2.5*0.5e-12/((1.71*0.1)**4.-0.1**4.)*4.*(1.71*0.1)**(4.-1.) # tes thermal conductivity [watts/kelvin]
		self.tes_heat_capacity = 33.e-3*2.5*0.5e-12/((1.71*0.1)**4.-0.1**4.)*4.*(1.71*0.1)**(4.-1.) # tes heat capacity [joule/kelvin]
		self.bias_current_amplitude = 33.e-6 * np.sqrt(2.) # bias current [ampere]
		self.ac_frequency = 1.e6 # AC bias driving frequency [hz]
		self.mux_frequency = 1.e6 # LC filter frequency [hz]
		self.mux_lc_inductor = 65.e-6 # Mux LC filter inductor [henry]
		self.mux_lc_capacitor = 1. / (65.e-6 * (2. * PI * 1.e6) * (2. * PI * 1.e6)) # Mux LC filter capacitor [farad]
	
	def print_info(self):
		print('\n')
		print('TES Model from Irwin and Hilton, 2005.')
		print('Default values for the DC biased TES:')
		print('\n')
		print('Normal resistance:%e ohm'%self.tes_normal_resistance)
		print('Shunt resistor:%e ohm'%self.shunt_resistor)
		print('SQUID input inductor:%e H'%self.squid_input_inductor)
		print('TES Log Sensitivity at 0.5 ohm (for arctan approx):%e'%self.tes_log_sensitivity_alpha)
		print('Focal plane temperature:%e K'%self.temperature_focal_plane)
		print('TES critical temperature:%e K (1.71 x bath temperature)'%self.tes_transition_temperature)
		print('TES leg thermal carrier:%e (n=b+1, b = 3 for phonons or 1 for electrons)'%self.tes_leg_thermal_carrier_exponent)
		print('TES normal time-constant:%e s'%self.tes_normal_time_constant)
		print('Total optical loading power expected:%e W'%self.optical_loading_power)
		print('Saturation power expected:%e W (Psat = 2.5 x optical expected)'%self.tes_saturation_power)
		print('TES leg thermal conductivity:%e W/K (P_sat*n*Tc^(n-1)/(Tc^n-Tb^n))'%self.tes_leg_thermal_conductivity)
		print('TES heat capacity:%e J/K (tau*G)'%self.tes_heat_capacity)
		print('Bias current amplitude:%e A'%self.bias_current_amplitude)
		print('AC bias frequency:%e Hz'%self.ac_frequency)
		print('LC filter frequency:%e Hz'%self.mux_frequency)
		print('LC filter inductance:%e H'%self.mux_lc_inductor)
		print('LC filter capacitance:%e F'%self.mux_lc_capacitor)
		print('\n')
		print('To change one or more parameter use the modify methods.')
		print('Example: to modify bias_current use modify_bias_current. Same sintax for other parameters.')
		print('\n')

	def modify_squid_input_inductor(self, squid_input_inductor):
		self.squid_input_inductor = squid_input_inductor  

	def modify_shunt_resistor(self, shunt_resistor):
		self.shunt_resistor = shunt_resistor

	def modify_temperature_focal_plane(self, temperature_focal_plane):
		self.temperature_focal_plane = temperature_focal_plane  

	def modify_tes_normal_resistance(self, tes_normal_resistance):
		self.tes_normal_resistance = tes_normal_resistance

	def modify_tes_log_sensitivity_alpha(self, tes_log_sensitivity_alpha):
		self.tes_log_sensitivity_alpha = tes_log_sensitivity_alpha

	def modify_tes_leg_thermal_carrier_exponent(self, tes_leg_thermal_carrier_exponent):
		self.tes_leg_thermal_carrier_exponent = tes_leg_thermal_carrier_exponent

	def modify_tes_normal_time_constant(self, tes_normal_time_constant):
		self.tes_normal_time_constant = tes_normal_time_constant

	def modify_optical_loading_power(self, optical_loading_power):
		self.optical_loading_power = optical_loading_power

	def modify_tes_saturation_power(self, tes_saturation_power):
		self.tes_saturation_power = tes_saturation_power

	def modify_tes_transition_temperature(self, tes_transition_temperature):
		self.tes_transition_temperature = tes_transition_temperature

	def modify_tes_leg_thermal_conductivity(self, tes_leg_thermal_conductivity):
		self.tes_leg_thermal_conductivity = tes_leg_thermal_conductivity

	def modify_tes_heat_capacity(self, tes_heat_capacity):
		self.tes_heat_capacity = tes_heat_capacity

	def modify_bias_current_amplitude(self, bias_current_amplitude):
		self.bias_current_amplitude = bias_current_amplitude

	def modify_ac_frequency(self, ac_frequency):
		self.ac_frequency = ac_frequency

	def modify_mux_frequency(self, mux_frequency):
		self.mux_frequency = mux_frequency

	def modify_mux_lc_inductor(self, mux_lc_inductor):
		self.mux_lc_inductor = mux_lc_inductor

	def modify_mux_lc_capacitor(self, mux_lc_capacitor):
		self.mux_lc_capacitor = mux_lc_capacitor


	def resistance_vs_temperature(self, temperature):

		alpha = self.tes_log_sensitivity_alpha
		Tc = self.tes_transition_temperature
		Rn = self.tes_normal_resistance

		steepness_rt = alpha * PI * (Tc * Tc + 1.) * Rn / 2. / Tc
	
		return Rn * (np.arctan((temperature - Tc) * steepness_rt) + PI / 2.) / PI


	def differential_equations(self, current_dcurrentdt_temperature, time, loading_power, bath_temperature):

		I = current_dcurrentdt_temperature[0]
		J = current_dcurrentdt_temperature[1]
		T = current_dcurrentdt_temperature[2]
		R = self.resistance_vs_temperature(T)
		Rs = self.shunt_resistor
		G = self.tes_leg_thermal_conductivity
		n = self.tes_leg_thermal_carrier_exponent
		L_s = self.squid_input_inductor
		L_mux = self.mux_lc_inductor
		C_mux = self.mux_lc_capacitor
		C = self.tes_heat_capacity

		dVdt = self.bias_current_amplitude * R * Rs * 2. * PI * self.ac_frequency * np.cos(2. * PI * self.ac_frequency * time) / (R + Rs)
		Pb = G * (T**n - bath_temperature**n) / (n * T**(n - 1.))
		Pr = I * I * R

		return [J, (dVdt - J * R - I / C_mux) / (L_mux + L_s), (-Pb + Pr + loading_power) / C]


	def dR_dT(self, temperature):
		
		h_t = 1.e-9
		
		return (self.resistance_vs_temperature(temperature + h_t) - self.resistance_vs_temperature(temperature - h_t)) / 2. / h_t 


	def alpha_calculation(self, temperature):
		
		return temperature * self.dR_dT(temperature) / self.resistance_vs_temperature(temperature)


	def loop_gain_calculation(self, temperature, current):
		
		return self.alpha_calculation(temperature) * current * current * self.resistance_vs_temperature(temperature) / self.tes_leg_thermal_conductivity / temperature

	
	def time_constant_calculation(self, temperature, current):
		
		return self.tes_normal_time_constant / (self.loop_gain_calculation(temperature, current) + 1.)


	def current_responsivity_calculation(self, temperature, current):
		
		return - 1. / current / self.resistance_vs_temperature(temperature) * self.loop_gain_calculation(temperature, current) / (self.loop_gain_calculation(temperature, current) + 1.)



# Solve & update Dc TES I & T with Runge-Kutta method
def TesDcRungeKuttaSolver(time_array, bias_current_array, loading_power_array, bath_temperature_array, tesdcmodel = TesDcModel()):

	if len(time_array) != len(bias_current_array):
		print('Error!!! Input arrays need to have same size!')
		return

	if len(time_array) != len(loading_power_array):
		print('Error!!! Input arrays need to have same size!')
		return

	if len(time_array) != len(bath_temperature_array):
		print('Error!!! Input arrays need to have same size!')
		return

	step = len(time_array)
	h = abs(time_array[1]-time_array[0])
	I = np.zeros(step)
	T = np.zeros(step)

	# Initial conditions
	T0 = tesdcmodel.tes_transition_temperature
	I0 = tesdcmodel.bias_current * tesdcmodel.shunt_resistor / (tesdcmodel.resistance_vs_temperature(T0) + tesdcmodel.shunt_resistor)
	IT0 = [I0, T0]

	for i in range(int(step)):
		kI1_tmp, kT1_tmp = tesdcmodel.differential_equations(IT0, bias_current_array[i], loading_power_array[i], bath_temperature_array[i])
		kI1 = h*kI1_tmp
		kT1 = h*kT1_tmp
	
		IT0_tmp = [IT0[0]+0.5*kI1, IT0[1]+0.5*kT1]
		kI2_tmp, kT2_tmp = tesdcmodel.differential_equations(IT0_tmp, bias_current_array[i], loading_power_array[i], bath_temperature_array[i])
		kI2 = h*kI2_tmp
		kT2 = h*kT2_tmp
		
		IT0_tmp = [IT0[0]+0.5*kI2, IT0[1]+0.5*kT2]
		kI3_tmp, kT3_tmp = tesdcmodel.differential_equations(IT0_tmp, bias_current_array[i], loading_power_array[i], bath_temperature_array[i])
		kI3 = h*kI3_tmp
		kT3 = h*kT3_tmp
	
		IT0_tmp = [IT0[0]+kI3, IT0[1]+kT3]
		kI4_tmp, kT4_tmp = tesdcmodel.differential_equations(IT0_tmp, bias_current_array[i], loading_power_array[i], bath_temperature_array[i])
		kI4 = h*kI4_tmp
		kT4 = h*kT4_tmp

		IT0 = [IT0[0]+(kI1+2*kI2+2*kI3+kI4)/6., IT0[1]+(kT1+2*kT2+2*kT3+kT4)/6.]
		I[i]=IT0[0]
		T[i]=IT0[1]
	
	return I, T

# Solve & update Ac TES I & T with Runge-Kutta method
def TesAcRungeKuttaSolver(time_array, loading_power_array, bath_temperature_array, tesacmodel = TesAcModel()):

	if len(time_array) != len(loading_power_array):
		print('Error!!! Input arrays need to have same size!')
		return

	if len(time_array) != len(bath_temperature_array):
		print('Error!!! Input arrays need to have same size!')
		return

	step = len(time_array)
	h = abs(time_array[1]-time_array[0])
	I = np.zeros(step)
	T = np.zeros(step)

	# Initial conditions
	T0 = tesacmodel.tes_transition_temperature 
	J0 = 0.
	I0 = tesacmodel.bias_current_amplitude * tesacmodel.shunt_resistor / (tesacmodel.resistance_vs_temperature(T0) + tesacmodel.shunt_resistor)
	IJT0 = [I0, J0, T0]

	for i in range(int(step)):
		kI1_tmp, kJ1_tmp, kT1_tmp = tesacmodel.differential_equations(IJT0, time_array[i], loading_power_array[i], bath_temperature_array[i])
		kI1 = h*kI1_tmp
		kJ1 = h*kJ1_tmp
		kT1 = h*kT1_tmp
	
		IJT0_tmp = [IJT0[0]+0.5*kI1, IJT0[1]+0.5*kJ1, IJT0[2]+0.5*kT1]
		kI2_tmp, kJ2_tmp, kT2_tmp = tesacmodel.differential_equations(IJT0_tmp, time_array[i], loading_power_array[i], bath_temperature_array[i])
		kI2 = h*kI2_tmp
		kJ2 = h*kJ2_tmp
		kT2 = h*kT2_tmp
		
		IJT0_tmp = [IJT0[0]+0.5*kI2, IJT0[1]+0.5*kJ2, IJT0[2]+0.5*kT2]
		kI3_tmp, kJ3_tmp, kT3_tmp = tesacmodel.differential_equations(IJT0_tmp, time_array[i], loading_power_array[i], bath_temperature_array[i])
		kI3 = h*kI3_tmp
		kJ3 = h*kJ3_tmp
		kT3 = h*kT3_tmp
	
		IJT0_tmp = [IJT0[0]+kI3, IJT0[1]+kJ3, IJT0[2]+kT3]
		kI4_tmp, kJ4_tmp, kT4_tmp = tesacmodel.differential_equations(IJT0_tmp, time_array[i], loading_power_array[i], bath_temperature_array[i])
		kI4 = h*kI4_tmp
		kJ4 = h*kJ4_tmp
		kT4 = h*kT4_tmp

		IJT0 = [IJT0[0]+(kI1+2*kI2+2*kI3+kI4)/6., IJT0[1]+(kJ1+2*kJ2+2*kJ3+kJ4)/6., IJT0[2]+(kT1+2*kT2+2*kT3+kT4)/6.]
		I[i]=IJT0[0]
		T[i]=IJT0[2]
	
	return I, T


