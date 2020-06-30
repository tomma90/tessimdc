# TES bolo model (see K.Irwin, G.C.Hilton, Transition Edge Sensors, 2005)

import numpy as np

PI = np.pi

class DcBiasedTesParameters:

    def __init__(self):
        self.squid_input_inductor = 1.e-9 # SQUID input inductor [henry]
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


class ResistanceTemperature(DcBiasedTesParameters):
	
    def RT(self, temperature):
    	steepness_rt = self.tes_log_sensitivity_alpha * PI * (self.tes_transition_temperature * self.tes_transition_temperature + 1.) * self.tes_normal_resistance / 2. / self.tes_transition_temperature
    	return self.tes_normal_resistance * (np.arctan((temperature - self.tes_transition_temperature) * steepness_rt) + PI / 2.) / PI


class CurrentTemperatureDiffEquations(ResistanceTemperature):

    def bolo(self, current_temperature, bias_current, loading_power, bath_temperature):
    	I = current_temperature[0]
    	T = current_temperature[1]
    	R = self.RT(T)
    	V = bias_current * R * self.shunt_resistor / (R + self.shunt_resistor)
    	Pb = self.tes_leg_thermal_conductivity * (T**self.tes_leg_thermal_carrier_exponent - bath_temperature**self.tes_leg_thermal_carrier_exponent) / (self.tes_leg_thermal_carrier_exponent * T**(self.tes_leg_thermal_carrier_exponent - 1.))
    	Pr = V * V / R
    	dITdt = [(V - I * self.shunt_resistor - I * R) / self.squid_input_inductor, (-Pb + Pr + loading_power) / self.tes_heat_capacity]
    	return dITdt


#solve & update TES I & T with Runge-Kutta method
def RungeKuttaSolver(time_array, bias_current, loading_power, bath_temperature, class_type = CurrentTemperatureDiffEquations):

    if len(time_array) != len(bias_current):
	print('Error!!! Input arrays need to have same size!')
	return

    if len(time_array) != len(loading_power):
	print('Error!!! Input arrays need to have same size!')
	return

    if len(time_array) != len(bath_temperature):
	print('Error!!! Input arrays need to have same size!')
	return

    step = len(time_array)
    h = abs(time_array[1]-time_array[0])
    I = np.zeros(step)
    T = np.zeros(step)

    detector_output = class_type()

    T0 = detector_output.tes_transition_temperature
    I0 = detector_output.bias_current * detector_output.shunt_resistor / (detector_output.RT(T0) + detector_output.shunt_resistor)
    IT0 = [I0, T0]

    for i in range(int(step)):
    	kI1_tmp, kT1_tmp = detector_output.bolo(IT0, bias_current[i], loading_power[i], bath_temperature[i])
    	kI1 = h*kI1_tmp
    	kT1 = h*kT1_tmp

    	IT0_tmp = [IT0[0]+0.5*kI1, IT0[1]+0.5*kT1]
    	kI2_tmp, kT2_tmp = detector_output.bolo(IT0_tmp, bias_current[i], loading_power[i], bath_temperature[i])
    	kI2 = h*kI2_tmp
    	kT2 = h*kT2_tmp
    
    	IT0_tmp = [IT0[0]+0.5*kI2, IT0[1]+0.5*kT2]
    	kI3_tmp, kT3_tmp = detector_output.bolo(IT0_tmp, bias_current[i], loading_power[i], bath_temperature[i])
    	kI3 = h*kI3_tmp
    	kT3 = h*kT3_tmp

    	IT0_tmp = [IT0[0]+kI3, IT0[1]+kT3]
    	kI4_tmp, kT4_tmp = detector_output.bolo(IT0_tmp, bias_current[i], loading_power[i], bath_temperature[i])
    	kI4 = h*kI4_tmp
    	kT4 = h*kT4_tmp

    	IT0 = [IT0[0]+(kI1+2*kI2+2*kI3+kI4)/6., IT0[1]+(kT1+2*kT2+2*kT3+kT4)/6.]
    	I[i]=IT0[0]
    	T[i]=IT0[1]
    
    return I, T



