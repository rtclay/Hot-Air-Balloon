import math
'''FLIGHT MODEL
Based on the 1976 USAA'''

temp_naught = 288.15
g_naught = 9.80665 
n_avogadro = 6.022169*math.pow(10, 36)
k_boltzmann = 1.380622*math.pow(10, -23)
r_gas_constant = 8.31432*math.pow(10, 3)

gas_molecules= ["N2", "O2", "Ar", "CO2", "Ne", "He", "Kr", "Xe", "CH4", "H2"]
gas_masses= [28.0134, 31.9988, 39.948, 44.00995, 20.183, 4.0026, 83.80, 131.30, 14.04303, 2.01594]
gas_frac_volume= [.78084, .209476, .00934, .000314, .00001818, .0000524, .0000114, .00000087, .000002, .0000005]

heights = 1000*[0, 11,20, 32, 47, 51, 71, 84.8520]
l_temp_gradients = [-6.5, 0, 1, 2.8, 0, -2.8, -2]

p_naught = 1.013250*math.pow(10, 5)
r_naught = 6356766 #meters

def get_grav_accel(gmet_height):
    '''Return an approximation of the gravitational acceleration for a given geometric altitude
    based on equation 17
    '''
    return g_naught*math.pow(r_naught/(r_naught*get_gpot_height(gmet_height)), 2)


def get_gpot_height(gmet_height):
    '''Return the Geopotential height for a given geometric altitude
    based on equation 19
    '''
    return (r_naught*gmet_height)/(r_naught-gmet_height)



def get_temp(gpot_height):
    '''Return the molecular-scale temperature in Kelvin, given a geopotential gpot_height in meters
    based on equation 23
    equal to kinetic temperature up to 80km
    '''
    if gpot_height <=heights[1]:
        return temp_naught+l_temp_gradients[0]*gpot_height
    elif gpot_height <=heights[2]:
        return temp_naught+l_temp_gradients[1]*(gpot_height-heights[1])
    elif gpot_height <=heights[3]:
        return temp_naught+l_temp_gradients[2]*(gpot_height-heights[2])
    elif gpot_height <=heights[4]:
        return temp_naught+l_temp_gradients[3]*(gpot_height-heights[3])
    elif gpot_height <=heights[5]:
        return temp_naught+l_temp_gradients[4]*(gpot_height-heights[4])
    elif gpot_height <=heights[6]:
        return temp_naught+l_temp_gradients[5]*(gpot_height-heights[5])
    elif gpot_height <=heights[7]:
        return temp_naught+l_temp_gradients[6]*(gpot_height-heights[6])

def get_pressure(gpot_height):
    '''Return the pressure given a geopotential height in meters
    based on equation 33a and 33b
    use 33A when temp gradient is not 0
    use 33B when temp gradient is 0 
    '''
    
    FIX ME PLEASE
    return math.exp((-get_grav_accel(gmet_height)*))
    
    if gpot_height <=heights[1]:
        return temp_naught+l_temp_gradients[0]*gpot_height
    elif gpot_height <=heights[2]:
        return temp_naught+l_temp_gradients[1]*(gpot_height-heights[1])
    elif gpot_height <=heights[3]:
        return temp_naught+l_temp_gradients[2]*(gpot_height-heights[2])
    elif gpot_height <=heights[4]:
        return temp_naught+l_temp_gradients[3]*(gpot_height-heights[3])
    elif gpot_height <=heights[5]:
        return temp_naught+l_temp_gradients[4]*(gpot_height-heights[4])
    elif gpot_height <=heights[6]:
        return temp_naught+l_temp_gradients[5]*(gpot_height-heights[5])
    elif gpot_height <=heights[7]:
        return temp_naught+l_temp_gradients[6]*(gpot_height-heights[6])



#Perfect gas is valid below 86km altitude
def perfect_gas(n_moles_of_gas, volume, temperature):
    '''Return P = NRT/V'''
    return (n_moles_of_gas*r_gas_constant*temperature)/volume

#dP = -g*

def get_deriv(func):
    '''Returns the derivative function of func'''
    