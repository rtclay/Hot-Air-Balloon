import math
'''FLIGHT MODEL
Based on the 1976 USAA
and also http://wahiduddin.net/calc/density_altitude.htm
'''

temp_naught = 288.15
g_naught = 9.80665
n_avogadro = 6.022169 * math.pow(10, 36)
k_boltzmann = 1.380622 * math.pow(10, -23)
r_gas_constant = 8.31432 * math.pow(10, 3)

gas_molecules = ["N2", "O2", "Ar", "CO2", "Ne", "He", "Kr", "Xe", "CH4", "H2"]
gas_masses = [28.0134, 31.9988, 39.948, 44.00995, 20.183, 4.0026, 83.80, 131.30, 14.04303, 2.01594]
gas_frac_volume = [.78084, .209476, .00934, .000314, .00001818, .0000524, .0000114, .00000087, .000002, .0000005]
M_dry = 28.964 #molecular weight of dry air gm/mol 

heights = [x*1000 for x in [0.0, 11.0, 20.0, 32.0, 47, 51, 71, 84.8520]]
l_temp_gradients = [-6.5, 0.0, 1.0, 2.8, 0, -2.8, -2]

p_naught = 1.013250 * math.pow(10, 5)
r_naught = 6356766.0 #meters

def get_grav_accel(gmet_height):
    '''Return an approximation of the gravitational acceleration for a given geometric altitude
    based on equation 17
    '''
#    print "\ngetting grav accel"
#    
#    print g_naught *  math.pow(r_naught / (r_naught + get_gpot_height(gmet_height)), 2)
    return g_naught * math.pow(r_naught / (r_naught + get_gpot_height(gmet_height)), 2)


def get_gpot_height(gmet_height):
    '''Return the Geopotential height for a given geometric altitude
    based on equation 19
    '''
    return (r_naught * gmet_height) / (r_naught - gmet_height)



def get_temp(gpot_height):
    '''Return the molecular-scale temperature in Kelvin, given a geopotential gpot_height in meters
    based on equation 23
    equal to kinetic temperature up to 80km
    '''
    if gpot_height <= heights[1]:
        return temp_naught + l_temp_gradients[0] * gpot_height
    elif gpot_height <= heights[2]:
        return temp_naught + l_temp_gradients[1] * (gpot_height - heights[1])
    elif gpot_height <= heights[3]:
        return temp_naught + l_temp_gradients[2] * (gpot_height - heights[2])
    elif gpot_height <= heights[4]:
        return temp_naught + l_temp_gradients[3] * (gpot_height - heights[3])
    elif gpot_height <= heights[5]:
        return temp_naught + l_temp_gradients[4] * (gpot_height - heights[4])
    elif gpot_height <= heights[6]:
        return temp_naught + l_temp_gradients[5] * (gpot_height - heights[5])
    elif gpot_height <= heights[7]:
        return temp_naught + l_temp_gradients[6] * (gpot_height - heights[6])

def get_pressure(gmet_height):
    '''Return the pressure given a geopotential height in meters
    based on equation 33a and 33b
    use 33A when temp gradient is not 0
    use 33B when temp gradient is 0 
    
    temp gradients defined in l_temp_gradients
    per page 9, the molecular mass is essentially the same for each molecule up to 80km
    '''

    gpot_height = get_gpot_height(gmet_height)
    grav_accel = get_grav_accel(gmet_height)
    temp = get_temp(gpot_height)


    if gpot_height <= heights[1]:
        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[0])
        return p_naught * math.pow(temp / (temp + l_temp_gradients[0] * (gmet_height - heights[0])), exponent)
    elif gpot_height <= heights[2]:
        return p_naught * math.exp((-grav_accel) * M_dry * (gmet_height - heights[1]) / (r_gas_constant * temp))
    elif gpot_height <= heights[3]:
        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[2])
        return p_naught * math.pow(temp / (temp + l_temp_gradients[2] * (gmet_height - heights[2])), exponent)
    elif gpot_height <= heights[4]:
        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[3])
        return p_naught * math.pow(temp / (temp + l_temp_gradients[3] * (gmet_height - heights[3])), exponent)
    elif gpot_height <= heights[5]:
        return p_naught * math.exp((-grav_accel) * M_dry * (gmet_height - heights[4]) / (r_gas_constant * temp))
    elif gpot_height <= heights[6]:
        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[5])
        return p_naught * math.pow(temp / (temp + l_temp_gradients[5] * (gmet_height - heights[5])), exponent)
    elif gpot_height <= heights[7]:
        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[6])
        return p_naught * math.pow(temp / (temp + l_temp_gradients[6] * (gmet_height - heights[6])), exponent)



#Perfect gas is valid below 86km altitude
def perfect_gas(n_moles_of_gas, volume, temperature):
    '''Return P = NRT/V'''
    return (n_moles_of_gas * r_gas_constant * temperature) / volume

#dP = -g*

def get_deriv(func):
    '''Returns the derivative function of func'''
    pass

for H in range(500, 20000, 500):
    "%10.2f %10.2f %6.2f %10.2f %10.2f"%(H, get_gpot_height(H), get_grav_accel(H), get_pressure(H), get_temp(get_gpot_height(H)))

#H=5000
#print "%10.2f %10.2f %6.2f %10.2f %10.2f"%(H, get_gpot_height(H), get_grav_accel(H), get_pressure(H), get_temp(get_gpot_height(H)))
