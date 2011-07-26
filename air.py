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

heights = [x * 1000 for x in [0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.8520]]
l_temp_gradients = [x / 1000 for x in [-6.5, 0.0, 1.0, 2.8, 0, -2.8, -2.0]]

p_naught = 1.013250 * math.pow(10, 5)
r_naught = 6356766.0 #meters
universal_grav = 6.67384 * math.pow(10, -11)
mass_earth = 5.9722 * math.pow(10, 24)



#Nomex http://www2.dupont.com/Energy_Solutions/en_US/assets/downloads/E56.pdf
#Nylon http://en.wikipedia.org/wiki/Nylon 

densities = dict({"steel": 7850,
                  "aluminum": 2600,
                  "silk": 1340,
                  "nylon": 1150,
                  "Nomex": 670,

                  })

def get_grav_accel(gmet_height):
    '''Return an approximation of the gravitational acceleration for a given geometric altitude
    based on equation 17
    '''
#    print "\ngetting grav accel"
#    
#    print g_naught *  math.pow(r_naught / (r_naught + get_gpot_height(gmet_height)), 2)
    return g_naught * math.pow(r_naught / (r_naught + gmet_height), 2)


def get_gpot_height(gmet_height):
    '''Return the Geopotential height for a given geometric altitude
    based on equation 19
    '''
    return (r_naught * gmet_height) / (r_naught + gmet_height)

def get_gmet_height(gpot_height):
    '''Return the Geometric height for a given Geopotential altitude
    based on equation 19
    '''
    return (r_naught * gpot_height) / (r_naught - gpot_height)

def get_temp(gpot_height):
    '''Return the molecular-scale temperature in Kelvin, given a geopotential gpot_height in meters
    based on equation 23
    equal to kinetic temperature up to 80km
    implemented recursively
    '''
    if gpot_height <= heights[1]:
        return temp_naught + l_temp_gradients[0] * gpot_height
    elif gpot_height <= heights[2]:
        return get_temp(heights[1]) + l_temp_gradients[1] * (gpot_height - heights[1])
    elif gpot_height <= heights[3]:
        return get_temp(heights[2]) + l_temp_gradients[2] * (gpot_height - heights[2])
    elif gpot_height <= heights[4]:
        return get_temp(heights[3]) + l_temp_gradients[3] * (gpot_height - heights[3])
    elif gpot_height <= heights[5]:
        return get_temp(heights[4]) + l_temp_gradients[4] * (gpot_height - heights[4])
    elif gpot_height <= heights[6]:
        return get_temp(heights[5]) + l_temp_gradients[5] * (gpot_height - heights[5])
    elif gpot_height <= heights[7]:
        return get_temp(heights[6]) + l_temp_gradients[6] * (gpot_height - heights[6])

#def old_get_pressure(gmet_height):
#    '''Return the pressure given a geometric height in meters
#    based on equation 33a and 33b
#    use 33A when temp gradient is not 0
#    use 33B when temp gradient is 0 
#    
#    temp gradients defined in l_temp_gradients
#    per page 9, the molecular mass is essentially the same for each molecule up to 80km
#    '''
#
#    gpot_height = get_gpot_height(gmet_height)
#    grav_accel = get_grav_accel(gmet_height)
#    temp = get_temp(gpot_height)
#
#
#    if gpot_height <= heights[1]:
#        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[0])
#        return p_naught * math.pow(temp / (temp + l_temp_gradients[0] * (gpot_height - heights[0])), exponent)
#    elif gpot_height <= heights[2]:
#        return get_pressure(heights[1]) * math.exp((-grav_accel) * M_dry * (gpot_height - heights[1]) / (r_gas_constant * temp))
#    elif gpot_height <= heights[3]:
#        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[2])
#        return get_pressure(heights[2]) * math.pow(temp / (temp + l_temp_gradients[2] * (gpot_height - heights[2])), exponent)
#    elif gpot_height <= heights[4]:
#        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[3])
#        return get_pressure(heights[3]) * math.pow(temp / (temp + l_temp_gradients[3] * (gpot_height - heights[3])), exponent)
#    elif gpot_height <= heights[5]:
#        return get_pressure(heights[4]) * math.exp((-grav_accel) * M_dry * (gpot_height - heights[4]) / (r_gas_constant * temp))
#    elif gpot_height <= heights[6]:
#        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[5])
#        return get_pressure(heights[5]) * math.pow(temp / (temp + l_temp_gradients[5] * (gpot_height - heights[5])), exponent)
#    elif gpot_height <= heights[7]:
#        exponent = grav_accel * M_dry / (r_gas_constant * l_temp_gradients[6])
#        return get_pressure(heights[6]) * math.pow(temp / (temp + l_temp_gradients[6] * (gpot_height - heights[6])), exponent)

def get_pressure(gmet_height):
    '''Return the pressure given a geometric height in meters
    based on NASA page http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
    '''
    gpot_height = get_gpot_height(gmet_height)
    grav_accel = get_grav_accel(gmet_height)
    temp = get_temp(gpot_height)

    sector = get_sector(gpot_height)

    if sector == 0:
        return 101.29 * math.pow((temp) / 288.08, 5.256)
    elif sector == 1:
        return 22.65 * math.exp(1.73 - .000157 * gpot_height)
    else:
        return 2.488 * math.pow((temp) / 216.6, -11.388)

def get_sector(gpot_height):
    '''Return the sector of the atmosphere the height falls into, according to Table 4
    see list heights'''
    if gpot_height <= heights[1]:
        return 0
    elif gpot_height <= heights[2]:
        return 1
    elif gpot_height <= heights[3]:
        return 2
    elif gpot_height <= heights[4]:
        return 3
    elif gpot_height <= heights[5]:
        return 4
    elif gpot_height <= heights[6]:
        return 5
    elif gpot_height <= heights[7]:
        return 6

def get_thermal_conductivity(temperature):
    '''Return a coefficient of thermal conductivity for air given a temperature
    based on equation 53
    '''
    return (2.64638 * math.pow(10, -3) * math.pow(temperature, 3.0 / 2.0)) / (temperature + (245.4 * math.pow(10, -12.0 / temperature)))

#Perfect gas is valid below 86km altitude
def perfect_gas(n_moles_of_gas, volume, temperature):
    '''Return P = NRT/V'''
    return (n_moles_of_gas * r_gas_constant * temperature) / volume

def get_density_of_air(gmet_height, temperature=None):
    '''Return density of air
    takes a value for temperature of air, or if blank calculates it from height
    based on NASA page http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
    '''
    gpot_height = get_gpot_height(gmet_height)
    if temperature is None:
        temperature = get_temp(gpot_height)

    return get_pressure(gmet_height) / (0.2869 * temperature)

def get_weight(density, volume, gmet_height):
    '''Return the weight in Newtons of the object
    F= G*(m1*m2)/r^2
    G= 6.67384 * 10^-11
    '''
    mass = density * volume
    weight = universal_grav * (mass_earth) / math.pow(r_naught + gmet_height, 2) * mass
    return weight

def x_squared(x):
    return x * x
def deriv(x):
    return 2 * x
def integ(x):
    return 1.0 / 3 * x * x * x

def num_int0(func, start, stop, step_size=.1):
    sum = 0
    for x in range(int((stop - start) / step_size)):
        sum = sum + (func(x * step_size + start) * step_size)
    return sum

def num_int1(func, start, stop, step_size=.1):
    sum = 0
    for x in range(int((stop - start) / step_size)):
        sum = sum + .5 * (func(x * step_size + start) * step_size) + .5 * (func((x + 1) * step_size + start) * step_size)
    return sum

def num_int2(func, start, stop, step_size=.1):
    sum = 0
    for x in range(int((stop - start) / step_size)):
        k1 = step_size * func(x * step_size + start)
        k2 = step_size * func((x + .5) * step_size + start)
#        k3= step_size*func((x+.5)*step_size+start)
        k3 = k2
        k4 = step_size * func((x + 1) * step_size + start)

        sum = sum + (k1 / 6.0) + (k2 / 3.0) + (k3 / 3.0) + (k4 / 6.0)
    return sum



envelope_radius = 500.0
envelope_thickness = 0.01

envelope_volume = math.pi * 4 / 3 * math.pow(envelope_radius, 3)

altitude = 1000 #this is the geometric elevation of the center of the envelope

def get_cylinder_weight(gmet_height, material, cyl_height, cyl_radius, temperature):
    '''Return the weight of a vertically oriented cylinder'''
#    print "height is ", str(gmet_height)
#    print "material is ", str(material)
#    print "cyl_height is ", str(cyl_height)
#    print "cyl_radius is ", str(cyl_radius)
#    print "temperature is ", str(temperature)
    
    volume = math.pi * cyl_radius * cyl_radius * cyl_height

    if material == "air":
        density = get_density_of_air(gmet_height, temperature)
    else:
        density = densities[material]

    return get_weight(density, volume, gmet_height)

def get_sphere_weight(gmet_height, radius, temperature_function, material=None):
    '''Return the weight of a sphere with center at gmet_height, of material, of radius (meters), with temperature being a function defined by gmet height
    
    '''
    radius = int(radius)

    if material is None:
        material = "air"

    step_size = .1 # step size in meters
    num_steps = int((2 * radius) / step_size)
    #starts integrating at bottom of sphere, goes to top of sphere
    #has a cylinder at each height
    #treats bottom of sphere as angle pi, goes to angle 0 at top of sphere


    #cos theta times radius is how far along the height of the sphere we are, so use acos to find theta
    #sin theta times radius is the cylinder radius of the slice at that height in the sphere
#    get_theta = lambda step: math.acos(-1+ (x / num_steps) *math.pi)
#    theta = get_theta(x)
#    get_cyl_radius= lambda step: math.sin(math.acos(-1+ (x / num_steps) *math.pi))* radius
#    temperature = temperature_function(gmet_height)
    
    
#    for step in range(1000):
#        print -1 + (float(step) / num_steps) * 2, math.acos(-1 + (float(step) / num_steps) * 2), math.sin(math.acos(-1 + (float(step) / num_steps) * 2)) * radius
#    get_cyl_weight = lambda step: get_cylinder_weight(gmet_height + (-1 + float(step) / num_steps) * 2 * radius, material, step_size, math.sin(math.acos(-1 + (float(step) / num_steps) * 2)) * radius, temperature_function(gmet_height + (-1 + float(step) / num_steps) * 2 * radius))

    
    def get_cyl_weight(iterative_height):
        '''Return the weight of a specific cylinder, given the cylinder's location in the sphere'''
        vertical_distance_from_center= iterative_height - gmet_height
        base = vertical_distance_from_center/radius
        
        return get_cylinder_weight(gmet_height + vertical_distance_from_center, material, step_size, math.sin(math.acos(base)) * radius, temperature_function(gmet_height + vertical_distance_from_center))
    
    

    return num_int2(get_cyl_weight, gmet_height - radius, gmet_height + radius, step_size)





#print "%10s, %10s, %6s, %10s, %10s" % ("GMet", "Gpot", "Gaccel", "Pressure", "temp")
#for H in range(0, 50000, 500):
#
#    print "%10.2f, %10.2f, %6.2f, %10.2f, %10.2f," % (H, get_gpot_height(H), get_grav_accel(H), get_pressure(H), get_temp(get_gpot_height(H)))

#print get_sphere_weight(10000, 500, get_temp, "air")
rad = 500
print "%10s %20.2f" % ( "air", get_sphere_weight(500, rad, get_temp, "air"))
vol = 4.0/3*math.pi*rad*rad*rad
for material in densities.keys():
    print "%10s, %20.2f, %20.2f" % ( material, get_sphere_weight(10000, rad, get_temp, material), get_weight(densities[material], vol, 10000))








#start = 0
#stop = 18
#
#print num_int0(x_squared, start,stop)
#print num_int1(x_squared, start,stop)
#print num_int2(x_squared, start,stop)
#print integ(stop)-integ(start)
