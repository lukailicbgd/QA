import us_std
import numpy as np

'''
This module provides some function to compute Rayleigh scattering
parameters in the atmosphere. It is based on

Anthony Bucholtz, "Rayleigh-scattering calculations for the terrestrial atmosphere", 
Applied Optics 34, no. 15 (May 20, 1995): 2765-2773.    
'''

def depolarization_factor(wavelength):
    '''Returns the depolarization factor for a specific wavelength.
    It uses the values given in Bucholtz (1995) and a linear interpolation 
    for values between them. Works for wavlenght between 200nm and 1100nm.
    
    Input: wavelength in nm
    Output: depolarization factor (no units) 
    
    Anthony Bucholtz, "Rayleigh-scattering calculations for the terrestrial atmosphere", 
    Applied Optics 34, no. 15 (May 20, 1995): 2765-2773.
    '''
    minvalue = np.min(wavelength) #In case the input is an array
    maxvalue = np.max(wavelength)
    if (minvalue < 200) or (maxvalue>1100):
        raise ValueError
    
    #Values from Bucholtz, with a point added in 1100nm to include Nd:YAG wl. in list.
    wl_original = np.array([0.200,0.205,0.210,0.215,0.220,0.225,0.230,0.240,0.250,0.260,0.270,
                   0.280,0.290,0.300,0.310,0.320,0.330,0.340,0.350,0.360,0.370,0.380,0.390,
                   0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.1])
    wl_nm = wl_original * 1000 #convert original data from micrometers to nanometers 
    
    pn= np.array([4.545,4.384,4.221,4.113,4.004,3.895,3.785,3.675,3.565,3.455,3.400,
                    3.289,3.233,3.178,3.178,3.122,3.066,3.066,3.010,3.010,3.010,2.955,2.955,
                    2.955,2.899,2.842,2.842,2.786,2.786,2.786,2.786,2.730,2.730,2.730,2.730,2.730, 2.730])*10**-2
    
    result = np.interp(wavelength, wl_nm, pn)
    return result

def kings_factor(wavelength):
    '''Caclulate the king's correction factor defined as 
    F_k = (6 + 3*rn)/(6 - 7*rn)
    where rn the depolarization factor.
    Input: Wavelength in nm
    Output: Kign correction factor F (no units)
    '''
    
    rn = depolarization_factor(wavelength)
    f_king = (6 + 3*rn)/(6 - 7*rn)
    return f_king

def n_air(wavelength):
    '''The refractive index of air at a specific wavelength. Calculated for 
    standard air at T = 15 degsC. 
    Input: Wavelenght in nm
    Output: Air refractive index (no untis)   
    wavelenght[nm]
    355            0.03010
    532            0.02842
    1064           0.02730
    '''
    
    minvalue = np.min(wavelength)
    if minvalue <230:
        print("Warning: The formula use for calculating air's refractive index is valid for wl>230nm")

    wl_micrometers = wavelength / 1000.0 # Convert nm to um 
    
    s = 1/wl_micrometers # the reciprocal of wavelength
    c1 = 5791817.0
    c2 = 238.0185
    c3 = 167909.0
    c4 =57.362
    ns = 1 + (c1/(c2 -s**2) + c3/(c4 -s**2))*10**-8
    return ns

def scattering_cross_section(wavelength):
    ''' Calculates the Rayleigh-scattering cross section per molecule.
    Input: Wavelength in nm.
    Output: Rayleigh-scattering cross section in  squared centimeters.
    '''
    
    #Just to be sure that the wavelength is a float
    wavelength = float(wavelength)
    #Molecular number density for standard air
    Ns = 2.54743*10**19 #in cm**-3
    
    #Wavelength of radiation
    wl_cm = wavelength *10**-7 #in cm
    
    #King's correction factor
    f_k = kings_factor(wavelength) # no units
    
    #Refractive index of air
    n = n_air(wavelength)
    
    #first part of the equation
    f1 = (24*np.pi**3)/(wl_cm**4*Ns**2)
    #second part of the equation
    f2 = (n**2 -1)**2/(n**2 +2)**2
    
    #result
    sigma = f1 * f2 *f_k
    return sigma

def phase_function(theta, wavelength):
    ''' Calculates the phase function at an angle theta for a specific wavelegth.
    Input:  theta in rads
            wavelength in nm
    
    The formula is derived from Bucholtz (1995). A different formula is given in 
    Miles(2001). 
    
    The use of this formula insetad of the wavelenght independent 3/4(1+cos(th)**2)
    improves the results for back and forward scatterring by ~1.5%
    
    Anthony Bucholtz, "Rayleigh-scattering calculations for the terrestrial atmosphere", 
    Applied Optics 34, no. 15 (May 20, 1995): 2765-2773.  

    R. B Miles, W. R Lempert, and J. N Forkey, "Laser Rayleigh scattering", 
    Measurement Science and Technology 12 (2001): R33-R51
    '''
    
    r = depolarization_factor(wavelength)
    gamma = r/(2-r)
    
    #first part of the equation
    f1 = 3/(4*(1 + 2 *gamma))
    #second part of the equation
    f2 = (1+3*gamma) + (1-gamma)* (np.cos(theta))**2
    #results
    p = f1 *f2
    
    return p

def scattering(numerical_density, wavelength):
#def scattering(h, wavelength):
    '''Calculates the molecular volume scattering coefficient (alpha) at 
    height h at a specific wavelength.Makes use of the standard atmosphere data.
    Formula: alpha = N * sigma
    Input: h in meters
           wavelength in nm
    Output: volume scattering coefficient in cm**-1
    '''
    
    alpha = numerical_density * scattering_cross_section(wavelength)
#   alpha = numerical_density(h) * scattering_cross_section(wavelength)
    return alpha

def angular_scattering(numerical_density,theta,wavelength):
#def angular_scattering(h,theta,wavelength):
    '''Calculates the angular molecular volume scattering coefficient (da/dW) at 
    height h, angle theta to the propagation direction, at a specific wavelength.
    Makes use of the standard atmosphere data.
    Formula: da/dW = alpha* Phase(theta)/4pi
    Input: h in meters
           theta in rad
           wavelength in nm
    Output: angular volume scattering coefficient in ???
    '''
    
    a = scattering(numerical_density,wavelength)
#   a = scattering(h,wavelength)
    phase = phase_function(theta, wavelength)
    da = a*phase/(4*np.pi)
    return da

def angular_scattering_cross_section(theta,wavelength):
    ''' Calculates the angular rayleigh scattering cross section per molecule.
    Input: theta in rad
           Wavelength in nm.
    
    Output: angular rayleigh-scattering cross section in ???
    '''
    
    phase = phase_function(theta, wavelength)/(4*np.pi)
    dsigma = scattering_cross_section(wavelength) * phase
    return dsigma
