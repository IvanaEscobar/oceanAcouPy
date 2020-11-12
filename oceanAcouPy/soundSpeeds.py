# EE348N: Ocean Acoustics 
# Sound Speeds: 
# Sound speed profiles used in class
# Primary author: Ivana Escobar
# Last Edited: 24 Oct 2020

### PACKAGES ###
from numpy import angle, log10, linspace, array, ndarray, arange,\
                sin, cos, arccos, arctan, sqrt, exp, pi, nan

### FUNCTIONS ###
def soundSpeed (t,s,z, lat=0, eqn='mackenzie81'):
# Inputs:
#   t : tempurature [degC]
#   s : salinity [ppt]
#   z : column depth [m]
#   lat : latitude [deg]
#   eqn : mackenzie81, leroy08
# Returns:
#   c : seawater sound speed [m/s]    
    if (eqn=='mackenzie81'):
        return  1448.96 + 4.59*t - 0.05304*t**2 + 2.374e-4*t**3 + \
                (1.340 - 0.01025*t)*(s - 35) + 0.01630*z + \
                1.675e-7*z**2 - 7.139e-13*t*z**3
    elif (eqn=='leroy08'):
        return  1402.5  +5*t - 5.44e-2*t**2 + 2.1e-4*t**3 + 1.33*s - \
                1.23e-2*t*s + 8.7e-5*t**2*s + 1.56e-2*z + \
                2.55e-7*z**2 - 7.3e-12*z**3 + 1.2e-6*z*(lat-45) - \
                9.5e-13*t*z**3 + 3e-7*t**2*z + 1.43e-5*s*z

def cMunk (z):
# Inputs:
#   z : depth [m]
# Returns:
#   Munk sound speed profile [m/s]
    c0 = 1500.
    eps = 0.00737
    zt=2*(z-1300)/1300
    return c0*( 1 + eps*(zt - 1 + exp(-zt)) )

def cDD (z): 
# Inputs:
#   z : depth [m]
# Returns:
#   double-ducted sound speed profile: 
    return 10*sin(f1(z)*(5000.-z)) +\
                        70*sin(f2(z)*(5000.-z)) +\
                        0.014*z + 1480.
f1 = lambda z : -3e-18*(5000-z)**4
f2 = lambda z :  3e-19*(5000-z)**4
