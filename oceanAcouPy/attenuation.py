# EE348N: Ocean Acoustics 
# Attenuation built from visuous, boric acid, and magnesium sulfate 
# contributions

from numpy import ndarray, sqrt, exp

### FUNCTIONS ###
def K (pH):
# Inputs:
#   pH : acid-base ratio [1]
# Returns:
#   K : logarithmic transformation of pH [1]
    return 10**(pH-8)
#-------------------------------------------------------------------------------
def viscAlpha (f,t,z):
# Inputs:
#   f : frequency [Hz]
#   t : tempurature [degC]
#   z : column depth [m]
# Returns:
#   alpha : attenuation contribution from viscosity [dB/m]
    return 4.9e-10*exp(z/1.7e4 - t/27) * f**2

#-------------------------------------------------------------------------------
def boricAlpha (f,t,s,pH):
# Inputs:
#   f : frequency [Hz]
#   t : tempurature [degC]
#   s : salinity [ppt]
#   pH : acid-base ratio [1]
# Returns:
#   alpha : attenuation contribution from boric acid [dB/m]
    fB = 7.8e2*sqrt(s/35)*exp(t/26)
    if (type(pH) != ndarray and pH > 2): # convert pH to K
        return 1.06e-4*fB*K(pH)**0.776 * f**2 / (f**2 + fB**2)
    else: # pH is in units of K already
        return 1.06e-4*fB*pH**0.776 * f**2 / (f**2 + fB**2)

#-------------------------------------------------------------------------------
def magnesiumAlpha (f,t,s,z):
# Inputs:
#   f : frequency [Hz]
#   t : tempurature [degC]
#   s : salinity [ppt]
#   z : column depth [m]
# Returns:
#   alpha : attenuation contribution from Mg [dB/m]
    fMg = 4.2e4*exp(t/17)
    return 5.2e-4*(1+t/43)*s/35*fMg*exp(-z/6e3) * f**2 / \
            (f**2 + fMg**2)

#-------------------------------------------------------------------------------
def alphaTotal (f,t,s,z,pH):
    return viscAlpha(f,t,z) + boricAlpha(f,t,s,pH) \
                + magnesiumAlpha(f,t,s,z)

#-------------------------------------------------------------------------------
alpha = lambda f : alphaTotal(f,T,S,D,pH)
av = lambda f : viscAlpha(f,T,D)
aB = lambda f : boricAlpha(f,T,S,pH)
aMg = lambda f : magnesiumAlpha(f,T,S,D)
