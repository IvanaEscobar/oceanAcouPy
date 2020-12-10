# EE348N: Ocean Acoustics 
# Bubbles lecture

from numpy import log10 

### FUNCTIONS ###
def tsBubble (sigma):
    return 10*log10(abs(sigma))
#-------------------------------------------------------------------------------

def sigmaBubble (f, z, a):
    fb = 3.25*sqrt(1 + z*0.1) / a
    deltar = 0.0025*f**(1/3)
    return a*a / (((f/fb)**2-1) + deltar*1j)**2
