# EE348N: Ocean Acoustics 
# Detection Threshold (DT) from Albersheim equation (1981)

from numpy import log, log10, sqrt

### FUNCTIONS ###
def DT (pd, pfa, m):
    a = log(0.62/pfa)
    b = log(pd/(1-pd))
    return -5*log10(m) + \
            6.2+ 4.54/sqrt(m+0.44)*log10(a+0.12*a*b+1.7*b) #[dB]
