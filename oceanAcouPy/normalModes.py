# EE348N: Ocean Acoustics 
# Normal Mode solutions

from numpy import log10, tan, sqrt, pi
from .reflCoeff import cComplex, RefC_FB, RefC_FB_shear
from .soundSpeed import cDD

### FUNCTIONS ###
k = lambda f: 2*pi*f / cw
krn = lambda f,n: sqrt(k(f)**2 - kzn(n)**2)

#-------------------------------------------------------------------------------
def phaseSpeed (f, n):
# Returns:
#   phase speed for ideal waveguide
    return 2*pi*f / krn(f,n)

#-------------------------------------------------------------------------------
def groupSpeed (f, n):
# Returns:
#   group speed for ideal waveguide
    return cw**2 * krn(f,n) / (2*pi*f)

#-------------------------------------------------------------------------------
def TL (r,z,M,f):
# Returns:
#   Transmission loss for ideal waveguide
    # M number of modes that fit in the waveguide
    psiSum = sum( [Psi(n,zs)*Psi(n,z)*\
            exp(1j*krn(f,n)*r)/sqrt(krn(f,n))  for n in M] )
    return -20*log10( abs( sqrt(2*pi/r)/rhow*psiSum ) )

#-------------------------------------------------------------------------------
kP = lambda c,f : 2*pi*f / c

#-------------------------------------------------------------------------------
# vertical wavenumbers for trapped modes ONLY
kzP = lambda kr,c,f : sqrt(kP(c,f)**2 - kr**2)

#-------------------------------------------------------------------------------
# Characteristic Equation for waveguide with lossless bottom
CE = lambda kr,f : tan(kzP(kr,cw,f)*D) +\
                complex(0, rhob/rhow*kzP(kr,cw,f)/kzP(kr,cb,f) )
