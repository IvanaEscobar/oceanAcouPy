# EE348N: Ocean Acoustics 
# Beamforming to find array gain

from numpy import log10, array, sin, exp, pi

### FUNCTIONS ###
def bPhase (theta, f, c, r):
    return 2*pi*f/c * r * sin(theta * pi / 180)

#-------------------------------------------------------------------------------
def beamPattern ( theta, thetaS, f, c, r, noise ):
    steerPhase = bPhase(thetaS, f, c, r)
    beamResponse = []
    for t in theta:
        lookPhase = bPhase(t, f, c, r)
        beamResponse.append( sum( exp(-1j*steerPhase)*\
                                 (exp(1j*lookPhase) + noise) ) )
    return 20*log10( abs(array(beamResponse))/len(r)  )

#-------------------------------------------------------------------------------
def beamPatternData ( theta, dataS, f, c, r, noise ):
# Same as beam pattern, but with steering data provided
    beamResponse = []
    for t in theta:
        lookPhase = phase(t, f, c, r)
        beamResponse.append( sum( dataS*\
                                 (exp(1j*lookPhase) + noise) ) )
    return 20*log10( abs(array(beamResponse))/len(r)  )

#-------------------------------------------------------------------------------
def beamPatternIntegrand ( theta, thetaS, f, c, r, noise ):
    steerPhase = bPhase(thetaS, f, c, r)
    lookPhase = bPhase(theta, f, c, r)
    return abs( sum( exp(-1j*steerPhase)*\
                        (exp(1j*lookPhase) + noise) ) / len(r) )**2

#-------------------------------------------------------------------------------
def beamPatternWindow ( thetaS, theta, f, c, r, noise, window ):
# window is generated using a SciPy signal function
    steerPhase = phase(thetaS, f, c, r)
    beamResponse = []
    for t in theta:
        lookPhase = phase(t, f, c, r)
        beamResponse.append( sum( window*exp(-1j*steerPhase)*\
                                 (exp(1j*lookPhase) + noise) ) )
    return 20*log10( abs(array(beamResponse))/ sum(window) )

#-------------------------------------------------------------------------------
def arrayGain (BFIntegral):
    return 10*log10(abs(2*pi / BFIntegral))
