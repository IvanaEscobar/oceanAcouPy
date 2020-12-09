# EE348N: Ocean Acoustics 
# Sediment roughness equations for finding backscattering strength with 
# perturbation theory (PT) or Kirchhoff approximation (KA)

from numpy import exp, arange, linspace, array, floor,\
                cos, sin, tan, arccos, arctan, sqrt, pi, \
                where, amax, exp, log10, append

### FUNCTIONS ###
def braggK(kw, theta):
    return 2*kw*cos(theta * pi /180)
#-------------------------------------------------------------------------------

def powerLaw (w1, k, gamma):
    return w1/(k**gamma)
#-------------------------------------------------------------------------------

theta2 = lambda c1,c2,theta1: arccos(c2/c1 * cos(theta1))
#-------------------------------------------------------------------------------

def ptSigma(df, W,kw,thetal, seaFloor=True):
    c1 = df.loc['water','c']
    rho1 = df.loc['water','rho']
    c2 = df.loc['sediment','c']
    rho2 = df.loc['sediment','rho']
    theta2l = oa.reflCoeff.Theta2(c1,c2,thetal)
    # convert theta to radians
    thetal = thetal*pi/180; theta2l = theta2l*pi/180

    U = lambda th,th2 : ( (rho2*c2*cos(th) - rho1*c1*cos(th2))**2 \
                    + rho2*rho2*c2*c2 - rho1*rho1*c1*c1 ) / \
                    (rho2*c2*sin(th) + rho1*c1*sin(th2))**2
    if (seaFloor):
        u = U(thetal,theta2l)
    else:
        u = -1
    return 4*(kw*sin(thetal))**4*abs(u)**2*W
#-------------------------------------------------------------------------------

def kaSigma (df, slopeVar, thetal, seaFloor=True):
    # convert theta to radians
    thetal = thetal*pi/180
    if (seaFloor):
        r90 = oa.reflCoeff.RefC (df, 'water', 'sediment', [90])
    else:
        r90 = -1
    return abs(r90)**2 * exp( -1/(2*slopeVar*tan(thetal)**2) ) / \
            (8*pi**slopeVar*sin(thetal)**4)
#-------------------------------------------------------------------------------

def lambertsSigma(thetal):
    # convert theta to radians
    thetal = thetal*pi/180
    return 10**(-2.7)*sin(thetal)**2
#-------------------------------------------------------------------------------

def roughnessPiersonMosko (w, k):
    g = 9.81 # [m/s]
    Omega = sqrt(k*g)
    return 8.1e-3*g**2 * exp(-0.74*g/w/Omega)/Omega**5
#-------------------------------------------------------------------------------

def sigmaGS ( a ):
    return a**4/4
#-------------------------------------------------------------------------------

def sigmaRayBS ( k, a ):
    return 25*a*a*(k*a)**4 / 36
#-------------------------------------------------------------------------------

def TS ( k,a ):
    if (k*a > 1):
        #geometric solution
        return 10*log10(abs(sigmaGS(a)))
    else:
        return 10*log10(abs(sigmaRayBS( k,a )))
#-------------------------------------------------------------------------------
