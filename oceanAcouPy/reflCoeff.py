# EE348N: Ocean Acoustics 
# Relfection coefficient for materials with and without shear components.

from numpy import log10, array, sin, cos, arccos, pi, nan

### FUNCTIONS ###
def Theta2 (c1, c2, theta1):
# Inputs:
#   c1 : sound speed of top fluid [m/s]
#   c2 : sound speed of bottom fluid [m/s]
#   theta1: grazing angle [deg] --> trig fxns use [rad]
# Returns:
#   theta2: transmission angle [deg]
    return arccos(c2 * cos(pi*theta1/180.) / c1) * 180./pi

#-------------------------------------------------------------------------------
def cComplex(cr, alpha):
# Inputs:
#   cr : real part of sound speed [m/s]
#   alpha : attenuation [dB/lambda]
# Returns:
#   c : lossy sound speed [m/s] (complex)
    return complex( cr, -alpha*cr/(17.372*pi) )

#-------------------------------------------------------------------------------
def Zimp (rho,c,theta):
# Inputs:
#   rho : density [kg/m3] or [g/cm3]
#   c : sound speed [m/s]
#   theta: angle [deg] --> trig fxns use [rad]
# Returns:
#   Z: effective impedance of fluid [kg/m2/s]
    if (abs(sin(pi*theta/180.)) == 0):
        return complex(nan)
    else:
        return rho*c/sin(pi*theta/180.)

#-------------------------------------------------------------------------------
def Ztot (rho,cs,thetas,cp,thetap):
# Inputs:
#   rho : density [kg/m3] or [g/cm3]
#   cs : shear sound speed [m/s]
#   thetas: shear angle [deg] --> trig fxns need [rad]
#   cp : compressional sound speed [m/s]
#   thetap: compressional angle [deg] --> trig fxns need [rad]
# Returns:
#   Ztot: effective impedance of solid [kg/m2/s]
    Zp = Zimp(rho,cp,thetap)
    Zs = Zimp(rho,cs,thetas)
    return Zp*cos(thetas*pi/90)**2 + Zs*sin(thetas*pi/90)**2

#-------------------------------------------------------------------------------
def RefC (df, mat1, mat2, thetaList):
# Inputs:
#   df : pandas DataFrame
#   mat1 : material name [string] from a Pandas DataFrame
#          c : sound speed [m/s]
#          rho: density [kg/m3]
#          alpha : attenuation [dB/lambda]
#   mat2 : material name [string] from a Pandas DataFrame
#          c : sound speed [m/s]
#          rho: density [kg/m3]
#          alpha : attenuation [dB/lambda]
#          cs : shear sound speed [m/s]
#          alphas : shear attenuation [dB/lambda]
#   thetaList : grazing angle [deg]
# Returns:
#   Rayleigh's relfection coefficient [1] (complex)
    c1 = complex(df.loc[mat1,'c'])
    if (mat1 != 'water'):
        c1 = cComplex(df.loc[mat1,'c'], df.loc[mat1,'alpha'])
    c2 = cComplex(df.loc[mat2,'c'], df.loc[mat2,'alpha'])
    df.loc[mat1,'cc']=c1
    df.loc[mat2,'cc']=c2

    rho1 = df.loc[mat1,'rho']
    rho2 = df.loc[mat2,'rho']

    cs = cComplex(df.loc[mat2,'cs'], df.loc[mat2,'alphas'])
    df.loc[mat2,'csc']=cs

    z1 = []; z2 = []
    for theta in thetaList:
        theta2 = Theta2(c1,c2,theta) 
        z1.append( Zimp(rho1,c1,theta) )
        if (cs == 0.):
            z2.append( Zimp(rho2,c2,theta2) )
        else:
            thetas = Theta2(c1,cs,theta)
            z2.append( Ztot(rho2,cs,thetas,\
                            c2,theta2) )
    z1 = array(z1)
    z2 = array(z2)
    return (z2 - z1) / (z1 + z2)

#-------------------------------------------------------------------------------
def RefC_FB(theta1,c1,rho1,c2,rho2):
    theta2 = arccos(c2/c1*cos(theta1))
    R = (rho2*c2/sin(theta2) - rho1*c1/sin(theta1))/\
        (rho2*c2/sin(theta2) + rho1*c1/sin(theta1))
    if (theta1 == 0):
        R=1
    return R

#-------------------------------------------------------------------------------
def RefC_FB_shear(theta1,c1,rho1,c2,rho2,cs):
    theta2 = arccos(c2/c1*cos(theta1))
    thetas = arccos(cs/c1*cos(theta1))

    Zp = rho2*c2/sin(theta2)
    Zs = rho2*c2/sin(thetas)
    Ztot = Zp*cos(2*thetas)**2 + Zs*sin(2*thetas)**2

    R = (Ztot - rho1*c1/sin(theta1))/(Ztot + rho1*c1/sin(theta1))
    if (theta1 == 0):
        R=1
    return R

#-------------------------------------------------------------------------------
def phi(f,c,h,theta):
# Inputs:
#   f : frequency [Hz]
#   c : sound speed [m/s]
#   h : thickness of middle layer [m]
#   theta : grazing angle [deg] (array)
# Returns:
#   phi2 : phase delay of middle layer [deg]
    return 2*pi*f*h*sin(theta*pi/180) / c

#-------------------------------------------------------------------------------
def Rlayered (data, mat1,mat2,mat3, h2, theta1, f):
# Inputs:
#   data : pandas dataframe containing ocean bottom properties
#         c : sound speed [m/s]
#         rho : density [g/cm3]
#   mat1 : material name [string]
#   mat2 : material name [string]
#   mat3 : material name [string]
#   h2 : thickness of middle layer [m]
#   theta1 : grazing angle [deg] (array)
#   f : frequency [Hz]
# Returns:
#   R : Rayleigh's relfection coefficient [1] (array)
    c1 = data.loc[mat1,'cc']
    rho1 = data.loc[mat1,'rho']
    c2 = data.loc[mat2,'cc']
    rho2 = data.loc[mat2,'rho']
    c3 = data.loc[mat3,'cc']
    rho3 = data.loc[mat3,'rho']

    RR = []
    for theta in theta1:
        theta2=( Theta2(c1,c2,theta) )
        theta3=( Theta2(c1,c3,theta) )
        z1=( Z(rho1,c1,theta ) )
        z2=( Z(rho2,c2,theta2) )
        z3=( Z(rho2,c3,theta3) )
        phi2=( phi(f,c2,h2,theta2) )
        RR.append( complex( z2*(z3-z1),\
                             (z1*z3 - z2**2)*tan(phi2) )/\
                     complex( z2*(z3+z1),\
                             (z1*z3 - z2**2)*tan(phi2) ) )
    return array(RR)

#-------------------------------------------------------------------------------
def layerFrequencies(m,c,h, wl='quarter'):
    frequency = c*(2.*m-1)/(4.*h)
    if (wl == 'half'):
        return frequency*2.
    return frequency
