# EE348N: Ocean Acoustics 
# functions designated for insruction of semester topics course 
# in ocean acoustics. 
# Author: Ivana Escobar
# Last Edited: 24 Oct 2020

### FUNCTIONS ###

def Theta2 (c1, c2, theta1):
# Inputs:
#   c1 : sound speed of top fluid [m/s]
#   c2 : sound speed of bottom fluid [m/s]
#   theta1: grazing angle [deg] --> trig fxns use [rad]
# Returns:
#   theta2: transmission angle [deg]
    return arccos(c2 * cos(pi*theta1/180.) / c1) * 180./pi

def cComplex(cr, alpha):
# Inputs:
#   cr : real part of sound speed [m/s]
#   alpha : attenuation [dB/lambda]
# Returns:
#   c : lossy sound speed [m/s] (complex)
    return complex( cr, -alpha*cr/(17.372*pi) )

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

def RefC (mat1, mat2, thetaList):
# Inputs:
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

    theta2 = []
    thetas = []
    z1 = []; z2 = []
    for theta in thetaList:
        theta2.append( Theta2(c1,c2,theta) )
        z1.append( Zimp(rho1,c1,theta) )
        if (cs == 0.):
            z2.append( Zimp(rho2,c2,theta2[-1]) )
        else:
            thetas.append( Theta2(c1,cs,theta) )
            z2.append( Ztot(rho2,cs,thetas[-1],\
                            c2,theta2[-1]) )
    z1 = array(z1)
    z2 = array(z2)
    return (z2 - z1) / (z1 + z2)

def phaseDeg (arr):
# Inputs:
#   arr: array of complex data type
# Returns:
#   phase: array of phase of num [deg]
#          ie. math.atan2(num.imag, num.real)
    ph = []

    for num in arr:
        ph.append( angle(num, deg=True) )
    ph = array(ph)
    return ph

def BL (refc):
# Inputs:
#   refc : Rayleigh's relfection coefficient [1]
# Returns:
#   BL : bottom loss [dB]
    return -20*log10( abs(refc) )

def intromissionCheck (mat1):
# Inputs:
#   mat1 : material name [string] from a Pandas DataFrame
#          cc : complex sound speed [m/s]
#          rho: density [kg/m3] or [g/cm3]
# Returns:
#   statement for angle of intromission
    c1 = df.loc[mat1,'cc']
    rho1 = df.loc[mat1,'rho']

    for rho2,c2,name in \
    zip(df.loc[:,'rho'], df.loc[:,'cc'], df.index):
        if (c1 > c2 and rho1*c1 < rho2*c2):
            print('%s has an angle of intromission' %name)
            print('\ttheta_I is %0.3f' %thetaI(c1,rho1,c2,rho2).real)
    return None

def thetaI (c1,r1, c2,r2):
# Inputs:
#   c1 : sound speed of top material [m/s]
#   r1 : density of top material [g/cm3]
#   c2 : sound speed of bottom material [m/s]
#   r2 : density of bottom material [g/cm3]
# Returns:
#   intromission angle : [deg]
    return arctan( sqrt( (1-(c2/c1)**2)/\
                        ((r2*c2/r1/c1)**2-1) ) ) * 180/pi

def criticalangleCheck (mat1):
# Inputs:
#   mat1 : material name [string] from a Pandas DataFrame
#          cc : complex sound speed [m/s]
#          rho: density [kg/m3] or [g/cm3]
# Returns:
#   statement with critical angle
    c1 = df.loc[mat1,'cc']
    rho1 = df.loc[mat1,'rho']

    for rho2,c2,name in \
    zip(df.loc[:,'rho'], df.loc[:,'cc'], df.index):
        if (c1 < c2):
            print('%s has a critical angle' %name)
            print('\ttheta_c is %0.3f'%thetac(c1,c2))
    return None

def thetac (c1,c2):
# Inputs:
#   c1 : sound speed of top material [m/s]
#   c2 : sound speed of bottom material [m/s]
# Returns:
#   critical angle : [deg]
    return arccos( c1/c2 ) * 180/pi
