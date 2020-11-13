# EE348N: Ocean Acoustics 
# Find bottom loss. determine critical angle, angle of intromission, and shear
# critical angle

from numpy import angle, log10, array, arccos, arctan, sqrt, pi

### FUNCTIONS ###
def phaseDeg (arr):
# Inputs:
#   arr: array of complex data type
# Returns:
#   phase: array of phase of num [deg]
#          ie. math.atan2(num.imag, num.real)
    return array( [angle(num, deg=True) for num in arr] )

#-------------------------------------------------------------------------------
def BL (refc):
# Inputs:
#   refc : Rayleigh's relfection coefficient [1]
# Returns:
#   BL : bottom loss [dB]
    return -20*log10( abs(refc) )

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
def thetac (c1,c2):
# Inputs:
#   c1 : sound speed of top material [m/s]
#   c2 : sound speed of bottom material [m/s]
# Returns:
#   critical angle : [deg]
    return arccos( c1/c2 ) * 180/pi

#-------------------------------------------------------------------------------
def intromissionCheck (df, mat1):
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

#-------------------------------------------------------------------------------
def criticalangleCheck (df, mat1):
# Inputs:
#   df : pandas DataFrame
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
#-------------------------------------------------------------------------------
def shearcriticalangleCheck (df, mat1):
# Inputs:
#   df : pandas DataFrame
#   mat1 : material name [string] from a Pandas DataFrame
#          cc : complex sound speed [m/s]
#          rho: density [kg/m3] or [g/cm3]
# Returns:
#   statement with shear critical angle
    c1 = df.loc[mat1,'cc']
    rho1 = df.loc[mat1,'rho']

    for rho2,cs2,name in \
    zip(df.loc[:,'rho'], df.loc[:,'csc'], df.index):
        if (c1 < cs2):
            print('%s has a shear critical angle' %name)
            print('\ttheta_c_s is %0.3f' %thetac(c1,cs2,'shear').real)
    return None
