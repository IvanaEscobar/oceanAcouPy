# EE348N: Ocean Acoustics 
# Method of Images for the pekeris waveguide

from numpy import log10, array, arange, arctan, sqrt, exp, pi
from .relfCoeff import cComplex, RefC_FB, RefC_FB_shear
from .soundSpeed import cDD

### FUNCTIONS ###
# wavenumber k:
#    c: sound speed, f: frequency
kP = lambda c,f : 2*pi*f / c

#-------------------------------------------------------------------------------
def Rnm (n,m,r,z,zs,D):
# Inputse
#   n : number of image group in sum (int)
#   m : number of image in image group (int)
#   r : receiver range [m]
#   z : receiver depth [m]
#   zs : source depth [m]
#   D : waveguide depth [m]
# Returns:
#   distance between image and receiver
    return sqrt(r**2 + Znm(n,m,z,zs,D)**2)

#-------------------------------------------------------------------------------
def Znm (n,m,z,zs,D):
# Inputs:
#   n : number of image group in sum (int)
#   m : number of image in image group (int)
#   z : receiver depth [m]
#   zs : source depth [m]
#   D : waveguide depth [m]
# Returns:
#   vertical distance between image source and receiver depth
    if (m==1):
        return z-zs-2*D*n
    if (m==2):
        return z+zs-2*D*(n+1)
    if (m==3):
        return z+zs+2*D*n
    if (m==4):
        return z-zs+2*D*(n+1)

#-------------------------------------------------------------------------------
def pressureMOIPekeris (r,z,df,mat1,mat2,f,zs,D,N):
# Inputs:
#   r : receiver range [m]
#   z : receiver depth [m]
#   df : pandas DataFrame
#   mat1 : material name [string] from a Pandas DataFrame
#   mat2 : material name [string] from a Pandas DataFrame
#   f : frequency [Hz]
#   zs : source depth [m]
#   D : waveguide depth [m]
#   N : arange of number of image groups (ints)
# Returns:
#   pressure derived by MOI for a Pekeris waveguide (complex)
    R1 = -1
    k1 = kP(df.loc[mat1,'c'], f)
    c1 = df.loc[mat1,'c']
    if (mat1 != 'water'):
        c1 = cComplex(df.loc[mat1,'c'], df.loc[mat1,'alpha'])
    c2 = cComplex(df.loc[mat2,'c'], df.loc[mat2,'alpha'])
    rho1 = df.loc[mat1, 'rho'] 
    rho2 = df.loc[mat2, 'rho'] 

    cs = cComplex(df.loc[mat2,'cs'], df.loc[mat2,'alphas'])

    p = []
    for n in N:
        theta2 = [arctan(abs(Znm(n,a,z,zs,D))/r ) for a in arange(4)+1]
        rnm = [Rnm(n,a,r,z,zs,D) for a in arange(4)+1]
        if (cs == 0.0):
            R2 = RefC_FB(theta2, c1, rho1, c2, rho2)
        else:
            R2 = RefC_FB_shear(theta2, c1, rho1, c2, rho2, cs)

        V1 = (R1*R2[0])**n
        V2 = R2[1]*(R1*R2[1])**n
        V3 = R1*(R1*R2[2])**n
        V4 = (R1*R2[3])**(n+1)

        p.append( V1*exp(1j*k1*rnm[0])/rnm[0] +\
                  V2*exp(1j*k1*rnm[1])/rnm[1] +\
                  V3*exp(1j*k1*rnm[2])/rnm[2] +\
                  V4*exp(1j*k1*rnm[3])/rnm[3] )
    return sum(p)

#-------------------------------------------------------------------------------
def TLMOIPekeris (r,z,df,mat1,mat2,f,zs,D,N):
# Inputs:
#   r : receiver range [m]
#   z : receiver depth [m]
#   df : pandas DataFrame
#   mat1 : material name [string] from a Pandas DataFrame
#   mat2 : material name [string] from a Pandas DataFrame
#   f : frequency [Hz]
#   zs : source depth [m]
#   D : waveguide depth [m]
#   N : arange of number of image groups (ints)
# Returns:
#   Transmission loss [dB]
    return -20*log10( abs(pressureMOIPekeris(r,z,df,mat1,mat2,f,zs,D,N)) )

#-------------------------------------------------------------------------------
def travelTimes (rayDf, IDlist):
# Inputs:
#   rayDf : pandas DataFrame
#   IDlist : list of indicies
# Returns:
#   travel times array [s]
    tau =[]
    for ii in IDlist:
        ds = array(rayDf[ii].loc[:,'s'][1:])-\
             array(rayDf[ii].loc[:,'s'][:-1])
        csi = cDD(rayDf[ii].loc[:,'ray'])
        tau.append(sum(ds/csi[:-1]))
    return array(tau)
