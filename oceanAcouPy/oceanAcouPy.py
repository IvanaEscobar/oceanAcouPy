# EE348N: Ocean Acoustics 
# functions designated for insruction of semester topics course 
# in ocean acoustics. 
# Primary author: Ivana Escobar
# Last Edited: 24 Oct 2020

### FUNCTIONS ###
def soundSpeed (t,s,z, lat=0, eqn='mackenzie81'):
# Inputs:
#   t : tempurature [degC]
#   s : salinity [ppt]
#   z : column depth [m]
#   lat : latitude [deg]
#   eqn : mackenzie81, leroy08
# Returns:
#   c : seawater sound speed [m/s]    
    if (eqn=='mackenzie81'):
        return  1448.96 + 4.59*t - 0.05304*t**2 + 2.374e-4*t**3 + \
                (1.340 - 0.01025*t)*(s - 35) + 0.01630*z + \
                1.675e-7*z**2 - 7.139e-13*t*z**3
    elif (eqn=='leroy08'):
        return  1402.5  +5*t - 5.44e-2*t**2 + 2.1e-4*t**3 + 1.33*s - \
                1.23e-2*t*s + 8.7e-5*t**2*s + 1.56e-2*z + \
                2.55e-7*z**2 - 7.3e-12*z**3 + 1.2e-6*z*(lat-45) - \
                9.5e-13*t*z**3 + 3e-7*t**2*z + 1.43e-5*s*z

def cMunk (z):
# Inputs:
#   z : depth [m]
# Returns:
#   Munk sound speed profile [m/s]
    c0 = 1500.
    eps = 0.00737
    zt=2*(z-1300)/1300
    return c0*( 1 + eps*(zt - 1 + exp(-zt)) )

# double-ducted sound speed profile: 
#     z : depth [m]
cDD = lambda z : 10*sin(f1(z)*(5000.-z)) +\
                        70*sin(f2(z)*(5000.-z)) +\
                        0.014*z + 1480.
f1 = lambda z : -3e-18*(5000-z)**4
f2 = lambda z :  3e-19*(5000-z)**4

def K (pH):
# Inputs:
#   pH : acid-base ratio [1]
# Returns:
#   K : logarithmic transformation of pH [1]
    return 10**(pH-8)

def viscAlpha (f,t,z):
# Inputs:
#   f : frequency [Hz]
#   t : tempurature [degC]
#   z : column depth [m]
# Returns:
#   alpha : attenuation contribution from viscosity [dB/m]
    return 4.9e-10*exp(z/1.7e4 - t/27) * f**2

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

def alphaTotal (f,t,s,z,pH):
    return viscAlpha(f,t,z) + boricAlpha(f,t,s,pH) \
                + magnesiumAlpha(f,t,s,z)

alpha = lambda f : alphaTotal(f,T,S,D,pH)
av = lambda f : viscAlpha(f,T,D)
aB = lambda f : boricAlpha(f,T,S,pH)
aMg = lambda f : magnesiumAlpha(f,T,S,D)

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

def RefC_FB(theta1,c1,rho1,c2,rho2):
    theta2 = arccos(c2/c1*cos(theta1))
    R = (rho2*c2/sin(theta2) - rho1*c1/sin(theta1))/(rho2*c2/sin(theta2) + rho1*c1/sin(theta1))
    if (theta1 == 0):
        R=1
    return R

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

def thetac (c1,c2):
# Inputs:
#   c1 : sound speed of top material [m/s]
#   c2 : sound speed of bottom material [m/s]
# Returns:
#   critical angle : [deg]
    return arccos( c1/c2 ) * 180/pi

# wavenumber k:
#    c: sound speed, f: frequency
kP = lambda c,f : 2*pi*f / c

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

def bPhase (theta, f, c, r):
    return 2*pi*f/c * r * sin(theta * pi / 180)

def beamPattern ( theta, thetaS, f, c, r, noise ):
    steerPhase = bPhase(thetaS, f, c, r)
    beamResponse = []
    for t in theta:
        lookPhase = bPhase(t, f, c, r)
        beamResponse.append( sum( exp(-1j*steerPhase)*\
                                 (exp(1j*lookPhase) + noise) ) )
    return 20*log10( abs(array(beamResponse))/len(r)  )

def beamPatternIntegrand ( theta, thetaS, f, c, r, noise ):
    steerPhase = bPhase(thetaS, f, c, r)
    lookPhase = bPhase(theta, f, c, r)
    return abs( sum( exp(-1j*steerPhase)*\
                        (exp(1j*lookPhase) + noise) ) / len(r) )**2

def arrayGain (BFIntegral):
    return 10*log10(2*pi / BFIntegral)
