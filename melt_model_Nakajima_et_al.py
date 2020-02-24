import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import sys
import seaborn as sns
from pylab import scatter
import pylab
import scipy.optimize as op
from figure3 import test_figure3
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rc
import matplotlib as mpl
from scipy.interpolate import interp1d
#from mpl_toolkits.axes_grid1 import make_axes_locatable


#coefficients

Mmar=6.39e23
R0=1.5717e6
M0=6.39e22
a0=0.3412
a1=-8.90e-3
a2=9.1442e-4
a3=-7.4332e-5
GG=6.67408e-11
impact_angle = [0, 30, 60, 90]

cm_data = np.loadtxt("vik/vik.txt")
vik_map = LinearSegmentedColormap.from_list('vik', cm_data)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


# ---  input data  -------

#vimp = vesc
Mtotal = 0.999827862
Mtotal = 0.50
gamma = 0.1
vel = 1.0

#vimp > vesc
#Mtotal =  1.87793434
#gamma = 0.5
#vel =  1.29999995 


Mt =   (1.0 - gamma)*Mtotal  # target mass in the Martian mass
Mi = gamma*Mtotal # impactor mass in the Martian mass
ang = impact_angle[0] # this has to be chosen from 0, 30, 60, 90
#vel = 1.0 # this is normalized by the escape velocity
Lp = 5.2e6 # specific energy for melting 
Latentheat= 10e5 # latent heat
outputfigurename = 'output.eps' #output figure name


# --- end of input data ---


def radius(mass):
    lnMM0=np.log(mass/M0)
    gamma=a0+a1*lnMM0+a2*lnMM0**2+a3*lnMM0**3
    return R0*(mass/M0)**gamma


def gamma_calc(mass):
    lnMM0=np.log(mass/M0)
    gamma=a0+a1*lnMM0+a2*lnMM0**2+a3*lnMM0**3
    return gamma

def legendre(n,x):
    if n==0:
        return 1
    elif n==1:
        return x
    elif n==2:
        return 0.5*(3*x**2.0 - 1.0)
    elif n==3:
        return 0.5*(5*x**3.0-3*x)
    elif n==4:
        return 1.0/8.0*(35*x**4.0 - 30* x**2.0 + 3)
    elif n==5:
        return 1.0/8.0*(63*x**5.0 - 70 * x**3.0 - 15 * x)
    elif n==6:
        return 1.0/8.0*(231*x**6.0  - 315*x**4.0 + 105*x**2.0 - 5)

def Mantle_mass_model_highV(gamma,x):
    y=0.0
    num=1.0
    if (x <= np.pi/6.0):
        y= num
    elif ( x < np.pi/3.0):
        y= -(num-gamma)/(np.pi/6.0)*(x-np.pi/6.0)+num
    elif (x <= np.pi/2.0):
        y = gamma

    return y


def create_model(theta, r, x): #r is the radius, x is the angle
    y_model= theta[0]*r**(-2)*legendre(0,np.cos(x)) +  theta[1]*r**(-2)*legendre(1,np.cos(x)) +  theta[2]*r**(-2)*legendre(2,np.cos(x))\
             +  theta[3]*r**(-1)*legendre(0,np.cos(x)) +  theta[4]*r**(-1)*legendre(1,np.cos(x)) +  theta[5]*r**(-1)*legendre(2,np.cos(x))\
             +  theta[6]*r**(0)*legendre(0,np.cos(x))  +  theta[7]*r**(0)*legendre(1,np.cos(x))  +  theta[8]*r**(0)*legendre(2,np.cos(x)) \
             +  theta[9]*r**(1)*legendre(0,np.cos(x)) +  theta[10]*r**(1)*legendre(1,np.cos(x)) +  theta[11]*r**(1)*legendre(2,np.cos(x)) \
             +  theta[12]*r**(2)*legendre(0,np.cos(x)) +  theta[13]*r**(2)*legendre(1,np.cos(x)) +  theta[14]*r**(2)*legendre(2,np.cos(x))

    return y_model


#  ---- computing the structure of a planet ----
rho_P = [line.split() for line in open('rho_p.dat')]
rho_input=np.zeros(shape=(0,0))
P_input=np.zeros(shape=(0,0))

for m in range(1,len(rho_P)):
    rho_input = np.append(rho_input,float(rho_P[m][0])*1e3)
    P_input = np.append(P_input,float(rho_P[m][1])*1e9)

rho_P_function=interp1d(P_input, rho_input)

def compute_density(P):

    P=abs(P)
    if P < P_input[0]:
        return rho_input[0]
    else:
        return rho_P_function(P)




def compute_pressure(Mt, Rmelt):
    Rt=radius(Mt)
    Rmelt=Rmelt*Rt

    if Rmelt==1.0:
        return 0.0

    dr=(Rt-Rmelt)/100.0 #random number

    P=0.0
    r=Rt
    Mass=Mt


    while (r > Rmelt):
        rho=compute_density(P)
        P=P-rho*GG*Mass/r**2.0*dr
        Mass=Mass-4*np.pi*rho*r**2.0*dr
        r=r-dr

    if Rmelt==1.0:
        return 0.0
    else:
        return -P*1e-9

def compute_coreradius(Mt):
    Rt=radius(Mt)

    dr=Rt/100.0 #random number

    P=0.0
    r=Rt
    Mass=Mt
    CoreMass=0.3*Mt


    while (Mass >  CoreMass):
        rho=compute_density(P)
        P=P-rho*GG*Mass/r**2.0*dr
        Mass=Mass-4*np.pi*rho*r**2.0*dr
        r=r-dr
    
    return r/Rt

# --- end of computing the structure of a planet


if vel == 1.0:
    vimp_vesc = 0
else:
    vimp_vesc = 1


Mt = Mt*Mmar
Mi = Mi*Mmar
Rt = radius(Mt) 
Ri = radius(Mi)
Rti = radius(Mt + Mi)
vesc = np.sqrt(2.0*GG*(Mt + Mi)/(Rt + Ri))
ratio = Mi/Mt
ang = ang/180.0*np.pi
targetmassfraction = Mt/(Mt + Mi)


if vimp_vesc==0:
    levels = np.arange(-2, 25, 2)
else:
    levels = np.arange(0, 25, 2)

# potential energy
dPE = (- 0.6 - 0.6*ratio**2.0/(Ri/Rt) - ratio/(1.0+Ri/Rt) + 0.6*(1.0 + ratio)**2.0/(Rti/Rt))*GG*Mt**2.0/Rt
# initial kinetic energy
dKE =  ratio/(1 + Ri/Rt)*GG*Mt**2.0/Rt*vel**2.0


dx=90.0
dy=1
ap=dx/dy*0.618
 

if vimp_vesc==0:


    m0= 0.9232925
    m1= 0.076443
    h0= 0.73550873
    h1= 0.12824856
    h2= -0.02226599

    e0=0.18432
    e1=0.06338
    e2=0.00353
    e3=0.06389
    e4=0.10604
    e5=-0.18243
    e6=0.0279
    
else:

    e0=0.01935
    e1=0.04506
    e2=0.11079
    e3=0.17159
    e4=0.14955
    e5=-0.11511
    e6=-0.01596
    
    h0=0.67129416
    h1=0.35726835
    h2=-0.24558039

# Mantle_mass_model: Eq 6, Eq 8
# h_model: Eq 5, Eq 7
# IE_model: Eq 3

if vimp_vesc==0:
    Mantle_mass_model=m0*legendre(0,np.cos(ang)) +  m1*legendre(1,np.cos(ang)) 
    h_model=h0*legendre(0,ang) +  h1*legendre(1,ang)  +  h2*legendre(2,ang)
            
else:
    Mantle_mass_model = Mantle_mass_model_highV(targetmassfraction, ang)
    h_model = h0*legendre(0,np.cos(ang)) +  h1*legendre(1,np.cos(ang))  +  h2*legendre(2,np.cos(ang))

IE_model=  e0*legendre(0,np.cos(ang)) + e1*legendre(1,np.cos(ang)) + e2*legendre(2,np.cos(ang)) + e3*legendre(3,np.cos(ang))+  e4*legendre(4,np.cos(ang)) + e5*legendre(5,np.cos(ang)) + e6*legendre(6,np.cos(ang))

# computing the internal energy (Equation 3)
u_ave=h_model*IE_model*(dPE + dKE)/(0.70*Mantle_mass_model*(Mt + Mi))       
f_model = h_model*IE_model*(dPE + dKE)/(0.70*Mantle_mass_model*(Mt+Mi))/Lp

if f_model>1:
    f_model = 1.0

# -- reading coefficients from coef.txt
coef_read = [line.split() for line in open('coef.txt')]
theta=np.zeros(shape=(len(coef_read),len(coef_read[1])))

for m in range(0,len(coef_read)):
    for n in range(0,len(coef_read[1])):
        theta[m][n]=float(coef_read[m][n])
# -------------


melt_model=np.zeros(4)

Mplanet = Mantle_mass_model*(Mt + Mi) #planetary mass
rcore = compute_coreradius(Mplanet) #core radius

# grid spacing for calculating the magma ocean geometry
rr=np.linspace(rcore,1.0,30)
theta_angle=np.linspace(-np.pi, np.pi, 60)
nt=int(len(theta_angle))
nr=int(len(rr))

drr = (rr[1]-rr[0])/rr[len(rr)-1]
dangle = theta_angle[1]-theta_angle[0] 
du = np.zeros(shape=(nr,nt))
number = np.zeros(shape=(nr,nt))
rmax_meltpool_model = 1.0 # magma ocean depth. 1.0: no magma ocean, 0.0: the entire mantle is molten 
            
for m in range(0,nr):
    for n in range(0,nt):
        du[m][n] = create_model(theta[impact_angle.index(ang/np.pi*180)], rr[m], theta_angle[n])
                
du=du*u_ave

meltV=0.0     
totalV=0.0

for m in range(0,nr):
    for n in range(0,nt):
        dV=np.abs(np.pi*rr[m]**2.0*np.sin(theta_angle[n])*drr*dangle)
        totalV=totalV+dV
                
        if du[m][n] > Latentheat:
            meltV=meltV+dV
            # magma ocean depth is measured at psi = 0
            if rmax_meltpool_model > rr[m] and n == int(0.5*nt):
                rmax_meltpool_model = rr[m]  
        
melt_model=meltV/totalV



# --- estimating the magma ocean depth and corresponding pressure



#rmax_meltpool_model = max(rcore, rmax_meltpool_model)
Pmax_meltpool_model = compute_pressure(Mplanet,rmax_meltpool_model)

# assuming the same melt volume as the melt pool
rmax_global_model = (1.0-meltV/totalV)**0.333
Pmax_global_model = compute_pressure(Mplanet, rmax_global_model)

# assuming the conventional melt model (Eq 4)
rmax_conventional_model = (1.0-f_model)**0.333
Pmax_conventional_model =  compute_pressure(Mplanet, rmax_conventional_model)

print ("magma ocean depth and pressure for a melt pool model: " + str(rmax_meltpool_model) + ", " + str(Pmax_meltpool_model) + " GPa")
print ("magma ocean depth and pressure for a global magma ocean model: " + str(rmax_global_model) + ", " + str(Pmax_global_model) + " GPa")
print ("magma ocean depth and pressure for a conventional model: " + str(rmax_conventional_model) + ", " + str(Pmax_conventional_model) + " GPa")

       
du=du*1e-5


#  --- plotting options ----
border_w=2
font_xylabels=10
ax=0
dx=360.0
dy=0.5
ap=dx/dy*0.618

fig1=plt.figure(figsize=(10,6.128*2))
#ax = fig1.add_subplot(111,adjustable='box',aspect=ap)
ax = fig1.add_subplot(111,adjustable='box', polar=True)

CS=ax.contourf(theta_angle,rr,du,cmap=vik_map,vmin=5,vmax=15,levels=levels)

ax.set_rmax(1.0); ax.set_rmin(0.0)
ax.set_thetamin(-179.9);
ax.set_thetamax(180)
ax.set_xticks(np.array([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]) / 180. * np.pi)
ax.set_yticks([0.6, 0.8, 1.0])
ax.set_theta_zero_location('N')

cNorm = mpl.colors.Normalize(vmin=5, vmax=20)
ax3 = fig1.add_axes([0.27, 0.05, 0.45, 0.015])  # left, bottom, width, height (range 0 to 1)
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=vik_map, norm=cNorm, orientation='horizontal')
cb1.set_label('Internal Energy Gain ($10^5$ J/kg)')

fig1.savefig(outputfigurename)

#--------


