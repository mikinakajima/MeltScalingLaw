import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import sys
import scipy.optimize as op
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rc
import matplotlib as mpl
from scipy.interpolate import interp1d


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

# color maps
cm_data = np.loadtxt("vik/vik.txt")
vik_map = LinearSegmentedColormap.from_list('vik', cm_data)
cm_data2 = np.loadtxt("turku/turku.txt")
turku_map = LinearSegmentedColormap.from_list('turku', cm_data2)

# font
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

# default values
Mtotal = 2.0
gamma = 0.5
vel = 2.0 #this means that vimp=vesc
entropy0 = 1100
outputfigurename = 'output.eps' #output figure name

# ---  input data  -- reading from input.txt ---- 
input_list = [line.split() for line in open('input.txt')]

for i in range(0,len(input_list)):
    if 'Mtotal' in input_list[i]:
        Mtotal = float(input_list[i][2])
    if 'gamma' in input_list[i]:
        gamma = float(input_list[i][2])
    if 'vel' in input_list[i]:
        vel = float(input_list[i][2])
    if 'outputfigurename' in input_list[i]:
        outputfigurename = input_list[i][2]
        
    if 'entropy0' in input_list[i]:
        entropy0 = float(input_list[i][2])
        if entropy0==1100:
            entropyfile = 'rho_u_S1100.dat'
            initial_S0 =  '$S_0=1100$ J/K/kg'
            
        elif entropy0==3160:
            entropyfile = 'rho_u_S3160.dat'
            initial_S0 =  '$S_0=3160$ J/K/kg'
    
        else:
            print('no entropy file exists with this entropy value -- choose either 1100 or 3160')
            sys.exit()
    if 'angle' in input_list[i]:
        ang = float(input_list[i][2])
        if ang not in impact_angle:
            print('Chosee impact angle among 0, 30, 60, and 90 degrees')
            sys.exit()            

Mt =   (1.0 - gamma)*Mtotal  # target mass in the Martian mass
Mi = gamma*Mtotal # impactor mass in the Martian mass
EM = 5.2e6 # specific energy for melting 
Latentheat= 7.18e5 # latent heat
rho_P = [line.split() for line in open(entropyfile)] #relationship between rho-P assuming S0=3160 J/K/kg. We also assume that rho-P structure is the same at S0=1100 J/K/kg.

levels = np.arange(-2, 100, 2)
vmin_value=5; vmax_value=40

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
    y_model = 0.0
    k = 0
    for i in range(0,5):
        for j in range(0,3):
            y_model = y_model + theta[k]*r**(i-2)*legendre(j, np.cos(x))
            k=k+1

    return y_model


#  ---- computing the structure of a planet ----

rho_input=np.zeros(shape=(0,0))
P_input=np.zeros(shape=(0,0))
U_input=np.zeros(shape=(0,0))

for m in range(1,len(rho_P)):
    rho_input = np.append(rho_input,float(rho_P[m][0]))
    P_input = np.append(P_input,float(rho_P[m][1])*1e9)
    U_input = np.append(U_input,float(rho_P[m][2]))

rho_P_function=interp1d(P_input, rho_input)
rho_U_function=interp1d(rho_input, U_input)

def compute_density(P):
    P=abs(P)
    if P < P_input[0]:
        return rho_input[0]
    else:
        return rho_P_function(P)

def compute_internal_energy(rho):
    U=rho_U_function(rho)
    if U < 0.0:
        return 0.0
    else:
        return U

def integrate_internal_energy(Mt):
    Rt=radius(Mt)
    Rc=compute_coreradius(Mt)

    press=0.0
    r=Rt
    Mass=Mt

    rr = np.linspace(1.0, Rc, 20)*r
    dr = np.abs(rr[1]-rr[0])

    u = np.zeros(shape=(0,0))
    P = np.zeros(shape=(0,0))

    for i in range(0, len(rr)):
        rho = compute_density(press)
        u = np.append(u, compute_internal_energy(rho))
        press = press + rho*GG*Mass/rr[i]**2.0*dr
        P = np.append(P, press)
        Mass=Mass-4*np.pi*rho*rr[i]**2.0*dr

    return u, P, rr/r, Rt



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
        P=P+rho*GG*Mass/r**2.0*dr
        Mass=Mass-4*np.pi*rho*r**2.0*dr
        r=r-dr

    if Rmelt==1.0:
        return 0.0
    else:
        return P*1e-9

def compute_coreradius(Mt):
    Rt=radius(Mt)

    dr=Rt/1000.0 #random number

    P=0.0
    r=Rt
    Mass=Mt
    CoreMass=0.3*Mt


    while (Mass > CoreMass):
        rho=compute_density(P)
        P=P+rho*GG*Mass/r**2.0*dr
        Mass=Mass-4*np.pi*rho*r**2.0*dr
        r=r-dr

    return r/Rt

# --- end of computing the structure of a planet


# merging criteria from Genda e tal 2012
def v_cr(GammaG, theta):
    theta_G=1-np.sin(theta)
    return 2.43*GammaG**2.0*theta_G**2.5-0.0408*GammaG+1.86*theta_G**2.50+1.08

Mt = Mt*Mmar
Mi = Mi*Mmar
Rt = radius(Mt) 
Ri = radius(Mi)
Rti = radius(Mt + Mi)
vesc = np.sqrt(2.0*GG*(Mt + Mi)/(Rt + Ri))
ratio = Mi/Mt
ang = ang/180.0*np.pi
targetmassfraction = Mt/(Mt + Mi)


# potential energy
dPE = (- 0.6 - 0.6*ratio**2.0/(Ri/Rt) - ratio/(1.0+Ri/Rt) + 0.6*(1.0 + ratio)**2.0/(Rti/Rt))*GG*Mt**2.0/Rt
# initial kinetic energy
dKE =  ratio/(1 + Ri/Rt)*GG*Mt**2.0/Rt*vel**2.0


dx=90.0
dy=1
ap=dx/dy*0.618


# reading parameter coefficients for melt model
parameter_list = [line.split() for line in open('parameter.txt')]
para0=np.array(parameter_list[0][:]).astype(np.float)
para1=np.array(parameter_list[1][:]).astype(np.float)

critical_velocity=v_cr((Mt-Mi)/(Mt+Mi), ang)

if vel <= critical_velocity: # merging
    Mantle_mass_model=para0[10]*legendre(0,np.cos(ang)) +  para0[11]*legendre(1,np.cos(ang))
    h_model=para0[0]*legendre(0,ang) +  para0[1]*legendre(1,ang)  +  para0[2]*legendre(2,ang)
else: #no merging
    Mantle_mass_model = Mantle_mass_model_highV(targetmassfraction, ang)    
    h_model = para1[0]*legendre(0,np.cos(ang)) +  para1[1]*legendre(1,np.cos(ang))  +  para1[2]*legendre(2,np.cos(ang))

if vel <=1.0:
    ee=para0[3:10]
elif vel <=1.1:
    ee=para0[3:10] + (vel-1.0)/0.1 * (para1[3:10]-para0[3:10])
else:
    ee=para1[3:10]

IE_model=  ee[0]*legendre(0,np.cos(ang)) + ee[1]*legendre(1,np.cos(ang)) + ee[2]*legendre(2,np.cos(ang)) + ee[3]*legendre(3,np.cos(ang)) +  ee[4]*legendre(4,np.cos(ang)) +  ee[5]*legendre(5,np.cos(ang)) +  ee[6]*legendre(6,np.cos(ang))


# computing the internal energy (Equation 3)
u_ave=h_model*IE_model*(dPE + dKE)/(0.70*Mantle_mass_model*(Mt + Mi))       
f_model = h_model*IE_model*(dPE + dKE)/(0.70*Mantle_mass_model*(Mt+Mi))/EM

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


u, P, r, rplanet = integrate_internal_energy(Mplanet)

r_U_function=interp1d(r, u)
r_P_function=interp1d(r, P)

# grid spacing for calculating the magma ocean geometry
rr=np.linspace(rcore,1.0,30)

theta_angle=np.linspace(-np.pi, np.pi, 60)
nt=int(len(theta_angle))
nr=int(len(rr))

drr = (rr[1]-rr[0])/rr[len(rr)-1]
dangle = theta_angle[1]-theta_angle[0] 
du = np.zeros(shape=(nr,nt))
du_gain = np.zeros(shape=(nr,nt))
number = np.zeros(shape=(nr,nt))
du_melt = np.zeros(shape=(nr,nt))
du_gain_melt = np.zeros(shape=(nr,nt))
rmax_meltpool_model = 1.0 # magma ocean depth. 1.0: no magma ocean, 0.0: the entire mantle is molten 


# make the internal energy as a function of r

for m in range(0,nr):
    for n in range(0,nt):
        du[m][n] = create_model(theta[impact_angle.index(ang/np.pi*180)], rr[m], theta_angle[n])
        du_gain[m][n] = du[m][n]

du=du*u_ave
du_gain=du_gain*u_ave


for m in range(0,nr):
    for n in range(0,nt):
        du_initial = float(r_U_function(rr[m]))
        du[m][n] = du[m][n] + du_initial


for m in range(0,nr):
    for n in range(0,nt):
        Press =  r_P_function(rr[m])
        Tmelt = (2500.0 + 26.0 * Press*1e-9 - 0.052 * (Press*1e-9) ** 2.0)*1000.0

        if du[m][n] > Tmelt:
            du_melt[m][n] = 1.0
        else:
            du_melt[m][n] = 0.0

meltV=0.0     
totalV=0.0

for m in range(0,nr):
    for n in range(0,nt):
        dV=np.abs(np.pi*rr[m]**2.0*np.sin(theta_angle[n])*drr*dangle)
        totalV=totalV+dV
        
        if du_gain[m][n] > Latentheat:
            meltV=meltV+dV
            du_gain_melt[m][n] = 1.0

            # magma ocean depth is measured at psi = 0
            if rmax_meltpool_model > rr[m] and np.abs(theta_angle[n])<np.pi/6.0:
                rmax_meltpool_model = rr[m]
        
melt_model=meltV/totalV


# --- estimating the magma ocean depth and corresponding pressure

#rmax_meltpool_model = max(rcore, rmax_meltpool_model)
Pmax_meltpool_model = compute_pressure(Mplanet,rmax_meltpool_model)

# assuming the same melt volume as the melt pool
rmax_global_model = (1.0-meltV/totalV*(1.0-rcore**3.0))**0.333
Pmax_global_model = compute_pressure(Mplanet, rmax_global_model)

# assuming the conventional melt model (Eq 4)
rmax_conventional_model = (1.0-f_model*(1.0-rcore**3.0))**0.333
Pmax_conventional_model =  compute_pressure(Mplanet, rmax_conventional_model)

print ("magma ocean depth and pressure for a melt pool model: " + str(rplanet*1e-3*(1.0-rmax_meltpool_model)) + " km, " + str(Pmax_meltpool_model) + " GPa")
print ("magma ocean depth and pressure for a global magma ocean model: " + str(rplanet*1e-3*(1.0-rmax_global_model)) + " km, " + str(Pmax_global_model) + " GPa")
print ("magma ocean depth and pressure for a conventional model: " + str(rplanet*1e-3*(1.0-rmax_conventional_model)) + " km, " + str(Pmax_conventional_model) + " GPa")

       
du=du*1e-5


#  --- plotting options ----
border_w=2
font_xylabels=10
ax=0
dx=360.0
dy=0.5
ap=dx/dy*0.618

fig1=plt.figure(figsize=(10,6.128*2))
ax= [0, 0, 0, 0]
ax[0] = fig1.add_subplot(221,adjustable='box', polar=True)
ax[1] = fig1.add_subplot(222,adjustable='box', polar=True)
ax[2] = fig1.add_subplot(223,adjustable='box', polar=True)
ax[3] = fig1.add_subplot(224,adjustable='box', polar=True)

level2=[int(Latentheat*1e-5)]
CS=ax[0].contourf(theta_angle,rr,du_gain*1e-5,cmap=vik_map,vmin=5,vmax=20,levels=levels)
CS2=ax[0].contour(CS,levels=level2, colors="white",linewidths=1, vmin=5, vmax=15)
#CS3=ax[1].contourf(theta_angle, rr, du_gain_melt, vmin=0, vmax=1.0,cmap='inferno')
CS3=ax[1].contourf(theta_angle, rr, du_gain_melt, vmin=0, vmax=1.0,cmap=plt.cm.inferno)

#you have to choose either of this to make the figure
CS4=ax[2].contourf(theta_angle,rr,du,cmap=turku_map, vmin=vmin_value,vmax=vmax_value,levels=levels)
CS6=ax[3].contourf(theta_angle, rr, du_melt, vmin=0, vmax=1.0, cmap=plt.cm.inferno)


for i in range(0,len(ax)):
    ax[i].set_rmax(1.0); ax[i].set_rmin(0.0)
    ax[i].set_thetamin(-179.9);
    ax[i].set_thetamax(180)
    ax[i].set_xticks(np.array([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]) / 180. * np.pi)
    ax[i].set_yticks([0.6, 0.8, 1.0])
    ax[i].tick_params(labelcolor='grey')
    ax[i].set_theta_zero_location('N')

   
    
ax[0].text(np.pi/5, 1.6, '(a) Internal Energy Gain', fontsize=15, color="black")
ax[1].text(np.pi/5, 1.6, '(b) Mantle Melt Mass Fraction', fontsize=15, color="black")
ax[2].text(np.pi/5, 1.6, '(c) Total Internal Energy', fontsize=15, color="black")
ax[3].text(np.pi/5, 1.6, '(d) Mantle Melt Mass Fraction', fontsize=15, color="black")
ax[2].text(np.pi/2, 0.4, initial_S0, fontsize=10, color="black")
ax[3].text(np.pi/2, 0.4, initial_S0, fontsize=10, color="black")


# color bars
cNorm = mpl.colors.Normalize(vmin=5, vmax=20)
ax3 = fig1.add_axes([0.05, 0.05, 0.25, 0.015])  # left, bottom, width, height (range 0 to 1)
cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=vik_map, norm=cNorm, orientation='horizontal')

cNorm = mpl.colors.Normalize(vmin=vmin_value, vmax=vmax_value)
ax4 = fig1.add_axes([0.375, 0.05, 0.25, 0.015])  # left, bottom, width, height (range 0 to 1)
cb2 = mpl.colorbar.ColorbarBase(ax4, cmap=turku_map, norm=cNorm, orientation='horizontal')

cNorm = mpl.colors.Normalize(vmin=0, vmax=1)
ax5 = fig1.add_axes([0.7, 0.05, 0.25, 0.015])  # left, bottom, width, height (range 0 to 1)
cb3 = mpl.colorbar.ColorbarBase(ax5, cmap=plt.cm.inferno, norm=cNorm, orientation='horizontal')


cb1.set_label('Internal Energy Gain ($10^5$ J/kg)')
cb2.set_label('Total Internal Energy ($10^5$ J/kg)')
cb3.set_label('Melt Fraction')

left = 0.05  # the left side of the subplots of the figure
right = 0.95  # the right side of the subplots of the figure
bottom = 0.1  # the bottom of the subplots of the figure
top = 0.9  # the top of the subplots of the figure
wspace = 0.1  # the amount of width reserved for blank space between subplots
hspace = 0.5  # the amount of height reserved for white space between subplots

fig1.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)


plt.show()
plt.close()

fig1.savefig(outputfigurename)

#--------


