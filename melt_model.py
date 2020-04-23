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

# TODO
# fix h model
# magma ocean depth is fixed at psi = 0

class Model:

    def __init__(self, Mtotal=2.0, gamma=0.5, vel=2.0, entropy0=1100, impact_angle=90,
                 outputfigurename="output.eps", use_tex=False):
        self.Mmar = 6.39e23  # mass of Mars
        self.R0 = 1.5717e6  # impactor radius
        self.M0 = 6.39e22  # scaling coefficient
        self.a0 = 0.3412  # planetary mass-radius relationship
        self.a1 = -8.90e-3  # planetary mass-radius relationship
        self.a2 = 9.1442e-4  # planetary mass-radius relationship
        self.a3 = -7.4332e-5  # planetary mass-radius relationship
        self.GG = 6.67408e-11  # gravitational constant
        self.impact_angle_choices = [0.0, 30.0, 60.0, 90.0]  # choice of impact angle

        self.impact_angle = float(impact_angle)  # impactor impact angle with target


        # color maps
        self.cm_data = np.loadtxt("vik/vik.txt")
        self.vik_map = LinearSegmentedColormap.from_list('vik', self.cm_data)
        self.cm_data2 = np.loadtxt("turku/turku.txt")
        self.turku_map = LinearSegmentedColormap.from_list('turku', self.cm_data2)

        # font
        rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
        if use_tex:
            rc('text', usetex=True)  # this will not work if there in no native LaTeX installation

        # default values
        self.Mtotal = Mtotal  # total mass
        self.gamma = gamma  # impactor-to-total-mass ratio
        self.vel = vel  # impact velocity normalized by the escape velocity (this means that vimp = vesc(i.e. v_inf = 0))
        self.entropy0 = entropy0  # initial entropy (assuming adiabatic, dS/dr = 0)
        # default value of 1100 J/K/kg represents an adiabatic mantle with a surface temperature of 300K
        # default value of 3160 J/K/kg represents an adiabatic mantle with a surface temperature of 2000K
        self.check_model_integrity()
        self.entropyfile = 'rho_u_S{}.dat'.format(self.entropy0)
        self.outputfigurename = outputfigurename  # output figure name

        self.Mt = (1.0 - self.gamma) * self.Mtotal  # target mass

        self.Mi = self.gamma * self.Mtotal  # impactor mass

        self.EM = 5.2e6  # specific energy needed in order for melting to occur
        self.latent_heat = 7.18e5  # latent heat
        self.rho_P = [line.split() for line in open(
            self.entropyfile)]  # relationship between rho-P assuming S0=3160 J/K/kg. We also assume that rho-P structure is the same at S0=1100 J/K/kg.

        self.levels = np.arange(-2, 100, 2)
        self.vmin_value = 5
        self.vmax_value = 40
        # --- end of input data ---

        # calculating rho-P relationship of planetary interior assuming that the mantle has a constant entropy.
        self.rho_input = np.zeros(shape=(0, 0))  # density model
        self.P_input = np.zeros(shape=(0, 0))  # pressure model
        self.U_input = np.zeros(shape=(0, 0))  # internal energy model

        for m in range(1, len(self.rho_P)):  # reading from rho_u_SXXXX.dat (XXXX is either 1100 or 3160).
            self.rho_input = np.append(self.rho_input, float(self.rho_P[m][0]))
            self.P_input = np.append(self.P_input, float(self.rho_P[m][1]) * 1e9)  # converting from GPa to Pa.
            self.U_input = np.append(self.U_input, float(self.rho_P[m][2]))

        self.rho_P_function = interp1d(self.P_input, self.rho_input)  # generating interpolation

        self.rho_U_function = interp1d(self.rho_input, self.U_input)

        self.theta_angle = None
        self.rr = None
        self.du = None
        self.du_gain = None
        self.du_melt = None
        self.du_gain_melt = None

        self.check_model_integrity()

    def check_model_integrity(self):
        if float(self.entropy0) != 1100.0 and float(self.entropy0) != 3160.0:
            print("Please choose an entropy (entropy0) value of 1100 or 3160.")
            sys.exit(1)
        if self.impact_angle not in self.impact_angle_choices:
            print("Please choose an impact angle of 0, 30, 60 or 90 degrees!")
            sys.exit(1)

    # mass-radius relationship (see Section S.1.1. in our paper)
    def __radius(self, mass):
        lnMM0 = np.log(mass / self.M0)
        gamma = self.a0 + self.a1 * lnMM0 + self.a2 * lnMM0 ** 2 + self.a3 * lnMM0 ** 3
        return self.R0 * (mass / self.M0) ** gamma

    # Calculating a large Gamma(M). See Section S1.1. in our paper)
    def __gamma_calc(self, mass):
        lnMM0 = np.log(mass / self.M0)
        gamma = self.a0 + self.a1 * lnMM0 + self.a2 * lnMM0 ** 2 + self.a3 * lnMM0 ** 3
        return gamma

    # legendre polynomial functions, solutions to a legendre DE
    def __legendre(self, n, x):
        if n == 0:
            return 1
        elif n == 1:
            return x
        elif n == 2:
            return 0.5 * (3 * x ** 2.0 - 1.0)
        elif n == 3:
            return 0.5 * (5 * x ** 3.0 - 3 * x)
        elif n == 4:
            return 1.0 / 8.0 * (35 * x ** 4.0 - 30 * x ** 2.0 + 3)
        elif n == 5:
            return 1.0 / 8.0 * (63 * x ** 5.0 - 70 * x ** 3.0 - 15 * x)
        elif n == 6:
            return 1.0 / 8.0 * (231 * x ** 6.0 - 315 * x ** 4.0 + 105 * x ** 2.0 - 5)

    # mantle mass model for a high velocity impact. See Equation (8)

    def __Mantle_mass_model_highV(self, gamma, x):
        y = 0.0
        num = 1.0
        if (x <= np.pi / 6.0):
            y = num
        elif (x < np.pi / 3.0):
            y = -(num - gamma) / (np.pi / 6.0) * (x - np.pi / 6.0) + num
        elif (x <= np.pi / 2.0):
            y = gamma
        return y

    # r is the radius, x is the angle for the Legendre formulation. See Equation (13)
    def __create_model(self, theta, r, x):

        y_model = 0.0
        k = 0
        for i in range(0, 5):
            for j in range(0, 3):
                y_model = y_model + theta[k] * r ** (i - 2) * self.__legendre(j, np.cos(x))
                k = k + 1
        return y_model

    # pressure calculation. If an input pressure is smaller than 0, the minium density is returned
    def __compute_density(self, P):

        P = abs(P)
        if P < self.P_input[0]:
            return self.rho_input[0]
        else:
            return self.rho_P_function(P)

    # internal energy calculation
    def __compute_internal_energy(self, rho):
        U = self.rho_U_function(rho)
        if U < 0.0:
            return 0.0
        else:
            return U

    # given a planetary mass, this calculate density, internal energy, pressure, and mass as a function of planetary radius
    def __integrate_internal_energy(self, Mt):
        Rt = self.__radius(Mt)
        Rc = self.__compute_coreradius(Mt)

        press = 0.0
        r = Rt
        Mass = Mt

        rr = np.linspace(1.0, Rc, 20) * r
        dr = np.abs(rr[1] - rr[0])

        u = np.zeros(shape=(0, 0))
        P = np.zeros(shape=(0, 0))

        for i in range(0, len(rr)):
            rho = self.__compute_density(press)
            u = np.append(u, self.__compute_internal_energy(rho))
            press = press + rho * self.GG * Mass / rr[i] ** 2.0 * dr
            P = np.append(P, press)
            Mass = Mass - 4 * np.pi * rho * rr[i] ** 2.0 * dr

        return u, P, rr / r, Rt

    # compute pressure at radius = Rmelt
    def __compute_pressure(self, Mt, Rmelt):
        Rt = self.__radius(Mt)
        Rmelt = Rmelt * Rt

        if Rmelt == 1.0:
            return 0.0

        dr = (Rt - Rmelt) / 100.0


        P = 0.0
        r = Rt
        Mass = Mt

        while (r > Rmelt):
            rho = self.__compute_density(P)
            P = P + rho * self.GG * Mass / r ** 2.0 * dr
            Mass = Mass - 4 * np.pi * rho * r ** 2.0 * dr
            r = r - dr

        if Rmelt == 1.0:
            return 0.0
        else:
            return P * 1e-9

    # compute core radius
    def __compute_coreradius(self, Mt):
        Rt = self.__radius(Mt)
        dr = Rt / 1000.0
        P = 0.0
        r = Rt
        Mass = Mt
        CoreMass = 0.3 * Mt

        while (Mass > CoreMass):
            rho = self.__compute_density(P)
            P = P + rho * self.GG * Mass / r ** 2.0 * dr
            Mass = Mass - 4 * np.pi * rho * r ** 2.0 * dr
            r = r - dr

        return r / Rt

    # --- end of computing the structure of a planet

    # merging criteria from Genda et al 2012 (equation 16)
    def __v_cr(self, GammaG, theta):
        theta_G = 1 - np.sin(theta)
        return 2.43 * GammaG ** 2.0 * theta_G ** 2.5 - 0.0408 * GammaG + 1.86 * theta_G ** 2.50 + 1.08

    # plotting results
    def plot_model(self, save=False):
        if self.theta_angle is None:
            print("Please run the model before plotting.")
            sys.exit(1)

        fig1 = plt.figure(figsize=(10, 6.128 * 2))
        ax = [0, 0, 0, 0]
        ax[0] = fig1.add_subplot(221, adjustable='box', polar=True)
        ax[1] = fig1.add_subplot(222, adjustable='box', polar=True)
        ax[2] = fig1.add_subplot(223, adjustable='box', polar=True)
        ax[3] = fig1.add_subplot(224, adjustable='box', polar=True)

        level2 = [int(self.latent_heat * 1e-5)]  # this shows a white contour line for panel a

        CS = ax[0].contourf(self.theta_angle, self.rr, self.du_gain, cmap=self.vik_map, vmin=5, vmax=20,
                            levels=self.levels)
        CS2 = ax[0].contour(CS, levels=level2, colors="white", linewidths=1, vmin=5, vmax=15)
        CS3 = ax[1].contourf(self.theta_angle, self.rr, self.du_gain_melt, vmin=0, vmax=1.0, cmap=plt.cm.inferno)
        CS4 = ax[2].contourf(self.theta_angle, self.rr, self.du, cmap=self.turku_map, vmin=self.vmin_value,
                             vmax=self.vmax_value, levels=self.levels)
        CS6 = ax[3].contourf(self.theta_angle, self.rr, self.du_melt, vmin=0, vmax=1.0, cmap=plt.cm.inferno)

        # setting for polar plots
        for i in range(0, len(ax)):
            ax[i].set_rmax(1.0)
            ax[i].set_rmin(0.0)
            ax[i].set_thetamin(-179.9)
            ax[i].set_thetamax(180)
            ax[i].set_xticks(np.array([-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180]) / 180. * np.pi)
            ax[i].set_yticks([0.6, 0.8, 1.0])
            ax[i].tick_params(labelcolor='grey')
            ax[i].set_theta_zero_location('N')

        ax[0].text(np.pi / 5, 1.6, '(a) Internal Energy Gain', fontsize=15, color="black")
        ax[1].text(np.pi / 5, 1.6, '(b) Mantle Melt Mass Fraction', fontsize=15, color="black")
        ax[2].text(np.pi / 5, 1.6, '(c) Total Internal Energy', fontsize=15, color="black")
        ax[3].text(np.pi / 5, 1.6, '(d) Mantle Melt Mass Fraction', fontsize=15, color="black")
        ax[2].text(np.pi / 2, 0.4, '$S_0=$' + str(self.entropy0) + ' J/K/kg', fontsize=10, color="black")
        ax[3].text(np.pi / 2, 0.4, '$S_0=$' + str(self.entropy0) + ' J/K/kg', fontsize=10, color="black")

        # color bars
        cNorm = mpl.colors.Normalize(vmin=5, vmax=20)

        ax3 = fig1.add_axes([0.05, 0.05, 0.25, 0.015])
        cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=self.vik_map, norm=cNorm, orientation='horizontal')

        cNorm = mpl.colors.Normalize(vmin=self.vmin_value, vmax=self.vmax_value)
        ax4 = fig1.add_axes([0.375, 0.05, 0.25, 0.015])
        cb2 = mpl.colorbar.ColorbarBase(ax4, cmap=self.turku_map, norm=cNorm, orientation='horizontal')

        cNorm = mpl.colors.Normalize(vmin=0, vmax=1)
        ax5 = fig1.add_axes([0.7, 0.05, 0.25, 0.015])

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

        if save:
            fig1.savefig(self.outputfigurename)

    def run_model(self):


        """

        :return:
        """

        self.check_model_integrity()

        Mt = self.Mt * self.Mmar  # recalculate target mass
        Mi = self.Mi * self.Mmar  # recalculate impactor mass
        Rt = self.__radius(Mt)  # calculate radius of the target using mass-radius relationships
        Ri = self.__radius(Mi)  # calculate radius of the impactor using mass-radius relationships
        Rti = self.__radius(
            Mt + Mi)  # calculate radius of the target + impactor (perfectly merged body) using mass-radius relationships

        vesc = np.sqrt(2.0 * self.GG * (Mt + Mi) / (Rt + Ri))  # calculate escape velocity
        ratio = self.Mi / self.Mt  # impactor-to-target mass ratio (not the same as gamma!)
        ang = self.impact_angle / 180.0 * np.pi  # convert impact angle degrees to radians
        targetmassfraction = Mt / (Mt + Mi)  # mass fraction of the target relative to total mass (target + impactor)

        # potential energy
        dPE = (- 0.6 - 0.6 * ratio ** 2.0 / (Ri / Rt) - ratio / (1.0 + Ri / Rt) + 0.6 * (1.0 + ratio) ** 2.0 / (
                Rti / Rt)) * self.GG * Mt ** 2.0 / Rt
        # where the first two terms are gravitational binding energies of the target and impactor bodies
        # and the third term is the gravitational energy of the impactor body in the gravity potential of the target body
        # and the fourth term is the gravitational binding energy of the post-impact body under the assumption that the target and impactor perfectly merge

        # initial kinetic energy
        dKE = ratio / (1 + Ri / Rt) * self.GG * Mt ** 2.0 / Rt * self.vel ** 2.0

        # reading parameter coefficients for melt model
        parameter_list = [line.split() for line in open('parameter.txt')]
        para0 = np.array(parameter_list[0][:]).astype(np.float)  # parameters for vimp=vesc cases. See Table S.5
        para1 = np.array(parameter_list[1][:]).astype(np.float)  # parameters for vimp>1.1vesc cases. See Table S.6

        critical_velocity = self.__v_cr((Mt - Mi) / (Mt + Mi),
                                        ang)  # critical velocity (Genda et al 2012). See equation 16 in our paper

        if self.vel <= critical_velocity:  # merging
            Mantle_mass_model = para0[10] * self.__legendre(0, np.cos(ang)) + para0[11] * self.__legendre(1,
                                                                                                          np.cos(
                                                                                                              ang))  # mantle mass fitting model at vimp=vesc. See Equation 7
            h_model = para0[0] * self.__legendre(0, ang) + para0[1] * self.__legendre(1, ang) + para0[
                # mantle heating partitioning model at vimp=vesc. See Equation 6 (I have to fix this TO DO)
                2] * self.__legendre(2, ang)  # fitting model
        else:  # no merging
            Mantle_mass_model = self.__Mantle_mass_model_highV(targetmassfraction,
                                                               ang)  # mantle mass fitting model at vimp>1.1vesc. See Equation 8
            h_model = para1[0] * self.__legendre(0, np.cos(ang)) + para1[1] * self.__legendre(1, np.cos(ang)) + para1[
                # mantle heating partitioning model at vimp>1.1vesc.  See Equation 6
                2] * self.__legendre(2, np.cos(ang))  # fitting model

        # recall that the impact velocity is normalized to the escape velocity
        if self.vel <= 1.0:
            ee = para0[3:10]  # internal energy fitting model at vimp=vesc. See Equation 4 and Table S.5
        elif self.vel <= 1.1:
            ee = para0[3:10] + (self.vel - 1.0) / 0.1 * (para1[3:10] - para0[
                                                                       3:10])  # internal energy fitting model at vesc<vimp<1.1vesc. See Section 4.1
        else:
            ee = para1[3:10]  # internal energy fitting model at vimp>1.1vesc. See Equation 4 and Table S.6


        IE_model = ee[0] * self.__legendre(0, np.cos(ang)) + ee[1] * self.__legendre(1, np.cos(ang)) + ee[
            2] * self.__legendre(2, np.cos(ang)) + ee[
                       3] * self.__legendre(3, np.cos(ang)) + ee[4] * self.__legendre(4, np.cos(ang)) + ee[
                       5] * self.__legendre(5, np.cos(ang)) + ee[
                       6] * self.__legendre(6, np.cos(ang))  # internal energy model

        # computing the internal energy (Equation 3)
        u_ave = h_model * IE_model * (dPE + dKE) / (
                    0.70 * Mantle_mass_model * (Mt + Mi))  # this is Delta U = Delta IE. See quation 10 and Section 3.2
        f_model = h_model * IE_model * (dPE + dKE) / (
                    0.70 * Mantle_mass_model * (Mt + Mi)) / self.EM  # Mantle melt mass fraction. See Equation 9.


        if f_model > 1:
            f_model = 1.0

        # -- reading coefficients from coef.txt
        # Heat distribution model within mantle. See equation 13 and Table S.7
        coef_read = [line.split() for line in open('coef.txt')]
        theta = np.zeros(shape=(len(coef_read), len(coef_read[1])))

        for m in range(0, len(coef_read)):
            for n in range(0, len(coef_read[1])):
                theta[m][n] = float(coef_read[m][n])
        # -------------

        melt_model = np.zeros(4)

        Mplanet = Mantle_mass_model * (Mt + Mi)  # planetary mass
        rcore = self.__compute_coreradius(Mplanet)  # core radius

        u, P, r, rplanet = self.__integrate_internal_energy(
            Mplanet)  # calculating internal energy, pressure as a function of planetary radial distance as well as planetary radius

        r_U_function = interp1d(r, u)  # making a function of the internal energy as a function of planetary radius
        r_P_function = interp1d(r, P)  # making a function of the pressure  as a function of planetary radius

        # grid spacing for calculating the magma ocean geometry
        self.rr = np.linspace(rcore, 1.0,
                              30)  # radial spacing - this value 30 can be changed to a different value depending on the radial resolution you need

        self.theta_angle = np.linspace(-np.pi, np.pi,
                                       60)  # angle spacing (psi) - this value 60 can be changed to a different value depending on the angle resoultion you need
        nt = int(len(self.theta_angle))  # size of angle (psi) array
        nr = int(len(self.rr))  # size of radius array

        drr = (self.rr[1] - self.rr[0]) / self.rr[len(self.rr) - 1]  # radial grid size
        dangle = self.theta_angle[1] - self.theta_angle[0]  # angle grid size
        self.du = np.zeros(shape=(nr, nt))  # internal energy
        self.du_gain = np.zeros(shape=(nr, nt))  # internal energy gain
        number = np.zeros(shape=(nr, nt))
        self.du_melt = np.zeros(shape=(nr,
                                       nt))  # melt model w considering the initial temperature profile. this is 0 or 1; if a given part of the mantle is molten, this value is 1 otherwise 0
        self.du_gain_melt = np.zeros(shape=(nr,
                                            nt))  # melt model w/o considering the initial temperature profile. this is 0 or 1; if a given part of the mantle is molten, this value is 1 otherwise 0

        rmax_meltpool_model = 1.0  # magma ocean depth. 1.0: no magma ocean, 0.0: the entire mantle is molten

        # make the internal energy as a function of r

        for m in range(0, nr):
            for n in range(0, nt):
                self.du[m][n] = self.__create_model(theta[self.impact_angle_choices.index(ang / np.pi * 180)],
                                                    self.rr[m], self.theta_angle[n])

                self.du_gain[m][n] = self.du[m][n]

        du = self.du * u_ave
        du_gain = self.du_gain * u_ave

        for m in range(0, nr):
            for n in range(0, nt):
                du_initial = float(r_U_function(self.rr[m]))  # initial internal energy profile
                du[m][n] = du[m][n] + du_initial  # total internal energy


        for m in range(0, nr):
            for n in range(0, nt):
                Press = r_P_function(self.rr[m])
                Tmelt = (2500.0 + 26.0 * Press * 1e-9 - 0.052 * (Press * 1e-9) ** 2.0) * 1000.0

                if du[m][n] > Tmelt:
                    self.du_melt[m][n] = 1.0  # this portion of the mantle is molten
                else:
                    self.du_melt[m][n] = 0.0  # this portion of the mantle is NOT molten


        meltV = 0.0  # melt volume
        totalV = 0.0  # total volume

        for m in range(0, nr):
            for n in range(0, nt):
                dV = np.abs(np.pi * self.rr[m] ** 2.0 * np.sin(
                    self.theta_angle[n]) * drr * dangle)  # calculating an incremental  volume at this radius and angle
                totalV = totalV + dV

                if du_gain[m][n] > self.latent_heat:  # if du_gain is larger than latent heat
                    meltV = meltV + dV
                    self.du_gain_melt[m][n] = 1.0  # this part is considered molten

                    # magma ocean depth is measured at psi = 0
                    if rmax_meltpool_model > self.rr[m] and np.abs(
                            self.theta_angle[n]) < np.pi / 6.0:  # actually this should be a bit smaller (TO DO)
                        rmax_meltpool_model = self.rr[m]

        melt_model = meltV / totalV  # calculating melt volume


        # --- estimating the magma ocean depth and corresponding pressure

        # rmax_meltpool_model = max(rcore, rmax_meltpool_model)
        Pmax_meltpool_model = self.__compute_pressure(Mplanet, rmax_meltpool_model)

        # assuming the same melt volume as the melt pool
        rmax_global_model = (1.0 - meltV / totalV * (1.0 - rcore ** 3.0)) ** 0.333
        Pmax_global_model = self.__compute_pressure(Mplanet, rmax_global_model)

        # assuming the conventional melt model (Eq 4)
        rmax_conventional_model = (1.0 - f_model * (1.0 - rcore ** 3.0)) ** 0.333
        Pmax_conventional_model = self.__compute_pressure(Mplanet, rmax_conventional_model)

        print("magma ocean depth and pressure for a melt pool model: " + str(
            rplanet * 1e-3 * (1.0 - rmax_meltpool_model)) + " km, " + str(Pmax_meltpool_model) + " GPa")
        print("magma ocean depth and pressure for a global magma ocean model: " + str(
            rplanet * 1e-3 * (1.0 - rmax_global_model)) + " km, " + str(Pmax_global_model) + " GPa")
        print("magma ocean depth and pressure for a conventional model: " + str(
            rplanet * 1e-3 * (1.0 - rmax_conventional_model)) + " km, " + str(Pmax_conventional_model) + " GPa")


        self.du = du * 1e-5  # normalized by 10^5 J/kg
        self.du_gain = du_gain * 1e-5

        d = {
            "impact velocity": self.vel,
            "impact angle": self.impact_angle,
            "critical velocity": critical_velocity,
            "planetary mass": Mplanet,
            "core radius": rcore,
            "max depth (global model) (km)": rplanet * 1e-3 * (1.0 - rmax_global_model),
            "max pressure (global model)": Pmax_global_model,
            "max depth (conventional model) (km)": rplanet * 1e-3 * (1.0 - rmax_conventional_model),
            "max pressure (conventional model)": Pmax_conventional_model,
            "max depth (melt pool model) (km)": rplanet * 1e-3 * (1.0 - rmax_meltpool_model),
            "max pressure (melt pool model)": Pmax_meltpool_model,
            "melt fraction": f_model,
            "melt volume": meltV,
            "total volume": totalV,
            "internal energy": self.du,
            "average internal energy": u_ave,
            "internal energy gain": self.du_gain,
            "internal energy of the melt (considering initial temperature profile)": self.du_melt,
            "internal energy of the melt (not considering initial temperature profile)": self.du_gain_melt
        }

        return d

