
# Melt Scaling law 

This is a scaling law described in 

Nakajima, M., Golabek, G. J., Wuennemann, K., Rubie, D. C., Burger, C., Melosh, H. J., Jacobson, S. A., Manske, L., Hull, S. D. 2021. Scaling laws for the geometry of an impact-induced magma ocean. Earth and Planetary Science Letters, 568, 116983 [(link)](https://www.sciencedirect.com/science/article/pii/S0012821X21002429).

The major difference between the publication and this code is that 45 degree impacts are now taken into account. Additionally, there were bugs in the Legendre polynomials (sign mistakes), which are corrected in this version.

To use the scaling law, please take the following steps.

## Getting Started


### Prerequisites

This code is based on Python 3 and may not be compatible with other versions of Python. Make sure you have python3 on your machine and the following libraries.

```
matplotlib
numpy
sys
scipy

```

Download ScientificColorMaps7.zip from [here](https://zenodo.org/record/5501399#.YlSGbzfMIqs). Make sure to have lapaz, vik, and turku folders in the directory where you have the rest of the codes.


### Running the script

The code is currently packaged as a single class, which can be easily imported within a different `.py` in the local directory.
A single function executes the model and a different function plots the outputs. While advanced users can access and manipulate data tracked in the class through object-oriented relationships, the code is packaged in this way so that the models can be build with minimal user input.
An example of such a file, `example.py`, is shown below:

```
from melt_model import Model

m = Model(Mtotal=8.9, gamma=0.09, vel=1.0, entropy0=1100, impact_angle=90, outputfigurename="output.eps", use_tex=False)
m.run_model()
m.plot_model(save=False)
```

Here, `m = Model()` instantiates the model class (default class input parameters are shown for redundancy), `m.run_model()` executes the full series of calculations, and `m.plot_model()` plots the outputs.
This will generate the following output as well as a plot window.

```
mantle depth: 3007.93 km
mantle volume fraction: 0.13 (+0.15, -0.09)
magma ocean depth and pressure for a melt pool model: 726.05(+311.17, -414.89) km, 23.41(+10.86, -13.93) GPa
magma ocean depth and pressure for a global magma ocean model: 249.06(+311.23, -180.89) km, 7.52(+10.18, -5.51) GPa
magma ocean depth and pressure for a conventional model: 471.77(+347.14, -309.89) km, 14.73(+11.97, -9.9) GPa
```

These are calculated depth (in km) and pressure (in GPa) at the base of a magma ocean depth for a melt pool model, global magma ocean model, and conventional model. <br />
Additionally, a figure is generated with the `Model.save_model()` function. <br />

![output.png](https://github.com/mikinakajima/MeltScalingLaw/blob/master/output.png)

If you want to produce output in data (impact velocity, for example), type (in the example file)

```
print(data["impact velocity"])
```

The list of output is described in Line 656-681 in melt_model.py.

### Arguments in the Model class

#### `Class Model(Mtotal=2.0, gamma=0.5, vel=2.0, entropy0=1100, impact_angle=90, outputfigurename="output.eps", use_tex=False)`
The main class of the melt model code. Default input parameters are shown.

The input parameters are:
- `Mtotal`: Total mass normalized by Mars mass.
- `gamma`: Imapctor-to-total mass ratio.
- `vel`: Imapct velocity normalized by mutual escape velocity.
- `entropy0`: Initial mantle entropy (before impact).
- `impact_angle`: Impact angle (0 is a head-on collision and 90 is the most oblique impact. Choose from 0, 30, 45, 60, 90 degrees).
- `outputfigurename`: Name of the output figure name (relevant if `save=True` is set in `Model.plot_model(save=True)`).
- `use_tex`: Will render the figure with LaTeX if your computer has a native installation of LaTeX.

Additional model parameters stored in `class Model` that are available for access and updated after running `Model.run_model()` are:
- `Model.Mmar`: the mass of Mars
- `Model.R0`: the impactor radius
- `Model.M0`: the scaling coefficient
- `Model.a0`: coefficient in planetary mass-radius relationship
- `Model.a1`: coefficient in planetary mass-radius relationship
- `Model.a2`: coefficient in planetary mass-radius relationship
- `Model.a3`: coefficient in planetary mass-radius relationship
- `Model.GG`: the gravitational constant
- `Model.impact_angle_choices`: the available choices for impact angle, in degrees, for which the model is calibrated
- `Model.Mtotal`: total mass (impactor + target)
- `Model.gamma`: impact velocity normalized by the escape velocity (this means that vimp = vesc(i.e. v_inf = 0))
- `Model.Mt`: target mass
- `Model.Mi`: impactor mass
- `Model.EM`: specific energy needed in order for melting to occur
- `Model.latent_heat`: latent heat
- `Model.rho_P`: relationship between rho-P assuming S0=3160 J/K/kg. We also assume that rho-P structure is the same at S0=1100 J/K/kg.
- `Model.rho_input`: density model
- `Model.P_input`: pressure model
- `Model.U_input`: internal energy model
- `Model.rho_P_function`: 1-D pressure-density interpolation model
- `Model.rho_U_function`:  1-D density-energy interpolation model
- `Model.theta_angle`: angle spacing (psi) (this value 60 can be changed to a different value depending on the angle resolution you need)
- `Model.rr`: radial spacing (this value of 30 can be changed to a different value depending on the radial resolution you need)
- `Model.du`: internal energy (normalized by 10^5 J/kg)
- `Model.du_gain`: internal energy gain (normalized by 10^5 J/kg)
- `Model.du_melt`: melt model considering the initial temperature profile (this is 0 or 1; if a given part of the mantle is molten, this value is 1 otherwise 0)
- `Model.du_gain_melt`: melt model not considering the initial temperature profile (this is 0 or 1; if a given part of the mantle is molten, this value is 1 otherwise 0)

#### `Model.run_model()`

No input parameters required. Runs the full series of impact-induced melting calculations with objects stored in the `Model` class.
This function returns a dictionary object with numerical and model results.

The calculation process is as follows:
1. A check for correct user inputs is executed.
2. The target and impactor masses are calculated.
3. The radii of the target and impactor are calculated using planetary mass-radius relationships.
4. The radius of the perfectly-merged target and impactor using planetary mass-radius relationships.
5. The potential energy of the target-impactor system is calculated.
6. The initial kinetic energy of the target-impactor system is calculated.
7. Parameter coefficients from external files are read into the model.
8. The critical velocity is calculated (See Equation 16 in Genda et al. 2012).
9. The merging scenario between the impactor and target is modeled with the fitting model.
10. The no-merging scenario between the impactor and target is modeled with the fitting model.
11. The internal energy fitting models are applied.
12. The mantle melt mass fraction is calculated.
13. The heat distribution within the mantle is calculated.
14. The planetary mass and core radius are calculated.
15. The internal energy and pressure within the body as a function of radius is calculated.
16. The geometry of the magma ocean is produced.
17. The melt models are recalculated.
18. The internal energy profile as a function of depth is calculated.
19. The magma ocean depth and corresponding pressure is calculated.

#### `Model.save_model(save=False)`
Produces the plot shown above.

- `save`: Saves the figure in the local directory under the name `outputfigurename` set in `class Model`.


### Questions?

If you have questions, please feel free to reach out to Miki Nakajima (mnakajima@rochester.edu). If you use this scaling law, please be sure to cite our paper, Nakajima et al., 2021. 

Nakajima, M., Golabek, G. J., Wuennemann, K., Rubie, D. C., Burger, C., Melosh, H. J., Jacobson, S. A., Manske, L., Hull, S. D. 2021. Scaling laws for the geometry of an impact-induced magma ocean. Earth and Planetary Science Letters, 568, 116983 [(link)](https://www.sciencedirect.com/science/article/pii/S0012821X21002429)


