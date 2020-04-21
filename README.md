# Melt Scaling law 

This is a scaling law described in Nakajima et al. in review. To use the scaling law, please take the following steps.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

This code is based on Python 3 and may not be compatible with other versions of Python. Make sure you have python3 on your machine and the following libraries.

```
matplotlib
numpy
sys
scipy

```

Download color maps from Fabio Crameri's website by typing

```
./DownloadColormaps.sh 
```


### Running the script

The code is currently packaged as a single class, which can be easily imported within a different `.py` in the local directory.
An example of such a file, `example.py`, is shown below:

```
from melt_model_Nakajima_et_al import Model

m = Model(Mtotal=8.9, gamma=0.09, vel=1.0, entropy0=1100, impact_angle=90, outputfigurename="output.eps", use_tex=False)
m.run_model()
m.plot_model(save=False)
```

Here, `m = Model()` instantiates the model class (default class input parameters are shown for redundancy), `m.run_model()` executes the full series of calculations, and `m.plot_model()` plots the outputs.
This will generate the following output as well as a plot window.

```
magma ocean depth and pressure for a melt pool model: 2278.4837243704346 km, 86.15126263895684 GPa
magma ocean depth and pressure for a global magma ocean model: 898.1255099313564 km, 29.474287387562157 GPa
magma ocean depth and pressure for a conventional model: 303.5506251133121 km, 9.21756896904488 GPa
```

These are calculated depth (in km) and pressure (in GPa) at the base of a magma ocean depth for a melt pool model, global magma ocean model, and conventional model. <br />
Additionally, a figure is generated with the `Model.save_model()` function. <br />

![output.png](https://github.com/mikinakajima/MeltScalingLaw/blob/master/output.png)

### Arguments in the Model class

#### `Class Model(Mtotal=2.0, gamma=0.5, vel=2.0, entropy0=1100, impact_angle=90, outputfigurename="output.eps", use_tex=False)`
The main class of the melt model code. Default parameters shown.

- `Mtotal`: Total mass normalized by Mars mass.
- `gamma`: Imapctor-to-total mass ratio.
- `vel`: Imapct velocity normalized by mutual escape velocity.
- `entropy0`: Initial mantle entropy (before impact).
- `impact_angle`: Impact angle (0 is a head-on collision and 90 is the most oblique impact. Choose from 0, 30, 60, 90 degrees).
- `outputfigurename`: Name of the output figure name (relevant if `save=True` is set in `Model.plot_model(save=True)`).
- `use_tex`: Will render the figure with LaTeX if your computer has a native installation of LaTeX.

#### `Model.run_model()`

No input parameters required. Runs the full series of impact-induced melting calculations with objects stored in the `Model` class.

#### `Model.save_model(save=False)`
Produces the plot shown above.

- `save`: Saves the figure in the local directory under the name `outputfigurename` set in `class Model`.


### Questions?

If you have questions, please feel free to reach out to Miki Nakajima (mnakajima@rochester.edu). If you use this scaling law, please be sure to cite our paper, Nakajima et al., in review.


