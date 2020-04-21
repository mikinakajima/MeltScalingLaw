# Melt Scaling law 

This is a scaling law described in Nakajima et al. in preview. The current version is found here:
https://arxiv.org/abs/2004.04269

To use the scaling law, please take the following steps.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Make sure you have python3 on your machine and the following libraries.

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

To run a script, simply type

```
python  melt_model_Nakajima_et_al.py
```
This will generate the following output

```
magma ocean depth and pressure for a melt pool model: 2278.4837243704346 km, 86.15126263895684 GPa
magma ocean depth and pressure for a global magma ocean model: 898.1255099313564 km, 29.474287387562157 GPa
magma ocean depth and pressure for a conventional model: 303.5506251133121 km, 9.21756896904488 GPa
```

These are calculated depth (in km) and pressure (in GPa) at the base of a magma ocean depth for a melt pool model, global magma ocean model, and conventional model. <br />
Additionally, a figure called output.eps is generated <br />

![output.png](https://github.com/mikinakajima/MeltScalingLaw/blob/master/output.png)

You can modify the input.txt to change the input values.
```
Mtotal: Total mass normalized by Mars mass 
gamma: Imapctor-to-total mass ratio 
vel: Imapct velocity normalized by mutual escape velocity 
entropy0: Initial mantle entropy (before impact)
angle: Impact angle (0 is a head-on collision and 90 is the most oblique impact. Choose from 0, 30, 60, 90 degrees) 
outputfigurename: Name of the output figure name 
```


### Questions?

If you have questions, please feel free to reach out to Miki Nakajima (mnakajima@rochester.edu). If you use this scaling law, please be sure to cite our paper, Nakajima et al., in review.
https://arxiv.org/abs/2004.04269

