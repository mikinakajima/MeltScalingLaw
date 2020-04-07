# MeltScalingLaw
Melt Scaling law from Nakajima et al. in rev.

This is a scaling law described in Nakajima et al. in review.

1. Download color maps from Fabio Crameri's website by typing  <br />
./DownloadColormaps.sh 

2. To change input parameters, open input.txt file <br />

Mtotal: Total mass normalized by Mars mass <br />
gamma: Imapctor-to-total mass ratio <br />
vel: Imapct velocity normalized by mutual escape velocity <br />
entropy0: Initial mantle entropy (before impact) <br />
angle: Impact angle (0 is a head-on collision and 90 is the most oblique impact. Choose from 0, 30, 60, 90 degrees) <br />
outputfigurename: Name of the output figure name <br />

3. To calculate the magma ocean pressures and magma ocean sizes, simply type  <br />

python melt_model_Nakajima_et_al.py  <br />

This will calculate three different pressures;  <br />
magma ocean depth and pressure for a melt pool model: 0.518000000000015, 54.70739524218776 GPa
magma ocean depth and pressure for a global magma ocean model: 0.0, 2.5580251604379935e+31 GPa
magma ocean depth and pressure for a conventional model: 0.8555730908432091, 14.211187973666862 GPa


