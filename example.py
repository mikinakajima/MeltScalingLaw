from melt_model import Model

m = Model(Mtotal=10.0, gamma=0.01, vel=1.0, entropy0=3160, impact_angle=30, outputfigurename="output.png", use_tex=False)
data = m.run_model()
m.plot_model(save=True)

# if you want to produce an output, for example impact velocity, write
print('fractional melt volume: ' + str(data["fractional melt volume"]))
print('fractional melt volume (max): ' + str(data["fractional melt volume max"]))
print('fractional melt volume (min): ' + str(data["fractional melt volume min"]))
print('fractional melt mass: '         + str(data["fractional melt mass"]))
print('fractional melt mass (max): '   + str(data["fractional melt mass max"]))
print('fractional melt mass (min): '   + str(data["fractional melt mass min"]))
print('total mantle mass: ' + str(data["total mantle mass"]))
print('total mantle mass: ' + str(data["total volume"]))
