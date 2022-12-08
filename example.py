from melt_model import Model

m = Model(Mtotal=8.9, gamma=0.09, vel=1.0, entropy0=1100, impact_angle=90, outputfigurename="output.png", use_tex=False)
data = m.run_model()
m.plot_model(save=True)

# if you want to produce an output, for example impact velocity, write
# print(data["impact velocity"])
