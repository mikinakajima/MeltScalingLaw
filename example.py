from melt_model_Nakajima_et_al import Model

m = Model(Mtotal=8.9, gamma=0.09, vel=1.0, entropy0=1100, impact_angle=90, outputfigurename="output.eps", use_tex=False)
m.run_model()
m.plot_model(save=False)
