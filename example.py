from melt_model_Nakajima_et_al import Model

m = Model(Mtotal=2.0, gamma=0.5, vel=2.0, entropy0=1100, outputfigurename="output.eps", use_tex=False)
m.run_model()
m.plot_model(save=False)
