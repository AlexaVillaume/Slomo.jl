using Slomo

model = JeansModel(Halos.NFWModel(), SersicModel(), ConstantBetaModel(0.5))
R = 10 .^ collect(-1:0.1:1)
betas = collect(-0.5:0.1:0.5)
params = [Dict(:beta => b) for b in betas]

s = sigma_los_parallel(model, R, params)

