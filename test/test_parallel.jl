using Distributed

@everywhere using Slomo

model = JeansModel(Halos.NFWModel(), SersicModel(), ConstantBetaModel(0.5))
R = 10 .^ collect(-1:0.1:1)
betas = [Dict(:beta => b) for b in collect(-0.5:0.1:0.5)]

s = sigma_los_parallel(model, R, betas)

