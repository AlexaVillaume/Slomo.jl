using Slomo

mass_model = NFWModel(50., 1e7)
tracer = Slomo.Tracer(Slomo.SersicModel(5.0, 4.0), Slomo.IsotropicModel())
R = collect(0.1:0.1:10)
parameters = Dict(:rhos => 1e7, :rs => 40.0, :Re => 10.0, :n => 4.0)
