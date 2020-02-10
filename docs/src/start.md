# Getting started

```julia
using Slomo
mass_model = Halos.NFWModel()
tracer_model = SersicModel()
anisotropy_model = ConstantBetaModel()
jean_model = JeansModel(mass_model, tracer_model, anisotropy_model)
```
